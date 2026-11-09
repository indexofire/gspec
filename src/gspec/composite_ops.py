"""Running a ``gspec.composite/1`` document: screen, type, fuse.

The pipeline is a second-stage orchestration over the existing contracts:
each stage screens its database through :func:`gspec.screening.run_screen`
(the same path ``gspec screen`` uses, so per-database ``screen.json``
defaults — including the ``gspec.screen/2`` coverage model — apply), then
each file is designated through the database's own typing document. Nothing
about thresholds or phenotype logic lives here: the pipeline only decides
*which stage's accepted call wins* and keeps every stage's evidence.

Grouping: stages that name the same database share one screen run, so a
locus database carrying several tier schemes (``tier_a``/``tier_b``/
``tier_c``) is aligned once and typed once — three stages for the price of
one. Cluster stages honour per-stage ``min_cluster_cov``/``min_gene_cov``
overrides (one locus database, several gates) and therefore screen per
distinct gate pair.
"""

import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from gspec import db
from gspec.bundled import resolve_screen_datadir
from gspec.cluster import require_typing
from gspec.composite_models import CompositeDocument, CompositeStage
from gspec.engines import OutputFormat
from gspec.errors import DatabaseError, DependencyError
from gspec.formats.json import utc_timestamp
from gspec.report import ScreeningParams
from gspec.screen_defaults import screen_defaults
from gspec.screening import find_database, run_screen
from gspec.typing_gene import evaluate_gene_calls
from gspec.typing_input import read_typing_input
from gspec.typing_models import TypingDocument

COMPOSITE_RESULT_SCHEMA = "gspec.composite_result/1"


def _primer_path(stage: CompositeStage, document: CompositeDocument) -> Path:
    """A stage's primer file, resolved next to the pipeline document."""
    path = Path(stage.primers or "")
    if path.is_absolute() or document.source_dir is None:
        return path
    return document.source_dir / path


def _load_primers(stage: CompositeStage, document: CompositeDocument) -> tuple[Any, ...]:
    """Load a PCR stage's primer pairs through gseek (lazy, typed error)."""
    try:
        import gseek
    except ImportError as exc:  # pragma: no cover - dependency gate
        raise DependencyError(
            "pcr stages need the gseek package with ispcr support (install gseek>=0.11)",
            code="MISSING_DEPENDENCY",
            context={"stage": stage.name},
        ) from exc
    if not hasattr(gseek, "read_primers"):
        raise DependencyError(
            "the installed gseek has no in-silico PCR support (need gseek>=0.11)",
            code="MISSING_DEPENDENCY",
            context={"stage": stage.name, "gseek": getattr(gseek, "__version__", "?")},
        )
    pairs = tuple(gseek.read_primers(_primer_path(stage, document)))
    if stage.report_all_products:
        return pairs
    needed = {product.name for rule in stage.rules for product in rule.products}
    known = {pair.name for pair in pairs}
    missing = sorted(needed - known)
    if missing:
        raise DatabaseError(
            f"pcr stage {stage.name!r} names product(s) the primer set lacks: {', '.join(missing)}",
            code="COMPOSITE_MALFORMED",
            context={"stage": stage.name, "products": ", ".join(missing)},
        )
    return tuple(pair for pair in pairs if pair.name in needed)


def _primer_worker(
    payload: tuple[str, tuple[Any, ...], int, int],
) -> tuple[str, dict[str, list[int]]]:
    """One file's in-silico PCR, in a worker process.

    Module-level (not a closure) so the payload pickles: primer matching is a
    CPU-bound byte scan that does not release the GIL, so threads bought
    nothing here and processes are what makes ``--jobs`` pay.
    """
    import gseek

    path, primers, max_mismatches, max_edits = payload
    if not primers:
        # Nothing this file's call could need: skip the scan entirely.
        return path, {}
    result = gseek.ispcr(Path(path), primers, max_mismatches=max_mismatches, max_edits=max_edits)
    products: dict[str, list[int]] = {}
    for amplicon in result.amplicons:
        products.setdefault(amplicon.pair, []).append(amplicon.size)
    return path, products


def _needed_products(stage: CompositeStage, phenotype: str | None) -> set[str]:
    """The products this call could possibly need.

    A rule with ``when`` can only fire on a matching base type, so a genome
    called HS2 does not need the HS4A/HS4B/HS5_31 scans. Unconditional rules
    (the assay's own overrides) are always scanned.
    """
    base = _base_type(phenotype) if phenotype else None
    needed: set[str] = set()
    for rule in stage.rules:
        if rule.when and (base is None or base not in {_base_type(item) for item in rule.when}):
            continue
        needed.update(product.name for product in rule.products)
    return needed


def _primer_products(
    stage: CompositeStage,
    paths: list[Path],
    primers: tuple[Any, ...],
    jobs: int,
    *,
    phenotypes: dict[str, str | None] | None = None,
) -> dict[str, dict[str, list[int]]]:
    """Positive products per file: ``{file: {pair: [sizes]}}``.

    ``jobs > 1`` fans the scan out over worker PROCESSES — the primer search is
    a pure-Python byte scan holding the GIL, so a thread pool would serialise
    it (measured 1.01x); processes give the expected ~Nx.
    """
    payloads = []
    for path in paths:
        subset = primers
        if not stage.report_all_products and phenotypes is not None:
            needed = _needed_products(stage, phenotypes.get(str(path)))
            subset = tuple(pair for pair in primers if pair.name in needed)
        payloads.append((str(path), subset, stage.max_mismatches, stage.max_edits))
    if jobs <= 1 or len(payloads) <= 1:
        return dict(_primer_worker(payload) for payload in payloads)
    try:
        with ProcessPoolExecutor(max_workers=min(jobs, len(payloads))) as executor:
            return dict(executor.map(_primer_worker, payloads))
    except (OSError, RuntimeError):  # pragma: no cover - restricted sandboxes
        # A platform that forbids extra processes must not break the pipeline:
        # fall back to the sequential scan (same results, slower).
        return dict(_primer_worker(payload) for payload in payloads)


def _base_type(phenotype: str) -> str:
    """The call without its parenthetical bracket (``HS4c (… complex)`` → ``HS4c``)."""
    return phenotype.split(" (", 1)[0].strip()


def _primer_refine(
    stage: CompositeStage,
    current: tuple[str, str, float],
    products: dict[str, list[int]],
) -> tuple[str, str, float] | None:
    """Apply the stage's rules to one file's call; None when nothing applies."""
    phenotype = current[0]
    base = _base_type(phenotype)
    for rule in stage.rules:
        if rule.when and base not in {_base_type(when) for when in rule.when}:
            continue
        if not all(
            any(
                product.min_size <= size <= product.max_size
                for size in products.get(product.name, ())
            )
            for product in rule.products
        ):
            continue
        size = next(
            size
            for size in products[rule.products[0].name]
            if rule.products[0].min_size <= size <= rule.products[0].max_size
        )
        return (rule.set, stage.name, float(size))
    return None


@dataclass(frozen=True)
class Call:
    """One stage's designation of one file."""

    phenotype: str | None
    score: float
    confidence: str
    scheme: str


@dataclass(frozen=True)
class FileCall:
    """The fused result: the winning stage, its call, and every stage's raw call."""

    file: str
    phenotype: str

    tier: str
    stage: str
    score: float
    evidence: tuple[tuple[str, Call | str | None], ...]


def bundled_panel_root(document: CompositeDocument) -> Path | None:
    """The pipeline's own database directory (``db_root``), when it declares one."""
    if not document.db_root or document.source_dir is None:
        return None
    root = document.source_dir / document.db_root
    return root if root.is_dir() else None


def _stage_database(document: CompositeDocument, name: str, datadir: Path | None) -> db.Database:
    """Resolve a stage's database: the bundle's copy first, then the datadir."""
    root = bundled_panel_root(document)
    if root is not None:
        for database in db.discover_databases(root):
            if database.name == name:
                return database
    return find_database(resolve_screen_datadir(datadir, name), name, quiet=True)


def _stage_scheme(document: TypingDocument, stage: CompositeStage) -> str:
    """The scheme a stage reads; defaults to the document's only scheme."""
    names = [scheme.name for scheme in document.schemes]
    if stage.scheme is None:
        if len(names) != 1:
            raise DatabaseError(
                f"pipeline stage {stage.name!r} names no scheme, but database"
                f" {stage.db!r} declares {len(names)} ({', '.join(names)})",
                code="COMPOSITE_MALFORMED",
                context={"stage": stage.name, "db": stage.db, "schemes": names},
            )
        return names[0]
    if stage.scheme not in names:
        raise DatabaseError(
            f"pipeline stage {stage.name!r} names scheme {stage.scheme!r}, which database"
            f" {stage.db!r} does not declare (has: {', '.join(names)})",
            code="COMPOSITE_MALFORMED",
            context={"stage": stage.name, "db": stage.db, "scheme": stage.scheme},
        )
    return stage.scheme


def _gene_params(database: db.Database, defaults: Any, threads: int) -> ScreeningParams:
    """The gene-path parameters: mode-wide ``screen.json`` defaults + pipeline threads."""
    contigs = defaults.contigs if defaults is not None else None
    return ScreeningParams(
        db=database.name,
        minid=(contigs.minid if contigs and contigs.minid is not None else 80.0),
        mincov=(contigs.mincov if contigs and contigs.mincov is not None else 80.0),
        threads=threads,
        merge_fragments=bool(contigs.merge_fragments) if contigs else False,
        coverage=(contigs.coverage or "best_hsp") if contigs else "best_hsp",
        dedup=(contigs.dedup or "query") if contigs else "query",
        min_aligned_len=contigs.min_aligned_len if contigs else None,
    )


def _inputs_digest(paths: list[Path]) -> str:
    """A short digest of the input list: identity, size and mtime per file.

    The cache key must change when the *inputs* change, not just when their
    count does — otherwise a second single-file run silently reuses the first
    file's screening table.
    """
    digest = hashlib.sha256()
    for path in sorted(str(item) for item in paths):
        digest.update(path.encode("utf-8"))
        try:
            stat = Path(path).stat()
        except OSError:  # pragma: no cover - unreadable input fails later
            digest.update(b"|missing")
            continue
        digest.update(f"|{stat.st_size}|{stat.st_mtime_ns}".encode())
        digest.update(b"\n")
    return digest.hexdigest()[:16]


def _cache_path(cache_dir: Path, key: str, suffix: str) -> Path:
    """One stage's screen result on disk (keyed by database + parameters)."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / f"{key}.{suffix}"


def _screen_gene_table(
    database: db.Database,
    paths: list[Path],
    params: ScreeningParams,
    *,
    datadir: Path | None,
    jobs: int,
    quiet: bool,
    cache_dir: Path | None = None,
) -> str:
    """The stage's gspec screen TSV (gene path), built in-process.

    With ``cache_dir`` the table is reused across runs: screening is the
    expensive half of a pipeline, and re-tuning only the typing layer (or the
    fusion) should not re-align thousands of genomes.
    """
    key = (
        f"{database.name}-{params.minid:g}-{params.mincov:g}"
        f"-{params.coverage}-{params.dedup}-{params.min_aligned_len}"
        f"-{'m' if params.merge_fragments else 'n'}-{_inputs_digest(paths)}"
    )
    if cache_dir is not None:
        cached = _cache_path(cache_dir, key, "tsv")
        if cached.is_file():
            return cached.read_text(encoding="utf-8")
    text = run_screen(
        paths,
        database.name,
        datadir,
        params.minid,
        params.mincov,
        params.threads,
        jobs,
        None,
        True,  # quiet: the pipeline owns progress reporting
        False,  # noheader
        False,  # nopath: FILE stays the caller's label so stages key by input
        False,
        OutputFormat.tsv,
        merge_fragments=params.merge_fragments,
        coverage=params.coverage,
        dedup=params.dedup,
        min_aligned_len=params.min_aligned_len,
        database=database,
    )
    if cache_dir is not None:
        _cache_path(cache_dir, key, "tsv").write_text(text, encoding="utf-8")
    return text


def _gene_calls(
    database: db.Database,
    table_text: str,
    document: TypingDocument,
    scheme: str,
) -> dict[str, Call]:
    """Per-file designation of the gene-path table under one scheme."""
    table = read_typing_input([], stdin_text=table_text)
    calls: dict[str, Call] = {}
    for entry in table.files:
        evaluated = evaluate_gene_calls(entry.calls, document)
        call = evaluated.get(scheme)
        if call is None:
            continue
        calls[entry.file] = Call(
            phenotype=call.phenotype,
            score=call.score if call.phenotype is not None else 0.0,
            confidence=call.confidence,
            scheme=scheme,
        )
    return calls


def _cluster_cache_key(db_name: str, stage: CompositeStage, paths: list[Path]) -> str:
    """The cluster stage's cache key: database, both gates, and the inputs."""
    cluster_cov = stage.min_cluster_cov if stage.min_cluster_cov is not None else 96.0
    gene_cov = stage.min_gene_cov if stage.min_gene_cov is not None else 90.0
    return f"{db_name}-C{cluster_cov:g}-G90-g{gene_cov:g}-{_inputs_digest(paths)}"


def _cluster_calls(
    database: db.Database,
    paths: list[Path],
    stage: CompositeStage,
    scheme: str,
    *,
    datadir: Path | None,
    threads: int,
    jobs: int,
    cache_dir: Path | None = None,
) -> dict[str, Call]:
    """Per-file designation of a cluster database (integrated typing)."""
    key = _cluster_cache_key(database.name, stage, paths)
    text: str | None = None
    if cache_dir is not None:
        cached = _cache_path(cache_dir, key, "json")
        if cached.is_file():
            text = cached.read_text(encoding="utf-8")
    if text is None:
        text = run_screen(
            paths,
            database.name,
            datadir,
            80.0,
            80.0,  # the cluster engine's default gene-engine pair (guard default)
            threads,
            jobs,
            None,
            True,
            False,
            True,
            False,
            OutputFormat.json,
            database=database,
            min_gene_cov=stage.min_gene_cov if stage.min_gene_cov is not None else 90.0,
            min_gene_id=90.0,
            min_cluster_cov=(stage.min_cluster_cov if stage.min_cluster_cov is not None else 96.0),
        )
        if cache_dir is not None:
            _cache_path(cache_dir, key, "json").write_text(text, encoding="utf-8")
    document = json.loads(text)
    calls: dict[str, Call] = {}
    for entry in document.get("files", []):
        phenotypes = entry.get("phenotypes")
        if phenotypes is None:
            # gspec.cluster/1 (single-scheme document): the best-locus call.
            best = entry.get("best") or {}
            detail = best.get("phenotype_detail") or {}
            calls[entry["file"]] = Call(
                phenotype=best.get("phenotype"),
                score=detail.get("score", 0.0) or 0.0,
                confidence=detail.get("confidence", "low"),
                scheme=scheme,
            )
            continue
        call = phenotypes.get(scheme)
        if call is None:
            continue
        calls[entry["file"]] = Call(
            phenotype=call.get("phenotype"),
            score=call.get("score", 0.0) or 0.0,
            confidence=call.get("confidence", "low"),
            scheme=scheme,
        )
    return calls


def _apply_conflicts(
    document: CompositeDocument,
    winner: tuple[CompositeStage, Call],
    calls: dict[str, Call | None],
    stages: dict[str, CompositeStage],
) -> tuple[CompositeStage, Call]:
    """Swap a winning call for a declared override stage's call.

    Published assays carry documented cross-reactivity exceptions ("the HS41
    marker also amplifies in HS55 strains: trust the locus there"). The pipeline
    applies them between the cascade and any PCR refinement, and only when both
    competing calls were accepted on their own gates.
    """
    stage, call = winner
    base = _base_type(call.phenotype or "")
    for rule in document.conflicts:
        when_stage = rule.when.stage
        if when_stage is not None and when_stage != stage.name:
            continue
        if _base_type(rule.when.phenotype) != base:
            continue
        target = _base_type(rule.override.phenotype)
        candidates = (
            [rule.override.stage]
            if rule.override.stage is not None
            else list(stages)  # document order (dicts insert in stage order)
        )
        for name in candidates:
            if name is None:
                continue
            override_stage = stages.get(name)
            override_call = calls.get(name)
            if override_stage is None or not _accept(override_call, override_stage):
                continue
            assert override_call is not None  # _accept guarantees it
            if _base_type(override_call.phenotype or "") != target:
                continue
            return override_stage, override_call
    return winner


#: Per-process k-mer set cache: ``(db dir, set name) -> frozenset[int]``.
_KMER_SETS: dict[tuple[str, str], frozenset[int]] = {}


def _kmer_reference_sets(
    db_dir: Path, stage: CompositeStage
) -> tuple[int, dict[str, frozenset[int]]]:
    """Load the stage's k size and every declared set (cached per process)."""
    from gspec.kmerdb import load_kmer_set, read_kmerdb

    document = read_kmerdb(db_dir)
    if document is None:
        raise DatabaseError(
            f"kmer stage {stage.name!r} needs a k-mer substrate, but"
            f" database {stage.db!r} carries no kmer.json",
            code="KMERDB_MALFORMED",
            context={"stage": stage.name, "db": stage.db},
        )
    missing = sorted({rule.set for rule in stage.sets} - set(document.sets))
    if missing:
        raise DatabaseError(
            f"kmer stage {stage.name!r} names set(s) the database lacks: {', '.join(missing)}",
            code="COMPOSITE_MALFORMED",
            context={"stage": stage.name, "sets": ", ".join(missing)},
        )
    loaded: dict[str, frozenset[int]] = {}
    for rule in stage.sets:
        key = (str(db_dir), rule.set)
        cached = _KMER_SETS.get(key)
        if cached is None:
            cached = load_kmer_set(db_dir, rule.set)
            _KMER_SETS[key] = cached
        loaded[rule.set] = cached
    return document.k, loaded


def _kmer_worker(
    payload: tuple[str, str, tuple[tuple[str, str, float], ...], int, int],
) -> tuple[str, dict[str, tuple[str, float, float]]]:
    """One file's k-mer designation: ``(set -> phenotype, fraction, min_fraction)``.

    Module-level for the process pool: k-mer counting is a pure-Python scan, so
    ``--jobs`` must use processes here too.
    """
    import gspec.kmerdb as kmerdb

    path, db_dir, rules, k, stage_min = payload
    query = kmerdb.count_fasta_kmers(Path(path), k)
    calls: dict[str, tuple[str, float, float]] = {}
    for set_name, phenotype, min_fraction in rules:
        key = (db_dir, set_name)
        reference = _KMER_SETS.get(key)
        if reference is None:
            reference = kmerdb.load_kmer_set(Path(db_dir), set_name)
            _KMER_SETS[key] = reference
        calls[set_name] = (
            phenotype,
            kmerdb.set_hit_fraction(query, reference),
            min_fraction,
        )
    _ = stage_min
    return path, calls


def _kmer_calls(
    database: db.Database,
    stage: CompositeStage,
    paths: list[Path],
    *,
    jobs: int,
) -> dict[str, Call]:
    """Designate every file from its best k-mer set above the set's floor."""
    k, _ = _kmer_reference_sets(database.path, stage)
    rules = tuple((rule.set, rule.phenotype, rule.min_fraction) for rule in stage.sets)
    payloads = [(str(path), str(database.path), rules, k, 0.0) for path in paths]
    if jobs <= 1 or len(payloads) <= 1:
        results = [_kmer_worker(payload) for payload in payloads]
    else:
        try:
            with ProcessPoolExecutor(max_workers=min(jobs, len(payloads))) as executor:
                results = list(executor.map(_kmer_worker, payloads))
        except (OSError, RuntimeError):  # pragma: no cover - restricted sandboxes
            results = [_kmer_worker(payload) for payload in payloads]

    calls: dict[str, Call] = {}
    for path, scores in results:
        accepted = [
            (phenotype, fraction)
            for phenotype, fraction, floor in scores.values()
            if fraction >= floor
        ]
        if not accepted:
            continue
        phenotype, fraction = max(accepted, key=lambda item: item[1])
        calls[path] = Call(
            phenotype=phenotype,
            score=round(fraction, 4),
            confidence="high" if fraction >= 0.5 else "low",
            scheme="kmer",
        )
    return calls


def _accept(call: Call | None, stage: CompositeStage) -> bool:
    """Whether a stage's call is usable: designated, and clearing the stage gates.

    A stage that fell through to its own fallback (``confidence == "low"``, e.g.
    a marker scheme answering ``marker-negative``) designates nothing: accepting
    it would stop the pipeline on a layer that explicitly found no evidence, so
    fallbacks always fall through to the next stage.
    """
    if call is None or call.phenotype is None:
        return False
    if call.confidence == "low":
        return False
    if call.score < stage.min_score:
        return False
    return not (stage.require_confidence == "high" and call.confidence != "high")


def run_composite(
    paths: list[Path],
    document: CompositeDocument,
    *,
    datadir: Path | None,
    threads: int = 1,
    jobs: int = 1,
    quiet: bool = True,
    cache_dir: Path | None = None,
) -> list[FileCall]:
    """Screen every stage's database, designate every file, then fuse in stage
    order (first accepted call wins; every stage's call stays as evidence).

    ``kind="pcr"`` stages run after the database cascade: they refine the
    winning call (override, or bracket a complex type) from the primer products
    that amplified, and their evidence lands next to the database calls.
    """
    db_stages = tuple(stage for stage in document.stages if stage.kind == "db")
    kmer_stages = tuple(stage for stage in document.stages if stage.kind == "kmer")
    evidence_stages = tuple(stage for stage in document.stages if stage.kind in ("db", "kmer"))
    primer_stages = tuple(stage for stage in document.stages if stage.kind == "pcr")
    databases: dict[str, db.Database] = {}
    typing_docs: dict[str, TypingDocument] = {}
    schemes: dict[str, str] = {}
    evidence: dict[str, dict[str, Call]] = {stage.name: {} for stage in document.stages}

    # One screen per (db, cluster gate pair); gene stages of the same db share
    # a single table, cluster stages re-screen only when the gates differ.
    for stage in db_stages:
        database = databases.get(stage.db)
        if database is None:
            database = _stage_database(document, stage.db, datadir)
            databases[stage.db] = database
            typing_docs[stage.db] = require_typing(database)
        schemes[stage.name] = _stage_scheme(typing_docs[stage.db], stage)

    gene_tables: dict[str, str] = {}
    for stage in db_stages:
        database = databases[stage.db]
        scheme = schemes[stage.name]
        if database.kind == "cluster":
            evidence[stage.name] = _cluster_calls(
                database,
                paths,
                stage,
                scheme,
                datadir=datadir,
                threads=threads,
                jobs=jobs,
                cache_dir=cache_dir,
            )
            continue
        table = gene_tables.get(stage.db)
        if table is None:
            defaults = screen_defaults(database)
            table = _screen_gene_table(
                database,
                paths,
                _gene_params(database, defaults, threads),
                datadir=datadir,
                jobs=jobs,
                quiet=quiet,
                cache_dir=cache_dir,
            )
            gene_tables[stage.db] = table
        evidence[stage.name] = _gene_calls(database, table, typing_docs[stage.db], scheme)

    for stage in kmer_stages:
        database = _stage_database(document, stage.db, datadir)
        evidence[stage.name] = _kmer_calls(database, stage, paths, jobs=jobs)

    staged: dict[str, dict[str, Call | None]] = {}
    for path in paths:
        label = str(path)
        staged[label] = {stage.name: evidence[stage.name].get(label) for stage in evidence_stages}

    stage_by_name = {item.name: item for item in evidence_stages}
    cascade: dict[str, tuple[CompositeStage, Call] | None] = {}
    for label in (str(path) for path in paths):
        call_map = staged[label]
        winner: tuple[CompositeStage, Call] | None = None
        for stage in evidence_stages:
            call = call_map[stage.name]
            if _accept(call, stage):
                assert call is not None
                winner = _apply_conflicts(document, (stage, call), call_map, stage_by_name)
                break
        cascade[label] = winner

    # PCR products per stage: {file: {pair: [sizes]}}, scanning only the pairs
    # this file's call could need.
    primer_products: dict[str, dict[str, dict[str, list[int]]]] = {}
    for stage in primer_stages:
        primers = _load_primers(stage, document)
        phenotypes = {
            label: (winner[1].phenotype if winner is not None else None)
            for label, winner in cascade.items()
        }
        primer_products[stage.name] = _primer_products(
            stage, paths, primers, jobs, phenotypes=phenotypes
        )

    def primer_evidence(stage: CompositeStage, label: str) -> tuple[tuple[str, str], ...]:
        products = primer_products.get(stage.name, {}).get(label, {})
        return tuple((name, f"{sizes[0]}bp") for name, sizes in sorted(products.items()))

    results: list[FileCall] = []
    for label in (str(path) for path in paths):
        call_map = staged[label]
        winner = cascade[label]
        if winner is None:
            phenotype, tier, stage_name, score = document.fallback, "-", "-", 0.0
        else:
            stage, call = winner
            phenotype = call.phenotype or document.fallback
            tier, stage_name, score = stage.tier, stage.name, call.score

        # PCR refinement: the assay is authoritative for the call it names.
        for pcr_stage in primer_stages:
            refined = _primer_refine(
                pcr_stage,
                (phenotype, stage_name, score),
                primer_products.get(pcr_stage.name, {}).get(label, {}),
            )
            if refined is None:
                continue
            phenotype, stage_name, score = refined
            tier = pcr_stage.tier
            break

        evidence_items: tuple[tuple[str, Call | str | None], ...] = tuple(
            (name, call_map[name]) for name in call_map
        )
        for pcr_stage in primer_stages:
            products = primer_evidence(pcr_stage, label)
            evidence_items += (
                (pcr_stage.name, "+".join(f"{name}({size})" for name, size in products) or None),
            )
        results.append(
            FileCall(
                file=label,
                phenotype=phenotype,
                tier=tier,
                stage=stage_name,
                score=score,
                evidence=evidence_items,
            )
        )
    return results


def _evidence_cell(value: Call | str | None) -> str:
    """One stage's evidence as a table cell: a call, a PCR product list, or ``-``."""
    if isinstance(value, str):
        return value
    if value is None or not value.phenotype:
        return "-"
    return f"{value.phenotype}({value.score:.4f})"


def render_composite_tsv(
    results: list[FileCall], document: CompositeDocument, *, noheader: bool = False
) -> str:
    """One row per file: the fused call plus every stage's raw call as evidence."""
    lines: list[str] = []
    if not noheader:
        lines.append("FILE\tPIPELINE\tPHENOTYPE\tTIER\tSTAGE\tSCORE\tEVIDENCE")
    for result in results:
        evidence = "|".join(f"{name}={_evidence_cell(value)}" for name, value in result.evidence)
        lines.append(
            "\t".join(
                [
                    result.file,
                    document.name,
                    result.phenotype,
                    result.tier,
                    result.stage,
                    f"{result.score:.4f}",
                    evidence,
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _evidence_document(value: Call | str | None) -> Any:
    """One stage's evidence in the JSON result: call object, PCR string, or null."""
    if isinstance(value, str):
        return {"products": value.split("+") if value else []}
    if value is None or not value.phenotype:
        return None
    return {
        "phenotype": value.phenotype,
        "score": round(value.score, 4),
        "confidence": value.confidence,
        "scheme": value.scheme,
    }


def render_composite_json(
    results: list[FileCall], document: CompositeDocument, *, now: datetime | None = None
) -> str:
    """``gspec.composite_result/1`` — the fused call plus per-stage evidence."""
    stamp = utc_timestamp(now or datetime.now(UTC))
    payload = {
        "schema": COMPOSITE_RESULT_SCHEMA,
        "pipeline": document.name,
        "created_at": stamp,
        "stages": [
            {
                "name": stage.name,
                "kind": stage.kind,
                "db": stage.db,
                "scheme": stage.scheme,
                "tier": stage.tier,
                "min_score": stage.min_score,
                "require_confidence": stage.require_confidence,
                "primers": stage.primers,
            }
            for stage in document.stages
        ],
        "files": [
            {
                "file": result.file,
                "phenotype": result.phenotype,
                "tier": result.tier,
                "stage": result.stage,
                "score": round(result.score, 4),
                "evidence": {name: _evidence_document(value) for name, value in result.evidence},
            }
            for result in results
        ],
    }
    return json.dumps(payload, indent=2) + "\n"


__all__ = [
    "COMPOSITE_RESULT_SCHEMA",
    "Call",
    "FileCall",
    "render_composite_json",
    "render_composite_tsv",
    "run_composite",
]
