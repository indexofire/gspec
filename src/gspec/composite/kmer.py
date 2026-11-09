"""The `kmer` composite stage: canonical-set counting mapped onto phenotypes."""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from gspec.composite.calls import Call
from gspec.composite.models import CompositeStage
from gspec.db import core as db
from gspec.errors import DatabaseError


def kmer_reference_sets(
    db_dir: Path, stage: CompositeStage
) -> tuple[int, dict[str, frozenset[int]]]:
    """Load the stage's k size and every declared set (cached per process)."""
    from gspec.cluster.kmerdb import load_kmer_set, read_kmerdb

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
        cached = KMER_SETS.get(key)
        if cached is None:
            cached = load_kmer_set(db_dir, rule.set)
            KMER_SETS[key] = cached
        loaded[rule.set] = cached
    return document.k, loaded


def kmer_worker(
    payload: tuple[str, str, tuple[tuple[str, str, float], ...], int, int],
) -> tuple[str, dict[str, tuple[str, float, float]]]:
    """One file's k-mer designation: ``(set -> phenotype, fraction, min_fraction)``.

    Module-level for the process pool: k-mer counting is a pure-Python scan, so
    ``--jobs`` must use processes here too.
    """
    import gspec.cluster.kmerdb as kmerdb

    path, db_dir, rules, k, stage_min = payload
    query = kmerdb.count_fasta_kmers(Path(path), k)
    calls: dict[str, tuple[str, float, float]] = {}
    for set_name, phenotype, min_fraction in rules:
        key = (db_dir, set_name)
        reference = KMER_SETS.get(key)
        if reference is None:
            reference = kmerdb.load_kmer_set(Path(db_dir), set_name)
            KMER_SETS[key] = reference
        calls[set_name] = (
            phenotype,
            kmerdb.set_hit_fraction(query, reference),
            min_fraction,
        )
    _ = stage_min
    return path, calls


def kmer_calls(
    database: db.Database,
    stage: CompositeStage,
    paths: list[Path],
    *,
    jobs: int,
) -> dict[str, Call]:
    """Designate every file from its best k-mer set above the set's floor."""
    k, _ = kmer_reference_sets(database.path, stage)
    rules = tuple((rule.set, rule.phenotype, rule.min_fraction) for rule in stage.sets)
    payloads = [(str(path), str(database.path), rules, k, 0.0) for path in paths]
    if jobs <= 1 or len(payloads) <= 1:
        results = [kmer_worker(payload) for payload in payloads]
    else:
        try:
            with ProcessPoolExecutor(max_workers=min(jobs, len(payloads))) as executor:
                results = list(executor.map(kmer_worker, payloads))
        except (OSError, RuntimeError):  # pragma: no cover - restricted sandboxes
            results = [kmer_worker(payload) for payload in payloads]

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


KMER_SETS: dict[tuple[str, str], frozenset[int]] = {}
