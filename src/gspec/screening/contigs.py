"""The blastn contig-screening use-case plus helpers shared by both engines.

The minimap2 engine use-cases (--r1/--r2 reads, --aligner minimap2 assemblies)
live in screening_reads.py; the cluster engine's orchestration and db-kind
guards live in screening_cluster.py (stage-3 split, both under the LOC
ceiling); this module keeps the abricate-parity contig pipeline, the shared
database lookup, and the run_screen dispatcher. OutputFormat/AlignerEnum are
imported from gspec.engines so every screening surface sees the same enums.
"""

from collections.abc import Iterator
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

import click

from gspec.align.blast import ensure_blast, screen_file
from gspec.composite.models import DOCUMENT_FILENAME
from gspec.db import core as db
from gspec.db.bundled import find_bundled, materialize_bundled, resolve_screen_datadir
from gspec.engines import AlignerEnum, Emit, OutputFormat
from gspec.errors import DatabaseError, InputError, ensure_input_file, usage_fail
from gspec.formats.json import render_json
from gspec.formats.md import md_report_chunk, md_report_preamble
from gspec.formats.tsv import tsv_file_chunk, tsv_preamble
from gspec.paths import manifest_path
from gspec.progress import (
    ChunkSink,
    map_inputs,
    screen_progress,
    terminal_rule,
    warn_oversubscription,
)
from gspec.report import Report, ScreeningParams
from gspec.screening.batch import screen_batch
from gspec.screening.cluster import (
    reject_cluster_engine_flags,
    reject_gene_engine_flags,
    resolve_cluster_params,
    run_cluster_screen,
)


def _resolve_inputs(files: list[Path] | None, fofn: Path | None) -> list[Path]:
    """Input files: --fofn (lines stripped, empties dropped) REPLACES positionals."""
    if fofn is not None:
        if not fofn.is_file():
            raise InputError(
                f"--fofn file not found: {fofn}",
                code="INPUT_NOT_FOUND",
                context={"file": str(fofn)},
            )
        try:
            fofn_text = fofn.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError) as exc:
            raise InputError(
                f"--fofn file could not be read: {fofn}",
                code="INPUT_NOT_FOUND",
                context={"file": str(fofn)},
            ) from exc
        inputs = [Path(line.strip()) for line in fofn_text.splitlines() if line.strip()]
    elif files:
        inputs = list(files)
    else:
        usage_fail("no input files given (positional FILEs or --fofn)")
    for path in inputs:
        ensure_input_file(path)
    return inputs


def run_composite_screen(
    inputs: list[Path],
    db_name: str,
    document_path: Path,
    *,
    datadir: Path | None,
    threads: int,
    jobs: int,
    quiet: bool,
    noheader: bool,
    output_format: OutputFormat,
    emit: Emit | None,
    stdout_emit: bool,
) -> str:
    """Run a composite database and render its fused typing result.

    A composite carries ``composite.json`` instead of a panel: `screen` runs its
    stages (panels, k-mer sets, in-silico PCR) and renders the fused call — the
    TSV table or ``gspec.composite_result/1`` JSON.
    """
    from gspec.composite.models import read_composite
    from gspec.composite.render import render_composite_json, render_composite_tsv
    from gspec.composite.run import run_composite

    if output_format is OutputFormat.md:
        usage_fail("--format md is not available for an analysis scheme (use tsv or json)")
    document = read_composite(document_path)
    results = run_composite(
        inputs, document, datadir=datadir, threads=threads, jobs=jobs, quiet=quiet
    )
    if output_format is OutputFormat.json:
        text = render_composite_json(results, document)
    else:
        text = render_composite_tsv(results, document, noheader=noheader)
    if emit is not None:
        emit(text)
    return "" if emit is not None and not stdout_emit else text


def find_database(datadir: Path, name: str, *, quiet: bool = False) -> db.Database:
    """Look up a database by name under the datadir; unknown names list what
    exists.

    A BUNDLED database (gspec.db.bundled) that is absent from the datadir is
    auto-materialized into it first — deterministic, zero network — so the
    first screen works with no ``db fetch``. Already-materialized bundles
    short-circuit above; a bundle directory left without its certifying
    manifest (interrupted build) is rebuilt by the same path.
    """
    databases = db.discover_databases(datadir)
    for database in databases:
        if database.name == name:
            return database
    if find_bundled(name) is not None:
        if not manifest_path(datadir / name).is_file():
            try:
                materialize_bundled(name, datadir, quiet=quiet)
            except DatabaseError as exc:
                # A concurrent first-use screen may have won the writer lock
                # and published the same bundle between our discovery and
                # materialize. That is success for us: re-discover and use it.
                if exc.code != "DB_ALREADY_EXISTS":
                    raise
        # Re-discover even when the manifest was already there: a concurrent
        # materializer may have published between the first discovery and this
        # point (the manifest check then skips materialization).
        for database in db.discover_databases(datadir):
            if database.name == name:
                return database
    from gspec.composite.build import composite_names

    available = ", ".join(entry.name for entry in databases) or "(none)"
    # ``datadir`` here is the resolved databases directory; schemes live under it.
    schemes = composite_names(datadir)
    scheme_hint = f". Analysis schemes: {', '.join(schemes)}" if schemes else ""
    raise DatabaseError(
        f"Database {name} is not in {datadir}. Available: {available}{scheme_hint}",
        code="DATABASE_NOT_FOUND",
        context={"db": name, "datadir": str(datadir)},
    )


def _reject_database_flags_for_composite(
    name: str,
    *,
    minid: float,
    mincov: float,
    min_gene_cov: float,
    min_gene_id: float,
    min_cluster_cov: float,
    merge_fragments: bool,
    coverage: Literal["best_hsp", "union"],
    dedup: Literal["query", "subject"],
    min_aligned_len: int | None,
    aligner: object,
    detail: bool,
) -> None:
    """A composite configures its stages in its own document, not through db flags."""
    non_default = (
        (minid != 80.0, "--minid"),
        (mincov != 80.0, "--mincov"),
        (min_gene_cov != 90.0, "--min-gene-cov"),
        (min_gene_id != 90.0, "--min-gene-id"),
        (min_cluster_cov != 96.0, "--min-cluster-cov"),
        (merge_fragments, "--merge-fragments"),
        (coverage != "best_hsp", "the coverage model"),
        (dedup != "query", "the dedup key"),
        (min_aligned_len is not None, "--min-aligned-len"),
        (aligner is not None, "--aligner"),
        (detail, "--detail"),
    )
    offending = [flag for set_, flag in non_default if set_]
    if offending:
        usage_fail(
            f"{', '.join(offending)} apply to single-panel screening;"
            f" {name!r} is a composite database — configure its stages in its own"
            f" {DOCUMENT_FILENAME}"
        )


def run_screen(
    files: list[Path] | None,
    db_name: str,
    datadir: Path | None,
    minid: float,
    mincov: float,
    threads: int,
    jobs: int,
    fofn: Path | None,
    quiet: bool,
    noheader: bool,
    nopath: bool,
    debug: bool,
    output_format: OutputFormat,
    merge_fragments: bool = False,
    database: db.Database | None = None,
    coverage: Literal["best_hsp", "union"] = "best_hsp",
    allele_sketch: bool = False,
    dedup: Literal["query", "subject"] = "query",
    min_aligned_len: int | None = None,
    aligner: AlignerEnum | None = None,
    min_gene_cov: float = 90.0,
    min_gene_id: float = 90.0,
    min_cluster_cov: float = 96.0,
    detail: bool = False,
    emit: Emit | None = None,
    stdout_emit: bool = False,
) -> str:
    """Screen each input file in order and render; returns the full output
    for the caller to echo (MCP contract), and when ``emit`` is given also
    hands each rendered chunk to it as soon as the chunk exists — tsv/csv/md
    stream per completed file (input order, head-of-line under
    ``--jobs``: file i waits for 1..i; md leads with its static frontmatter,
    counts live in the JSON document only), while json is a single document
    emitted once at the end. Emission order is document order, so chunk
    concatenation is always byte-identical to the buffered return value.

    The per-run gates (blastn presence via ``ensure_blast``, dbtype via one
    ``blastdbcmd -info``) fire once up front, so MISSING_DEPENDENCY and
    DATABASE_NOT_INDEXED errors precede any "Processing:" stderr lines.
    ``merge_fragments`` is the opt-in cross-contig gene-fragment merge
    (gspec extension; default off keeps abricate parity). ``coverage``/
    ``dedup``/``min_aligned_len`` are the ``gspec.screen/2`` evidence-model
    knobs (defaults reproduce the abricate-parity single-HSP fold).

    A cluster-kind database dispatches to the minimap2 cluster engine
    (gspec.cluster) after the kind guards: gene-engine flags on a cluster
    db and cluster flags on a gene db are usage errors; the gene path below
    the dispatch is byte-identical to the pre-cluster behavior.
    """
    inputs = _resolve_inputs(files, fofn)
    # A caller that already resolved the database (composite stages) keeps it:
    # re-resolving by name would look for a panel under the datadir.
    if database is None:
        database = find_database(resolve_screen_datadir(datadir, db_name), db_name, quiet=quiet)
    if database.kind == "composite":
        _reject_database_flags_for_composite(
            db_name,
            minid=minid,
            mincov=mincov,
            min_gene_cov=min_gene_cov,
            min_gene_id=min_gene_id,
            min_cluster_cov=min_cluster_cov,
            merge_fragments=merge_fragments,
            coverage=coverage,
            dedup=dedup,
            min_aligned_len=min_aligned_len,
            aligner=aligner,
            detail=detail,
        )
        return run_composite_screen(
            inputs,
            db_name,
            database.document_path or (database.path / DOCUMENT_FILENAME),
            datadir=datadir,
            threads=threads,
            jobs=jobs,
            quiet=quiet,
            noheader=noheader,
            output_format=output_format,
            emit=emit,
            stdout_emit=stdout_emit,
        )
    if not 0.0 < minid <= 100.0:
        usage_fail(f"--minid must be in (0, 100]: got {minid}")
    if not 0.0 <= mincov <= 100.0:
        usage_fail(f"--mincov must be in [0, 100]: got {mincov}")
    if threads < 1:
        usage_fail(f"--threads must be >= 1: got {threads}")
    if jobs < 1:
        usage_fail(f"--jobs must be >= 1: got {jobs}")
    if database.kind == "cluster":
        reject_gene_engine_flags(minid, mincov, merge_fragments, aligner)
        return run_cluster_screen(
            inputs,
            database,
            resolve_cluster_params(database, min_gene_cov, min_gene_id, min_cluster_cov, threads),
            output_format=output_format,
            noheader=noheader,
            nopath=nopath,
            quiet=quiet,
            debug=debug,
            detail=detail,
            emit=emit,
            stdout_emit=stdout_emit,
            jobs=jobs,
        )
    if detail:
        usage_fail("--detail applies to cluster databases only (use --format json)")
    reject_cluster_engine_flags(min_gene_cov, min_gene_id, min_cluster_cov)
    params = ScreeningParams(
        db=db_name,
        minid=minid,
        mincov=mincov,
        threads=threads,
        merge_fragments=merge_fragments,
        coverage=coverage,
        dedup=dedup,
        min_aligned_len=min_aligned_len,
        allele_sketch=allele_sketch,
    )
    ensure_blast()
    dbtype = db.blast_db_info(database.sequences_path).dbtype
    warn_oversubscription(jobs, threads, quiet=quiet)

    def screen_one(path: Path) -> Report:
        if task_id is not None:
            # Bar active: its description carries the current file (notes suppressed).
            bar.describe(path.name)
        elif not quiet:
            click.echo(f"Processing: {path}", err=True)
        report = screen_file(
            path, database, params, dbtype=dbtype, debug=debug, merge_fragments=merge_fragments
        )
        if task_id is None and not quiet:
            click.echo(f"Found {len(report.hits)} genes in {path}", err=True)
        return report

    batch_reports: list[Report] | None = None
    batch_checked = False

    def iter_reports() -> Iterator[Report]:
        nonlocal batch_checked, batch_reports

        if not batch_checked:
            batch_checked = True
            batch_reports = screen_batch(
                inputs,
                database,
                params,
                dbtype=dbtype,
                debug=debug,
                merge_fragments=merge_fragments,
                jobs=jobs,
            )
        if batch_reports is not None:
            for report in batch_reports:
                bar.advance()
                yield report
            return
        yield from map_inputs(inputs, screen_one, jobs=jobs, on_done=bar.advance)

    with screen_progress(inputs, quiet=quiet) as (bar, task_id):
        sink = ChunkSink(emit, stdout=stdout_emit, bar=bar)

        if output_format is OutputFormat.json:
            sink(render_json(list(iter_reports()), params, now=datetime.now(UTC)))
        elif output_format is OutputFormat.md:
            # Static frontmatter leads (knowable before file 1), then one
            # section per file the moment it completes — exactly like tsv
            # chunks, head-of-line in input order under --jobs.
            sink(md_report_preamble(params, now=datetime.now(UTC)))
            for report in iter_reports():
                sink(md_report_chunk(report))
        else:
            as_csv = output_format is OutputFormat.csv
            sink(tsv_preamble(csv=as_csv, noheader=noheader))
            separator = terminal_rule() if stdout_emit and not as_csv else ""
            if separator and not noheader:
                sink(separator)
            first = True
            for report in iter_reports():
                if not first and separator:
                    sink(separator)
                first = False
                sink(tsv_file_chunk(report, csv=as_csv, nopath=nopath))
            if separator:
                sink(separator)
    return sink.text()
