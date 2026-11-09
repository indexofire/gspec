"""The `db` composite stage on a GENE database: one screen table, then designation."""

from pathlib import Path

from gspec.composite.calls import Call
from gspec.composite.stages import cache_path, inputs_digest
from gspec.db import core as db
from gspec.engines import OutputFormat
from gspec.report import ScreeningParams
from gspec.screening.contigs import run_screen
from gspec.typing.gene import evaluate_gene_calls
from gspec.typing.input import read_typing_input
from gspec.typing.models import TypingDocument


def screen_gene_table(
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
        f"-{'m' if params.merge_fragments else 'n'}-{inputs_digest(paths)}"
    )
    if cache_dir is not None:
        cached = cache_path(cache_dir, key, "tsv")
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
        cache_path(cache_dir, key, "tsv").write_text(text, encoding="utf-8")
    return text


def gene_calls(
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
