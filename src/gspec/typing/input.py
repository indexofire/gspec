"""`gspec typing` input side: screen report tables → per-FILE gene calls.

Reads gspec/abricate contig screen tables and gspec reads tables (the
``#SAMPLE``/``BREADTH%`` shape), TSV or CSV, auto-detected per file, and
folds them into the evaluation input: one present
:class:`gspec.cluster.GeneCall` per (FILE, GENE) — the best row by
(%IDENTITY, %COVERAGE) descending, the first row winning ties, mirroring
the per-gene best-hit fold the inline engine applied to report hits
(gspec.typing.gene). Rows key on their FILE/#SAMPLE column, so one
designation can span several tables. Tables arrive as paths or as
already-read text (the piped-stdin table the CLI hands in, labeled ``-``).
Validation is typed at this boundary: a single DATABASE value across every
row of every file (mixed values are DATABASE_MISMATCH), at least one data
row (TYPING_NO_DATA — the database to evaluate is only knowable from rows),
and well-formed numeric columns (SCREEN_TABLE_MALFORMED). Table-shape
parsing lives in :mod:`gspec.typing.table`.
"""

from pathlib import Path

from pydantic import BaseModel

from gspec.cluster.core import GeneCall
from gspec.errors import InputError
from gspec.tables import STDIN_LABEL, read_table_text
from gspec.typing.table import ScreenRow, parse_table


class FileGeneCalls(BaseModel, frozen=True):
    """One screened FILE and its per-gene folded calls (fold order)."""

    file: str
    calls: dict[str, GeneCall]


class TypingInput(BaseModel, frozen=True):
    """The parsed designation input: the database every row screened
    against plus the per-FILE gene calls, in first-appearance order."""

    database: str
    files: tuple[FileGeneCalls, ...]


def _read_table(path: Path) -> str:
    """One result table's text; missing/unreadable/non-UTF-8 is an InputError."""
    return read_table_text(path, code="SCREEN_TABLE_MALFORMED", noun="screen result")


def _reject_cluster_table(label: str, text: str) -> None:
    """A cluster/2 TSV looks like a table but is not the gene-table contract;
    point the user at the JSON result that carries the scheme metadata."""
    header = text.splitlines()[0] if text.splitlines() else ""
    columns = {column.strip().upper() for column in header.split("\t")}
    if "SEROTYPE" in columns or any(column.endswith("_GROUP") for column in columns):
        raise InputError(
            "cluster screen tables cannot be re-typed (no scheme metadata);"
            " rerun screen with --format json and pass the JSON to typing",
            code="TYPING_CLUSTER_TABLE",
            context={"file": label},
        )


def _rows_from_screen_json(label: str, text: str) -> list[ScreenRow]:
    """Rows from a ``gspec.report/1`` screen document.

    The document carries the anchored allele-sketch score per hit, which the
    15-column table cannot: typing a screen JSON is therefore how the sketch
    tiebreak reaches the decision stage. The document is parsed with its own
    model, so a malformed file fails with the schema's own message.
    """
    from pydantic import ValidationError

    from gspec.formats.json import ReportDocument

    try:
        document = ReportDocument.model_validate_json(text)
    except ValidationError as exc:
        raise InputError(
            f"screen result JSON does not match gspec.report/1: {exc.error_count()} problem(s)",
            code="SCREEN_TABLE_MALFORMED",
            context={"file": label},
        ) from exc
    rows: list[ScreenRow] = []
    for entry in document.files:
        for hit in entry.hits:
            subject_len: int | None = None
            _, _, length_field = hit.coverage.rpartition("/")
            if length_field.isdigit():
                subject_len = int(length_field)
            rows.append(
                ScreenRow(
                    file=entry.file,
                    gene=hit.gene,
                    database=hit.database,
                    identity_pct=hit.identity_pct,
                    coverage_pct=hit.coverage_pct,
                    alignment_start_pct=None,
                    subject_len=subject_len,
                    sketch_pct=hit.sketch_pct,
                    phase_pct=hit.phase_pct,
                    contig=hit.sequence,
                    start=hit.start,
                    end=hit.end,
                )
            )
    return rows


def _looks_like_json(text: str) -> bool:
    """A screen JSON document starts with ``{`` once whitespace is stripped."""
    return text.lstrip().startswith("{")


def read_typing_input(paths: list[Path], stdin_text: str | None = None) -> TypingInput:
    """Parse and validate the designation input: every table's rows (the
    listed files plus, when ``stdin_text`` is given, the already-read stdin
    table) folded into per-FILE gene calls, with the single DATABASE they
    share."""
    rows: list[ScreenRow] = []
    labels: list[str] = []
    for path in paths:
        text = _read_table(path)
        labels.append(str(path))
        if _looks_like_json(text):
            rows.extend(_rows_from_screen_json(str(path), text))
            continue
        _reject_cluster_table(str(path), text)
        rows.extend(parse_table(str(path), text))
    if stdin_text is not None:
        labels.append(STDIN_LABEL)
        if _looks_like_json(stdin_text):
            rows.extend(_rows_from_screen_json(STDIN_LABEL, stdin_text))
        else:
            _reject_cluster_table(STDIN_LABEL, stdin_text)
            rows.extend(parse_table(STDIN_LABEL, stdin_text))
    if not rows:
        listing = ", ".join(labels)
        raise InputError(
            f"no screen result rows to type in: {listing}",
            code="TYPING_NO_DATA",
            context={"files": listing},
        )
    databases: list[str] = []
    for row in rows:
        if row.database not in databases:
            databases.append(row.database)
    if len(databases) > 1:
        raise InputError(
            "screen results mix databases; type each database's results separately:"
            f" {', '.join(databases)}",
            code="DATABASE_MISMATCH",
            context={"databases": ", ".join(databases)},
        )
    best: dict[str, dict[str, ScreenRow]] = {}
    for row in rows:
        current = best.setdefault(row.file, {}).get(row.gene)
        if current is None or (row.identity_pct, row.coverage_pct) > (
            current.identity_pct,
            current.coverage_pct,
        ):
            best[row.file][row.gene] = row
    return TypingInput(
        database=databases[0],
        files=tuple(
            FileGeneCalls(
                file=file,
                calls={
                    gene: GeneCall(
                        gene_id=gene,
                        # Coordinates when the source knows them (screen JSON):
                        # the phase guard needs the locus, not the offsets.
                        start=row.start or 0,
                        end=row.end or 0,
                        strand="+",
                        coverage_pct=row.coverage_pct,
                        identity_pct=row.identity_pct,
                        verdict="present",
                        alignment_start_pct=row.alignment_start_pct,
                        subject_len=row.subject_len,
                        sketch_pct=row.sketch_pct,
                        phase_pct=row.phase_pct,
                        contig=row.contig,
                    )
                    for gene, row in per_gene.items()
                },
            )
            for file, per_gene in best.items()
        ),
    )
