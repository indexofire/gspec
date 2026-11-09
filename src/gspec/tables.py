"""Abricate-format table text plumbing shared by the table consumers.

Two surfaces read the same report tables — the summary matrix
(:mod:`gspec.summary`) and the typing input boundary
(:mod:`gspec.typing_table`) — and differ only in the error code and the
noun they attach to a malformed row. Separator detection, row iteration and
the malformed-row envelope live here once; each consumer keeps a thin
``_malformed`` wrapper naming its own document.
"""

import csv
import io
from collections.abc import Iterator
from pathlib import Path

from gspec.errors import InputError

#: The table label for stdin input — the conventional ``-`` marker, also how
#: stdin renders in a document's ``source`` array and in error listings.
STDIN_LABEL = "-"


def table_lines(text: str) -> list[str]:
    """Split like Perl's ``while (<$fh>)``: on \\n, no phantom final line."""
    if not text:
        return []
    parts = text.split("\n")
    if parts[-1] == "":
        parts.pop()
    return parts


def detect_separator(text: str) -> str:
    """[gspec-extension] auto-detect the table separator: tab if the first
    line has any, else comma, else tab (single-column degenerate)."""
    lines = table_lines(text)
    first = lines[0] if lines else ""
    if "\t" in first:
        return "\t"
    if "," in first:
        return ","
    return "\t"


def iter_rows(text: str, separator: str) -> Iterator[tuple[int, list[str]]]:
    """Physical line number + cells, using the csv module for comma tables.

    TSV keeps the abricate-parity raw split (no quoting convention); CSV uses
    ``csv.reader`` so a quoted comma/newline inside a field round-trips.
    """
    if separator == ",":
        reader = csv.reader(io.StringIO(text))
        for columns in reader:
            if columns:
                yield reader.line_num, columns
        return
    for line_number, line in enumerate(table_lines(text), start=1):
        yield line_number, line.split(separator)


def malformed(code: str, noun: str, label: str, line_number: int, detail: str) -> InputError:
    """One malformed data row: ``noun`` names the document, ``code`` its
    error family; the file and physical line ride the context."""
    return InputError(
        f"malformed {noun} row: {detail} ({label}:{line_number})",
        code=code,
        context={"file": label, "line": str(line_number)},
    )


def read_table_text(path: Path, *, code: str, noun: str) -> str:
    """One table's text; missing/unreadable/non-UTF-8 is a typed input error."""
    if not path.is_file():
        raise InputError(
            f"{noun} file not found or unreadable: {path}",
            code="INPUT_NOT_FOUND",
            context={"file": str(path)},
        )
    try:
        return path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise InputError(
            f"{noun} is not valid UTF-8: {path}",
            code=code,
            context={"file": str(path)},
        ) from exc


__all__ = [
    "STDIN_LABEL",
    "detect_separator",
    "iter_rows",
    "malformed",
    "read_table_text",
    "table_lines",
]
