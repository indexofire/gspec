"""Typed JSON sidecar I/O: the one read/validate and write path.

Every database sidecar (manifest, typing, floors, features, profile, kmer,
info, pipeline) is a pydantic document read back through the same three
failure modes — a missing file, unreadable bytes, and a document that does
not validate — and written as indented JSON with a trailing newline. The
sidecar loaders previously restated all three per module; centralizing them
keeps the error code, the message wording, and the file-in-context shape
identical across documents (AGENTS.md §5 error envelope).

The caller names the document (``noun``) and the malformed code; the
missing-file code is always INPUT_NOT_FOUND (exit 5) while a malformed
document is a DatabaseError by default (exit 4) — ``error`` overrides that
for surfaces that classify input documents as user input.
"""

from pathlib import Path

from pydantic import BaseModel, ValidationError

from gspec.errors import DatabaseError, GspecError, ensure_input_file


def validation_reason(exc: ValidationError) -> str:
    """The validation failure rendered as one line (str(exc) is multiline)."""
    return " ".join(str(exc).split())


def read_json[M: BaseModel](
    model: type[M],
    path: Path,
    *,
    code: str,
    noun: str,
    error: type[GspecError] = DatabaseError,
) -> M:
    """Parse and validate one document.

    A missing/unreadable path is the typed INPUT_NOT_FOUND input error;
    unreadable bytes or a document that fails validation raise ``error``
    with ``code`` and the file in context.
    """
    ensure_input_file(path, noun)
    try:
        return model.model_validate_json(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError) as exc:
        raise error(
            f"could not read {noun} {path}: {exc}",
            code=code,
            context={"file": str(path)},
        ) from exc
    except ValidationError as exc:
        raise error(
            f"malformed {noun} {path}: {validation_reason(exc)}",
            code=code,
            context={"file": str(path)},
        ) from exc


def read_json_optional[M: BaseModel](
    model: type[M], path: Path, *, code: str, noun: str
) -> M | None:
    """:func:`read_json` for an optional sidecar: None when the file is absent."""
    return read_json(model, path, code=code, noun=noun) if path.is_file() else None


def write_json(document: BaseModel, path: Path, *, exclude_none: bool = False) -> None:
    """Write one document deterministically: indent 2, by alias (``schema``
    first), LF endings, UTF-8, trailing newline."""
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(document.model_dump_json(indent=2, by_alias=True, exclude_none=exclude_none))
        handle.write("\n")


__all__ = [
    "read_json",
    "read_json_optional",
    "validation_reason",
    "write_json",
]
