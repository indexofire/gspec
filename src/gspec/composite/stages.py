"""Stage resolution and shared screening parameters: which database a stage reads,
which scheme it designates with, and the cache keys that keep re-runs and
parallel workers cheap.
"""

import hashlib
from pathlib import Path
from typing import Any

from gspec.composite.models import CompositeDocument, CompositeStage
from gspec.db import core as db
from gspec.db.bundled import resolve_screen_datadir
from gspec.errors import DatabaseError
from gspec.report import ScreeningParams
from gspec.screening.contigs import find_database
from gspec.typing.models import TypingDocument


def bundled_panel_root(document: CompositeDocument) -> Path | None:
    """The pipeline's own database directory (``db_root``), when it declares one."""
    if not document.db_root or document.source_dir is None:
        return None
    root = document.source_dir / document.db_root
    return root if root.is_dir() else None


def stage_database(document: CompositeDocument, name: str, datadir: Path | None) -> db.Database:
    """Resolve a stage's database: the bundle's copy first, then the datadir."""
    root = bundled_panel_root(document)
    if root is not None:
        for database in db.discover_databases(root):
            if database.name == name:
                return database
    return find_database(resolve_screen_datadir(datadir, name), name, quiet=True)


def stage_scheme(document: TypingDocument, stage: CompositeStage) -> str:
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


def gene_params(database: db.Database, defaults: Any, threads: int) -> ScreeningParams:
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


def inputs_digest(paths: list[Path]) -> str:
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


def cache_path(cache_dir: Path, key: str, suffix: str) -> Path:
    """One stage's screen result on disk (keyed by database + parameters)."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / f"{key}.{suffix}"
