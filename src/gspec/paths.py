"""Where gspec keeps things on disk: the datadir and the database-directory layout.

Two related halves, one module. The datadir resolvers answer *where the
database root is* — ``--datadir`` > ``$GSPEC_DATADIR`` >
``$XDG_DATA_HOME/gspec/db`` > ``~/.local/share/gspec/db`` — while
:data:`MANIFEST_NAME`, :func:`manifest_path` and :data:`BUNDLED_ROOT` name the
files *inside* a database directory (and the wheel's bundled-database root).
Both are dependency-light path arithmetic (pathlib + :mod:`gspec.errors` only),
so every layer imports this module without a cycle.
"""

import os
from pathlib import Path
from typing import Literal

from gspec.errors import DatabaseError

#: Manifest filename written by 0.9 and later.
MANIFEST_NAME = "gspec-manifest.json"

#: The wheel's bundled-database root (``data/dbs``). Anchored here, not in the
#: reading module, so moving that module cannot silently relocate it.
BUNDLED_ROOT = Path(__file__).resolve().parent / "data" / "dbs"

ENV_DATADIR = "GSPEC_DATADIR"
ENV_XDG_DATA_HOME = "XDG_DATA_HOME"
DEFAULT_DATADIR = Path("~/.local/share/gspec/db")
DatadirSource = Literal["cli", "env", "xdg", "default"]


def manifest_path(db_dir: Path) -> Path:
    """The manifest file of a database directory."""
    return db_dir / MANIFEST_NAME


def datadir_source(cli_value: Path | None) -> tuple[Path, DatadirSource]:
    """Resolve the datadir location and where the choice came from.

    Precedence: ``--datadir`` > ``$GSPEC_DATADIR`` > ``$XDG_DATA_HOME/gspec/db``
    > ``~/.local/share/gspec/db``. An empty/relative
    ``$XDG_DATA_HOME`` is ignored per the XDG Base Directory specification
    (it must be an absolute path). The returned path is expanded and resolved
    but not required to exist.
    """
    if cli_value is not None:
        return Path(cli_value).expanduser().resolve(), "cli"
    env_value = os.environ.get(ENV_DATADIR)
    if env_value:
        return Path(env_value).expanduser().resolve(), "env"
    xdg_value = os.environ.get(ENV_XDG_DATA_HOME)
    if xdg_value:
        xdg_path = Path(xdg_value).expanduser()
        if xdg_path.is_absolute():
            return (xdg_path / "gspec" / "db").resolve(), "xdg"
    return DEFAULT_DATADIR.expanduser().resolve(), "default"


def datadir_path(cli_value: Path | None) -> Path:
    """Resolve the datadir LOCATION without requiring it to exist.

    Precedence and XDG handling live in :func:`datadir_source`; catalog/listing
    commands use this so a fresh install can be inspected before any database
    exists. Commands that need installed databases still call
    :func:`resolve_datadir`.
    """
    return datadir_source(cli_value)[0]


def resolve_datadir(cli_value: Path | None) -> Path:
    """Resolve the datadir and require it to exist.

    Expands ``~`` and resolves to an absolute path; raises DatabaseError when
    the resolved directory does not exist. Use :func:`datadir_path` for a
    non-validating resolution (catalog listing, bootstrap diagnostics).
    """
    datadir = datadir_path(cli_value)
    if not datadir.is_dir():
        raise DatabaseError(
            f"datadir does not exist: {datadir}",
            code="DATADIR_NOT_FOUND",
            context={"datadir": str(datadir)},
        )
    return datadir


def ensure_datadir(cli_value: Path | None) -> Path:
    """Resolve the datadir and mkdir it when absent (fresh-machine bootstrap).

    Write paths (``db fetch``, ``db build``, ``db list`` catalog bootstrap,
    ``db setup``) create a missing root; every path that actually consumes an
    installed database still requires it via :func:`resolve_datadir`.
    """
    root = datadir_path(cli_value)
    try:
        # New datadirs are private by default: databases affect screening
        # results, and shared multi-user access should be an explicit ACL/chmod
        # decision rather than an accidental world-readable directory.
        root.mkdir(mode=0o700, parents=True, exist_ok=True)
    except OSError as exc:
        raise DatabaseError(
            f"cannot create datadir: {root}",
            code="DATADIR_CREATE_FAILED",
            context={"datadir": str(root)},
        ) from exc
    return root
