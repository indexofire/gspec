"""Composite databases: ``gspec db build NAME --composite FILE``.

A composite is a database directory whose ``composite.json`` plans several
panels and methods instead of describing one panel. It is discovered like any
database (``db list`` shows ``KIND=composite``), addressed like any database
(``screen -d NAME``) and can carry its panels so it runs on a machine that has
none of them installed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from gspec.composite_build import (
    BUNDLE_DIRNAME,
    build_composite,
    composite_document_path,
    composite_names,
    composite_summary,
)
from gspec.composite_models import DOCUMENT_FILENAME, read_composite
from gspec.db import discover_databases
from gspec.errors import DatabaseError, InputError


def _database(datadir: Path, name: str = "panel_a") -> None:
    directory = datadir / name
    directory.mkdir(parents=True)
    (directory / "sequences").write_text(
        f">gspec|db={name}|gene=g|acc=a|func= x\nACGT\n", encoding="utf-8"
    )
    (directory / "gspec-manifest.json").write_text(
        json.dumps(
            {
                "schema": "gspec.manifest/1",
                "name": name,
                "source_urls": ["local"],
                "fetched_at": "2026-01-01T00:00:00Z",
                "sha256": "0" * 64,
                "n_records": 1,
                "dbtype": "nucl",
                "kind": "gene",
            }
        ),
        encoding="utf-8",
    )


def _document(tmp_path: Path, *, panel: str = "panel_a", primers: bool = False) -> Path:
    stages: list[dict[str, object]] = [{"name": "s1", "db": panel, "tier": "A"}]
    if primers:
        (tmp_path / "primers.tsv").write_text(
            "name\tforward\treverse\np\tACGT\tACGT\n", encoding="utf-8"
        )
        stages.append(
            {
                "name": "pcr",
                "kind": "pcr",
                "db": "primers",  # a label: PCR stages screen no database
                "primers": "primers.tsv",
                "rules": [{"products": [{"name": "p"}], "set": "X"}],
            }
        )
    payload = {
        "schema": "gspec.composite/1",
        "name": "demo",
        "fallback": "untypeable",
        "stages": stages,
    }
    path = tmp_path / "demo.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    root = tmp_path / "datadir"
    _database(root)
    return root


def test_build_creates_a_database_directory(tmp_path: Path, datadir: Path) -> None:
    installed = build_composite(_document(tmp_path), datadir)
    assert installed == composite_document_path("demo", datadir)
    assert installed.is_file()
    assert (installed.parent / "gspec-manifest.json").is_file()
    assert composite_names(datadir) == ("demo",)


def test_a_composite_is_discovered_as_a_database(tmp_path: Path, datadir: Path) -> None:
    build_composite(_document(tmp_path), datadir)
    kinds = {database.name: database.kind for database in discover_databases(datadir)}
    assert kinds == {"demo": "composite", "panel_a": "gene"}
    composite = next(d for d in discover_databases(datadir) if d.name == "demo")
    assert composite.document_path == composite.document_path


def test_build_rejects_a_missing_panel(tmp_path: Path, datadir: Path) -> None:
    with pytest.raises(DatabaseError) as exc_info:
        build_composite(_document(tmp_path, panel="not_installed"), datadir)
    assert exc_info.value.code == "DATABASE_NOT_FOUND"
    assert "not_installed" in str(exc_info.value)


def test_build_copies_relative_primer_sidecars(tmp_path: Path, datadir: Path) -> None:
    installed = build_composite(_document(tmp_path, primers=True), datadir)
    assert (installed.parent / "primers.tsv").is_file()


def test_rebuild_needs_force(tmp_path: Path, datadir: Path) -> None:
    build_composite(_document(tmp_path), datadir)
    with pytest.raises(DatabaseError) as exc_info:
        build_composite(_document(tmp_path), datadir)
    assert exc_info.value.code == "DB_ALREADY_EXISTS"
    build_composite(_document(tmp_path), datadir, force=True)


def test_bundle_copies_panels_and_pins_db_root(tmp_path: Path, datadir: Path) -> None:
    installed = build_composite(_document(tmp_path, primers=True), datadir, bundle=True)
    document = json.loads(installed.read_text(encoding="utf-8"))
    assert document["db_root"] == BUNDLE_DIRNAME
    assert (installed.parent / BUNDLE_DIRNAME / "panel_a" / "sequences").is_file()
    assert (installed.parent / "primers.tsv").is_file()

    moved = tmp_path / "elsewhere"
    moved.mkdir()
    import shutil

    shutil.copytree(installed.parent, moved / "demo")
    copied = read_composite(moved / "demo" / DOCUMENT_FILENAME)
    root = moved / "demo" / (copied.db_root or "")
    assert [database.name for database in discover_databases(root)] == ["panel_a"]


def test_bundle_rejects_a_missing_panel(tmp_path: Path, datadir: Path) -> None:
    with pytest.raises(DatabaseError) as exc_info:
        build_composite(_document(tmp_path, panel="absent"), datadir, bundle=True)
    assert exc_info.value.code == "DATABASE_NOT_FOUND"


def test_missing_primer_sidecar_is_a_typed_error(tmp_path: Path, datadir: Path) -> None:
    source = _document(tmp_path, primers=True)
    (tmp_path / "primers.tsv").unlink()
    with pytest.raises(InputError) as exc_info:
        build_composite(source, datadir)
    assert exc_info.value.code == "INPUT_NOT_FOUND"


def test_summary_lists_stages_and_panels(tmp_path: Path, datadir: Path) -> None:
    installed = build_composite(_document(tmp_path, primers=True), datadir)
    summary = composite_summary(read_composite(installed))
    assert summary["kind"] == "composite"
    assert summary["databases"] == ["panel_a"]
    assert [stage["kind"] for stage in summary["stages"]] == ["db", "pcr"]
