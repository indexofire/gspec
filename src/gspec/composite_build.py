"""Building a composite database: several panels and methods in one artifact.

A composite database is a database directory whose ``composite.json`` describes
an analysis plan instead of a single panel:

    <datadir>/cje_penner/composite.json          # the document
    <datadir>/cje_penner/cje_penner_primers.json # relative sidecars travel along
    <datadir>/cje_penner/dbs/cje_hs_marker/      # optional bundled panels
    <datadir>/cje_penner/dbs/cje_cps_locus/

It is discovered like any other database (``gspec db list`` shows it with
``KIND=composite``), addressed like any other database (``gspec screen -d
cje_penner``), and runs all of its stages — blastn panels, minimap2 loci, k-mer
sets, in-silico PCR — before emitting the fused typing result. There is no
separate command and no separate namespace: composite-ness is a property of the
database directory.

``--bundle`` copies every panel the stages name into ``dbs/<panel>`` and pins
``db_root``, so the directory runs on a machine where none of those databases is
installed. Without it the stages resolve against the datadir as usual.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

from gspec import config, db
from gspec.composite_models import DOCUMENT_FILENAME, CompositeDocument, read_composite
from gspec.errors import DatabaseError, InputError

BUNDLE_DIRNAME = "dbs"
MANIFEST_FILENAME = "gspec-manifest.json"


def composite_document_path(name: str, datadir: Path | None) -> Path:
    """``<datadir>/<name>/composite.json`` (may not exist)."""
    return config.datadir_path(datadir) / name / DOCUMENT_FILENAME


def composite_names(datadir: Path | None) -> tuple[str, ...]:
    """Installed composite database names, sorted."""
    root = config.datadir_path(datadir)
    if not root.is_dir():
        return ()
    return tuple(
        sorted(
            child.name
            for child in root.iterdir()
            if child.is_dir() and (child / DOCUMENT_FILENAME).is_file()
        )
    )


def _panel_names(document: CompositeDocument) -> set[str]:
    """The databases a composite's stages screen (PCR stages name no database)."""
    return {stage.db for stage in document.stages if stage.kind in ("db", "kmer")}


def _installed_databases(datadir: Path | None) -> dict[str, Path]:
    """Installed panel databases: name → directory."""
    root = config.datadir_path(datadir)
    return {database.name: database.path for database in db.discover_databases(root)}


def _write_manifest(destination: Path, name: str, document: CompositeDocument) -> None:
    """A manifest so ``db list``/``db info``/rename treat this as a database."""
    manifest = {
        "schema": "gspec.manifest/1",
        "name": name,
        "source_urls": ["local"],
        "fetched_at": "2026-01-01T00:00:00Z",
        "sha256": "0" * 64,
        "n_records": len(document.stages),
        "dbtype": "nucl",
        "kind": "composite",
    }
    (destination / MANIFEST_FILENAME).write_text(
        json.dumps(manifest, indent=1) + "\n", encoding="utf-8"
    )


def build_composite(
    source: Path,
    datadir: Path | None,
    *,
    name: str | None = None,
    force: bool = False,
    bundle: bool = False,
) -> Path:
    """Create ``<datadir>/<name>/`` from a ``gspec.composite/1`` document.

    Returns the installed document path. ``force`` replaces an existing
    database of that name. ``bundle`` copies the panels into ``dbs/`` and pins
    ``db_root`` so the directory is self-contained.
    """
    document = read_composite(source)
    composite_name = name or document.name
    if not composite_name or "/" in composite_name or composite_name.startswith("."):
        raise DatabaseError(
            f"invalid composite name {composite_name!r}",
            code="COMPOSITE_MALFORMED",
            context={"composite": composite_name},
        )
    installed = _installed_databases(datadir)
    unknown = sorted(_panel_names(document) - set(installed))
    if unknown:
        raise DatabaseError(
            f"composite {composite_name!r} names database(s) that are not installed:"
            f" {', '.join(unknown)} (install them first with `gspec db setup`/`gspec db fetch`)",
            code="DATABASE_NOT_FOUND",
            context={"composite": composite_name, "databases": ", ".join(unknown)},
        )

    destination = config.datadir_path(datadir) / composite_name
    if destination.exists():
        if not force:
            raise DatabaseError(
                f"database {composite_name!r} already exists at {destination}"
                " (pass --force to replace)",
                code="DB_ALREADY_EXISTS",
                context={"name": composite_name, "destination": str(destination)},
            )
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)

    installed_document = destination / DOCUMENT_FILENAME
    if bundle:
        bundled = document.model_copy(update={"db_root": BUNDLE_DIRNAME})
        installed_document.write_text(
            json.dumps(bundled.model_dump(by_alias=True, exclude_none=True), indent=1) + "\n",
            encoding="utf-8",
        )
        # Stages share panels (three tiers of one database): copy each once.
        for panel in sorted(_panel_names(document)):
            shutil.copytree(installed[panel], destination / BUNDLE_DIRNAME / panel)
    else:
        shutil.copy2(source, installed_document)

    # Relative sidecars (primer sets) travel with the document: the composite's
    # own resolution rule is "next to the document", so copying keeps it valid.
    source_dir = source.resolve().parent
    for stage in document.stages:
        if stage.kind != "pcr" or not stage.primers:
            continue
        primer = Path(stage.primers)
        if primer.is_absolute():
            continue
        origin = source_dir / primer
        if not origin.is_file():
            raise InputError(
                f"composite {composite_name!r} references primer file {primer} that does not"
                f" exist next to {source}",
                code="INPUT_NOT_FOUND",
                context={"file": str(origin), "composite": composite_name},
            )
        shutil.copy2(origin, destination / primer.name)

    _write_manifest(destination, composite_name, document)
    return installed_document


def composite_summary(document: CompositeDocument) -> dict[str, Any]:
    """The list view of one composite: stages, panels and tiers."""
    return {
        "name": document.name,
        "kind": "composite",
        "fallback": document.fallback,
        "databases": sorted(_panel_names(document)),
        "stages": [
            {
                "name": stage.name,
                "kind": stage.kind,
                "db": stage.db,
                "scheme": stage.scheme,
                "tier": stage.tier,
            }
            for stage in document.stages
        ],
        "conflicts": len(document.conflicts),
    }


__all__ = [
    "BUNDLE_DIRNAME",
    "build_composite",
    "composite_document_path",
    "composite_names",
    "composite_summary",
]
