"""``gspec screen -d <composite>``: a composite database through the screen surface.

A composite carries ``composite.json`` instead of a panel: `screen` runs its
stages (panels, k-mer sets, in-silico PCR), refuses single-panel flags (a
composite configures its own stages) and emits the fused typing result.
"""

from __future__ import annotations

import json
import random
from pathlib import Path

import pytest
from click.testing import CliRunner

from gspec.cli import app
from gspec.composite_build import build_composite
from gspec.dbbuild import build_database
from gspec.records import Record, write_records

DB = "scheme_gene"


def _synthetic_dna(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


SEQ_DEMO = _synthetic_dna(240, seed=11)
SEQ_OTHER = _synthetic_dna(240, seed=12)

RECORDS = (
    Record(db=DB, gene="demo", sequence=SEQ_DEMO, accession="SYN-1", product="demo"),
    Record(db=DB, gene="other", sequence=SEQ_OTHER, accession="SYN-2", product="other"),
)


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    db_dir = tmp_path / "datadir" / DB
    db_dir.mkdir(parents=True)
    write_records(RECORDS, db_dir / "records.jsonl")
    build_database(
        db_dir, name=DB, dbtype="nucl", source_urls=(), fetched_at="2026-09-17T00:00:00Z"
    )
    typing = tmp_path / "typing.json"
    typing.write_text(
        json.dumps(
            {
                "schema": "gspec.typing/2",
                "schemes": [
                    {
                        "name": "only",
                        "cutoff": 0.5,
                        "ambiguity_margin": 0.0,
                        "fallback": "untypeable",
                        "rules": [
                            {
                                "model": "best_hit",
                                "phenotype": "TypeA",
                                "gene": "demo",
                                "identity_floor": 90.0,
                                "coverage_floor": 90.0,
                            }
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    built = CliRunner().invoke(
        app,
        [
            "db",
            "build",
            DB,
            str(db_dir / "sequences"),
            "-T",
            str(typing),
            "-D",
            str(tmp_path / "datadir"),
            "-f",
        ],
    )
    assert built.exit_code == 0, built.output
    return tmp_path / "datadir"


@pytest.fixture()
def scheme(tmp_path: Path, datadir: Path) -> str:
    document = tmp_path / "scheme.json"
    document.write_text(
        json.dumps(
            {
                "schema": "gspec.composite/1",
                "name": "demo_scheme",
                "fallback": "untypeable",
                "stages": [{"name": "panel", "db": DB, "tier": "A"}],
            }
        ),
        encoding="utf-8",
    )
    build_composite(document, datadir, bundle=True)
    return "demo_scheme"


def _sample(tmp_path: Path) -> Path:
    path = tmp_path / "sample.fna"
    path.write_text(f">ctg1\n{SEQ_DEMO}\n", encoding="utf-8")
    return path


def test_screen_runs_an_installed_scheme(tmp_path: Path, datadir: Path, scheme: str) -> None:
    result = CliRunner().invoke(
        app, ["screen", "-d", scheme, str(_sample(tmp_path)), "-D", str(datadir)]
    )
    assert result.exit_code == 0, result.output
    header = result.stdout.splitlines()[0]
    assert header.split("\t")[:5] == ["FILE", "PIPELINE", "PHENOTYPE", "TIER", "STAGE"]
    row = result.stdout.splitlines()[1].split("\t")
    assert row[2] == "TypeA"
    assert row[4] == "panel"


def test_screen_scheme_json_is_the_pipeline_result(
    tmp_path: Path, datadir: Path, scheme: str
) -> None:
    result = CliRunner().invoke(
        app, ["screen", "-d", scheme, str(_sample(tmp_path)), "-D", str(datadir), "-f", "json"]
    )
    assert result.exit_code == 0, result.output
    document = json.loads(result.stdout)
    assert document["schema"] == "gspec.composite_result/1"
    assert document["files"][0]["phenotype"] == "TypeA"


def test_screen_rejects_database_flags_for_a_scheme(
    tmp_path: Path, datadir: Path, scheme: str
) -> None:
    result = CliRunner().invoke(
        app,
        ["screen", "-d", scheme, str(_sample(tmp_path)), "-D", str(datadir), "-i", "95"],
    )
    assert result.exit_code == 2
    assert "--minid" in result.stderr
    assert "composite database" in result.stderr


def test_unknown_name_is_a_database_not_found(tmp_path: Path, datadir: Path) -> None:
    result = CliRunner().invoke(
        app, ["screen", "-d", "nope", str(_sample(tmp_path)), "-D", str(datadir)]
    )
    assert result.exit_code != 0
    assert "DATABASE_NOT_FOUND" in result.stderr


def test_a_bundled_composite_runs_from_a_bare_datadir(
    tmp_path: Path, datadir: Path, scheme: str
) -> None:
    """The composite directory goes anywhere; only it needs to be in the datadir."""
    import shutil

    bare = tmp_path / "bare"
    bare.mkdir()
    shutil.copytree(datadir / scheme, bare / scheme)
    result = CliRunner().invoke(
        app, ["screen", "-D", str(bare), "-d", scheme, str(_sample(tmp_path))]
    )
    assert result.exit_code == 0, result.output
    assert result.stdout.splitlines()[1].split("\t")[2] == "TypeA"


def test_composite_build_pins_db_root_and_copies_panels(
    tmp_path: Path, datadir: Path, scheme: str
) -> None:
    """`db build -C --bundle` is how a distributable composite is made."""
    document = json.loads((datadir / scheme / "composite.json").read_text(encoding="utf-8"))
    assert document["db_root"] == "dbs"
    assert (datadir / scheme / "dbs" / DB / "sequences").is_file()
