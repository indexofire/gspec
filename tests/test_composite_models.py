"""``gspec.composite/1`` — the multi-database fusion document."""

import json
from pathlib import Path

import pytest

from gspec.composite.models import read_composite
from gspec.errors import DatabaseError, InputError


def _doc(**overrides: object) -> str:
    payload: dict[str, object] = {
        "schema": "gspec.composite/1",
        "name": "cje_penner",
        "stages": [
            {"name": "marker", "db": "cje_hs_marker", "tier": "B", "min_score": 0.3},
            {"name": "locus_a", "db": "cje_cps_locus", "scheme": "tier_a", "tier": "A"},
        ],
    }
    payload.update(overrides)
    return json.dumps(payload)


def test_valid_document_parses(tmp_path: Path) -> None:
    path = tmp_path / "p.json"
    path.write_text(_doc(), encoding="utf-8")
    document = read_composite(path)
    assert [stage.name for stage in document.stages] == ["marker", "locus_a"]
    assert document.stages[0].min_score == 0.3
    assert document.stages[1].scheme == "tier_a"
    # Acceptance defaults: any confidence, no score floor.
    assert document.stages[1].require_confidence == "any"
    assert document.stages[1].min_score == 0.0


def test_missing_file_is_input_not_found(tmp_path: Path) -> None:
    with pytest.raises(InputError) as exc_info:
        read_composite(tmp_path / "nope.json")
    assert exc_info.value.code == "INPUT_NOT_FOUND"


@pytest.mark.parametrize(
    "text",
    [
        "{nope",
        '{"schema": "gspec.pipeline/2", "name": "x", "stages": []}',
        '{"schema": "gspec.composite/1", "name": "x", "stages": []}',
        '{"schema": "gspec.composite/1", "name": "x", "stages": [{"name": "a", "db": "d", '
        '"min_score": 1.5}]}',
    ],
)
def test_malformed_documents_are_typed(tmp_path: Path, text: str) -> None:
    path = tmp_path / "p.json"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(DatabaseError) as exc_info:
        read_composite(path)
    assert exc_info.value.code == "COMPOSITE_MALFORMED"


def test_duplicate_stage_names_are_rejected(tmp_path: Path) -> None:
    """Stage names label evidence columns, so they must be unique."""
    path = tmp_path / "p.json"
    path.write_text(
        _doc(
            stages=[
                {"name": "locus", "db": "a"},
                {"name": "locus", "db": "b"},
            ]
        ),
        encoding="utf-8",
    )
    with pytest.raises(DatabaseError) as exc_info:
        read_composite(path)
    assert exc_info.value.code == "COMPOSITE_MALFORMED"
    assert "duplicate stage names" in str(exc_info.value)


def test_fallback_defaults_to_untypeable(tmp_path: Path) -> None:
    path = tmp_path / "p.json"
    path.write_text(_doc(), encoding="utf-8")
    assert read_composite(path).fallback == "untypeable"
