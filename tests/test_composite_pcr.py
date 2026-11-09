"""Pipeline PCR stages: refine a cascade call with in-silico PCR.

A PCR stage is the *refinement* half of a pipeline: the database stages
designate, then a primer set either overrides the call (the assay is
authoritative — HS29/HS55/HS58 style) or brackets a complex call into a
subtype (HS4c(HS16/64) style). Everything is declared in the pipeline document:
primer set, mismatch budget, and the ordered rules that map "this product
amplified inside this size window, given this call" to a new call.

Runs through ``gseek.ispcr`` (pure stdlib), so the same pipeline works on a
machine without exonerate/ipcress.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

gseek = pytest.importorskip("gseek")

from gspec.composite.models import read_composite  # noqa: E402
from gspec.composite.pcr import needed_products  # noqa: E402
from gspec.composite.run import apply_conflicts  # noqa: E402
from gspec.composite.stages import inputs_digest  # noqa: E402
from gspec.errors import DatabaseError  # noqa: E402

FWD = "ACGTACGTACGTACGTACGT"
REV = "TTTTGGGGCCCCAAAATTTT"


def _rc(sequence: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return sequence.translate(table)[::-1]


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _pipeline_doc(tmp_path: Path, *, product_size: tuple[int, int] = (100, 400)) -> Path:
    primers = _write(
        tmp_path / "primers.tsv",
        "name\tforward\treverse\tmin_size\tmax_size\n"
        f"HS29_target\t{FWD}\t{REV}\t{product_size[0]}\t{product_size[1]}\n",
    )
    document = {
        "schema": "gspec.composite/1",
        "name": "pcr_demo",
        "fallback": "untypeable",
        "stages": [
            {
                "name": "pcr",
                "kind": "pcr",
                "db": "primers",
                "primers": str(primers),
                "rules": [
                    {
                        "products": [
                            {
                                "name": "HS29_target",
                                "min_size": product_size[0],
                                "max_size": product_size[1],
                            }
                        ],
                        "set": "HS29",
                    }
                ],
            }
        ],
    }
    return _write(tmp_path / "composite.json", json.dumps(document))


def test_screen_cache_does_not_leak_between_input_sets(tmp_path: Path) -> None:
    """The cache key must cover the input list, not just its length."""

    first = _write(tmp_path / "a.fna", ">a\nACGT\n")
    second = _write(tmp_path / "b.fna", ">b\nACGTACGT\n")
    assert inputs_digest([first]) != inputs_digest([second])
    assert inputs_digest([first, second]) == inputs_digest([second, first])


def _conflict_document(**conflict: object) -> object:
    from gspec.composite.models import CompositeDocument

    payload = {
        "schema": "gspec.composite/1",
        "name": "conflict_demo",
        "fallback": "untypeable",
        "stages": [
            {"name": "marker", "db": "db_marker", "tier": "B"},
            {"name": "locus", "db": "db_locus", "tier": "A"},
        ],
    }
    if conflict:
        payload["conflicts"] = [conflict]
    return CompositeDocument.model_validate(payload)


def _call(phenotype: str, score: float = 0.95, confidence: str = "high") -> object:
    from gspec.composite.calls import Call

    return Call(phenotype=phenotype, score=score, confidence=confidence, scheme="only")


def test_conflict_rule_swaps_a_marker_call_for_the_locus_call() -> None:
    """CJT's documented HS41/HS55 exception, expressed in configuration."""
    from gspec.composite.models import CompositeStage

    document = _conflict_document(
        when={"stage": "marker", "phenotype": "HS41"},
        override={"stage": "locus", "phenotype": "HS55"},
    )
    stages = {
        "marker": CompositeStage(name="marker", db="db_marker", tier="B"),
        "locus": CompositeStage(name="locus", db="db_locus", tier="A"),
    }
    winner = (stages["marker"], _call("HS41", 1.0))
    calls = {"marker": winner[1], "locus": _call("HS55", 0.93)}
    stage, call = apply_conflicts(document, winner, calls, stages)
    assert (stage.name, call.phenotype, stage.tier) == ("locus", "HS55", "A")


def test_conflict_override_may_name_any_stage() -> None:
    """``override.stage`` omitted = whichever database stage called it."""
    from gspec.composite.models import CompositeStage

    document = _conflict_document(
        when={"stage": "marker", "phenotype": "HS41"},
        override={"phenotype": "HS55"},
    )
    stages = {
        "marker": CompositeStage(name="marker", db="db_marker", tier="B"),
        "locus_a": CompositeStage(name="locus_a", db="db_locus", tier="A"),
        "locus_b": CompositeStage(name="locus_b", db="db_locus", tier="B"),
    }
    winner = (stages["marker"], _call("HS41", 1.0))
    calls = {
        "marker": winner[1],
        "locus_a": _call("HS2", 0.99),
        "locus_b": _call("HS55", 0.88),
    }
    stage, call = apply_conflicts(document, winner, calls, stages)
    assert (stage.name, call.phenotype) == ("locus_b", "HS55")


def test_conflict_does_not_fire_on_unrelated_or_gated_calls() -> None:
    from gspec.composite.models import CompositeStage

    stages = {
        "marker": CompositeStage(name="marker", db="db_marker", tier="B"),
        "locus": CompositeStage(name="locus", db="db_locus", tier="A"),
    }
    document = _conflict_document(
        when={"stage": "marker", "phenotype": "HS41"},
        override={"stage": "locus", "phenotype": "HS55"},
    )
    marker_call = _call("HS41", 1.0)

    # 1) the override stage called something else
    unmatching = {"marker": marker_call, "locus": _call("HS53")}
    assert (
        apply_conflicts(document, (stages["marker"], marker_call), unmatching, stages)[1].phenotype
        == "HS41"
    )

    # 2) the override call is the stage's own fallback (nothing designated)
    gated = {"marker": marker_call, "locus": _call("HS55", confidence="low")}
    assert (
        apply_conflicts(document, (stages["marker"], marker_call), gated, stages)[1].phenotype
        == "HS41"
    )

    # 3) no declaration at all
    assert (
        apply_conflicts(
            _conflict_document(), (stages["marker"], marker_call), {"marker": marker_call}, stages
        )[1].phenotype
        == "HS41"
    )


def test_parallel_pcr_matches_sequential(tmp_path: Path) -> None:
    """``jobs>1`` runs the PCR scan in worker processes; results must not move."""
    from gspec.composite.run import run_composite

    genomes = [
        _write(tmp_path / f"asm{i}.fna", f">ctg1\n{FWD}{'A' * (140 + i * 10)}{_rc(REV)}\n")
        for i in range(4)
    ]
    document = read_composite(_pipeline_doc(tmp_path))
    sequential = run_composite(genomes, document, datadir=tmp_path, jobs=1)
    parallel = run_composite(genomes, document, datadir=tmp_path, jobs=2)
    assert [result.phenotype for result in parallel] == [result.phenotype for result in sequential]
    assert [result.file for result in parallel] == [str(path) for path in genomes]
    assert all(result.phenotype == "HS29" for result in parallel)


def test_only_the_pairs_the_current_call_needs_are_scanned() -> None:
    """A ``when``-gated rule's products are skipped for unrelated calls."""
    from gspec.composite.models import CompositeStage

    stage = CompositeStage.model_validate(
        {
            "name": "pcr",
            "kind": "pcr",
            "db": "primers",
            "primers": "p.tsv",
            "rules": [
                {"products": [{"name": "ALWAYS"}], "set": "X"},
                {
                    "when": ["HS4c"],
                    "products": [{"name": "HS4A"}, {"name": "HS4B"}],
                    "set": "HS4c (HS4A)",
                },
            ],
        }
    )
    assert needed_products(stage, "HS2") == {"ALWAYS"}
    assert needed_products(stage, "HS4c (HS4/13/43/50/62/65 complex)") == {
        "ALWAYS",
        "HS4A",
        "HS4B",
    }
    # No call at all: only the unconditional override can still apply.
    assert needed_products(stage, None) == {"ALWAYS"}


def test_indel_budget_is_configurable_on_the_stage(tmp_path: Path) -> None:
    """``max_edits`` enables the indel pass for a stage's primers."""
    from gspec.composite.run import run_composite

    forward = "ACGTACGTACGTACGTACGTACGTACGTAC"
    reverse = "TTTTGGGGCCCCAAAATTTTGGGGCCCCAA"
    deleted = forward[:15] + forward[16:]
    primers = _write(
        tmp_path / "primers.tsv",
        f"name\tforward\treverse\tmin_size\tmax_size\n"
        f"HS29_target\t{forward}\t{reverse}\t100\t400\n",
    )
    genome = _write(tmp_path / "asm.fna", f">ctg1\n{deleted}{'A' * 150}{_rc(reverse)}\n")
    document = json.loads(_pipeline_doc(tmp_path).read_text(encoding="utf-8"))
    document["stages"][0]["primers"] = str(primers)
    document["stages"][0]["rules"][0]["products"][0]["name"] = "HS29_target"
    path = _write(tmp_path / "indel.json", json.dumps(document))
    strict = run_composite([genome], read_composite(path), datadir=tmp_path)
    assert strict[0].phenotype == "untypeable"

    document["stages"][0]["max_edits"] = 1
    lenient = read_composite(_write(tmp_path / "indel1.json", json.dumps(document)))
    assert run_composite([genome], lenient, datadir=tmp_path)[0].phenotype == "HS29"


def test_pcr_stage_parses_with_its_rules(tmp_path: Path) -> None:
    document = read_composite(_pipeline_doc(tmp_path))
    stage = document.stages[0]
    assert stage.kind == "pcr"
    assert stage.primers is not None
    assert stage.rules[0].set == "HS29"
    assert stage.rules[0].products[0].name == "HS29_target"


def test_pcr_stage_requires_a_primer_source(tmp_path: Path) -> None:
    document = json.loads(_pipeline_doc(tmp_path).read_text(encoding="utf-8"))
    document["stages"][0].pop("primers")
    path = _write(tmp_path / "bad.json", json.dumps(document))
    with pytest.raises(DatabaseError) as exc_info:
        read_composite(path)
    assert exc_info.value.code == "COMPOSITE_MALFORMED"


def test_pcr_stage_requires_rules(tmp_path: Path) -> None:
    document = json.loads(_pipeline_doc(tmp_path).read_text(encoding="utf-8"))
    document["stages"][0]["rules"] = []
    path = _write(tmp_path / "bad.json", json.dumps(document))
    with pytest.raises(DatabaseError) as exc_info:
        read_composite(path)
    assert exc_info.value.code == "COMPOSITE_MALFORMED"


def test_pcr_rule_needs_a_product(tmp_path: Path) -> None:
    document = json.loads(_pipeline_doc(tmp_path).read_text(encoding="utf-8"))
    document["stages"][0]["rules"][0]["products"] = []
    path = _write(tmp_path / "bad.json", json.dumps(document))
    with pytest.raises(DatabaseError) as exc_info:
        read_composite(path)
    assert exc_info.value.code == "COMPOSITE_MALFORMED"


def test_amplified_product_overrides_the_call(tmp_path: Path) -> None:
    """The full run: a primer pair that amplifies rewrites the database call."""
    from gspec.composite.run import run_composite

    genome = _write(tmp_path / "asm.fna", f">ctg1\n{FWD}{'A' * 150}{_rc(REV)}\n")
    document = read_composite(_pipeline_doc(tmp_path))
    results = run_composite([genome], document, datadir=tmp_path, cache_dir=tmp_path / "cache")
    assert len(results) == 1
    assert results[0].phenotype == "HS29"
    assert results[0].stage == "pcr"


def test_absent_product_leaves_the_call_alone(tmp_path: Path) -> None:
    from gspec.composite.run import run_composite

    genome = _write(tmp_path / "asm.fna", ">ctg1\n" + "A" * 400 + "\n")
    result = run_composite([genome], read_composite(_pipeline_doc(tmp_path)), datadir=tmp_path)[0]
    assert result.phenotype == "untypeable"
    assert result.stage == "-"


def test_when_condition_gates_the_refinement(tmp_path: Path) -> None:
    """``when`` restricts a rule to the calls it may bracket."""
    from gspec.composite.run import run_composite

    document = json.loads(_pipeline_doc(tmp_path).read_text(encoding="utf-8"))
    rule = document["stages"][0]["rules"][0]
    rule["when"] = ["HS4c"]
    rule["set"] = "HS4c (HS16/64)"
    path = _write(tmp_path / "gated.json", json.dumps(document))
    genome = _write(tmp_path / "asm.fna", f">ctg1\n{FWD}{'A' * 150}{_rc(REV)}\n")
    result = run_composite([genome], read_composite(path), datadir=tmp_path)[0]
    # No database stage produced HS4c, so the bracket cannot apply.
    assert result.phenotype == "untypeable"


def test_primer_path_resolves_next_to_the_pipeline(tmp_path: Path) -> None:
    """A relative primer path is read from the pipeline document's directory."""
    from gspec.composite.run import run_composite

    genome = _write(tmp_path / "asm.fna", f">ctg1\n{FWD}{'A' * 150}{_rc(REV)}\n")
    document_json = json.loads(_pipeline_doc(tmp_path).read_text(encoding="utf-8"))
    document_json["stages"][0]["primers"] = "primers.tsv"
    path = _write(tmp_path / "relative.json", json.dumps(document_json))
    result = run_composite([genome], read_composite(path), datadir=tmp_path)[0]
    assert result.phenotype == "HS29"


def test_unknown_primer_file_is_a_typed_error(tmp_path: Path) -> None:
    from gspec.composite.run import run_composite

    genome = _write(tmp_path / "asm.fna", ">ctg1\nACGT\n")
    document = json.loads(_pipeline_doc(tmp_path).read_text(encoding="utf-8"))
    document["stages"][0]["primers"] = str(tmp_path / "nope.tsv")
    path = _write(tmp_path / "bad.json", json.dumps(document))
    with pytest.raises(Exception) as exc_info:
        run_composite([genome], read_composite(path), datadir=tmp_path)
    assert getattr(exc_info.value, "code", "") == "INPUT_NOT_FOUND"


def test_missing_rule_product_is_a_typed_error(tmp_path: Path) -> None:
    """A rule naming a pair the primer set lacks fails loudly at load time."""
    from gspec.composite.run import run_composite

    genome = _write(tmp_path / "asm.fna", ">ctg1\nACGT\n")
    document = json.loads(_pipeline_doc(tmp_path).read_text(encoding="utf-8"))
    document["stages"][0]["rules"][0]["products"][0]["name"] = "nope_target"
    path = _write(tmp_path / "bad.json", json.dumps(document))
    with pytest.raises(DatabaseError) as exc_info:
        run_composite([genome], read_composite(path), datadir=tmp_path)
    assert exc_info.value.code == "COMPOSITE_MALFORMED"


def test_only_the_pairs_the_rules_name_are_scanned(tmp_path: Path) -> None:
    """Evidence by default lists the rule-winning product, not every pair."""
    from gspec.composite.run import run_composite

    genome = _write(tmp_path / "asm.fna", f">ctg1\n{FWD}{'A' * 150}{_rc(REV)}\n")
    result = run_composite([genome], read_composite(_pipeline_doc(tmp_path)), datadir=tmp_path)[0]
    primer_evidence = dict(result.evidence)["pcr"]
    assert isinstance(primer_evidence, str)
    assert "HS29_target" in primer_evidence


def test_primer_mismatch_budget_is_configurable(tmp_path: Path) -> None:
    from gspec.composite.run import run_composite

    target = "T" + FWD[1:]
    genome = _write(tmp_path / "asm.fna", f">ctg1\n{target}{'A' * 150}{_rc(REV)}\n")
    document = json.loads(
        _pipeline_doc(tmp_path, product_size=(100, 400)).read_text(encoding="utf-8")
    )
    document["stages"][0]["max_mismatches"] = 0
    strict = read_composite(_write(tmp_path / "strict.json", json.dumps(document)))
    assert run_composite([genome], strict, datadir=tmp_path)[0].phenotype == "untypeable"

    document["stages"][0]["max_mismatches"] = 1
    lenient = read_composite(_write(tmp_path / "lenient.json", json.dumps(document)))
    assert run_composite([genome], lenient, datadir=tmp_path)[0].phenotype == "HS29"
