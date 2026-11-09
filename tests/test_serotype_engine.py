"""gspec.serotype/1 rule evaluation tests."""

from __future__ import annotations

import pytest

from gspec.cluster.core import GeneCall
from gspec.typing.serotype_engine import evaluate_serotype
from gspec.typing.serotype_models import SerotypeDocument, SerotypeRule


def _call(
    gene: str,
    *,
    coverage: float = 100.0,
    identity: float = 100.0,
    alignment_start: float | None = None,
) -> GeneCall:
    return GeneCall(
        gene_id=gene,
        start=1,
        end=100,
        strand="+",
        coverage_pct=coverage,
        identity_pct=identity,
        verdict="present",
        alignment_start_pct=alignment_start,
    )


def _document() -> SerotypeDocument:
    return SerotypeDocument(
        name="test",
        control_gene="control",
        rules=(
            SerotypeRule(phenotype="type A", requires=("a",)),
            SerotypeRule(phenotype="type B", requires=("b",)),
            SerotypeRule(
                phenotype="exempt type",
                requires=("x",),
                control_exempt=True,
            ),
        ),
    )


def test_serotype_rule_call_is_high_confidence() -> None:
    call = evaluate_serotype({"control": _call("control"), "a": _call("a")}, _document())
    assert (call.phenotype, call.confidence) == ("type A", "high")


def test_multiple_matching_rules_are_ambiguous() -> None:
    call = evaluate_serotype(
        {"control": _call("control"), "a": _call("a"), "b": _call("b")},
        _document(),
    )
    assert call.phenotype is None
    assert call.confidence == "ambiguous"
    assert {entry.phenotype for entry in call.ambiguous} == {"type A", "type B"}


def _nested_document(decision: str) -> SerotypeDocument:
    """A full-formula rule and the generic group rule it implies."""
    return SerotypeDocument(
        name="nested",
        control_gene="control",
        decision=decision,  # type: ignore[arg-type]
        rules=(
            SerotypeRule(phenotype="serovar X (4:i:1,2)", requires=("o4", "h_i", "h_12")),
            SerotypeRule(phenotype="group O:4, H not determined", requires=("o4",)),
        ),
    )


def test_default_decision_reports_nested_rules_as_ambiguous() -> None:
    """Given the default policy, When a specific rule and the generic rule it
    implies both match, Then the call stays ambiguous — existing profiles
    (shi_sero) keep their behaviour."""
    calls = {
        "control": _call("control"),
        "o4": _call("o4"),
        "h_i": _call("h_i"),
        "h_12": _call("h_12"),
    }
    call = evaluate_serotype(calls, _nested_document("ambiguous"))
    assert call.confidence == "ambiguous"
    assert {entry.phenotype for entry in call.ambiguous} == {
        "serovar X (4:i:1,2)",
        "group O:4, H not determined",
    }


def test_most_specific_decision_prefers_the_full_formula() -> None:
    """Under most_specific, the 3-marker rule outranks the 1-marker rule it
    implies, so the complete antigenic formula wins."""
    calls = {
        "control": _call("control"),
        "o4": _call("o4"),
        "h_i": _call("h_i"),
        "h_12": _call("h_12"),
    }
    call = evaluate_serotype(calls, _nested_document("most_specific"))
    assert (call.phenotype, call.confidence) == ("serovar X (4:i:1,2)", "high")


def test_most_specific_keeps_the_generic_rule_for_partial_formulas() -> None:
    """Given only the O antigen, Then the generic rule answers instead of
    falling back: "O group known, H undetermined" is a real result."""
    calls = {"control": _call("control"), "o4": _call("o4")}
    call = evaluate_serotype(calls, _nested_document("most_specific"))
    assert call.phenotype == "group O:4, H not determined"
    assert call.confidence == "high"


def test_most_specific_keeps_incomparable_matches_ambiguous() -> None:
    """Two rules that are not in a subset relation stay ambiguous: the policy
    resolves specificity, never genuine conflict."""
    document = SerotypeDocument(
        name="conflict",
        control_gene="control",
        decision="most_specific",
        rules=(
            SerotypeRule(phenotype="type A", requires=("a",)),
            SerotypeRule(phenotype="type B", requires=("b",)),
        ),
    )
    call = evaluate_serotype(
        {"control": _call("control"), "a": _call("a"), "b": _call("b")}, document
    )
    assert call.confidence == "ambiguous"


def _ranked_document() -> SerotypeDocument:
    """O:2 vs O:9: two near-identical alleles, decided by ranking.

    Both directions declare the gate — the DB pattern — so whichever allele
    scores higher wins and the other is blocked by its own gate.
    """
    return SerotypeDocument(
        name="ranked",
        control_gene="control",
        unresolved="no allele ranked",
        rules=(
            SerotypeRule(phenotype="O2", requires=("tyr_A",), outranks={"tyr_A": ("rfbE_D1",)}),
            SerotypeRule(phenotype="O9", requires=("rfbE_D1",), outranks={"rfbE_D1": ("tyr_A",)}),
        ),
    )


def _ranked_call(gene: str, identity: float, coverage: float = 100.0) -> dict[str, GeneCall]:
    return {gene: _call(gene, identity=identity, coverage=coverage)}


def test_outranks_prefers_the_higher_scoring_allele() -> None:
    """Given a group-A copy (100%) against a group-D1 copy (98.9%), Then the A
    rule fires and the D1 rule is blocked by its own gate."""
    calls = {
        "control": _call("control"),
        **_ranked_call("tyr_A", 100.0),
        **_ranked_call("rfbE_D1", 98.9),
    }
    call = evaluate_serotype(calls, _ranked_document())
    assert (call.phenotype, call.confidence) == ("O2", "high")


def test_outranks_blocks_the_lower_scoring_allele() -> None:
    """Given the D1 copy on top, Then the A rule's gate fails and D1 wins."""
    calls = {
        "control": _call("control"),
        **_ranked_call("tyr_A", 98.9),
        **_ranked_call("rfbE_D1", 100.0),
    }
    call = evaluate_serotype(calls, _ranked_document())
    assert (call.phenotype, call.confidence) == ("O9", "high")


def test_outranks_tie_reports_unresolved_not_a_guess() -> None:
    """A tie fails "strictly greater" on both sides, so neither allele is
    called. SeqSero2 defaults to O:9 here; gspec reports the missing decision
    instead of inheriting a preference from the tool it replaced."""
    calls = {
        "control": _call("control"),
        **_ranked_call("tyr_A", 100.0),
        **_ranked_call("rfbE_D1", 100.0),
    }
    call = evaluate_serotype(calls, _ranked_document())
    assert (call.phenotype, call.confidence) == ("no allele ranked", "low")


def test_outranks_ignores_an_absent_competitor() -> None:
    calls = {"control": _call("control"), **_ranked_call("tyr_A", 100.0)}
    call = evaluate_serotype(calls, _ranked_document())
    assert call.phenotype == "O2"


def test_outranks_competitor_below_the_floor_does_not_block() -> None:
    """A competitor that cannot pass its own floor is not evidence to beat."""
    calls = {
        "control": _call("control"),
        **_ranked_call("tyr_A", 100.0),
        **_ranked_call("rfbE_D1", 60.0),  # below the 80% rule floor
    }
    call = evaluate_serotype(calls, _ranked_document())
    assert call.phenotype == "O2"


def test_outranks_requires_a_positive_marker() -> None:
    """The gate must name a marker the rule itself requires."""
    with pytest.raises(ValueError, match="not a positive marker"):
        SerotypeRule(phenotype="bad", requires=("a",), outranks={"b": ("c",)})


def test_outranks_rejects_a_self_reference() -> None:
    with pytest.raises(ValueError, match="names itself"):
        SerotypeRule(phenotype="bad", requires=("a",), outranks={"a": ("a",)})


def test_missing_control_falls_back_without_control_exempt() -> None:
    call = evaluate_serotype({"a": _call("a")}, _document())
    assert (call.phenotype, call.confidence) == ("Not Shigella/EIEC", "low")


def test_control_exempt_rule_fires_without_control() -> None:
    call = evaluate_serotype({"x": _call("x")}, _document())
    assert (call.phenotype, call.confidence) == ("exempt type", "high")


def test_rule_floors_filter_weak_markers() -> None:
    document = SerotypeDocument(
        name="test",
        rules=(
            SerotypeRule(
                phenotype="strict",
                requires=("a",),
                coverage_floor=90.0,
                identity_floor=95.0,
            ),
        ),
    )
    call = evaluate_serotype({"a": _call("a", coverage=80.0)}, document)
    assert call.phenotype == "Shigella spp. (unresolved)"


def _sb6_sb10_document() -> SerotypeDocument:
    """The shipped SB6/SB10 rule shape: a 5'-intact wbaM is type 10, a
    5'-truncated (IS629) wbaM is type 6, and missing wbaM stays ambiguous."""
    return SerotypeDocument(
        name="test",
        rules=(
            SerotypeRule(
                phenotype="S. boydii 10",
                requires=("Sb6_wzx", "wbaM"),
                max_alignment_start_pct={"wbaM": 20.0},
            ),
            SerotypeRule(
                phenotype="S. boydii 6",
                requires=("Sb6_wzx", "wbaM"),
                max_coverage={"wbaM": 95.0},
                min_alignment_start_pct={"wbaM": 20.0},
            ),
            SerotypeRule(
                phenotype="S. boydii 6/10",
                requires=("Sb6_wzx",),
                excludes=("wbaM",),
            ),
        ),
    )


def test_wbam_alignment_start_separates_partial_assemblies() -> None:
    """Given wbaM alignments, When the 5' end is covered, Then the call is
    SB10 even with a 3'-truncated assembly; a 5'-truncated alignment (the
    IS629 insertion) is SB6."""
    document = _sb6_sb10_document()
    calls = {
        "Sb6_wzx": _call("Sb6_wzx"),
        "wbaM": _call("wbaM", alignment_start=0.0),
    }
    intact = evaluate_serotype(calls, document)
    three_prime_partial = evaluate_serotype(
        {"Sb6_wzx": _call("Sb6_wzx"), "wbaM": _call("wbaM", coverage=86.83, alignment_start=0.0)},
        document,
    )
    five_prime_partial = evaluate_serotype(
        {
            "Sb6_wzx": _call("Sb6_wzx"),
            "wbaM": _call("wbaM", coverage=73.67, alignment_start=26.33),
        },
        document,
    )
    assert intact.phenotype == "S. boydii 10"
    assert three_prime_partial.phenotype == "S. boydii 10"
    assert five_prime_partial.phenotype == "S. boydii 6"


def test_wbam_start_gate_requires_alignment_evidence() -> None:
    """Given reads tables (no alignment start), When wbaM is present, Then
    the SB6/SB10 rules do not match: breadth cannot locate the IS629
    junction, and the resolver stays unresolved instead of guessing."""
    document = _sb6_sb10_document()
    for coverage in (73.67, 100.0):
        call = evaluate_serotype(
            {"Sb6_wzx": _call("Sb6_wzx"), "wbaM": _call("wbaM", coverage=coverage)},
            document,
        )
        assert call.phenotype == "Shigella spp. (unresolved)"


def test_missing_wbam_keeps_sb6_sb10_ambiguous() -> None:
    """Given Sb6_wzx without wbaM, When the resolver runs, Then it keeps
    the explicit 6/10 ambiguity rather than guessing."""
    call = evaluate_serotype({"Sb6_wzx": _call("Sb6_wzx")}, _sb6_sb10_document())
    assert call.phenotype == "S. boydii 6/10"


def test_alignment_start_gates_validate_their_markers() -> None:
    """Given a gate on a non-positive marker or a non-percentage value,
    When the rule is constructed, Then validation rejects it."""
    with pytest.raises(ValueError, match="not a positive marker"):
        SerotypeRule(
            phenotype="bad",
            requires=("a",),
            max_alignment_start_pct={"other": 10.0},
        )
    with pytest.raises(ValueError, match="must be a percentage"):
        SerotypeRule(
            phenotype="bad",
            requires=("a",),
            min_alignment_start_pct={"a": 150.0},
        )
