"""``best_hit`` — the gene-path evidence-ranking rule model.

``exact_set``/``weighted_genes`` are boolean in practice: a satisfied rule
scores 1.0, so homologous references tie and declaration order decides.
``best_hit`` scores the gene's own coverage/identity instead, which is what a
full-locus (or allele) reference set needs to rank its candidates.
"""

from gspec.cluster.core import GeneCall
from gspec.typing.gene import evaluate_gene_calls, score_gene_scheme
from gspec.typing.models import TypingDocument


def _call(gene: str, coverage: float, identity: float) -> GeneCall:
    return GeneCall(
        gene_id=gene,
        start=1,
        end=100,
        strand="+",
        coverage_pct=coverage,
        identity_pct=identity,
        verdict="present",
    )


def _document(*rules: dict[str, object]) -> TypingDocument:
    return TypingDocument.model_validate(
        {
            "schema": "gspec.typing/2",
            "schemes": [
                {
                    "name": "tier",
                    "cutoff": 0.5,
                    "ambiguity_margin": 0.0,
                    "fallback": "untypeable",
                    "rules": list(rules),
                }
            ],
        }
    )


def _rule(phenotype: str, gene: str, *, coverage_floor: float = 40.0) -> dict[str, object]:
    return {
        "model": "best_hit",
        "phenotype": phenotype,
        "gene": gene,
        "identity_floor": 80.0,
        "coverage_floor": coverage_floor,
    }


def test_ranking_prefers_the_better_evidence() -> None:
    """The higher identity*coverage wins, not the declaration order."""
    scheme = _document(_rule("HS2", "HS2"), _rule("HS4c", "HS4c")).schemes[0]
    calls = {"HS2": _call("HS2", 92.9, 99.7), "HS4c": _call("HS4c", 97.9, 99.6)}
    scored = {rule.phenotype: rule.score for rule in score_gene_scheme(scheme, calls)}
    assert scored["HS4c"] > scored["HS2"]
    phenotype = evaluate_gene_calls(calls, _document(_rule("HS2", "HS2"), _rule("HS4c", "HS4c")))[
        "tier"
    ].phenotype
    assert phenotype == "HS4c"


def _call_with_len(gene: str, coverage: float, identity: float, length: int) -> GeneCall:
    return GeneCall(
        gene_id=gene,
        start=1,
        end=100,
        strand="+",
        coverage_pct=coverage,
        identity_pct=identity,
        verdict="present",
        subject_len=length,
    )


def test_rank_prefers_the_longer_reference() -> None:
    """``rank: identity_coverage_length`` orders on alignment strength, so a
    41 kb locus at 95% beats a 16 kb one at 99% even though its score is lower."""
    long_rule = {
        **_rule("HS33_35", "HS33_35"),
        "identity_floor": 90.0,
        "coverage_floor": 40.0,
        "rank": "identity_coverage_length",
    }
    short_rule = {
        **_rule("HS19", "HS19"),
        "identity_floor": 90.0,
        "coverage_floor": 40.0,
        "rank": "identity_coverage_length",
    }
    document = TypingDocument.model_validate(
        {
            "schema": "gspec.typing/2",
            "schemes": [
                {
                    "name": "tier",
                    "cutoff": 0.3,
                    "ambiguity_margin": 0.0,
                    "fallback": "untypeable",
                    "rules": [long_rule, short_rule],
                }
            ],
        }
    )
    calls = {
        "HS33_35": _call_with_len("HS33_35", 60.0, 95.0, 28_524),
        "HS19": _call_with_len("HS19", 95.0, 99.0, 16_751),
    }
    scored = {rule.phenotype: rule for rule in score_gene_scheme(document.schemes[0], calls)}
    assert scored["HS19"].score > scored["HS33_35"].score
    assert scored["HS33_35"].order_key > scored["HS19"].order_key
    assert evaluate_gene_calls(calls, document)["tier"].phenotype == "HS33_35"


def test_floors_gate_the_score_to_zero() -> None:
    scheme = _document(_rule("HS2", "HS2", coverage_floor=60.0)).schemes[0]
    below = {"HS2": _call("HS2", 45.0, 99.9)}
    above = {"HS2": _call("HS2", 60.0, 99.9)}
    assert score_gene_scheme(scheme, below)[0].score == 0.0
    assert score_gene_scheme(scheme, above)[0].score > 0.0


def test_score_metric_is_configurable() -> None:
    rule = {**_rule("HS2", "HS2"), "score": "coverage"}
    scheme = _document(rule).schemes[0]
    scored = score_gene_scheme(scheme, {"HS2": _call("HS2", 50.0, 100.0)})
    assert abs(scored[0].score - 0.5) < 1e-9


def test_absent_gene_scores_zero_and_falls_back() -> None:
    document = _document(_rule("HS2", "HS2"))
    result = evaluate_gene_calls({}, document)["tier"]
    assert result.phenotype == "untypeable"
    assert result.score == 0.0


def test_unknown_gene_is_rejected_at_build_validation() -> None:
    """``best_hit`` references go through the same UNKNOWN_GENE guard."""
    import pytest

    from gspec.typing.gene import validate_gene_typing

    document = _document(_rule("HS2", "nope"))
    with pytest.raises(Exception) as exc_info:
        validate_gene_typing(document, frozenset({"HS2"}))
    assert getattr(exc_info.value, "code", "") == "TYPING_UNKNOWN_GENE"
    assert "nope" in str(exc_info.value)
