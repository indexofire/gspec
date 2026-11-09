"""Evaluate gspec.serotype/1 rules over one file's gene calls."""

from __future__ import annotations

from collections.abc import Mapping

from gspec.cluster.core import GeneCall
from gspec.typing.results import PhenotypeScore, SchemeCall, ScoreComponent
from gspec.typing.serotype_models import (
    AntigenFamily,
    AnySerotypeDocument,
    SerotypeDocument,
    SerotypeDocumentV2,
    SerotypeRule,
)


def _present(call: GeneCall | None, coverage: float, identity: float) -> bool:
    return call is not None and call.coverage_pct >= coverage and call.identity_pct >= identity


def _evidence(call: GeneCall) -> float:
    """The ranking score: identity x coverage, in [0, 1].

    The same metric the `best_hit` typing rule ranks by, so a serotype profile
    and a typing scheme order competing alleles identically.
    """
    return (call.identity_pct / 100.0) * (call.coverage_pct / 100.0)


def _outranks(
    rule: SerotypeRule, calls: Mapping[str, GeneCall], *, coverage: float, identity: float
) -> bool:
    """Every ranked allele must beat each of its competitors outright.

    A tie fails the gate (strictly greater), so the tie-break belongs to rule
    order: the group the reference tools default to (SeqSero2 keeps O:9 on a
    tie) simply declares no competitor, and its rule wins the ordered pass.
    """
    for gene, competitors in rule.outranks.items():
        own = calls.get(gene)
        if not _present(own, coverage, identity):
            return False
        assert own is not None  # _present proved it
        own_score = _evidence(own)
        for competitor in competitors:
            other = calls.get(competitor)
            if _present(other, coverage, identity):
                assert other is not None
                if own_score <= _evidence(other):
                    return False
    return True


def _rule_matches(rule: SerotypeRule, calls: Mapping[str, GeneCall]) -> tuple[str, ...] | None:
    """The ordered positive markers when the rule matches, else None."""
    positives: list[str] = []
    for gene in rule.requires:
        if not _present(calls.get(gene), rule.coverage_floor, rule.identity_floor):
            return None
        positives.append(gene)
    if rule.requires_any:
        any_hit = next(
            (
                gene
                for gene in rule.requires_any
                if _present(calls.get(gene), rule.coverage_floor, rule.identity_floor)
            ),
            None,
        )
        if any_hit is None:
            return None
        positives.append(any_hit)
    for gene in rule.excludes:
        if _present(calls.get(gene), rule.coverage_floor, rule.identity_floor):
            return None
    if rule.outranks and not _outranks(
        rule, calls, coverage=rule.coverage_floor, identity=rule.identity_floor
    ):
        return None
    for gene, ceiling in rule.max_coverage.items():
        call = calls.get(gene)
        if call is not None and call.coverage_pct >= ceiling:
            return None
    for gene, minimum in rule.min_alignment_start_pct.items():
        call = calls.get(gene)
        if call is not None and (
            call.alignment_start_pct is None or call.alignment_start_pct < minimum
        ):
            return None
    for gene, maximum in rule.max_alignment_start_pct.items():
        call = calls.get(gene)
        if call is not None and (
            call.alignment_start_pct is None or call.alignment_start_pct > maximum
        ):
            return None
    return tuple(positives)


def _control_present(calls: Mapping[str, GeneCall], document: AnySerotypeDocument) -> bool:
    if document.control_gene is None:
        return True
    return _present(
        calls.get(document.control_gene),
        document.control_coverage_floor,
        document.control_identity_floor,
    )


def _most_specific(
    matched: list[tuple[SerotypeRule, tuple[str, ...]]],
) -> list[tuple[SerotypeRule, tuple[str, ...]]]:
    """Drop matches whose positive markers are a strict subset of another match.

    A serovar rule (O + H1 + H2) therefore outranks the generic group rule
    (O only) that must sit beneath it for partial formulas; matches that are
    not comparable — different markers, neither a subset — are kept, so real
    ambiguity is still reported.
    """
    positives = {id(rule): set(genes) for rule, genes in matched}
    kept: list[tuple[SerotypeRule, tuple[str, ...]]] = []
    for rule, genes in matched:
        own = positives[id(rule)]
        if any(own < other and other != own for key, other in positives.items() if key != id(rule)):
            continue
        kept.append((rule, genes))
    return kept


def _family_call(
    family: AntigenFamily, calls: Mapping[str, GeneCall], document: SerotypeDocumentV2
) -> tuple[str | None, float, tuple[str, ...]]:
    """The winning antigen of one family: (name|None, score, winning gene ids).

    Every reference allele is scored by identity x coverage (the `best_hit`
    metric); alleles within the document's ``margin`` of the winner make the
    family ambiguous (None), which is how near-identical alleles are pinned
    down — SeqSero2 compares two such scores by hand, gspec by the margin.
    """
    scored: list[tuple[float, str, tuple[str, ...]]] = []
    for entry in family.entries:
        best = 0.0
        winners: tuple[str, ...] = ()
        for gene in entry.genes:
            call = calls.get(gene)
            if not _present(call, document.coverage_floor, document.identity_floor):
                continue
            assert call is not None
            score = _evidence(call)
            if score > best:
                best, winners = score, (gene,)
            elif score == best:
                winners += (gene,)
        if best > 0.0:
            scored.append((best, entry.name, winners))
    if not scored:
        return None, 0.0, ()
    scored.sort(key=lambda item: item[0], reverse=True)
    best_score, name, winners = scored[0]
    if len(scored) > 1 and best_score - scored[1][0] <= document.margin:
        return None, best_score, winners
    return name, best_score, winners


def _lookup_formula(
    document: SerotypeDocumentV2, o: str | None, h1: str | None, h2: str | None
) -> str | None:
    """The serovar for a formula, trying the phase combinations the table holds.

    Monophasic serovars are tabulated with ``"-"``, so an unresolved phase is
    looked up as ``"-"`` first and only then reported as a partial formula.
    """
    if o is None:
        return None
    index = {(row.o, row.h1, row.h2): row.serovar for row in document.formulas}
    first, second = h1 or "-", h2 or "-"
    for key in ((o, first, second), (o, first, "-"), (o, "-", second), (o, "-", "-")):
        found = index.get(key)
        if found is not None:
            return found
    return None


def _evaluate_v2(calls: Mapping[str, GeneCall], document: SerotypeDocumentV2) -> SchemeCall:
    """Ranked antigen families, then the formula table."""
    control = _control_present(calls, document)
    resolved: dict[str, tuple[str | None, float, tuple[str, ...]]] = {
        family.key: _family_call(family, calls, document) for family in document.antigens
    }
    o, o_score, _o_genes = resolved.get("o", (None, 0.0, ()))
    h1 = resolved.get("h1", (None, 0.0, ()))
    h2 = resolved.get("h2", (None, 0.0, ()))
    components = tuple(
        ScoreComponent(name=f"antigen:{gene}", score=score)
        for _, score, genes in resolved.values()
        for gene in genes
    )
    notes = [
        f"formula {o or '-'}:{h1[0] or '-'}:{h2[0] or '-'}",
        *(
            f"{family_key.upper()}: {'/'.join(sorted(genes)) or '-'}"
            for family_key, (_, _, genes) in sorted(resolved.items())
        ),
    ]
    if not control:
        return SchemeCall(
            phenotype=document.fallback,
            score=0.0,
            confidence="low",
            notes=("control gene absent",),
        )
    serovar = _lookup_formula(document, o, h1[0], h2[0]) if control else None
    if serovar is not None:
        score = max(o_score, h1[1], h2[1])
        return SchemeCall(
            phenotype=serovar,
            score=round(score, 4),
            confidence="high",
            components=components,
            notes=tuple(notes),
        )
    if o is None:
        return SchemeCall(
            phenotype=document.unresolved,
            score=0.0,
            confidence="low",
            components=components,
            notes=tuple(notes),
        )
    formula = f"{o}:{h1[0] or '-'}:{h2[0] or '-'}"
    return SchemeCall(
        phenotype=document.partial.format(formula=formula),
        score=round(max(o_score, h1[1], h2[1]), 4),
        confidence="low",
        components=components,
        notes=tuple(notes),
    )


def evaluate_serotype(calls: Mapping[str, GeneCall], document: AnySerotypeDocument) -> SchemeCall:
    """One serotype profile decision, for either document shape."""
    if isinstance(document, SerotypeDocumentV2):
        return _evaluate_v2(calls, document)
    return _evaluate_v1(calls, document)


def _evaluate_v1(calls: Mapping[str, GeneCall], document: SerotypeDocument) -> SchemeCall:
    """One ordered marker-set decision, surfacing every matching rule.

    A rule whose ``control_exempt`` is set can still fire when the scheme
    control gene is absent (the ShigaTyper S. boydii 13 precedent). Multiple
    distinct phenotypes are an ambiguous call, never declaration-order luck.
    """
    control = _control_present(calls, document)
    matched: list[tuple[SerotypeRule, tuple[str, ...]]] = []
    for rule in document.rules:
        if not control and not rule.control_exempt:
            continue
        positive = _rule_matches(rule, calls)
        if positive is not None:
            matched.append((rule, positive))
    if matched and document.decision == "most_specific":
        matched = _most_specific(matched)
    if matched:
        phenotypes = {rule.phenotype for rule, _ in matched}
        positives = sorted({gene for _, genes in matched for gene in genes})
        components = tuple(ScoreComponent(name=f"marker:{gene}", score=1.0) for gene in positives)
        if len(phenotypes) == 1:
            rule, genes = matched[0]
            notes = list(rule.notes)
            notes.append("markers: " + ",".join(genes))
            if not control:
                notes.append("control gene absent; control_exempt rule used")
            return SchemeCall(
                phenotype=rule.phenotype,
                score=1.0,
                confidence="high",
                components=components,
                notes=tuple(notes),
            )
        candidates = tuple(
            PhenotypeScore(phenotype=phenotype, score=1.0) for phenotype in sorted(phenotypes)
        )
        return SchemeCall(
            phenotype=None,
            score=1.0,
            confidence="ambiguous",
            components=components,
            ambiguous=candidates,
            notes=(document.ambiguous, "candidates: " + ",".join(sorted(phenotypes))),
        )
    if not control:
        return SchemeCall(
            phenotype=document.fallback,
            score=0.0,
            confidence="low",
            notes=("control gene absent",),
        )
    return SchemeCall(
        phenotype=document.unresolved,
        score=0.0,
        confidence="low",
        notes=("no marker set matched",),
    )


__all__ = ["evaluate_serotype"]
