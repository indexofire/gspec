"""Evaluate gspec.serotype/1 rules over one file's gene calls."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Literal

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


@dataclass(frozen=True)
class _FamilyCall:
    """One family's decision — the mirror of SISTR's prediction record.

    ``perfect`` reproduces ``BlastReader.top_result()``: the call is perfect only
    when a 100 %-identity/100 %-coverage hit exists **and every hit for that gene
    names the same antigen**. ``ambiguous`` is the document's own ``margin`` rule —
    another antigen scored within it — which is what a cross-labelled reference
    looks like (`1,2` vs `1,2,7` are 99.54 % identical in the panel).
    """

    name: str | None
    score: float
    genes: tuple[str, ...]
    perfect: bool
    ambiguous: bool
    #: Every name within the document's ``margin`` of the winner, best first. The
    #: lookup widens all of them: a score gap inside the margin is the measured
    #: noise between sibling antigens (Newport's `e,h` scored 0.9103 against
    #: `e,n,x`'s 0.9152 on a real genome), so committing to one discards the name
    #: the upstream tool reports.
    ties: tuple[str, ...] = ()


def _family_call(
    family: AntigenFamily, calls: Mapping[str, GeneCall], document: SerotypeDocumentV2
) -> _FamilyCall:
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
        return _FamilyCall(name=None, score=0.0, genes=(), perfect=False, ambiguous=False)
    scored.sort(key=lambda item: item[0], reverse=True)
    best_score, name, winners = scored[0]
    names = {entry_name for _, entry_name, _ in scored}
    cutoff = best_score - document.margin
    ties: list[str] = []
    for score, entry_name, _ in scored:
        if score >= cutoff and entry_name not in ties:
            ties.append(entry_name)
    return _FamilyCall(
        name=name,
        score=best_score,
        genes=winners,
        perfect=best_score >= 1.0 and len(names) == 1,
        ambiguous=len(ties) > 1,
        ties=tuple(ties),
    )


def _formula_keys(o: str, h1: str | None, h2: str | None) -> list[tuple[str, str, str]]:
    """The phase combinations a table may hold for this call.

    Monophasic serovars are tabulated with ``"-"``, so an unresolved phase is
    looked up as ``"-"`` before the call is called partial.
    """
    first, second = h1 or "-", h2 or "-"
    return [(o, first, second), (o, first, "-"), (o, "-", second), (o, "-", "-")]


FALLBACK_PHASE = "-"


def _lookup_formula(
    document: SerotypeDocumentV2,
    o: tuple[str, ...],
    h1: tuple[str, ...],
    h2: tuple[str, ...],
    *,
    o_names: tuple[str, ...] = (),
) -> tuple[str, ...]:
    """The serovar(s) for a formula: exact first, then similarity widening.

    A single serovar is a call; several are genuine ambiguity (the upstream tools
    pipe-join them, gspec reports candidates); none is a partial formula.
    """
    index: dict[tuple[str, str, str], list[str]] = {}
    for row in document.formulas:
        index.setdefault((row.o, row.h1, row.h2), []).append(row.serovar)

    # An undecided O group widens to every serogroup the panel declares, which is
    # the upstream behaviour (SISTR replaces an unknown serogroup with the table's
    # whole column and reads the result as a pipe-joined candidate list). Without
    # it a genome whose O group is absent from the assembly cannot be named at all
    # — which is the common case on real data (`-:i:1,2` in SISTR's own output).
    # The H phases are NOT widened: their absence is a real monophasic signal and
    # the `"-"` alias rows already encode it.
    def tokens(
        family: Literal["o", "h1", "h2"], names: tuple[str, ...], *, missing_all: bool
    ) -> tuple[str, ...]:
        """The names a family's call stands for.

        A resolved name is widened to its similarity class *before* the lookup
        (SISTR replaces the predicted antigen with its `*_SIMILARITY_GROUPS` class
        and matches with `isin`), which is what reconciles cross-labelled
        references: `1,2` and `1,2,7` are 99.54 % identical. A missing O group
        widens to every serogroup instead; a missing H phase stays "-" because the
        `can_h2_be_missing` alias rows already encode that case.
        """
        if not names:
            if missing_all:
                return o_names or (FALLBACK_PHASE,)
            return (FALLBACK_PHASE,)
        members: list[str] = []
        for name in names:
            if name not in members:
                members.append(name)
            for group in document.similarity_groups.get(family, ()):
                if name in group:
                    members.extend(member for member in group if member not in members)
                    break
        return tuple(members)

    found: set[str] = set()
    for o_token in tokens("o", o, missing_all=True):
        for h1_token in tokens("h1", h1, missing_all=False):
            for h2_token in tokens("h2", h2, missing_all=False):
                for key in _formula_keys(o_token, h1_token, h2_token):
                    found.update(index.get(key, ()))
    return tuple(sorted(found))


def _evaluate_v2(calls: Mapping[str, GeneCall], document: SerotypeDocumentV2) -> SchemeCall:
    """Ranked antigen families, then the formula table."""
    control = _control_present(calls, document)
    empty = _FamilyCall(name=None, score=0.0, genes=(), perfect=False, ambiguous=False)
    resolved: dict[str, _FamilyCall] = {
        family.key: _family_call(family, calls, document) for family in document.antigens
    }
    o = resolved.get("o", empty)
    h1 = resolved.get("h1", empty)
    h2 = resolved.get("h2", empty)
    components = tuple(
        ScoreComponent(name=f"antigen:{gene}", score=call.score)
        for _, call in sorted(resolved.items())
        for gene in call.genes
    )
    notes = [
        f"formula {o.name or '-'}:{h1.name or '-'}:{h2.name or '-'}",
        *(
            f"{family_key.upper()}: {call.name or '-'}"
            + ("" if call.perfect else " (not a perfect match)")
            for family_key, call in sorted(resolved.items())
        ),
    ]
    # A family whose hits disagreed on the antigen name cannot be a perfect call.
    exact = o.perfect and h1.perfect and h2.perfect
    if not control:
        return SchemeCall(
            phenotype=document.fallback,
            score=0.0,
            confidence="low",
            notes=("control gene absent",),
        )
    o_names = tuple(
        entry.name for family in document.antigens if family.key == "o" for entry in family.entries
    )
    serovars = (
        _lookup_formula(document, o.ties, h1.ties, h2.ties, o_names=o_names) if control else ()
    )
    observed = (o.name, h1.name or FALLBACK_PHASE, h2.name or FALLBACK_PHASE)
    exact_row = any((row.o, row.h1, row.h2) == observed for row in document.formulas)
    score = max(o.score, h1.score, h2.score)
    if not exact_row and serovars:
        # The observed formula is not in the table: the name came from widening a
        # sibling's class (a cross-labelled reference) or from another reference
        # covering the phase, which is weaker evidence than reading the formula.
        notes.append("matched after widening the antigen names to their class")
    if len(serovars) == 1:
        return SchemeCall(
            phenotype=serovars[0],
            score=round(score, 4),
            confidence="high" if exact and exact_row else "low",
            components=components,
            notes=tuple(notes),
        )
    if len(serovars) > 1:
        return SchemeCall(
            phenotype=None,
            score=round(score, 4),
            confidence="ambiguous",
            components=components,
            ambiguous=tuple(
                PhenotypeScore(phenotype=serovar, score=round(score, 4)) for serovar in serovars
            ),
            notes=(*notes, "formula matched several serovars after widening"),
        )
    if o.name is None:
        return SchemeCall(
            phenotype=document.unresolved,
            score=0.0,
            confidence="low",
            components=components,
            notes=tuple(notes),
        )
    formula = f"{o.name}:{h1.name or '-'}:{h2.name or '-'}"
    return SchemeCall(
        phenotype=document.partial.format(formula=formula),
        score=round(max(o.score, h1.score, h2.score), 4),
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
