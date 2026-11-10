"""Evaluate gspec.serotype/1 rules over one file's gene calls."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Literal

from gspec.align.sketch import SKETCH_MARGIN
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


def _partial_evidence(call: GeneCall, document: SerotypeDocumentV2) -> bool:
    """A near-identical but truncated hit: the assembly broke the locus.

    Only the second (weaker) acceptance rule qualifies, and only when the
    document enables it; such a hit still ranks below a full-length one because
    the score is identity x coverage, and it can never be perfect.
    """
    if document.partial_coverage_floor <= 0:
        return False
    return (
        call.coverage_pct >= document.partial_coverage_floor
        and call.identity_pct >= document.partial_identity_floor
    )


def _evidence(call: GeneCall) -> float:
    """The ranking score: identity x coverage, in [0, 1].

    The same metric the `best_hit` typing rule ranks by, so a serotype profile
    and a typing scheme order competing alleles identically.
    """
    return (call.identity_pct / 100.0) * (call.coverage_pct / 100.0)


def _antigen_score(call: GeneCall, document: SerotypeDocumentV2) -> float:
    """How strongly one allele supports its antigen.

    Usually identity x coverage. A hit that only qualifies through the partial
    rule is truncated by the assembly, not by nature, so it is ranked by the
    identity it does show — otherwise a 99.85 %-identical 52 %-covered `wzx`
    loses to a 95 %-identical full-length cross-reactive sibling and the
    serogroup flips (measured on a real D1 genome). It still cannot be perfect.
    """
    if _partial_evidence(call, document) and not _present(
        call, document.coverage_floor, document.identity_floor
    ):
        return call.identity_pct / 100.0
    return _evidence(call)


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
class _Ranked:
    """One antigen's claim on a family: score, supporting alleles, sketch, length."""

    name: str
    score: float
    genes: tuple[str, ...]
    sketch: float | None
    #: True when the supporting hit cleared the strict floors (a partial hit is
    #: truncated by the assembly, and must not be trusted to the sketch).
    full: bool


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
    scored: list[_Ranked] = []
    for entry in family.entries:
        best = 0.0
        winners: tuple[str, ...] = ()
        sketch: float | None = None
        full_length = False
        for gene in entry.genes:
            call = calls.get(gene)
            if call is None:
                continue
            complete = _present(call, document.coverage_floor, document.identity_floor)
            if not (complete or _partial_evidence(call, document)):
                continue
            score = _antigen_score(call, document)
            if call.sketch_pct is not None and (sketch is None or call.sketch_pct > sketch):
                sketch = call.sketch_pct
            if score > best:
                best, winners, full_length = score, (gene,), complete
            elif score == best:
                winners += (gene,)
                full_length = full_length or complete
        if best > 0.0:
            scored.append(
                _Ranked(name=entry.name, score=best, genes=winners, sketch=sketch, full=full_length)
            )
    if not scored:
        return _FamilyCall(name=None, score=0.0, genes=(), perfect=False, ambiguous=False)
    scored.sort(key=lambda item: item.score, reverse=True)
    names = {item.name for item in scored}
    cutoff = scored[0].score - document.margin
    tied = _sketch_break([item for item in scored if item.score >= cutoff])
    best = tied[0]
    ties = tuple(dict.fromkeys(item.name for item in tied))
    return _FamilyCall(
        name=best.name,
        score=best.score,
        genes=best.genes,
        perfect=best.score >= 1.0 and len(names) == 1,
        ambiguous=len(ties) > 1,
        ties=ties,
    )


def _sketch_break(tied: list[_Ranked]) -> list[_Ranked]:
    """Narrow a margin tie with the anchored allele sketch, when it decides.

    Alignment cannot separate alleles that differ by a SNP under a gapped HSP
    (measured: `fliC` `e,h` vs `e,n,x` scored 0.9103 vs 0.9152), which is why the
    tie is normally widened and reported as candidates. When every tied antigen
    carries a sketch score, the sketch decides the order of magnitude that
    matters — whole-allele k-mer containment — so a clear sketch winner (a gap
    wider than ``SKETCH_MARGIN``) collapses the tie to that antigen. A tie the
    sketch cannot separate stays a tie: the candidates remain visible.
    """
    if len(tied) < 2:
        return tied
    if any(not item.full for item in tied):
        # A truncated hit's window covers only part of the allele, so its
        # containment tracks coverage, not allele identity: the sketch would
        # reward whichever antigen the assembly happened to leave intact. Real
        # case: `wzx` A and D1 both truncated at ~52 % of the same window, where
        # a sketch decided A on a 0.15-point identity difference.
        return tied
    sketches = [item.sketch for item in tied]
    if any(sketch is None for sketch in sketches):
        return tied
    top = max(sketch for sketch in sketches if sketch is not None)
    decided = [item for item in tied if (item.sketch or 0.0) >= top - 100.0 * SKETCH_MARGIN]
    if not decided:
        return tied
    # The sketch winner leads, so it names the family even when alignment ranked
    # it second: that is the whole point of the tiebreak.
    return sorted(decided, key=lambda item: (item.sketch or 0.0, item.score), reverse=True)


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


#: Two hits on one contig whose intervals overlap by more than this fraction of
#: the shorter one are the same physical gene.
SAME_LOCUS = 0.5

#: A family's other locus counts as an alternative only when it scores at least
#: this fraction of the contested one. Alleles of one antigen within a phase sit
#: inside ~1 %, so 2 % is the band: a wider one lets a weaker stray hit stand in
#: for a real gene (measured: a 10 % band cost Braenderup its `e,h` call, whose
#: own locus scored 0.86 against a 0.913 contested one).
_PHASE_ALTERNATIVE = 0.98


_PHASE_FAMILIES = ("h1", "h2")


def _separate_phases(
    calls: Mapping[str, GeneCall], document: SerotypeDocumentV2
) -> dict[str, GeneCall]:
    """Keep one physical flagellin gene out of both phase families.

    SISTR's `fliC` and `fljB` pools both carry an `l,v` allele and the two are
    ~99 % identical, so a genome's phase-1 gene is scored by the phase-2 pool as
    well (measured: one contig hit `h1 l,v` and `h2 l,v` at the same coordinates).
    That spurious phase-2 hit can outrank the real fljB and the serotype flips —
    SISTR reports `l,v` for the same genome.

    A gene is one phase or the other, so each *locus* may support one family. The
    locus is scored by its best allele in each family, and a shared locus goes to
    the family that would lose more by losing it (score here minus the best it
    reaches anywhere else). Real case — a Mendoza genome scores `l,v` 0.997 on its
    phase-1 gene and `1,2` 0.968 on its real `fljB`: phase 1 keeps `l,v` (losing
    it costs 0.1) and phase 2 moves to its own locus (giving up the spurious hit
    costs 0.03).
    """
    if not any(family.key in _PHASE_FAMILIES for family in document.antigens):
        return dict(calls)
    loci, _spans = _loci(calls, document)
    kept = dict(calls)
    for key, families in loci.items():
        if family_a := families.get("h1"):
            family_a_score, family_a_genes = family_a
        else:
            continue
        family_b = families.get("h2")
        if family_b is None:
            continue
        family_b_score, family_b_genes = family_b
        if not _same_antigen(calls, family_a_genes, family_b_genes):
            # Different antigens matched at one locus is ordinary cross-reactivity
            # between flagellin alleles, not a duplicated reference: leave it.
            continue
        loser = _phase_loser(loci, key, ("h1", family_a_score), ("h2", family_b_score))
        if loser is None:
            continue
        for gene in family_a_genes if loser == "h1" else family_b_genes:
            kept.pop(gene, None)
    return kept


def _loci(
    calls: Mapping[str, GeneCall], document: SerotypeDocumentV2
) -> tuple[
    dict[tuple[str, int], dict[str, tuple[float, list[str]]]],
    dict[tuple[str, int], tuple[int, int]],
]:
    """Group each family's hits into loci.

    Returns ``(groups, spans)`` where a group is
    ``(contig, first_start) -> family -> (best score, gene ids)`` and ``spans``
    maps the same key to the locus' interval. Two intervals on one contig that
    share more than ``SAME_LOCUS`` of the shorter one are the same physical gene.
    Hits without a locus (a frozen table source) are ignored.
    """
    groups: dict[tuple[str, int], dict[str, tuple[float, list[str]]]] = {}
    spans: dict[tuple[str, int], tuple[int, int]] = {}
    for gene, call in sorted(
        calls.items(), key=lambda item: (item[1].contig or "", item[1].start, item[1].end)
    ):
        family = gene.partition("_")[0]
        if family not in _PHASE_FAMILIES or call.contig is None:
            continue
        key = _locus_key(call, spans)
        span_start, span_end = spans.get(key, (call.start, call.end))
        spans[key] = (min(span_start, call.start), max(span_end, call.end))
        score = _antigen_score(call, document)
        bucket = groups.setdefault(key, {})
        best, genes = bucket.get(family, (0.0, []))
        bucket[family] = (max(best, score), [*genes, gene])
    return groups, spans


def _locus_key(call: GeneCall, spans: Mapping[tuple[str, int], tuple[int, int]]) -> tuple[str, int]:
    """The existing locus this hit overlaps, or a new key for it."""
    assert call.contig is not None
    width = call.end - call.start or 1
    for key, (start, end) in spans.items():
        if key[0] != call.contig:
            continue
        overlap = min(end, call.end) - max(start, call.start)
        if overlap > 0 and overlap / min(end - start or 1, width) > SAME_LOCUS:
            return key
    return (call.contig, call.start)


def _phase_loser(
    loci: Mapping[tuple[str, int], dict[str, tuple[float, list[str]]]],
    key: tuple[str, int],
    mine: tuple[str, float],
    theirs: tuple[str, float],
) -> str | None:
    """Which family gives up a contested locus, or None to leave both alone.

    Only a family with a *competitive* locus elsewhere yields (``_PHASE_ALTERNATIVE``
    of the contested score): without that guard the correct family's evidence is
    deleted whenever its best other hit is junk — measured on the 304-genome
    corpus, an unguarded version cost 0.6 points. When both families have one, the
    higher score keeps the locus; when neither does, both keep it (a monophasic or
    fragmented genome, where the evidence cannot be separated).
    """
    my_family, my_score = mine
    their_family, their_score = theirs
    mine_covered = _best_elsewhere(loci, key, my_family) >= _PHASE_ALTERNATIVE * my_score
    theirs_covered = _best_elsewhere(loci, key, their_family) >= _PHASE_ALTERNATIVE * their_score
    if mine_covered == theirs_covered:
        if not mine_covered:
            return None
        return my_family if my_score < their_score else their_family
    return my_family if mine_covered else their_family


def _same_antigen(calls: Mapping[str, GeneCall], first: list[str], second: list[str]) -> bool:
    """Do both families name the SAME antigen at this locus?

    The guard exists for one duplication: SISTR's `fliC` and `fljB` pools both
    carry an `l,v` allele, so one physical gene is scored as `h1 l,v` and
    `h2 l,v`. Two *different* names at one locus (`h1 e,h` against
    `h2 e,n,z15`) are genuine cross-reactivity between flagellin alleles — the
    ranker handles that, and dropping either would lose the right name.
    """
    mine = {_antigen_name(calls, gene) for gene in first}
    theirs = {_antigen_name(calls, gene) for gene in second}
    return bool(mine & theirs)


def _antigen_name(calls: Mapping[str, GeneCall], gene: str) -> str:
    """``h1_l_v__3`` -> ``l,v``."""
    del calls  # the name lives in the gene id
    antigen = gene.rsplit("__", 1)[0]
    return antigen.partition("_")[2].replace("_", ",")


def _best_elsewhere(
    loci: Mapping[tuple[str, int], dict[str, tuple[float, list[str]]]],
    key: tuple[str, int],
    family: str,
) -> float:
    """A family's best locus score other than the contested one."""
    return max(
        (
            families[family][0]
            for other, families in loci.items()
            if other != key and family in families
        ),
        default=0.0,
    )


def _evaluate_v2(calls: Mapping[str, GeneCall], document: SerotypeDocumentV2) -> SchemeCall:
    """Ranked antigen families, then the formula table."""
    control = _control_present(calls, document)
    calls = _separate_phases(calls, document)
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
