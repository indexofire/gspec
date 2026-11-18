"""Evaluate gspec.serotype/1 rules over one file's gene calls."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
from typing import Literal, NamedTuple

from gspec.align.sketch import SKETCH_MARGIN
from gspec.cluster.core import GeneCall
from gspec.typing.results import PhenotypeScore, SchemeCall, ScoreComponent
from gspec.typing.serotype_models import (
    AntigenEntry,
    AntigenFamily,
    AnySerotypeDocument,
    ODiscriminator,
    SerotypeDocument,
    SerotypeDocumentV2,
    SerotypeRule,
    o_discriminator_specs,
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
    #: The winning allele's own evidence, surfaced in the call's notes so the
    #: alignment detail a lab audits is in the answer, not only in the screen
    #: table (identity of the best allele; coverage of the same call).
    identity_pct: float = 0.0
    coverage_pct: float = 0.0
    #: The winning allele's reference length (subject length), for the O
    #: family's length-aware tie-break: different serogroups have different-
    #: length wzx/wzy genes, so a longer full gene beats a shorter cross-reactor.
    length: float = 0.0
    #: The winning hit's genomic locus, so the O reconciliation can tell a
    #: cross-reacting paralog at ANOTHER locus (length is evidence there) from a
    #: same-locus reference indel (length is not).
    contig: str | None = None
    start: int = 0
    end: int = 0


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
    #: The names sitting at the winner's EXACT score that are not similarity
    #: class mates (a class mate is a cross-label of the same antigen): every
    #: measured channel tied, so no mechanism chose between them. A singleton
    #: means the family is decided; a longer tuple feeds the equal-tie lookup
    #: pass instead of letting declaration order pick a winner.
    equal_ties: tuple[str, ...] = ()
    #: The winning allele's alignment evidence (0.0 on an absent family), for
    #: the per-family note: which allele named the antigen, at what identity
    #: and coverage.
    identity_pct: float = 0.0
    coverage_pct: float = 0.0
    #: The winning allele's reference length, for the O-family wzx/wzy consensus
    #: (a longer gene is the stronger evidence, SISTR's bitscore ordering).
    length: float = 0.0
    #: The winning hit's locus (``contig``/``start``/``end``), so the O
    #: reconciliation can tell a cross-reacting paralog at another locus from a
    #: same-locus reference indel.
    contig: str | None = None
    start: int = 0
    end: int = 0
    #: The winning hit cleared the STRICT floors (not the truncated-assembly
    #: rule): a locus the assembly failed to close is real evidence but is not
    #: whole, so it caps confidence below ``high``.
    full: bool = True
    #: No antigen name holds the winner's EXACT score besides the winner itself
    #: (``equal_ties``): a labelling conflict resolved by declaration order, which
    #: is again not ``high``-grade evidence. Cross-reactions that merely *score*
    #: something are not a conflict — requiring "only one antigen scored at all"
    #: is what made the old ``perfect`` flag unreachable (measured: 8/100).
    single_name: bool = True


def _same_class(
    document: SerotypeDocumentV2, family: Literal["o", "h1", "h2"], first: str, second: str
) -> bool:
    """Do two antigen names of one family share a similarity class?

    A class mate is a spelling or cross-label variant of the same antigen
    (SISTR merges `g,m` with `g,m,q`; the reference audit merged `g,p` with
    `g,p,s`/`g,p,u`), so an exact score tie between class mates is a labelling
    artefact, not two antigens at identical evidence.
    """
    for group in document.similarity_groups.get(family, ()):
        if first in group and second in group:
            return True
    return False


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
        evidence: GeneCall | None = None
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
                best, winners, full_length, evidence = score, (gene,), complete, call
            elif score == best:
                winners += (gene,)
                full_length = full_length or complete
        if best > 0.0:
            assert evidence is not None  # best > 0 implies a call was recorded
            scored.append(
                _Ranked(
                    name=entry.name,
                    score=best,
                    genes=winners,
                    sketch=sketch,
                    full=full_length,
                    identity_pct=evidence.identity_pct,
                    coverage_pct=evidence.coverage_pct,
                    length=float(evidence.subject_len or 0.0),
                    contig=evidence.contig,
                    start=evidence.start,
                    end=evidence.end,
                )
            )
    if not scored:
        return _FamilyCall(name=None, score=0.0, genes=(), perfect=False, ambiguous=False)
    scored.sort(key=lambda item: item.score, reverse=True)
    if family.key in _PHASE_FAMILIES and scored[0].score < document.phase_floor:
        # Below the naming bar a phase is absent, not weak: monophasic genomes
        # carry cross-reactions that would otherwise invent a second phase.
        return _FamilyCall(
            name=None, score=scored[0].score, genes=(), perfect=False, ambiguous=False
        )
    names = {item.name for item in scored}
    cutoff = scored[0].score - document.margin
    tied = _sketch_break([item for item in scored if item.score >= cutoff])
    if family.key == "o" and len(tied) > 1:
        # Length-aware tie-break for the O family: within the margin, a full
        # longer gene (the real wzx/wzy) beats a shorter cross-reacting paralog
        # (measured: a 1539 bp E1 wzx vs a 1293 bp B-like paralog at 99.93 %).
        tied = sorted(tied, key=lambda item: (item.length, item.score), reverse=True)
    best = tied[0]
    ties = tuple(dict.fromkeys(item.name for item in tied))
    equal_ties = tuple(
        dict.fromkeys(
            item.name
            for item in tied
            if item.score == best.score
            and (
                item.name == best.name
                or not _same_class(document, family.key, best.name, item.name)
            )
        )
    )
    return _FamilyCall(
        name=best.name,
        score=best.score,
        genes=best.genes,
        perfect=best.score >= 1.0 and len(names) == 1,
        ambiguous=len(ties) > 1,
        ties=ties,
        equal_ties=equal_ties,
        identity_pct=best.identity_pct,
        coverage_pct=best.coverage_pct,
        length=best.length,
        contig=best.contig,
        start=best.start,
        end=best.end,
        full=best.full,
        single_name=len(equal_ties) <= 1,
    )


def _same_locus(a: _FamilyCall, b: _FamilyCall) -> bool:
    """Do two family calls rest on the same physical gene?

    Same contig and overlapping intervals. That is what separates a
    cross-reacting paralog sitting elsewhere in the genome (length IS evidence:
    the real gene is the longer, fully-covered one) from two references of the
    same gene that differ by an indel (length is NOT evidence: a 1437 bp subset
    reference is covered perfectly by a 1539 bp superset allele).
    """
    if a.contig is None or b.contig is None or a.contig != b.contig:
        return False
    overlap = min(a.end, b.end) - max(a.start, b.start)
    if overlap <= 0:
        return False
    shorter = min(a.end - a.start or 1, b.end - b.start or 1)
    return overlap / shorter > SAME_LOCUS


def _numeric_o(
    document: SerotypeDocumentV2, serogroup: str | None, subspecies: str = "enterica"
) -> str | None:
    """The numeric O antigen for a serogroup letter, when exactly one applies.

    The build maps a letter only when every row of that serogroup reduces to one O
    formula (38 of 48 do: ``P`` is O:38, ``D1`` is 9,12). The rest — ``C2-C3`` is
    6,8 *or* 8,20 *or* 8 *or* 6,8,20, ``G`` is 13,22 *or* 13,23 — stay letters,
    because picking one would be a guess.
    """
    if serogroup is None:
        return None
    per_subspecies = document.o_numeric.get(subspecies) or document.o_numeric.get("enterica", {})
    return per_subspecies.get(serogroup)


def _drop_duplicate_phase(
    h1: _FamilyCall, h2: _FamilyCall
) -> tuple[_FamilyCall, _FamilyCall, str | None, bool]:
    """H1 and H2 cannot carry the same antigen — drop the weaker phase.

    The two phases are different flagellin genes (`fliC` / `fljB`), so an identical
    antigen name in both is a cross-reaction, not a serovar: **no WHO-2007 row has
    H1 == H2** (0 of 2541 rows, every subspecies). The panel used to report such a
    formula anyway — measured on a real genome, `P:k:k` came from an `h2 k` hit at
    only **87.15 %** identity while `h1 k` was 100 %.

    The weaker phase loses and the survivor is returned unchanged. A tie cannot
    say which one is the cross-reaction, so the drop is then by convention and the
    caller must cap confidence: the third element is the note to add, the fourth
    says "this needed a resolution".
    """
    if h1.name is None or h2.name is None or h1.name != h2.name:
        return h1, h2, None, False
    h2_loses = h1.score >= h2.score
    dropped, kept = (h2, h1) if h2_loses else (h1, h2)
    label, other = ("H2", "H1") if h2_loses else ("H1", "H2")
    tie = h1.score == h2.score
    note = (
        f"{label} {dropped.name!r} repeats {other} (no WHO row has H1 == H2):"
        f" phase dropped as a cross-reaction ({dropped.score:.4f} vs {kept.score:.4f})"
    )
    blank = replace(
        dropped,
        name=None,
        score=0.0,
        genes=(),
        ties=(),
        equal_ties=(),
        ambiguous=False,
        identity_pct=0.0,
        coverage_pct=0.0,
    )
    return (kept, blank, note, tie) if h2_loses else (blank, kept, note, tie)


def _reconcile_o(wzx: _FamilyCall, wzy: _FamilyCall) -> _FamilyCall:
    """SISTR's wzx/wzy reconciliation: agreement wins, one missing defers, and
    a disagreement goes to the longer (bitscore-like) evidence — but only when
    the two hits sit at different loci.

    At one locus the length comparison is meaningless (two references of the
    same gene, one a truncated subset of the other), so the disagreement is kept
    as a reported ambiguity instead of being decided by gene length.
    """
    if wzx.name is None:
        return wzy
    if wzy.name is None:
        return wzx
    if wzx.name == wzy.name:
        merged_ties = tuple(dict.fromkeys((*wzx.ties, *wzy.ties)))
        return replace(wzx, ties=merged_ties, ambiguous=len(merged_ties) > 1)
    if _same_locus(wzx, wzy):
        # Same physical gene: keep both names live so a discriminator (or the
        # honest "ambiguous" answer) decides, never the gene's length.
        merged_ties = tuple(dict.fromkeys((*wzx.ties, *wzy.ties)))
        return replace(
            wzx if wzx.score >= wzy.score else wzy,
            ties=merged_ties,
            ambiguous=len(merged_ties) > 1,
        )
    wzx_evidence = wzx.score * wzx.length
    wzy_evidence = wzy.score * wzy.length
    return wzx if wzx_evidence >= wzy_evidence else wzy


def _family_call_o(
    family: AntigenFamily, calls: Mapping[str, GeneCall], document: SerotypeDocumentV2
) -> _FamilyCall:
    """The O family, ranked per locus (wzx/wzy) and reconciled SISTR-style.

    wzx and wzy are independent O-antigen genes: a cross-reacting wzy must not
    override a clear wzx (and vice versa). Each locus ranks on its own; the two
    calls then reconcile (see :func:`_reconcile_o`). An O family without locus
    tags (a legacy panel) ranks as one mixed family, the pre-split behaviour.
    """
    by_locus: dict[str, list[AntigenEntry]] = {}
    for entry in family.entries:
        by_locus.setdefault(entry.locus or "", []).append(entry)
    if "wzx" not in by_locus or "wzy" not in by_locus:
        return _family_call(family, calls, document)
    wzx = _family_call(AntigenFamily(key="o", entries=tuple(by_locus["wzx"])), calls, document)
    wzy = _family_call(AntigenFamily(key="o", entries=tuple(by_locus["wzy"])), calls, document)
    return _reconcile_o(wzx, wzy)


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


def _formula_keys(
    o: str, h1: str | None, h2: str | None, *, alias_h2_absent: bool = False
) -> list[tuple[str, str, str]]:
    """The one formula a call stands for.

    Monophasic serovars are tabulated with ``"-"``, and an unresolved phase is
    already passed in as ``"-"``, so the exact key covers both cases. Substituting
    ``"-"`` for a phase that WAS resolved is wrong and expensive: a monophasic
    `i:-` genome then also matched every `B:-:-` serovar (measured on the in-house
    corpus, `I 1,4,[5],12:i:-` arrived with 15 unrelated candidates).
    """
    keys = [(o, h1 or "-", h2 or "-")]
    if alias_h2_absent and h2 and h2 != "-":
        # A weak h2 may be cross-reaction rather than a real phase: monophasic
        # serovars are tabulated with `-`, so the h2-absent row is the fallback.
        keys.append((o, h1 or "-", "-"))
    return keys


FALLBACK_PHASE = "-"


def _lookup_formula(
    document: SerotypeDocumentV2,
    o: tuple[str, ...],
    h1: tuple[str, ...],
    h2: tuple[str, ...],
    *,
    o_names: tuple[str, ...] = (),
    alias_h2_absent: bool = False,
    subspecies: str = "enterica",
    o_exact: bool = False,
) -> tuple[str, ...]:
    """The serovar(s) for a formula: exact first, then similarity widening.

    A single serovar is a call; several are genuine ambiguity (the upstream tools
    pipe-join them, gspec reports candidates); none is a partial formula.
    """
    index: dict[tuple[str, str, str], list[str]] = {}
    for row in document.formulas:
        if row.subspecies != subspecies:
            continue
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
            if family == "o" and o_exact:
                # A discriminator already settled this O group; widening it to
                # its similarity class would re-open the very ambiguity it
                # resolved (A/D1/D2 share one class, so a decided D2 would reach
                # every D1 row again).
                continue
            for group in document.similarity_groups.get(family, ()):
                if name in group:
                    members.extend(member for member in group if member not in members)
                    break
        return tuple(members)

    def search(keys_alias: bool) -> set[str]:
        found: set[str] = set()
        for o_token in tokens("o", o, missing_all=True):
            for h1_token in tokens("h1", h1, missing_all=False):
                for h2_token in tokens("h2", h2, missing_all=False):
                    keys = _formula_keys(o_token, h1_token, h2_token, alias_h2_absent=keys_alias)
                    for key in keys:
                        found.update(index.get(key, ()))
        return found

    found = search(keys_alias=False)
    if not found and alias_h2_absent:
        # The exact formula has no row: only then does a weak h2 fall back to the
        # h2-absent (monophasic) row. Merging both keys unconditionally made every
        # imperfect h2 double the candidate list.
        found = search(keys_alias=True)
    return tuple(sorted(found))


def _exact_serovars(
    document: SerotypeDocumentV2,
    observed: tuple[str | None, str, str],
    subspecies: str = "enterica",
) -> list[str]:
    """The serovars whose row IS the observed formula, declaration order first.

    The lookup's similarity widening exists for formulas the table does not
    carry; when the exact row exists, widening only adds sibling serovars the
    antigens did not name (measured: a Muenchen genome arrived as four
    candidates because the ``1,2`` class also reaches ``1,5`` rows) — the
    opposite of helpful. The first-declared name of a shared formula is the
    primary: the WHO table lists the common serovar first (``C2-C3:d:1,2`` is
    Muenchen, then Virginia), a convention the table itself carries.
    """
    o, h1, h2 = observed
    if o is None:
        return []
    return [
        row.serovar
        for row in document.formulas
        if (row.o, row.h1, row.h2) == (o, h1, h2) and row.subspecies == subspecies
    ]


def _formula_rows(
    document: SerotypeDocumentV2,
    o_tokens: tuple[str, ...],
    h1_tokens: tuple[str, ...],
    h2_tokens: tuple[str, ...],
    subspecies: str = "enterica",
) -> tuple[str, ...]:
    """Serovars whose (O, H1, H2) key sits in the three token sets, declaration
    order — a pure key lookup with no similarity widening."""
    found: list[str] = []
    for row in document.formulas:
        if (
            row.o in o_tokens
            and row.h1 in h1_tokens
            and row.h2 in h2_tokens
            and row.subspecies == subspecies
            and row.serovar not in found
        ):
            found.append(row.serovar)
    return tuple(found)


def _lookup_subspecies_serovar(
    document: SerotypeDocumentV2,
    subspecies: str,
    h1: str | None,
    h2: str | None,
    o_group: str | None = None,
) -> tuple[str, ...]:
    """The non-enterica serovar(s) whose row matches (subspecies, H1, H2), and
    optionally the O serogroup. The O-antigen references are enterica-specific
    and must not resolve a non-enterica O group on their own (a VI O:6,7 wzx is
    100 % identical to enterica O:11), so the subspecies' own formula rows —
    which carry the full ``NUMERAL O:H1:H2`` name — are the primary key; the O
    group only disambiguates when several rows share the H phases.
    """
    h1_key = h1 or FALLBACK_PHASE
    h2_key = h2 or FALLBACK_PHASE
    found: list[str] = []
    for row in document.formulas:
        if (
            row.subspecies == subspecies
            and row.h1 == h1_key
            and row.h2 == h2_key
            and (o_group is None or row.o == o_group)
            and row.serovar not in found
        ):
            found.append(row.serovar)
    return tuple(found)


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
) -> tuple[dict[str, GeneCall], dict[str, tuple[str, ...]]]:
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
        return dict(calls), {}
    loci, _spans = _loci(calls, document)
    kept = dict(calls)
    demoted: dict[str, list[str]] = {}
    for key, families in loci.items():
        if "h1" not in families or "h2" not in families:
            continue
        family_a_score, family_a_genes, phase_a = families["h1"]
        family_b_score, family_b_genes, phase_b = families["h2"]
        if not _same_antigen(family_a_genes, family_b_genes):
            # Different antigens matched at one locus is ordinary cross-reactivity
            # between flagellin alleles, not a duplicated reference: leave it.
            continue
        loser = _phase_by_ends(phase_a, phase_b) or _phase_loser(
            loci, key, ("h1", family_a_score), ("h2", family_b_score)
        )
        if loser is None:
            continue
        for gene in family_a_genes if loser == "h1" else family_b_genes:
            kept.pop(gene, None)
            # The locus may only be one phase, but the antigen it named is still
            # a lead: keep it as a low-confidence candidate instead of erasing it
            # (measured: Braenderup's `e,h` references fit the `fljB` locus better
            # than its own `fliC`, and deleting that name lost the serovar).
            demoted.setdefault(loser, []).append(_antigen_name(gene))
    return kept, {family: tuple(dict.fromkeys(names)) for family, names in demoted.items()}


def _loci(
    calls: Mapping[str, GeneCall], document: SerotypeDocumentV2
) -> tuple[
    dict[tuple[str, int], dict[str, tuple[float, list[str], float | None]]],
    dict[tuple[str, int], tuple[int, int]],
]:
    """Group each family's hits into loci.

    Returns ``(groups, spans)`` where a group is
    ``(contig, first_start) -> family -> (best score, gene ids, best phase)`` and
    ``spans``
    maps the same key to the locus' interval. Two intervals on one contig that
    share more than ``SAME_LOCUS`` of the shorter one are the same physical gene.
    Hits without a locus (a frozen table source) are ignored.
    """
    groups: dict[tuple[str, int], dict[str, tuple[float, list[str], float | None]]] = {}
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
        best, genes, phase = bucket.get(family, (0.0, [], None))
        if call.phase_pct is not None and (phase is None or call.phase_pct > phase):
            phase = call.phase_pct
        bucket[family] = (max(best, score), [*genes, gene], phase)
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


#: Below this score a family's surviving locus is too weak to justify erasing a
#: phase-demoted antigen name (see `_evaluate_v2`).
DEMOTED_SCORE_FLOOR = 0.9

#: Confidence floors on the weakest PRESENT antigen's score (identity x coverage).
#: Calibrated on the 100-genome corpus: with ``HIGH_SCORE_FLOOR = 0.99`` the
#: resulting ``high`` tier holds 23/100 genomes and every one of them agrees with
#: both reference tools; ``MEDIUM_SCORE_FLOOR = 0.90`` is where the calls that
#: still name something stop being defensible (the `low` tier is 3/5 correct).
HIGH_SCORE_FLOOR = 0.99
MEDIUM_SCORE_FLOOR = 0.90

#: The six S. enterica subspecies' Kauffmann-White roman numerals. The
#: subspecies itself is a whole-genome property decided by the k-mer signature
#: (``gspec.subspecies``); the serotype engine only renders it.
Confidence = Literal["high", "medium", "ambiguous", "low"]

SUBSPECIES_NUMERAL = {
    "enterica": "I",
    "salamae": "II",
    "arizonae": "IIIa",
    "diarizonae": "IIIb",
    "houtenae": "IV",
    "indica": "VI",
}

#: The synthetic gene id prefix the screen pipeline uses to hand the subspecies
#: into the typing engine (``subspecies_indica``, ``subspecies_diarizonae``...).
SUBSPECIES_GENE_PREFIX = "subspecies_"


def _detect_subspecies(calls: Mapping[str, GeneCall]) -> str | None:
    """The subspecies name threaded through a synthetic gene call, if present."""
    for gene in calls:
        if gene.startswith(SUBSPECIES_GENE_PREFIX):
            name = gene[len(SUBSPECIES_GENE_PREFIX) :]
            if name in SUBSPECIES_NUMERAL:
                return name
    return None


#: The end-based phase call needs this much separation (percentage points) to
#: overrule the score ordering; below it the ends are as ambiguous as the middle.
PHASE_END_GAP = 5.0


def _phase_by_ends(phase_a: float | None, phase_b: float | None) -> str | None:
    """Which family a contested locus really is, from the ends of the window.

    `fliC` and `fljB` `l,v` alleles share ~99 % of their middle — which is why one
    gene scores in both families — but differ at the 5'/3' ends (the head/tail
    markers SeqSero2 uses). Each family's hits carry how well the locus' ends match
    that family's own alleles, so the comparison is direct evidence rather than the
    score ordering the marginal-gain rule has to infer from.
    """
    if phase_a is None or phase_b is None:
        return None
    if abs(phase_a - phase_b) < PHASE_END_GAP:
        return None
    return "h2" if phase_a > phase_b else "h1"


def _phase_loser(
    loci: Mapping[tuple[str, int], dict[str, tuple[float, list[str], float | None]]],
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


def _same_antigen(first: list[str], second: list[str]) -> bool:
    """Do both families name the SAME antigen at this locus?

    The guard exists for one duplication: SISTR's `fliC` and `fljB` pools both
    carry an `l,v` allele, so one physical gene is scored as `h1 l,v` and
    `h2 l,v`. Two *different* names at one locus (`h1 e,h` against
    `h2 e,n,z15`) are genuine cross-reactivity between flagellin alleles — the
    ranker handles that, and dropping either would lose the right name.
    """
    mine = {_antigen_name(gene) for gene in first}
    theirs = {_antigen_name(gene) for gene in second}
    return bool(mine & theirs)


def _antigen_name(gene: str) -> str:
    """``h1_l_v__3`` -> ``l,v``."""
    antigen = gene.rsplit("__", 1)[0]
    return antigen.partition("_")[2].replace("_", ",")


def _best_elsewhere(
    loci: Mapping[tuple[str, int], dict[str, tuple[float, list[str], float | None]]],
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


class _DiscriminatorResult(NamedTuple):
    """The O call after (or without) the SNP-level tie-break."""

    family: _FamilyCall
    notes: list[str]


def _discriminate_o(
    family_call: _FamilyCall,
    calls: Mapping[str, GeneCall],
    document: SerotypeDocumentV2,
) -> _DiscriminatorResult | None:
    """Break an O-family tie with SNP-level markers, when a gap decides.

    Several pairs need a marker (A/D1 on `rfbS`, D1/D2 on `wbaV`), so the
    document carries a list; each is consulted in order and the first whose
    contested names are actually tied speaks. Returns None when no discriminator
    applies, or when the markers disagree too little to decide — the last case is
    common on short-read assemblies, where a partial locus shows a few tenths of
    a point and the honest answer stays "ambiguous".
    """
    specs = o_discriminator_specs(document)
    if not specs or len(family_call.ties) < 2:
        return None
    result: _DiscriminatorResult | None = None
    for spec in specs:
        result = _apply_discriminator(family_call, calls, spec)
        if result is None:
            continue
        # A spec that actually decided ends the search; a "cannot decide" note
        # is kept so the next spec still gets a chance to speak.
        if result.family.name != family_call.name or result.family.ties != family_call.ties:
            return result
        family_call = result.family
    return result


def _apply_discriminator(
    family_call: _FamilyCall,
    calls: Mapping[str, GeneCall],
    spec: ODiscriminator,
) -> _DiscriminatorResult | None:
    """One discriminator's verdict, or None when the tie does not span its names."""
    contested = [name for name in family_call.ties if name in spec.favours]
    if len(contested) < 2:
        return None
    scored: dict[str, float] = {}
    thin: list[str] = []
    for name in contested:
        # A group scores by the BEST of its alleles: a group's own diversity can
        # rival the between-group gap, so one reference per group would decide on
        # which strain happened to be curated.
        best: float | None = None
        for gene in spec.favours[name]:
            call = calls.get(gene)
            if call is None or call.coverage_pct < spec.min_coverage:
                continue
            if best is None or call.identity_pct > best:
                best = call.identity_pct
        if best is None:
            thin.append(name)
            continue
        scored[name] = best
    label = "/".join(contested)
    marker = spec.name or spec.note or "discriminator"
    if len(scored) < 2:
        reason = (
            f"{label} not separated: the discriminator is present for"
            f" {', '.join(sorted(scored)) or 'neither group'} only"
        )
        if thin:
            reason += f" ({', '.join(thin)} marker under {spec.min_coverage:.0f} % coverage)"
        return _DiscriminatorResult(family_call, [reason])
    ranked: list[tuple[str, float]] = sorted(scored.items(), key=lambda item: item[1], reverse=True)
    winner = ranked[0][0]
    loser = ranked[1][0]
    gap = ranked[0][1] - ranked[1][1]
    if gap < spec.min_gap_pct:
        return _DiscriminatorResult(
            family_call,
            [
                f"{winner}/{loser} not separated: {marker}"
                f" gap {gap:.2f} pp below {spec.min_gap_pct:.2f}"
            ],
        )
    return _DiscriminatorResult(
        # The marker settled the pair, so the family is decided again — the
        # exact-formula call may fire with the marker's name.
        replace(family_call, name=winner, ties=(winner,), equal_ties=(winner,)),
        [f"{winner} over {loser}: {marker} gap {gap:.2f} pp"],
    )


def _sum_gene_evidence(
    genes: tuple[str, ...], calls: Mapping[str, GeneCall], document: SerotypeDocumentV2
) -> float:
    """Sum of identity x coverage over the genes that clear the strict floors.

    The O:22/O:23 sub-factor is a three-gene verdict (SeqSero2 sums the scores
    of its O:22-gene1/2/3 and O:23-gene1/2/3), so a single weak hit cannot
    outweigh two strong ones the way a best-score maximum would.
    """
    total = 0.0
    for gene in genes:
        call = calls.get(gene)
        if call is not None and _present(call, document.coverage_floor, document.identity_floor):
            total += _evidence(call)
    return total


def _best_gene_evidence(
    genes: tuple[str, ...], calls: Mapping[str, GeneCall], document: SerotypeDocumentV2
) -> float:
    """Best identity x coverage among the genes that clear the strict floors."""
    return max(
        (
            _evidence(calls[gene])
            for gene in genes
            if gene in calls
            and _present(calls[gene], document.coverage_floor, document.identity_floor)
        ),
        default=0.0,
    )


def _break_o_subfactor(
    serovars: tuple[str, ...] | list[str],
    calls: Mapping[str, GeneCall],
    document: SerotypeDocumentV2,
) -> tuple[str | None, str] | None:
    """Break an O:22/O:23 serovar tie with the dedicated sub-factor genes.

    Returns ``(winner, note)`` when the pair is recognised and the genes decide,
    ``(None, note)`` when the pair is recognised but the evidence is too thin or
    too close to call, and ``None`` when the tie does not span a known pair.
    """
    spec = document.o_subfactor
    if spec is None or len(serovars) < 2:
        return None
    pair = next((p for p in spec.pairs if set(p) == set(serovars)), None)
    if pair is None:
        return None
    o22 = _sum_gene_evidence(spec.o22_genes, calls, document)
    o23 = _sum_gene_evidence(spec.o23_genes, calls, document)
    loser = pair[1] if o22 > o23 else pair[0]
    winner = pair[0] if o22 > o23 else pair[1]
    gap = abs(o22 - o23)
    if not spec.o22_genes:
        # One-sided: only the O:23 sub-factor genes are curated, so their
        # presence proves O:23 and their absence leaves O:22 (a documented
        # limitation until a curated 13,22 reference supplies the O:22 side).
        if o23 <= 0.0:
            return (
                pair[0],
                f"{pair[0]} over {pair[1]}: no O:23 sub-factor gene cleared the"
                f" floors (one-sided O:23 marker, absence assumed O:22)",
            )
        return (
            pair[1],
            f"{pair[1]} over {pair[0]}: O:23 sub-factor genes {o23:.3f}"
            f" ({spec.note or 'dedicated O:23 markers'})",
        )
    if o22 <= 0.0 and o23 <= 0.0:
        return (
            None,
            f"{pair[0]}/{pair[1]} not separated: no O:22/O:23 sub-factor gene cleared the floors",
        )
    if gap < spec.min_gap:
        return (
            None,
            f"{pair[0]}/{pair[1]} not separated: O:22/O:23 gene evidence"
            f" {o22:.3f} vs {o23:.3f}, gap {gap:.3f} below {spec.min_gap:.3f}",
        )
    return (
        winner,
        f"{winner} over {loser}: O:22/O:23 sub-factor genes {o22:.3f} vs {o23:.3f}"
        f" ({spec.note or 'dedicated O:22/O:23 markers'})",
    )


def _special_notes(
    phenotype: str, calls: Mapping[str, GeneCall], document: SerotypeDocumentV2
) -> tuple[str, ...]:
    """Post-decision serovar refinements (oafA, gntR/STM3356, sdf...)."""
    out: list[str] = []
    for check in document.special_checks:
        if check.serovar not in ("*", phenotype):
            continue
        mine = _best_gene_evidence(check.markers, calls, document)
        if check.when == "absent":
            # A deletion: the note fires only when the marker is entirely below
            # the floors (a truncated assembly of a present gene still scores).
            if mine <= 0.0:
                out.append(check.note)
            continue
        if mine <= 0.0:
            continue
        theirs = _best_gene_evidence(check.versus, calls, document) if check.versus else 0.0
        if mine > theirs:
            out.append(check.note)
    return tuple(out)


def _class_exact_serovars(
    document: SerotypeDocumentV2,
    o_tokens: tuple[str, ...],
    h1: _FamilyCall,
    h2: _FamilyCall,
    *,
    subspecies: str = "enterica",
) -> tuple[str, ...]:
    """Rows for the observed formula read through each antigen's own class.

    A bracket spelling is the same antigen (`g,[s],t` == `g,s,t`), so the exact
    row lookup must try every class mate of the H-phase winner before falling
    back to the wide pass. That is what separates two cases which both look like
    a margin tie:

    * Gateshead — the winner `g,[s],t` has the exact row `g,s,t`, so the answer
      is read off the table and the wide pass (which also reaches `g,t` and `i`,
      both in SISTR's large `g`-family class) never runs;
    * Newport — the winner `e,n,x` has no row at all, so the wide pass must run
      to reach the `e,h` row the upstream tool reports.

    The O tokens are the caller's whole tie set (a genuine E1/E4 margin tie must
    keep both groups reachable); a discriminator-decided O group arrives already
    narrowed to one name.
    """
    if not o_tokens:
        return ()

    def names(family: Literal["o", "h1", "h2"], call: _FamilyCall) -> tuple[str, ...]:
        if call.name is None:
            return (FALLBACK_PHASE,)
        members = [call.name]
        for group in document.similarity_groups.get(family, ()):
            if call.name in group:
                members.extend(member for member in group if member not in members)
                break
        return tuple(members)

    h1_tokens = names("h1", h1)
    h2_tokens = names("h2", h2)
    found: list[str] = []
    for row in document.formulas:
        if row.subspecies != subspecies:
            continue
        if (
            row.o in o_tokens
            and row.h1 in h1_tokens
            and row.h2 in h2_tokens
            and row.serovar not in found
        ):
            found.append(row.serovar)
    return tuple(found)


def _evaluate_v2(calls: Mapping[str, GeneCall], document: SerotypeDocumentV2) -> SchemeCall:
    """Ranked antigen families, then the formula table."""
    control = _control_present(calls, document)
    calls, demoted = _separate_phases(calls, document)
    empty = _FamilyCall(name=None, score=0.0, genes=(), perfect=False, ambiguous=False)
    resolved: dict[str, _FamilyCall] = {
        family.key: (
            _family_call_o(family, calls, document)
            if family.key == "o"
            else _family_call(family, calls, document)
        )
        for family in document.antigens
    }
    for family_key, extra in demoted.items():
        call = resolved.get(family_key)
        if call is None or not extra or call.score >= DEMOTED_SCORE_FLOOR:
            continue
        if call.name is None:
            # The phase floor already said this family is absent (a monophasic
            # genome): resurrecting the demoted name would invent a second phase
            # and pull in every serovar sharing the first one.
            continue
        # The phase assignment is solid, but the family's own locus is weak: this
        # is where erasing the demoted name loses the serovar (Braenderup's `e,h`
        # references fit its fljB locus better than its own fliC, at 0.86). Keep the
        # name as a candidate only while the remaining evidence is poor.
        merged = tuple(dict.fromkeys((*call.ties, *extra)))
        resolved[family_key] = replace(call, ties=merged, ambiguous=len(merged) > 1)
    o = resolved.get("o", empty)
    o_before = o
    discriminator = _discriminate_o(o, calls, document)
    if discriminator is not None:
        # The marker's verdict replaces the family call everywhere downstream
        # (notes, lookup), not just the local binding — otherwise the per-family
        # note could still name the antigen the marker overruled.
        o = discriminator.family
        resolved["o"] = o
    # A discriminator that actually narrowed the tie makes the O group exact: the
    # lookup must not widen it back to its similarity class (A/D1/D2 share one).
    o_exact = o.ties != o_before.ties
    h1 = resolved.get("h1", empty)
    h2 = resolved.get("h2", empty)
    # An identical antigen in both phases is a cross-reaction, not a serovar.
    h1, h2, phase_note, phase_tie = _drop_duplicate_phase(h1, h2)
    if phase_note is not None:
        resolved["h1"], resolved["h2"] = h1, h2
    components = tuple(
        ScoreComponent(name=f"antigen:{gene}", score=call.score)
        for _, call in sorted(resolved.items())
        for gene in call.genes
    )
    # --- The confidence ladder (calibrated on the 100-genome corpus) ---------
    # `perfect` above asks for a 100 %-identity/100 %-coverage hit on every
    # antigen, which real assemblies almost never produce (measured: 8/100
    # genomes), so a ladder keyed on it labelled everything `low` and told a lab
    # nothing. The three signals below are the ones that actually separate a
    # trustworthy call from a lucky one, and they were measured rather than
    # guessed: on the corpus the resulting `high` tier is 23/100 genomes and
    # **100 % of them agree with both reference tools**.
    present_scores = [call.score for call in (o, h1, h2) if call.name is not None]
    weakest = min(present_scores) if present_scores else 0.0
    any_ambiguous = o.ambiguous or h1.ambiguous or h2.ambiguous or phase_tie
    # Whole, unambiguous evidence on every antigen is what `high` means; a
    # truncated locus or two names at the same score is named but not certain.
    all_full = all(call.full for call in (o, h1, h2) if call.name is not None)
    all_single = all(call.single_name for call in (o, h1, h2) if call.name is not None)

    def _confidence(*, widened: bool, shared: bool, evidence_only: bool = False) -> Confidence:
        """How much the evidence supports the name that was produced.

        ``high`` = a single-row formula read off the table, every family decided
        on whole-length evidence with one name, no class widening, and the
        weakest antigen at or above ``HIGH_SCORE_FLOOR`` (≈99 % identity at full
        coverage). ``medium`` covers every case where the name is real but
        something had to be resolved for it: a family tie broken by a
        discriminator, a one-to-many formula row (the name is declaration order),
        a widened token, truncated coverage, or softer identity. ``low`` is
        reserved for evidence that barely names anything.
        """
        if weakest < MEDIUM_SCORE_FLOOR:
            return "low"
        if evidence_only:
            return "medium"
        if (
            not any_ambiguous
            and all_full
            and all_single
            and not widened
            and not shared
            and weakest >= HIGH_SCORE_FLOOR
        ):
            return "high"
        return "medium"

    if not control:
        return SchemeCall(
            phenotype=document.fallback,
            score=0.0,
            confidence="low",
            notes=("control gene absent",),
        )
    notes = [
        f"formula {o.name or '-'}:{h1.name or '-'}:{h2.name or '-'}",
        *(
            f"{family_key.upper()}: {call.name or '-'}"
            + (
                f" via {call.genes[0]} ({call.identity_pct:.2f}% id, {call.coverage_pct:.2f}% cov)"
                if call.name is not None and call.genes
                else ""
            )
            + ("" if call.perfect else " (not a perfect match)")
            for family_key, call in sorted(resolved.items())
        ),
        # A near-tie the ranker did not collapse: the runner-up names stay
        # visible next to the definitive call, so a lab can audit the margin
        # (measured: `g,m` at 0.9993 against `g,q` at 0.9987 on an Enteritidis
        # genome, one SNP of difference inside the 0.01 margin). Names the
        # evidence could not separate at all report through the candidate list.
        *(
            f"{family_key.upper()} near-tie within {document.margin}: {', '.join(call.ties[1:])}"
            for family_key, call in sorted(resolved.items())
            if call.name is not None and len(call.ties) > 1 and len(call.equal_ties) == 1
        ),
    ]
    if discriminator is not None:
        notes.extend(discriminator.notes)
    if phase_note is not None:
        notes.append(phase_note)
    # The numeric O for a serogroup the WHO table pins to one formula: the bridge
    # to SeqSero2's notation (`P` is what it prints as `38`). Ambiguous serogroups
    # get no note — there is no single answer to print.
    numeric_group = _numeric_o(document, o.name)
    if numeric_group is not None:
        notes.append(f"O group {o.name} = O:{numeric_group}")
    observed = (o.name, h1.name or FALLBACK_PHASE, h2.name or FALLBACK_PHASE)
    score = max(o.score, h1.score, h2.score)
    subspecies = _detect_subspecies(calls)
    if subspecies is not None and subspecies != "enterica":
        # Non-enterica subspecies are reported by their own formula rows (the
        # Kauffmann-White ``NUMERAL O:H1:H2`` names), looked up by the H phases
        # ONLY: the O-antigen references are enterica-specific and must not
        # resolve a non-enterica O group (a VI O:6,7 wzx is 100 % identical to
        # enterica O:11, which would flip the group to F).
        numeral = SUBSPECIES_NUMERAL[subspecies]
        candidates = _lookup_subspecies_serovar(document, subspecies, h1.name, h2.name)
        if len(candidates) != 1:
            # Several rows share the H phases: disambiguate by the resolved O
            # serogroup, which maps consistently for most non-enterica rows.
            candidates = _lookup_subspecies_serovar(document, subspecies, h1.name, h2.name, o.name)
        if candidates:
            formula = candidates[0]
        else:
            numeric = document.o_numeric.get(subspecies) or document.o_numeric.get("enterica", {})
            numeric_o = numeric.get(o.name, o.name) if o.name else "-"
            formula = f"{numeral} {numeric_o}:{h1.name or '-'}:{h2.name or '-'}"
        note = (
            f"subspecies {subspecies} ({numeral}) outside the enterica panel; reported as {formula}"
        )
        return SchemeCall(
            phenotype=formula,
            score=round(score, 4),
            confidence=_confidence(widened=False, shared=len(candidates) > 1, evidence_only=True),
            components=components,
            notes=(*notes, note),
        )
    o_names = tuple(
        entry.name for family in document.antigens if family.key == "o" for entry in family.entries
    )
    decided = len(o.equal_ties) <= 1 and len(h1.equal_ties) <= 1 and len(h2.equal_ties) <= 1
    exact_serovars = _exact_serovars(document, observed) if decided else []
    class_exact = False
    if not exact_serovars and decided:
        # No row spells the formula the way the panel does: read it through each
        # antigen's own class first (a bracket spelling is the same antigen), and
        # only widen over every tied name when even that finds nothing.
        o_tokens: tuple[str, ...] = (o.name,) if o_exact and o.name else o.ties
        exact_serovars = list(_class_exact_serovars(document, o_tokens, h1, h2))
        class_exact = bool(exact_serovars)
    if exact_serovars:
        # The table carries the observed formula: that row IS the call. The
        # widening pass exists for formulas it does not carry, and running it
        # here only imports sibling serovars the antigens never named.
        if len(exact_serovars) > 1 and class_exact:
            # Two DIFFERENT rows matched through the class: genuine ambiguity, not
            # one formula with two names, so no primary may be picked.
            notes.append("formula matched several serovars after widening")
            return SchemeCall(
                phenotype=None,
                score=round(score, 4),
                confidence="ambiguous",
                components=components,
                ambiguous=tuple(
                    PhenotypeScore(phenotype=serovar, score=round(score, 4))
                    for serovar in exact_serovars
                ),
                notes=tuple(notes),
            )
        if len(exact_serovars) > 1:
            broken = _break_o_subfactor(exact_serovars, calls, document)
            if broken is not None:
                winner, discriminator_note = broken
                if winner is not None:
                    notes.append(discriminator_note)
                    notes.extend(_special_notes(winner, calls, document))
                    return SchemeCall(
                        phenotype=winner,
                        score=round(score, 4),
                        # A serovar-level tie broken by dedicated markers: real
                        # evidence, but a resolved ambiguity rather than a
                        # formula read off the table.
                        confidence=_confidence(widened=class_exact, shared=True),
                        components=components,
                        notes=tuple(notes),
                    )
                notes.append(discriminator_note)
            notes.append(f"formula shared with {', '.join(exact_serovars[1:])}")
        if class_exact:
            notes.append("matched after widening the antigen names to their class")
        notes.extend(_special_notes(exact_serovars[0], calls, document))
        return SchemeCall(
            phenotype=exact_serovars[0],
            score=round(score, 4),
            confidence=_confidence(widened=class_exact, shared=len(exact_serovars) > 1),
            components=components,
            runner_up=(
                PhenotypeScore(phenotype=exact_serovars[1], score=round(score, 4))
                if len(exact_serovars) > 1
                else None
            ),
            notes=tuple(notes),
        )
    if not decided:
        # Some family has two antigens at identical evidence (`g,m` against
        # `g,q` at 99.93 %/100 % on a real genome — every channel tied): no
        # mechanism chose, so the rows those names key are the candidates —
        # a tight set, never the margin siblings' rows or the class
        # cross-product the widening pass would import.
        tied_serovars = _formula_rows(
            document,
            o.equal_ties if o.name is not None else o_names,
            h1.equal_ties if h1.name is not None else (FALLBACK_PHASE,),
            h2.equal_ties if h2.name is not None else (FALLBACK_PHASE,),
        )
        if tied_serovars:
            notes.append(
                "antigens tied at identical evidence: "
                + "; ".join(
                    f"{family_key} {{{', '.join(call.equal_ties)}}}"
                    for family_key, call in sorted(resolved.items())
                    if len(call.equal_ties) > 1
                )
            )
            if len(tied_serovars) == 1:
                # Both tied names key the same row: the tie cannot change the
                # serovar, so the call stands (with the tie on record).
                return SchemeCall(
                    phenotype=tied_serovars[0],
                    score=round(score, 4),
                    confidence=_confidence(widened=False, shared=False),
                    components=components,
                    notes=tuple(notes),
                )
            return SchemeCall(
                phenotype=None,
                score=round(score, 4),
                confidence="ambiguous",
                components=components,
                ambiguous=tuple(
                    PhenotypeScore(phenotype=serovar, score=round(score, 4))
                    for serovar in tied_serovars
                ),
                notes=tuple(notes),
            )
    serovars = _lookup_formula(
        document,
        o.ties,
        h1.ties,
        h2.ties,
        o_names=o_names,
        alias_h2_absent=not h2.perfect,
        o_exact=o_exact,
    )
    if serovars:
        # The observed formula is not in the table: the name came from widening a
        # sibling's class (a cross-labelled reference) or from another reference
        # covering the phase, which is weaker evidence than reading the formula.
        notes.append("matched after widening the antigen names to their class")
    if len(serovars) == 1:
        # The widened pass found exactly one row: a name, but evidence weaker
        # than reading the formula off the table (hence never `high`).
        notes.extend(_special_notes(serovars[0], calls, document))
        return SchemeCall(
            phenotype=serovars[0],
            score=round(score, 4),
            confidence=_confidence(widened=True, shared=False),
            components=components,
            notes=tuple(notes),
        )
    if len(serovars) > 1:
        broken = _break_o_subfactor(serovars, calls, document)
        if broken is not None:
            winner, discriminator_note = broken
            if winner is not None:
                notes.append(discriminator_note)
                notes.extend(_special_notes(winner, calls, document))
                return SchemeCall(
                    phenotype=winner,
                    score=round(score, 4),
                    confidence=_confidence(widened=True, shared=True),
                    components=components,
                    notes=tuple(notes),
                )
            notes.append(discriminator_note)
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
    numeric_o = _numeric_o(document, o.name) or o.name
    formula = f"{numeric_o}:{h1.name or '-'}:{h2.name or '-'}"
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
