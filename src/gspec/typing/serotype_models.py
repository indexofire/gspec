"""Serotype profile documents: the marker-set profile and the formula table.

Two shapes, one contract file (``serotype.json``):

- ``gspec.serotype/1`` — an ordered marker-set profile with per-marker floors
  and ambiguity policy (``shi_sero``, ``sen_sero``).
- ``gspec.serotype/2`` — ranked antigen families plus a formula table, the
  shape a full Salmonella panel needs (tens of O groups, hundreds of H alleles,
  thousands of O:H1:H2 rows). See the block below.

Deliberately separate from ``gspec.typing``: a serotype call combines several
antigen markers and needs its own fallback and ambiguity policy, while the
rule-scheme engine answers "which rule wins" over one gene table.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from gspec.errors import DatabaseError
from gspec.jsonio import read_json

SEROTYPE_FILENAME = "serotype.json"


class SerotypeRule(BaseModel, frozen=True):
    """One ordered marker-set rule; every matched rule is surfaced.

    ``max_coverage`` is a per-marker ceiling (strictly below the value); it
    expresses the S. boydii 6 IS629-interrupted ``wbaM`` case without turning
    the whole rule into a low-coverage call.  ``min_alignment_start_pct`` and
    ``max_alignment_start_pct`` gate a marker's alignment start (0 = the 5'
    end is covered), separating a 5'-interrupted allele from a 3'-truncated
    assembly.  A rule that declares a start gate requires alignment evidence;
    reads tables carry none, so such rules do not match there.
    """

    phenotype: str
    requires: tuple[str, ...] = ()
    requires_any: tuple[str, ...] = ()
    excludes: tuple[str, ...] = ()
    max_coverage: dict[str, float] = Field(
        default_factory=dict,
        description="Per-marker coverage ceiling; the rule fails at or above this percentage.",
    )
    outranks: dict[str, tuple[str, ...]] = Field(
        default_factory=dict,
        description=(
            "Ranked allele gate: each named marker must be present AND score"
            " strictly higher (identity x coverage) than every listed"
            " competitor. This is how near-identical alleles are separated"
            " (O:2 vs O:9: the two tyr copies are 98.5-99.2% identical, so a"
            " presence floor cannot decide them)."
        ),
    )
    min_alignment_start_pct: dict[str, float] = Field(
        default_factory=dict,
        description=(
            "Per-marker minimum alignment start as percent of marker length"
            " (0 = the 5' end is covered)."
        ),
    )
    max_alignment_start_pct: dict[str, float] = Field(
        default_factory=dict,
        description=(
            "Per-marker maximum alignment start as percent of marker length"
            " (0 = the 5' end is covered)."
        ),
    )
    coverage_floor: float = Field(default=50.0, ge=0.0, le=100.0)
    identity_floor: float = Field(default=80.0, ge=0.0, le=100.0)
    control_exempt: bool = False
    notes: tuple[str, ...] = ()

    @model_validator(mode="after")
    def _has_positive_marker(self) -> SerotypeRule:
        if not (self.requires or self.requires_any):
            raise ValueError("serotype rule needs requires or requires_any")
        if set(self.requires) & set(self.requires_any):
            raise ValueError("requires and requires_any overlap")
        if set(self.requires) & set(self.excludes):
            raise ValueError("requires and excludes overlap")
        if set(self.requires_any) & set(self.excludes):
            raise ValueError("requires_any and excludes overlap")
        positives = {*self.requires, *self.requires_any}
        for name, gates in (
            ("max_coverage", self.max_coverage),
            ("min_alignment_start_pct", self.min_alignment_start_pct),
            ("max_alignment_start_pct", self.max_alignment_start_pct),
        ):
            for gene, value in gates.items():
                if gene not in positives:
                    raise ValueError(f"{name} gene {gene!r} is not a positive marker")
                if not 0.0 <= value <= 100.0:
                    raise ValueError(f"{name}[{gene!r}] must be a percentage, got {value}")
        for gene, competitors in self.outranks.items():
            if gene not in positives:
                raise ValueError(f"outranks gene {gene!r} is not a positive marker")
            if gene in competitors:
                raise ValueError(f"outranks[{gene!r}] names itself as a competitor")
            if len(set(competitors)) != len(competitors):
                raise ValueError(f"outranks[{gene!r}] lists a competitor twice")
        return self


class SerotypeDocument(BaseModel, frozen=True):
    """gspec.serotype/1 — one database's ordered serotype rules."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gspec.serotype/1"] = Field(default="gspec.serotype/1", alias="schema")
    name: str
    control_gene: str | None = None
    control_coverage_floor: float = Field(default=50.0, ge=0.0, le=100.0)
    control_identity_floor: float = Field(default=80.0, ge=0.0, le=100.0)
    fallback: str = "Not Shigella/EIEC"
    ambiguous: str = "No prediction (multiple marker sets)"
    unresolved: str = "Shigella spp. (unresolved)"
    #: How several matching rules are resolved. ``ambiguous`` (the default)
    #: reports the ambiguity — a marker set that matches two phenotypes is never
    #: decided by declaration order. ``most_specific`` drops any match whose
    #: positive markers are a strict subset of another match, so a full
    #: antigenic-formula rule beats the generic "O group only, H undetermined"
    #: rule beneath it; incomparable matches stay ambiguous.
    decision: Literal["ambiguous", "most_specific"] = "ambiguous"
    rules: tuple[SerotypeRule, ...]


# ---------------------------------------------------------------------------
# gspec.serotype/2 — antigen families + the antigenic-formula table.
#
# The /1 marker-set profile scales to panels of tens of rules. A full
# Salmonella panel is a different shape: ~46 O groups, ~300 H alleles and a
# reference table of >2,000 O:H1:H2 formulas, which cannot be written as one
# rule per serovar. /2 models that shape directly: each antigen family is a
# RANKED list of allele references (the allele that scores highest wins — the
# SeqSero2/SISTR "top hit" step), and the resulting formula is looked up in a
# table. The control-gene gate and the floors are the same as /1.
# ---------------------------------------------------------------------------


class AntigenEntry(BaseModel, frozen=True):
    """One named antigen variant and the reference sequences that prove it.

    ``genes`` holds one or more marker ids; the entry's score is the best
    (identity x coverage) among them that clears the document floors, so a
    variant covered by several reference alleles needs one entry, not several.
    """

    name: str = Field(min_length=1)
    genes: tuple[str, ...] = Field(min_length=1)

    @model_validator(mode="after")
    def _no_duplicate_genes(self) -> AntigenEntry:
        if len(set(self.genes)) != len(self.genes):
            raise ValueError(f"antigen {self.name!r} lists a gene twice")
        return self


class AntigenFamily(BaseModel, frozen=True):
    """One ranked antigen family (``o``, ``h1``, ``h2``)."""

    key: Literal["o", "h1", "h2"]
    entries: tuple[AntigenEntry, ...] = Field(min_length=1)

    @model_validator(mode="after")
    def _unique_names(self) -> AntigenFamily:
        names = [entry.name for entry in self.entries]
        if len(set(names)) != len(names):
            raise ValueError(f"family {self.key!r} declares an antigen name twice")
        return self


class SerotypeFormula(BaseModel, frozen=True):
    """One row of the antigenic-formula table: O:H1:H2 -> serovar.

    ``"-"`` marks a phase the formula does not carry (monophasic serovars) and
    is also what a genome without that phase resolves to, so the lookup treats
    the two as the same token.
    """

    o: str
    h1: str
    h2: str
    serovar: str = Field(min_length=1)


class ODiscriminator(BaseModel, frozen=True):
    """Optional SNP-level tie-breaker for the O family (``gspec.serotype/2``).

    Some serogroups cannot be separated by their antigen references at all: A
    (O:2) and D1 (O:9) share the rfb cluster and differ only in `rfbS`, whose
    alleles are 99.3 % identical over the gene (a single SNP in the 5' 609 bp
    SeqSero2 uses). Such a marker only speaks when the gap is decisive, so the
    document carries the floor and the engine stays with the ambiguity below it.
    """

    #: O-family antigen name -> the gene id that supports it.
    favours: dict[str, str]
    #: Minimum identity gap, in percentage points, between the two markers.
    min_gap_pct: float = Field(default=0.5, ge=0.0, le=100.0)
    #: Both markers must be present at least this well for the gap to mean anything.
    min_coverage: float = Field(default=60.0, ge=0.0, le=100.0)
    note: str = ""


class SerotypeDocumentV2(BaseModel, frozen=True):
    """gspec.serotype/2 — ranked antigen families plus the formula table."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gspec.serotype/2"] = Field(default="gspec.serotype/2", alias="schema")
    name: str
    control_gene: str | None = None
    control_coverage_floor: float = Field(default=50.0, ge=0.0, le=100.0)
    control_identity_floor: float = Field(default=80.0, ge=0.0, le=100.0)
    #: A reference allele counts only when it clears both floors.
    identity_floor: float = Field(default=80.0, ge=0.0, le=100.0)
    coverage_floor: float = Field(default=60.0, ge=0.0, le=100.0)
    #: A SECOND, weaker acceptance rule for a locus the assembly failed to close.
    #: Antigen genes sit in repeat-rich regions (the rfb cluster especially), so
    #: short-read assemblies break them across contigs: a real D1 genome showed
    #: `wzx` at 99.85 % identity over only 52 % union coverage, which the
    #: coverage floor then discarded as absent. A near-identical partial hit is
    #: evidence of fragmentation, not of absence, so it is accepted here — it
    #: still scores below a full-length hit (the ranker multiplies identity by
    #: coverage) and can never be perfect. Set ``partial_coverage_floor`` to 0 to
    #: disable the second rule.
    partial_coverage_floor: float = Field(default=45.0, ge=0.0, le=100.0)
    partial_identity_floor: float = Field(default=95.0, ge=0.0, le=100.0)
    #: Two antigen candidates closer than this (on identity x coverage) are
    #: ambiguous. It must sit below the measured gap between near-identical
    #: alleles (O:2 vs O:9 differ by 0.8-1.5 %): the calibration knob.
    margin: float = Field(default=0.01, ge=0.0, le=1.0)
    antigens: tuple[AntigenFamily, ...] = Field(min_length=1)
    formulas: tuple[SerotypeFormula, ...] = Field(min_length=1)
    #: Equivalent token classes per family, e.g. ``{"h1": [["z13,z28", "l,z13,[z28]"]]}``.
    #: The lookup is exact first; only when no row matches is a token widened to
    #: its class, which is how the upstream tools treat WKL bracket variants
    #: (SISTR's ``*_SIMILARITY_GROUPS``). Widening can match several serovars,
    #: which is reported as an ambiguous call, never a silent pick.
    similarity_groups: dict[Literal["o", "h1", "h2"], tuple[tuple[str, ...], ...]] = Field(
        default_factory=dict
    )
    #: SNP-level tie-breaker between serogroups whose references cannot separate
    #: them (see :class:`ODiscriminator`).
    o_discriminator: ODiscriminator | None = None
    fallback: str = "Not Salmonella"
    unresolved: str = "Salmonella subsp. enterica (no antigen set resolved)"
    partial: str = "Salmonella {formula} (serovar not in the formula table)"

    @model_validator(mode="after")
    def _families_are_unique(self) -> SerotypeDocumentV2:
        keys = [family.key for family in self.antigens]
        if len(set(keys)) != len(keys):
            raise ValueError("an antigen family is declared twice")
        return self


class _SchemaTag(BaseModel, frozen=True):
    """The ``schema`` tag alone, for loader dispatch (extra keys ignored)."""

    model_config = ConfigDict(populate_by_name=True, extra="ignore")

    schema_name: str | None = Field(default=None, alias="schema")


AnySerotypeDocument = SerotypeDocument | SerotypeDocumentV2


def read_serotype_document(path: Path) -> AnySerotypeDocument:
    """Read a serotype profile, dispatching on its schema tag.

    ``gspec.serotype/1`` is the marker-set profile; ``/2`` is the ranked-antigen
    + formula-table shape. An unknown tag is a typed database error, never a
    silent fallback to the older shape.
    """
    tag = read_json(_SchemaTag, path, code="SEROTYPE_MALFORMED", noun="serotype profile")
    if tag.schema_name == "gspec.serotype/2":
        return read_json(
            SerotypeDocumentV2, path, code="SEROTYPE_MALFORMED", noun="serotype profile"
        )
    if tag.schema_name in (None, "gspec.serotype/1"):
        return read_json(SerotypeDocument, path, code="SEROTYPE_MALFORMED", noun="serotype profile")
    raise DatabaseError(
        f"unsupported serotype schema tag: {tag.schema_name!r}",
        code="SEROTYPE_MALFORMED",
        context={"file": str(path), "schema": str(tag.schema_name)},
    )


def load_serotype(db_dir: Path) -> AnySerotypeDocument | None:
    """The optional serotype profile next to a gene database."""
    path = db_dir / SEROTYPE_FILENAME
    return read_serotype_document(path) if path.is_file() else None


def validate_serotype_references(document: AnySerotypeDocument, gene_ids: frozenset[str]) -> None:
    """Every marker a rule or antigen family references must exist."""
    if isinstance(document, SerotypeDocumentV2):
        for family in document.antigens:
            for entry in family.entries:
                for gene in entry.genes:
                    if gene not in gene_ids:
                        raise DatabaseError(
                            f"antigen {family.key}:{entry.name!r} references unknown marker"
                            f" {gene!r}",
                            code="SEROTYPE_UNKNOWN_GENE",
                            context={"family": family.key, "antigen": entry.name, "gene": gene},
                        )
        return
    for rule in document.rules:
        for gene in (
            *rule.requires,
            *rule.requires_any,
            *rule.excludes,
            *rule.max_coverage,
            *rule.min_alignment_start_pct,
            *rule.max_alignment_start_pct,
            *rule.outranks,
            *(c for v in rule.outranks.values() for c in v),
        ):
            if gene not in gene_ids:
                raise DatabaseError(
                    f"serotype rule {rule.phenotype!r} references unknown marker {gene!r}",
                    code="SEROTYPE_UNKNOWN_GENE",
                    context={"phenotype": rule.phenotype, "gene": gene},
                )


__all__ = [
    "SEROTYPE_FILENAME",
    "ODiscriminator",
    "SerotypeDocument",
    "SerotypeRule",
    "load_serotype",
    "read_serotype_document",
    "validate_serotype_references",
]
