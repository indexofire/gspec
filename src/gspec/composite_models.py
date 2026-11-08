"""``gspec.composite/1`` — declarative multi-database (multi-tool) fusion.

One typing document describes ONE database. Real serotyping workflows need
several evidence layers whose tools, thresholds and databases differ: a
marker panel (blastn, strict), a full-locus alignment (blastn + union
coverage, or minimap2), an optional PCR/SNP layer. Today those layers are
stitched together by scripts, which puts the *fusion policy* — not just its
parameters — outside the configuration.

A pipeline document makes the fusion declarative: an ordered list of
``stages``, each naming a database (and optionally one of its schemes) plus
the acceptance rule for its call. Stages evaluate in order and the first
accepted call wins; every stage's raw call is kept as evidence, so the
result carries the audit trail a two-layer call needs.

Schema::

    {
      "schema": "gspec.composite/1",
      "name": "cje_penner",
      "fallback": "untypeable",
      "stages": [
        {"name": "marker",  "db": "cje_hs_marker", "scheme": "hs_marker",
         "tier": "B", "min_score": 0.90},
        {"name": "locus_a", "db": "cje_cps_locus", "scheme": "tier_a",
         "tier": "A", "min_score": 0.85},
        {"name": "locus_c", "db": "cje_cps_locus", "scheme": "tier_c",
         "tier": "C", "min_score": 0.30}
      ]
    }

Stage fields: ``name`` (evidence label, unique), ``db``, ``scheme`` (None =
the database's only scheme), ``tier`` (label flowing to the result),
``min_score`` (0..1 acceptance floor), ``require_confidence`` (``high`` =
accept only high-confidence calls, ``any`` = accept ambiguous/fallback-free
calls as well).
"""

from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator

from gspec.errors import DatabaseError, InputError

COMPOSITE_SCHEMA = "gspec.composite/1"
#: The document's filename inside a composite database directory.
DOCUMENT_FILENAME = "composite.json"


class KmerSetRule(BaseModel, frozen=True):
    """One k-mer set → phenotype mapping for a ``kind: "kmer"`` stage.

    ``set`` names a set in the database's ``gspec.kmerdb/1`` sidecar (the
    differential sets cluster databases use for SNP-level pairs); the stage
    counts the query's canonical k-mers once and scores every set by the
    fraction of its k-mers the query carries. ``min_fraction`` is the floor a
    set must clear to designate.
    """

    set: str = Field(min_length=1)
    phenotype: str = Field(min_length=1)
    min_fraction: float = Field(default=0.5, ge=0.0, le=1.0)
    notes: tuple[str, ...] = ()


class PrimerProduct(BaseModel, frozen=True):
    """One primer pair's acceptance window inside a PCR rule."""

    name: str = Field(min_length=1)
    min_size: int = Field(default=20, ge=0)
    max_size: int = Field(default=20_000, ge=0)

    @model_validator(mode="after")
    def _window(self) -> "PrimerProduct":
        if self.min_size > self.max_size:
            raise ValueError(
                f"product {self.name!r}: min_size {self.min_size} exceeds max_size {self.max_size}"
            )
        return self


class PrimerRule(BaseModel, frozen=True):
    """One refinement rule: when ``products`` amplify inside their windows (and
    the current call matches ``when``, if given), the call becomes ``set``.

    Rules are evaluated in declaration order and the first applicable one wins,
    mirroring the if/elif cascade a published assay describes. ``when`` matches
    the call's base type, so ``HS4c (HS4/13/43/50/62/65 complex)`` matches
    ``"HS4c"``; ``set`` is the literal new phenotype, brackets included.
    """

    products: tuple[PrimerProduct, ...] = Field(min_length=1)
    when: tuple[str, ...] = ()
    set: str = Field(min_length=1)
    notes: tuple[str, ...] = ()


class CompositeStage(BaseModel, frozen=True):
    """One evidence stage: a database (+scheme) and its acceptance rule."""

    # "db" designates from a database (alignment), "kmer" from a database's
    # k-mer substrate (no alignment), "pcr" refines the cascade's call with an
    # in-silico PCR primer set (gseek.ispcr). One scheme can mix all three.
    kind: Literal["db", "kmer", "pcr"] = "db"
    name: str = Field(min_length=1)
    db: str = Field(min_length=1)
    scheme: str | None = None
    tier: str = "A"
    min_score: float = Field(default=0.0, ge=0.0, le=1.0)
    require_confidence: Literal["high", "any"] = "any"
    # Cluster databases: the locus decision gate (gspec.screen -C) and the
    # per-gene floors (-G), so one locus db can serve several tiers.
    min_cluster_cov: float | None = Field(default=None, ge=0.0, le=100.0)
    min_gene_cov: float | None = Field(default=None, ge=0.0, le=100.0)
    # PCR stages (kind="pcr"): the primer set and its ordered refinement rules.
    primers: str | None = None
    max_mismatches: int = Field(default=2, ge=0)
    #: Indel-tolerant pass budget (edits per primer site); 0 = substitutions only.
    max_edits: int = Field(default=0, ge=0)
    rules: tuple[PrimerRule, ...] = ()
    # kmer stages: which sets designate which phenotype.
    sets: tuple[KmerSetRule, ...] = ()
    # Scan the whole primer set for evidence (default: only the pairs the
    # rules name — the rest only cost a genome scan each).
    report_all_products: bool = False

    @model_validator(mode="after")
    def _kind_contract(self) -> "CompositeStage":
        if self.kind == "pcr":
            if not self.primers:
                raise ValueError(f"pcr stage {self.name!r} needs a primers path")
            if not self.rules:
                raise ValueError(f"pcr stage {self.name!r} needs at least one rule")
        elif self.rules:
            raise ValueError(
                f"stage {self.name!r} declares rules but kind is {self.kind!r} (rules are pcr-only)"
            )
        if self.kind == "kmer" and not self.sets:
            raise ValueError(f"kmer stage {self.name!r} needs at least one set rule")
        if self.kind != "kmer" and self.sets:
            raise ValueError(
                f"stage {self.name!r} declares sets but kind is {self.kind!r} (sets are kmer-only)"
            )
        return self


class StageCallSelector(BaseModel, frozen=True):
    """A (stage, phenotype) pair a conflict rule matches on.

    ``stage=None`` on an ``override`` means "any database stage" — the assay's
    exception is about the phenotype, not about which layer happens to have
    emitted it (the locus tier that wins varies with alignment strength).
    """

    stage: str | None = None
    phenotype: str = Field(min_length=1)


class ConflictRule(BaseModel, frozen=True):
    """A declared cross-stage contradiction: when the cascade picked ``when``
    and another stage independently called ``override``, the override wins.

    Published assays carry exactly this kind of exception ("the HS41 marker
    also amplifies in HS55 strains, so trust the locus there"). Declaring it
    keeps the exception in configuration instead of in code, and pins the two
    competing calls so an unrelated change cannot silently trigger it.
    """

    when: StageCallSelector
    override: StageCallSelector
    notes: tuple[str, ...] = ()


class CompositeDocument(BaseModel, frozen=True):
    """gspec.composite/1 — ordered evidence stages fused into one call."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gspec.composite/1"] = Field(default=COMPOSITE_SCHEMA, alias="schema")
    name: str = Field(min_length=1)
    stages: tuple[CompositeStage, ...] = Field(min_length=1)
    fallback: str = "untypeable"
    # Declared marker-vs-locus exceptions, applied after the cascade and
    # before any PCR refinement.
    conflicts: tuple[ConflictRule, ...] = ()
    #: Directory (relative to this document) holding the databases the stages
    #: name — the self-contained bundle form. Bundled copies win over installed
    #: databases of the same name, so a bundle runs on a machine that has none
    #: of them installed.
    db_root: str | None = None

    @model_validator(mode="after")
    def _db_root_is_relative(self) -> "CompositeDocument":
        if self.db_root is None:
            return self
        candidate = Path(self.db_root)
        if not self.db_root or candidate.is_absolute():
            raise ValueError("db_root must be a relative path inside the pipeline directory")
        if ".." in candidate.parts:
            raise ValueError(f"db_root must not traverse upwards: {self.db_root!r}")
        return self

    # The document's own directory, so relative primer paths resolve next to
    # the pipeline that names them (never serialized: it is where the file was
    # read from, not what the file says).
    source_dir: Path | None = Field(default=None, exclude=True)

    @property
    def schema_tag(self) -> str:
        return self.schema_name

    def validate_stages(self) -> None:
        """Stage names label evidence columns, so they must be unique."""
        names = [stage.name for stage in self.stages]
        duplicates = sorted({name for name in names if names.count(name) > 1})
        if duplicates:
            raise DatabaseError(
                f"pipeline {self.name!r} has duplicate stage names: {', '.join(duplicates)}",
                code="COMPOSITE_MALFORMED",
                context={"pipeline": self.name, "stages": duplicates},
            )


def read_composite(path: Path) -> CompositeDocument:
    """Parse one composite document; missing file = INPUT_NOT_FOUND,
    malformed content = COMPOSITE_MALFORMED (typing/floors precedent)."""
    if not path.is_file():
        raise InputError(
            f"pipeline file not found or unreadable: {path}",
            code="INPUT_NOT_FOUND",
            context={"file": str(path)},
        )
    try:
        document = CompositeDocument.model_validate_json(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError) as exc:
        raise DatabaseError(
            f"could not read pipeline {path}: {exc}",
            code="COMPOSITE_MALFORMED",
            context={"file": str(path)},
        ) from exc
    except ValidationError as exc:
        raise DatabaseError(
            f"malformed pipeline {path}: {' '.join(str(exc).split())}",
            code="COMPOSITE_MALFORMED",
            context={"file": str(path)},
        ) from exc
    document.validate_stages()
    return document.model_copy(update={"source_dir": path.resolve().parent})


__all__ = [
    "COMPOSITE_SCHEMA",
    "DOCUMENT_FILENAME",
    "CompositeDocument",
    "CompositeStage",
    "read_composite",
]
