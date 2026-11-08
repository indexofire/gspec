"""Screening parameters and the canonical Report model."""

from typing import Literal

from pydantic import BaseModel, Field

from gspec.hits import Hit


class ScreeningParams(BaseModel, frozen=True):
    """Parameters for one screening run (SPEC.md §3 defaults and bounds).

    ``coverage``/``dedup``/``min_aligned_len`` are the gspec extension that
    makes the blastn path's *evidence model* configurable (``gspec.screen/2``);
    their defaults reproduce the abricate-parity behaviour byte for byte.

    - ``coverage="best_hsp"`` — one row per dedup survivor, coverage is that
      single HSP's span (parity). ``coverage="union"`` — all of a subject's
      HSPs are aggregated into one hit: coverage is the *union* of the subject
      intervals (capped at 100), identity is aligned-length weighted, and the
      rendered span is the union's bounding box.
    - ``dedup`` — the dedup key: ``query`` (parity, one survival per query
      interval) or ``subject`` (per-subject HSPs survive).
    - ``min_aligned_len`` — drop aggregated hits whose total ungapped aligned
      length is below this many bases (a spurious-HSP guard for union mode).
    """

    db: str
    minid: float = Field(default=80.0, gt=0.0, le=100.0)
    mincov: float = Field(default=80.0, ge=0.0, le=100.0)
    threads: int = Field(default=1, ge=1)
    # Additive opt-in flag; serialized only when true so default reports stay
    # byte-identical to the pre-merge-fragments contract.
    merge_fragments: bool = False
    # Evidence model (gspec.screen/2). Defaults = SPEC.md §4 parity.
    coverage: Literal["best_hsp", "union"] = "best_hsp"
    dedup: Literal["query", "subject"] = "query"
    min_aligned_len: int | None = Field(default=None, ge=1)


class Report(BaseModel, frozen=True):
    """Canonical in-memory result of screening one input file.

    Hits are sorted (sequence, start) before construction; the tuple is the
    final, stable order. Screening is pure gene detection: typed and
    untyped gene databases produce byte-identical reports (designation is
    the ``gspec typing`` command's job).
    """

    file: str
    hits: tuple[Hit, ...]
