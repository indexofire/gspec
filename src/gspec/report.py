"""The screening result model: parameters, hits, reports, SPEC.md §4 processing.

One module owns what a gene call *is* and how a blastn row becomes one:
:class:`ScreeningParams` and :class:`Report` are the document models, while
:class:`Hit`/:class:`Fragment` plus the SPEC.md §4 steps (decode, dedupe,
COVERAGE_MAP) turn aligner output into them. ``--merge-fragments`` folding
lives in :mod:`gspec.align.fragments`.

No interval merging of any kind — upstream reports overlapping genes at
different query spans, and so do we. Subject ids decode through
:mod:`gspec.dbcodec`: legacy ``~~~`` ids parse by the frozen db.py rules,
``gspec|``-tagged ids by the strict native codec — a malformed native header
raises DatabaseError ``HEADER_MALFORMED`` (exit 4), never a silent fallback.
"""

import re
from collections.abc import Iterable, Iterator
from typing import TYPE_CHECKING, Literal

from pydantic import BaseModel, Field

from gspec.db.core import IDSEP
from gspec.dbcodec import decode_seqid, is_gspec_header

# Perl: $product =~ s/^\S+\s+// if $product =~ m/~~~/  — strips the leading
# makeblastdb id token only when followed by whitespace; a bare ~~~id stays.
# Native gspec| ids get the same substitution: makeblastdb prefixes stitle
# with the id, and there the prefix is the whitespace-free tagged seqid.
_LEADING_TOKEN_RE = re.compile(r"^\S+\s+")

if TYPE_CHECKING:
    # TYPE_CHECKING-only: blast -> report would otherwise be a cycle.
    from gspec.align.blast import BlastRow


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
    - ``allele_sketch`` — annotate each hit with its anchored allele-sketch score
      (``align/sketch.py``), the alignment-free evidence the serotype ranker uses
      to break near-ties. Off by default: it costs one k-mer count per hit.
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
    allele_sketch: bool = False


class Fragment(BaseModel, frozen=True):
    """One contributing fragment of a merged hit (``--merge-fragments``)."""

    contig: str
    start: int
    end: int
    strand: Literal["+", "-"]
    identity_pct: float
    coverage_pct: float


class Hit(BaseModel, frozen=True):
    """One surviving hit — raw semantic values only (display formatting is Phase 3)."""

    sequence: str
    start: int
    end: int
    strand: Literal["+", "-"]
    gene: str
    database: str
    accession: str
    product: str
    function: str
    s_start: int
    s_end: int
    s_len: int
    coverage_map: str
    gap_openings: int
    gaps: int
    identity_pct: float
    coverage_pct: float
    # Ungapped aligned columns (BLAST length - gaps); weight for fragment
    # identity means and anchor selection (internal, never rendered).
    aligned_len: int
    # --merge-fragments extension: merged hits carry their contributing
    # fragments; default-path hits never do (parity output is unchanged).
    merged: bool = False
    fragments: tuple[Fragment, ...] = ()
    # anchored allele sketch (align/sketch.py): the best containment of this
    # hit's genomic window in one reference allele of its antigen. Present only
    # when the database ships k-mer sketches and the caller asked for it, so
    # every existing report stays byte-identical.
    sketch_pct: float | None = None
    # phase affinity: how well this hit's window ENDS match its own family's
    # alleles of the same antigen. `fliC` and `fljB` `l,v` alleles are ~99 %
    # identical in the middle and differ at the ends, so comparing the two
    # families' values tells which phase a contested locus is.
    phase_pct: float | None = None


def swap_dedup(
    rows: Iterable["BlastRow"], *, dedup: Literal["query", "subject"] = "query"
) -> Iterator[tuple["BlastRow", int, int]]:
    """SPEC.md §4 steps 1-2: minus-strand swap (subject coords only) and dedup —
    first row wins, the key ignores strand, and the key is claimed even if the
    row is later coverage-filtered.

    ``dedup="query"`` (the abricate-parity default) keys on
    ``(qseqid, qstart, qend)``: one survival per query interval, so when
    homologous references share a query block only the first subject keeps it.
    ``dedup="subject"`` keys on ``(sseqid, qseqid, qstart, qend)`` instead, so
    every reference keeps its own HSPs — required by ``coverage="union"``,
    which aggregates all of a subject's blocks into one breadth figure.
    """
    seen: set[tuple[str | int, ...]] = set()
    for row in rows:
        if row.sstrand == "minus":
            s_start, s_end = row.send, row.sstart
        else:
            s_start, s_end = row.sstart, row.send
        key: tuple[str | int, ...] = (
            (row.qseqid, row.qstart, row.qend)
            if dedup == "query"
            else (row.sseqid, row.qseqid, row.qstart, row.qend)
        )
        if key in seen:
            continue
        seen.add(key)
        yield row, s_start, s_end


def hit_from_row(
    row: "BlastRow", s_start: int, s_end: int, coverage_pct: float, default_db: str
) -> Hit:
    """SPEC.md §4 steps 4-5: subject-id decode (native ``gspec|`` codec or
    legacy ``~~~`` rules) and product cleanup, assembled into a Hit."""
    header = decode_seqid(row.sseqid, default_db)
    product = row.stitle or "n/a"
    product = product.replace(",", "").replace("\t", "")
    if IDSEP in product or is_gspec_header(row.sseqid):
        product = _LEADING_TOKEN_RE.sub("", product, count=1)
    return Hit(
        sequence=row.qseqid,
        start=row.qstart,
        end=row.qend,
        strand="-" if row.sstrand == "minus" else "+",
        gene=header.gene,
        database=header.database,
        accession=header.accession,
        function=header.function,
        product=product,
        s_start=s_start,
        s_end=s_end,
        s_len=row.slen,
        coverage_map=minimap(s_start, s_end, row.slen, row.gapopen),
        gap_openings=row.gapopen,
        gaps=row.gaps,
        identity_pct=row.pident,
        coverage_pct=coverage_pct,
        aligned_len=row.length - row.gaps,
    )


def process_rows(
    rows: Iterable["BlastRow"],
    *,
    mincov: float,
    default_db: str,
    dedup: Literal["query", "subject"] = "query",
) -> list[Hit]:
    """Turn BLAST rows into hits, strictly in SPEC.md §4 order:

    1. minus-strand swap (subject coords only), 2. dedup on
    (qseqid, qstart, qend) — first row wins, key ignores strand, and the key
    is claimed even if the row is later coverage-filtered, 3. coverage filter
    on the unrounded float, 4. header decode (native ``gspec|`` codec or
    legacy ``~~~`` rules), 5. product cleanup.

    ``dedup="subject"`` widens the dedup key to include the subject id; see
    :func:`swap_dedup`.
    """
    hits: list[Hit] = []
    for row, s_start, s_end in swap_dedup(rows, dedup=dedup):
        # 3. coverage filter on the unrounded float (decode stays after the
        # filter: a malformed native header on a filtered row never raises)
        coverage_pct = 100.0 * (row.length - row.gaps) / row.slen
        if coverage_pct < mincov:
            continue
        hits.append(hit_from_row(row, s_start, s_end, coverage_pct, default_db))
    return hits


def process_rows_with_partials(
    rows: Iterable["BlastRow"],
    *,
    mincov: float,
    default_db: str,
    dedup: Literal["query", "subject"] = "query",
) -> tuple[list[Hit], list[Hit]]:
    """Merge-mode variant of :func:`process_rows`: returns ``(hits, partials)``
    where ``partials`` are the sub-mincov dedup survivors, decoded so the
    fragment-merge pass (:mod:`gspec.align.fragments`) can reason about them.

    Opt-in ``--merge-fragments`` path only — every dedup survivor is decoded
    (step 4 runs before the split), unlike the parity path. ``process_rows``
    stays the default and byte-identical to abricate.
    """
    hits: list[Hit] = []
    partials: list[Hit] = []
    for row, s_start, s_end in swap_dedup(rows, dedup=dedup):
        coverage_pct = 100.0 * (row.length - row.gaps) / row.slen
        hit = hit_from_row(row, s_start, s_end, coverage_pct, default_db)
        (hits if coverage_pct >= mincov else partials).append(hit)
    return hits, partials


class Report(BaseModel, frozen=True):
    """Canonical in-memory result of screening one input file.

    Hits are sorted (sequence, start) before construction; the tuple is the
    final, stable order. Screening is pure gene detection: typed and
    untyped gene databases produce byte-identical reports (designation is
    the ``gspec typing`` command's job).
    """

    file: str
    hits: tuple[Hit, ...]


def minimap(x: int, y: int, length: int, broken: int) -> str:
    """15-char coverage bar: '=' over [x, y] scaled onto the map width.

    Replicates the Perl int() truncation quirks verbatim — do not "fix" the
    1-based coordinate artifacts (yi may exceed width-1). ``broken`` is the
    gap-open count; a broken map is 14 boxes with '/' after box width//2.
    """
    width = 15 - (1 if broken > 0 else 0)
    scale = length / width
    xi = int(x / scale)
    yi = int(y / scale)
    chars: list[str] = []
    for i in range(width):
        chars.append("=" if xi <= i <= yi else ".")
        if broken > 0 and i == width // 2:
            chars.append("/")
    return "".join(chars)
