"""BLAST invocation, tabular parsing, and the screening pipeline (SPEC.md §3)."""

import re
import shlex
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import Literal

from pydantic import BaseModel

from gspec.align.fragments import merge_fragment_hits, union_coverage_map, union_coverage_pct
from gspec.align.gseek import engine, run_blast_gseek
from gspec.db.core import Database
from gspec.errors import DependencyError, GspecError, InputError
from gspec.proctools import run_tool
from gspec.report import (
    Hit,
    Report,
    ScreeningParams,
    hit_from_row,
    process_rows,
    process_rows_with_partials,
    swap_dedup,
)
from gspec.screening.floors import apply_hit_floors, gene_floors
from gspec.seqconvert import SeqFormat, detect_format, to_fasta_lines

BLAST_FIELDS = [
    "qseqid",
    "qstart",
    "qend",
    "qlen",
    "sseqid",
    "sstart",
    "send",
    "slen",
    "sstrand",
    "evalue",
    "length",
    "pident",
    "gaps",
    "gapopen",
    "stitle",
]

_OUTFMT = "6 " + " ".join(BLAST_FIELDS)
_MIN_BLAST_VERSION = (2, 2, 30)


class BlastRow(BaseModel, frozen=True):
    """One outfmt-6 row, typed at the boundary (SPEC.md §3 field order)."""

    qseqid: str
    qstart: int
    qend: int
    qlen: int
    sseqid: str
    sstart: int
    send: int
    slen: int
    sstrand: str
    evalue: float
    length: int
    pident: float
    gaps: int
    gapopen: int
    stitle: str


def parse_blast_row(line: str) -> BlastRow:
    """Parse one tab-delimited outfmt-6 line; a row with != 15 columns is a
    hard error (upstream wording)."""
    fields = line.split("\t")
    if len(fields) != 15:
        raise GspecError("can not find sequence data", code="BLAST_PARSE_FAILED")
    try:
        return BlastRow(
            qseqid=fields[0],
            qstart=int(fields[1]),
            qend=int(fields[2]),
            qlen=int(fields[3]),
            sseqid=fields[4],
            sstart=int(fields[5]),
            send=int(fields[6]),
            slen=int(fields[7]),
            sstrand=fields[8],
            evalue=float(fields[9]),
            length=int(fields[10]),
            pident=float(fields[11]),
            gaps=int(fields[12]),
            gapopen=int(fields[13]),
            stitle=fields[14],
        )
    except (TypeError, ValueError) as exc:
        raise GspecError(
            f"can not parse BLAST row (numeric field): {line!r}",
            code="BLAST_PARSE_FAILED",
        ) from exc


def ensure_blast() -> None:
    """Require blastn >= 2.2.30 on PATH (SPEC.md §3 version gate)."""
    result = run_tool(["blastn", "-version"])
    if result.returncode != 0:
        raise DependencyError(
            f"blastn -version failed: {result.stderr.strip()}",
            code="BLAST_VERSION_CHECK_FAILED",
        )
    match = re.search(r"blastn:\s*(\d+)\.(\d+)\.(\d+)", result.stdout)
    if match is None:
        raise DependencyError(
            f"could not parse blastn version from: {result.stdout!r}",
            code="BLAST_VERSION_CHECK_FAILED",
        )
    major, minor, revision = (int(part) for part in match.groups())
    if (major, minor, revision) < _MIN_BLAST_VERSION:
        raise DependencyError(
            f"blastn {major}.{minor}.{revision} is older than 2.2.30",
            code="BLAST_VERSION_TOO_OLD",
        )


def _normalize(query: Path) -> tuple[str, SeqFormat]:
    """Convert one input file to FASTA text in-process (seqconvert), wrapping
    any parser failure in the blast-path error shape
    ``invalid input file {query}: <reason>``.

    The whole normalized FASTA is held in memory and fed to blast via
    ``input=`` — genome-scale inputs are a few MB, a deliberate trade-off:
    ``subprocess.run``'s communicate() owns stdin/stdout concurrency, so the
    retired ``any2fasta | blastn`` Popen chain needs no stderr-drain thread
    here anymore.
    """
    try:
        fmt = detect_format(query)
        text = "".join(to_fasta_lines(query, fmt))
    except InputError as exc:
        raise InputError(
            f"invalid input file {query}: {exc}",
            code="INVALID_INPUT",
            context={"file": str(query)},
        ) from exc
    return text, fmt


def _pipeline(fasta_text: str, argv: list[str]) -> str:
    """``<argv>`` reading the normalized FASTA on stdin; returns stdout as text."""
    result = run_tool(argv, stdin_text=fasta_text)
    if result.returncode != 0:
        raise GspecError(
            f"{argv[0]} failed: {result.stderr.strip()}",
            code="BLAST_FAILED",
            context={"binary": argv[0]},
        )
    return result.stdout


def _hsps(params: ScreeningParams) -> list[str]:
    """HSP-selection flags: parity culling, or none for union coverage.

    ``-culling_limit 1`` keeps one HSP per subject (abricate parity), which is
    exactly what makes single-HSP coverage measure only the largest block of a
    split locus. ``coverage="union"`` needs every block, so culling is dropped
    there and BLAST's default (unlimited) HSP count applies.
    """
    return [] if params.coverage == "union" else ["-culling_limit", "1"]


def run_blastn(
    query: Path,
    database: Database,
    params: ScreeningParams,
    *,
    dbtype: Literal["nucl", "prot"],
    debug: bool = False,
) -> list[BlastRow]:
    """Run the normalize (native seqconvert) -> blastn/blastx pipeline for one
    query file.

    The database must already be indexed; protein databases switch to blastx
    without -perc_identity (upstream quirk, minid silently ignored). With
    ``debug``, echo the normalization step and the exact blast argv to stderr
    (abricate --debug parity). ``dbtype`` comes from the caller resolving it
    once per run, so dependency/index errors surface before any per-file work.
    """
    if dbtype == "prot":
        argv = [
            "blastx",
            "-task",
            "blastx-fast",
            "-seg",
            "no",
            "-db",
            str(database.sequences_path),
            "-outfmt",
            _OUTFMT,
            "-num_threads",
            str(params.threads),
            "-evalue",
            "1E-20",
            *_hsps(params),
            "-max_target_seqs",
            "10000",
        ]
        sys.stderr.write("--minid is not applied to protein databases (abricate parity)\n")
    else:
        argv = [
            "blastn",
            "-task",
            "blastn",
            "-dust",
            "no",
            "-perc_identity",
            str(params.minid),
            "-db",
            str(database.sequences_path),
            "-outfmt",
            _OUTFMT,
            "-num_threads",
            str(params.threads),
            "-evalue",
            "1E-20",
            *_hsps(params),
            "-max_target_seqs",
            "10000",
        ]
    fasta_text, fmt = _normalize(query)
    switch = engine()
    mode = switch.mode.value
    if debug:
        sys.stderr.write(f"gspec: normalize: {query} ({fmt.value})\n")
        if mode in {"legacy", "shadow"}:
            sys.stderr.write(f"gspec: run: {shlex.join(argv)}\n")

    def legacy_rows() -> list[BlastRow]:
        output = _pipeline(fasta_text, argv)
        return [parse_blast_row(line) for line in output.splitlines() if line.strip()]

    def gseek_rows() -> list[BlastRow]:
        return run_blast_gseek(fasta_text, database, params, dbtype=dbtype, debug=debug)

    return switch.run(legacy_rows, gseek_rows)


#: Two blocks of one reference whose reference spans overlap by more than this
#: fraction of the shorter span are alternate placements of the same locus, so
#: only the better one is kept. Disjoint spans (a locus split across contigs)
#: stay additive.
_SAME_LOCUS_OVERLAP = 0.5


def _distinct_loci(blocks: list[tuple["BlastRow", int, int]]) -> list[tuple["BlastRow", int, int]]:
    """Drop blocks that re-cover an already accepted block of the same reference."""
    ranked = sorted(
        blocks,
        key=lambda block: (block[0].length - block[0].gaps) * block[0].pident,
        reverse=True,
    )
    kept: list[tuple[BlastRow, int, int]] = []
    for block in ranked:
        _, start, end = block
        width = end - start + 1
        covered = False
        for _, other_start, other_end in kept:
            other_width = other_end - other_start + 1
            overlap = min(end, other_end) - max(start, other_start) + 1
            if overlap > 0 and overlap >= _SAME_LOCUS_OVERLAP * min(width, other_width):
                covered = True
                break
        if not covered:
            kept.append(block)
    return kept


def union_hits(
    rows: Iterable[BlastRow],
    *,
    mincov: float,
    min_aligned_len: int | None,
    default_db: str,
) -> list[Hit]:
    """Fold every reference's HSPs into ONE union-breadth hit (``coverage="union"``).

    The evidence model this implements is the one a divergent or repeat-laden
    locus needs: a 97%-identity locus aligns as several blocks, so single-HSP
    coverage measures the largest block while the true breadth (the union of
    the subject intervals) is near 100%. Per subject:

    - coverage = union of the subject intervals / subject length (capped at 100)
    - identity = aligned-length weighted mean of the block identities
    - sequence = comma-joined contigs (a locus split across contigs counts once)
    - anchor = the longest block (query span, product, strand)

    Blocks are grouped across contigs on purpose: union breadth is additive,
    so a split assembly reaches the same figure as a finished one.
    ``min_aligned_len`` drops subjects whose total ungapped aligned length is
    too small to be a real locus (spurious short-HSP guard).

    Two blocks that cover the SAME part of the reference are alternative loci,
    not fragments of one locus: a subject that is a perfect full-length gene at
    one site and a paralogous remnant at another would otherwise average to a
    mediocre identity (measured: a 98.9 % fljB block merged with an 80.6 %
    remnant reported 90.6 % and lost the antigen). Such overlaps are resolved by
    keeping the better block, so the union only ever sums disjoint spans.
    """
    groups: dict[str, list[tuple[BlastRow, int, int]]] = {}
    for row, s_start, s_end in swap_dedup(rows, dedup="subject"):
        groups.setdefault(row.sseqid, []).append((row, s_start, s_end))
    hits: list[Hit] = []
    for blocks in (_distinct_loci(group) for group in groups.values()):
        slen = blocks[0][0].slen
        if slen <= 0:
            continue
        aligned = sum(row.length - row.gaps for row, _, _ in blocks)
        if min_aligned_len is not None and aligned < min_aligned_len:
            continue
        spans = [(s_start, s_end) for _, s_start, s_end in blocks]
        coverage_pct = min(100.0, union_coverage_pct(spans, slen))
        if coverage_pct < mincov:
            continue
        anchor = max(blocks, key=lambda block: block[0].length - block[0].gaps)[0]
        identity_pct = (
            sum((row.length - row.gaps) * row.pident for row, _, _ in blocks) / aligned
            if aligned
            else 0.0
        )
        gap_openings = sum(row.gapopen for row, _, _ in blocks)
        hit = hit_from_row(
            anchor,
            min(start for _, start, _ in blocks),
            max(end for _, _, end in blocks),
            coverage_pct,
            default_db,
        ).model_copy(
            update={
                "sequence": ",".join(sorted({row.qseqid for row, _, _ in blocks})),
                "coverage_map": union_coverage_map(spans, slen, broken=gap_openings > 0),
                "gap_openings": gap_openings,
                "gaps": sum(row.gaps for row, _, _ in blocks),
                "identity_pct": identity_pct,
                "aligned_len": aligned,
            }
        )
        hits.append(hit)
    return hits


def report_from_rows(
    rows: Iterable[BlastRow],
    query: Path,
    database: Database,
    params: ScreeningParams,
    *,
    merge_fragments: bool = False,
) -> Report:
    """Fold typed BLAST rows into one sorted Report.

    Shared by the per-file and batch gseek paths; the caller owns the BLAST
    run and the per-run ``ensure_blast``/dbtype gates.
    """
    floors = gene_floors(database)
    if params.coverage == "union":
        hits = union_hits(
            rows,
            mincov=params.mincov,
            min_aligned_len=params.min_aligned_len,
            default_db=params.db,
        )
        if floors is not None:
            hits = apply_hit_floors(hits, floors)
    elif merge_fragments:
        hits, partials = process_rows_with_partials(
            rows, mincov=params.mincov, default_db=params.db, dedup=params.dedup
        )
        if floors is not None:
            hits = apply_hit_floors(hits, floors)
            partials = apply_hit_floors(partials, floors)
        hits = merge_fragment_hits(hits, partials, mincov=params.mincov)
    else:
        hits = process_rows(rows, mincov=params.mincov, default_db=params.db, dedup=params.dedup)
        if floors is not None:
            hits = apply_hit_floors(hits, floors)
    ordered = sorted(hits, key=lambda hit: (hit.sequence, hit.start))
    return Report(file=str(query), hits=tuple(ordered))


def screen_file(
    query: Path,
    database: Database,
    params: ScreeningParams,
    *,
    dbtype: Literal["nucl", "prot"],
    debug: bool = False,
    merge_fragments: bool = False,
) -> Report:
    """Screen one input file against one database into a sorted Report.

    The caller owns the per-run gates (``ensure_blast`` and the one-shot
    ``dbtype`` resolution) so a multi-file run pays each probe once. With
    ``merge_fragments`` (opt-in gspec extension), sub-mincov dedup survivors
    are folded into one hit per gene when their union coverage reaches mincov.
    """
    rows = run_blastn(query, database, params, dbtype=dbtype, debug=debug)
    return report_from_rows(rows, query, database, params, merge_fragments=merge_fragments)
