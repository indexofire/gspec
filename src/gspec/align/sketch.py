"""Anchored allele sketches: which reference allele sits in THIS window?

Antigen panels are full of near-identical references, and a BLAST HSP cannot
always tell them apart: gapped extension and partial coverage move an
identity x coverage score by more than the true difference between two alleles
(measured on real genomes: `fljB` `1,6` lost to `1,7`, `fliC` `e,h` lost to
`e,n,x`). Counting k-mers only inside the window the alignment already found asks
a cleaner question — "whose k-mers are here" — and it separated near-ties that
alignment could not (78.7 % -> 94.7 % agreement with SISTR's antigen tokens over
94 H1/H2 calls, 60 genomes).

Two shapes were measured and only the anchored one works:

- a whole *genome* against an antigen's k-mer set saturates (a 4.8 Mb query
  against a 1.5 kb set) and degenerates into "which antigen has the biggest
  reference pool" — 53 %, worse than alignment;
- the anchored window against each reference *allele* is the signal.

This module is deliberately free of BLAST types: it takes a window and a set of
allele sketches.
"""

from __future__ import annotations

import struct
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from gspec.cluster.kmerdb import count_sequence_kmers, read_kmerdb
from gspec.fasta import iter_fasta

#: Bases of flanking sequence included either side of the aligned window. The
#: window edges sit on the alignment boundary, which for a partial HSP can fall
#: inside the gene; the padding keeps the k-mers that straddle it.
WINDOW_PAD = 60

#: How far apart two sketch scores must be before the sketch is allowed to
#: decide a margin tie. 100 % containment of one allele against 95 % of another
#: is a real call; two alleles both at ~100 % are indistinguishable by k-mers.
SKETCH_MARGIN = 0.05

#: Lowest sketch score that names an antigen. Below this the window simply does
#: not carry the allele (a divergent locus, a partial hit, or the wrong family).
MIN_SKETCH = 0.5


@dataclass(frozen=True)
class AlleleSketch:
    """One reference allele's canonical k-mer set."""

    family: str
    antigen: str
    codes: frozenset[int]


def _gene_of(header: str) -> str | None:
    """The ``gene=`` field of a gspec/v1 header (None for other formats)."""
    for part in header.split("|"):
        if part.startswith("gene="):
            return part.split("=", 1)[1]
    return None


def _split_gene(gene: str) -> tuple[str, str]:
    """``h1_e_h__3`` -> ``("h1", "e,h")``.

    Allele suffixes (``__3``) are stripped and the antigen's underscores are the
    table's commas, so the name matches what ``serotype.json`` declares.
    """
    antigen = gene.rsplit("__", 1)[0]
    family, _, name = antigen.partition("_")
    return family, name.replace("_", ",")


@dataclass(frozen=True)
class SketchIndex:
    """Every reference allele of one database, grouped by antigen family."""

    k: int
    families: tuple[tuple[str, tuple[AlleleSketch, ...]], ...]

    def of(self, family: str) -> tuple[AlleleSketch, ...]:
        """The sketches of one family (empty when the family is not a panel)."""
        for name, sketches in self.families:
            if name == family:
                return sketches
        return ()


@lru_cache(maxsize=8)
def allele_sketches(db_dir: Path) -> SketchIndex | None:
    """Per-family allele sketches for one database, plus the k they use.

    Reads the packed ``.full.u64`` arrays ``db build --kmerdb`` materializes: one
    entry per reference allele, which is what makes an anchored comparison an
    allele-level question instead of a pool-level one.
    """
    document = read_kmerdb(db_dir)
    if document is None:
        return None
    k = document.k
    antigen_of: dict[str, str] = {}
    for record in iter_fasta(db_dir / "sequences"):
        gene = _gene_of(record.id)
        if gene is not None:
            antigen_of[record.id] = gene
    by_family: dict[str, list[AlleleSketch]] = {}
    for path in sorted((db_dir / "kmers").glob("*.full.u64")):
        gene = antigen_of.get(path.name[: -len(".full.u64")])
        if gene is None:
            continue
        raw = path.read_bytes()
        codes = frozenset(struct.unpack(f"<{len(raw) // 8}Q", raw))
        if not codes:
            continue
        family, antigen = _split_gene(gene)
        by_family.setdefault(family, []).append(AlleleSketch(family, antigen, codes))
    return SketchIndex(
        k=k, families=tuple((family, tuple(items)) for family, items in by_family.items())
    )


def window_kmers(
    contigs: Mapping[str, str],
    sequence_field: str,
    start: int,
    end: int,
    k: int,
    *,
    pad: int = WINDOW_PAD,
) -> frozenset[int]:
    """The canonical k-mers of a hit's genomic window.

    ``sequence_field`` is the hit's SEQUENCE cell, which union coverage renders as
    a comma-joined contig list (a locus split across contigs), and query
    coordinates are shared by every contig in it.
    """
    codes: set[int] = set()
    for contig in sequence_field.split(","):
        sequence = contigs.get(contig)
        if sequence:
            codes |= count_sequence_kmers(sequence[max(0, start - pad) : end + pad], k)
    return frozenset(codes)


def antigen_scores(
    window: frozenset[int],
    sketches: Iterable[AlleleSketch],
) -> dict[str, float]:
    """``antigen -> best allele containment`` for one window.

    The score is ``|window ∩ allele| / |allele|`` over the best allele of that
    antigen: an antigen is "present" when one of its alleles is fully covered,
    which is the reads-path definition applied to a locus instead of a read set.
    """
    if not window:
        return {}
    scores: dict[str, float] = {}
    for sketch in sketches:
        shared = len(window & sketch.codes)
        if not shared:
            continue
        fraction = shared / len(sketch.codes)
        if fraction > scores.get(sketch.antigen, 0.0):
            scores[sketch.antigen] = fraction
    return scores


__all__ = [
    "MIN_SKETCH",
    "SKETCH_MARGIN",
    "WINDOW_PAD",
    "AlleleSketch",
    "SketchIndex",
    "allele_sketches",
    "antigen_scores",
    "window_kmers",
]
