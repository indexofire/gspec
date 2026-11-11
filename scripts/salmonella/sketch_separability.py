#!/usr/bin/env python3
"""Can whole-genome similarity separate serovars at all? Measure before building.

Twice now the answer decided a design: a per-serovar reference-genome sketch (the
mash/cgMLST idea) would resolve the antigen ties the ranker cannot — but only if a
Jaccard threshold exists that puts every same-serovar pair above it and every
different-serovar pair below. This script answers that question from a corpus and a
truth table, so the next person who proposes a similarity channel gets a number
instead of an opinion.

Measured on 125 in-house Salmonella assemblies (k=21, bottom-1000 MinHash):

* within-serovar Jaccard: 0.630 / 0.864 / 0.979 (min/median/max, 135 pairs)
* between-serovar:        0.016 / 0.519 / 0.997 (7615 pairs)
* 5.7 % of the between pairs sit *above* the closest within pair — and the
  most-similar different-serovar pair (0.997) is closer than the most-similar
  same-serovar pair (0.979): no threshold exists, and it is not a tuning problem.

Sketch parameters are not the limit. A bottom-k sketch estimates Jaccard with a
standard error of about ``sqrt(J*(1-J)/s)`` (1.9 % at J=0.9, s=1000), and the
k-mer/ANI relation ``J ~ p**k / (2 - p**k)`` makes k=21 resolve a 0.5 % ANI
difference by 8 Jaccard points. The overlap above is biological (accessory
genome, phage, recombination, and clonal variants such as Typhimurium versus its
monophasic `1,4,[5],12:i:-` form), so no choice of k or s removes it.

Usage::

    python scripts/salmonella/sketch_separability.py \\
        --assemblies ~/assemblies --sistr-tab sistr_all.tab [--size 1000 --k 21]

Exit status is 0 when a threshold exists and 1 when it does not, so a pipeline can
branch on the answer.
"""

from __future__ import annotations

import argparse
import csv
import heapq
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from gspec.cluster.kmerdb import iter_canonical_kmer_codes
from gspec.fasta import iter_fasta

FASTA_SUFFIXES = (".fna", ".fa", ".fasta", ".fas")


def sketch(path: Path, *, size: int, k: int) -> frozenset[int]:
    """Bottom-k MinHash sketch: the ``size`` smallest canonical k-mer codes."""
    heap: list[int] = []
    for record in iter_fasta(path):
        for code in iter_canonical_kmer_codes(record.sequence, k):
            if len(heap) < size:
                heapq.heappush(heap, -code)
            elif code < -heap[0]:
                heapq.heapreplace(heap, -code)
    return frozenset(-code for code in heap)


def jaccard(first: frozenset[int], second: frozenset[int]) -> float:
    """MinHash Jaccard estimate (exact for the sketches themselves)."""
    if not first or not second:
        return 0.0
    return len(first & second) / len(first | second)


def _truth(path: Path) -> dict[str, tuple[str, ...]]:
    """``stem -> serovar set`` from a SISTR ``-f tab`` dump."""
    truth: dict[str, tuple[str, ...]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            stem = Path(row.get("fasta_filepath") or "").stem
            if not stem:
                continue
            names = tuple(name for name in (row.get("serovar") or "").split("|") if name)
            truth[stem] = names or (row.get("antigenic_formula") or "-",)
    return truth


def separable(within: list[float], between: list[float]) -> bool:
    """Does a threshold put every same-serovar pair above every different one?"""
    return bool(within and between) and min(within) > max(between)


def _single(label: tuple[str, ...]) -> bool:
    """Only unambiguous truth rows are usable for a separability test."""
    return len(label) == 1 and not label[0].startswith("-") and ":" not in label[0]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assemblies", type=Path, required=True)
    parser.add_argument("--sistr-tab", type=Path, required=True)
    parser.add_argument("--size", type=int, default=1000, help="bottom-k sketch size")
    parser.add_argument("--k", type=int, default=21)
    parser.add_argument("--out", type=Path, help="write the pairwise distances here (JSON)")
    args = parser.parse_args(argv)

    truth = _truth(args.sistr_tab)
    files = [
        path
        for path in sorted(args.assemblies.glob("*"))
        if path.suffix.lower() in FASTA_SUFFIXES and path.stem in truth
    ]
    usable = [path for path in files if _single(truth[path.stem])]
    if len(usable) < 10:
        raise SystemExit(
            f"only {len(usable)} assembly/assemblies carry a single serovar in the truth table;"
            " a separability test needs unambiguous labels"
        )
    print(f"sketching {len(usable)} genome(s) with k={args.k}, bottom-{args.size}")

    sketches: dict[str, frozenset[int]] = {}
    for index, path in enumerate(usable):
        sketches[path.stem] = sketch(path, size=args.size, k=args.k)
        if (index + 1) % 25 == 0:
            print(f"  {index + 1}/{len(usable)}", flush=True)

    within: list[float] = []
    between: list[float] = []
    pairs: list[dict[str, object]] = []
    stems = sorted(sketches)
    for position, first in enumerate(stems):
        for second in stems[position + 1 :]:
            distance = jaccard(sketches[first], sketches[second])
            same = truth[first] == truth[second]
            (within if same else between).append(distance)
            pairs.append(
                {"first": first, "second": second, "jaccard": round(distance, 4), "same": same}
            )

    def summary(values: list[float]) -> str:
        if not values:
            return "(none)"
        return (
            f"{len(values)} pair(s): min {min(values):.3f}"
            f" median {statistics.median(values):.3f} max {max(values):.3f}"
        )

    print(f"\nwithin a serovar : {summary(within)}")
    print(f"between serovars : {summary(between)}")
    threshold_exists = separable(within, between)
    if not threshold_exists and within and between:
        intrusions = sum(1 for value in between if value >= min(within))
        print(
            f"\nNOT separable: {intrusions} of {len(between)} between-serovar pair(s)"
            f" ({100 * intrusions / len(between):.1f} %) sit at or above the closest"
            f" within-serovar pair ({min(within):.3f})."
        )
        print("  the most similar different-serovar pairs:")
        for entry in sorted(pairs, key=lambda entry: -float(entry["jaccard"]))[:5]:
            if not entry["same"]:
                print(
                    f"    {entry['first']} / {entry['second']}"
                    f"  jaccard {entry['jaccard']}  ({truth[str(entry['first'])]}"
                    f" vs {truth[str(entry['second'])]})"
                )
        print(
            "\n  k and size are not the limit: bottom-k error is"
            " sqrt(J*(1-J)/size) and k rescales ANI differences; this overlap is"
            " biological. A similarity channel cannot designate serovars here."
        )
    else:
        print(f"\nseparable at {min(within):.3f} < J <= {max(between):.3f}")

    if args.out is not None:
        import json

        args.out.write_text(
            json.dumps(
                {
                    "k": args.k,
                    "size": args.size,
                    "separable": threshold_exists,
                    "within": [round(value, 4) for value in within],
                    "between": [round(value, 4) for value in between],
                    "pairs": pairs,
                },
                indent=1,
            )
            + "\n",
            encoding="utf-8",
        )
        print(f"\nwrote {args.out}")
    return 0 if threshold_exists else 1


if __name__ == "__main__":
    raise SystemExit(main())
