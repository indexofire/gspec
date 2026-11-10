#!/usr/bin/env python3
"""Audit cross-labelled reference alleles in an antigen panel.

Upstream antigen databases label sequences by antigen name, but the same sequence
(or a near-copy) sometimes appears under two names: SISTR's `1,2` and `1,2,7`
references are 99.5 % identical, `r` and `r,[i]` likewise. Because the ranker picks
the best-scoring allele by name, such a pair can bury the right antigen — a genome
carrying `1,2` matched a `1,2,7` reference at 100 % and lost the `1,2` rows entirely
before the similarity-class widening was added.

This script finds those pairs from the panel itself, using the k-mer substrate the
reads path already builds: two alleles are one class when the shorter one's k-mers
are (almost) all contained in the longer one. The output is a
``gspec.cross_labels/1`` document that ``import_upstream.py --cross-labels`` folds
into the document's ``similarity_groups``, so both names' rows stay reachable and
the ambiguity is reported instead of being resolved by bitscore noise.

Usage::

    python scripts/salmonella/audit_cross_labels.py \\
        --db-dir ~/.local/share/gspec/db/sen_sero_sistr --out cross_labels.json
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from gspec.align.sketch import allele_sketches

#: Containment of the shorter allele's k-mers in the longer one, above which the
#: two references cannot be told apart by a k-mer sketch.
CONTAINMENT = 0.98

#: Unrelated antigens share short conserved stretches; require a real overlap.
MIN_SHARED = 200


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--containment", type=float, default=CONTAINMENT)
    args = parser.parse_args(argv)

    index = allele_sketches(args.db_dir)
    if index is None:
        raise SystemExit(f"{args.db_dir}: no k-mer substrate (build with --kmerdb)")

    groups: dict[str, list[list[str]]] = {}
    pairs: list[dict[str, object]] = []
    for family, sketches in index.families:
        # One representative per antigen: the longest allele (the same rule the
        # k-mer substrate uses), so the comparison stays antigen-to-antigen.
        by_antigen: dict[str, frozenset[int]] = {}
        for sketch in sketches:
            if len(sketch.codes) > len(by_antigen.get(sketch.antigen, frozenset())):
                by_antigen[sketch.antigen] = sketch.codes
        names = sorted(by_antigen)
        adjacency: dict[str, set[str]] = defaultdict(set)
        for i, first in enumerate(names):
            for second in names[i + 1 :]:
                a, b = by_antigen[first], by_antigen[second]
                shared = len(a & b)
                if shared < MIN_SHARED:
                    continue
                containment = shared / min(len(a), len(b))
                if containment < args.containment:
                    continue
                adjacency[first].add(second)
                adjacency[second].add(first)
                pairs.append(
                    {
                        "family": family,
                        "first": first,
                        "second": second,
                        "containment": round(containment, 4),
                        "shared_kmers": shared,
                    }
                )
        # Connected components become one class: `1,2` -> `1,2,7` -> `1,[2],7`
        # is a chain of the same sequence under three labels.
        seen: set[str] = set()
        family_groups: list[list[str]] = []
        for name in names:
            if name in seen or name not in adjacency:
                continue
            stack, component = [name], []
            while stack:
                current = stack.pop()
                if current in seen:
                    continue
                seen.add(current)
                component.append(current)
                stack.extend(sorted(adjacency[current]))
            family_groups.append(sorted(component))
        if family_groups:
            groups[family] = family_groups

    document = {
        "schema": "gspec.cross_labels/1",
        "containment_floor": args.containment,
        "min_shared_kmers": MIN_SHARED,
        "groups": groups,
        "pairs": sorted(pairs, key=lambda pair: -float(pair["containment"])),
    }
    args.out.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")

    for family, family_groups in sorted(groups.items()):
        print(f"{family}: {len(family_groups)} cross-labelled class(es)")
        for component in family_groups:
            print("  " + " == ".join(component))
    worst = document["pairs"][:5]
    if worst:
        print("\ntightest pairs:")
        for pair in worst:
            print(
                f"  {pair['family']} {pair['first']} / {pair['second']}"
                f"  containment {pair['containment']}  shared {pair['shared_kmers']}"
            )
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
