#!/usr/bin/env python3
"""Call antigens from reads by k-mer containment (no alignment).

The panel ships one differential k-mer set per antigen (`kmer.json`,
`gspec.kmerdb/1`): the antigen's representative reference minus every other
antigen of the same family. Counting the reads once and scoring each set by
``set_hit_fraction`` answers "which O group / which H1 allele / which H2 allele"
without mapping — which is what the alignment path struggles with when an antigen
has several near-identical references (the reads split between them).

Usage::

    python scripts/salmonella/kmers_call_reads.py \\
        --db-dir src/gspec/data/dbs/sen_sero_sistr --fastq reads.fq [--k 27]

The sets are materialized on first use under the database directory
(`materialize_kmer_sets`), so the command works on a copied bundle too.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from gspec.cluster.kmerdb import (
    containment_to_identity,
    count_fastq_kmers,
    load_kmer_set,
    materialize_kmer_sets,
    read_kmerdb,
    set_hit_fraction,
)

#: A set is called when at least this fraction of its differential k-mers is
#: present in the reads. Real 30x data clears it by a wide margin; a genome of
#: another serovar leaves the set mostly absent.
MIN_FRACTION = 0.2


def _ensure_sets(db_dir: Path, *, k: int) -> tuple[dict[str, float], int]:
    """Materialize the declared sets (once) and return ``{set: fraction}``."""
    document = read_kmerdb(db_dir)
    if document is None:
        raise SystemExit(f"{db_dir}: no kmer.json (run import_upstream.py --kmerdb)")
    if not (db_dir / "kmers").is_dir():
        materialize_kmer_sets(db_dir, document, db_dir / "sequences")
        document = read_kmerdb(db_dir) or document
    return {}, document.k


def call(db_dir: Path, fastq: list[Path], *, min_fraction: float = MIN_FRACTION) -> dict[str, str]:
    """Best antigen per family, from the reads' k-mer containment."""
    _, k = _ensure_sets(db_dir, k=0)
    query = count_fastq_kmers(*fastq, k=k)
    families: dict[str, list[tuple[float, str]]] = {}
    document = read_kmerdb(db_dir)
    assert document is not None
    for name in document.sets:
        family = name.split("_", 1)[0]
        fraction = set_hit_fraction(query, load_kmer_set(db_dir, name))
        families.setdefault(family, []).append((fraction, name))
    calls: dict[str, str] = {}
    for family, scored in sorted(families.items()):
        scored.sort(reverse=True)
        best_fraction, best_name = scored[0]
        runner_up = scored[1][0] if len(scored) > 1 else 0.0
        if best_fraction < min_fraction:
            calls[family] = ""
            continue
        antigen = best_name.split("_", 1)[1] if "_" in best_name else best_name
        calls[family] = antigen
        print(
            f"  {family}: {antigen} ({best_fraction:.2f} of its differential k-mers,"
            f" runner-up {runner_up:.2f}, ~{containment_to_identity(best_fraction, k):.1f}% ANI)"
        )
    return calls


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-dir", type=Path, required=True)
    parser.add_argument("--fastq", type=Path, nargs="+", required=True)
    parser.add_argument("--min-fraction", type=float, default=MIN_FRACTION)
    args = parser.parse_args(argv)
    calls = call(args.db_dir, args.fastq, min_fraction=args.min_fraction)
    print("formula: " + ":".join(calls.get(family) or "-" for family in ("o", "h1", "h2")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
