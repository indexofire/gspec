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
import shutil
import sys
import tempfile
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


#: Per-source scratch dirs holding materialized sets. A bundle is read-only data:
#: writing ``kmers/`` into it would ship 8 MiB of packed arrays inside the wheel
#: (and made the source tree fail the bundled-size budget test).
_SCRATCH: dict[Path, Path] = {}


def _substrate(db_dir: Path) -> Path:
    """The directory that actually carries the packed sets.

    An installed database has them already. A source bundle does not (the bundle
    ships only ``kmer.json``, so ``db setup`` computes them at install time), so
    materialize a scratch copy instead of mutating the source.
    """
    if (db_dir / "kmers").is_dir():
        return db_dir
    cached = _SCRATCH.get(db_dir)
    if cached is not None:
        return cached
    scratch = Path(tempfile.mkdtemp(prefix="gspec-kmers-"))
    shutil.copyfile(db_dir / "kmer.json", scratch / "kmer.json")
    shutil.copyfile(db_dir / "sequences", scratch / "sequences")
    document = read_kmerdb(scratch)
    if document is None:
        raise SystemExit(f"{db_dir}: no kmer.json (run import_upstream.py --kmerdb)")
    materialize_kmer_sets(scratch, document, scratch / "sequences")
    _SCRATCH[db_dir] = scratch
    return scratch


def _ensure_sets(db_dir: Path) -> int:
    """Materialize the declared sets (once) and return ``k``."""
    document = read_kmerdb(_substrate(db_dir))
    if document is None:
        raise SystemExit(f"{db_dir}: no kmer.json (run import_upstream.py --kmerdb)")
    return document.k


def call(db_dir: Path, fastq: list[Path], *, min_fraction: float = MIN_FRACTION) -> dict[str, str]:
    """Best antigen per family, from the reads' k-mer containment."""
    k = _ensure_sets(db_dir)
    substrate = _substrate(db_dir)
    query = count_fastq_kmers(*fastq, k=k)
    families: dict[str, list[tuple[float, str]]] = {}
    document = read_kmerdb(substrate)
    assert document is not None
    for name in document.sets:
        family = name.split("_", 1)[0]
        fraction = set_hit_fraction(query, load_kmer_set(substrate, name))
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
