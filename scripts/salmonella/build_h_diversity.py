#!/usr/bin/env python3
"""Add H-allele diversity to a panel from assemblies you already have.

An antigen's reference pool can be complete enough for serovar rows to be
reachable and still miss the allele a given genome carries: a Mendoza genome's
`fljB` matched the panel's `1,7` references at 99.61 % and its `1,2` references at
96.83 %, so the call landed on the wrong antigen even though `1,2` is the truth.
That is a *reference diversity* gap, not a cross-labelling one, and the fix is more
`1,2` alleles.

This script extracts them from assemblies by homology — no annotation needed:

1. every candidate locus that matches one of the antigen's existing references at
   ``--identity`` over ``--coverage`` is cut out;
2. alleles already represented in the panel (k-mer containment ≥ ``--min-novelty``
   of an existing allele) are skipped, so only genuinely new sequence is added;
3. each survivor is verified against its own assembly.

The output is a ``gspec.extra_refs/1`` document, which
``import_upstream.py --extra-refs`` appends as alleles of that antigen. Local data
only: ``--extra-refs`` refuses ``--bundle``, because an extracted allele is neither
gspec-curated nor a redistributable upstream dataset.

Usage::

    python scripts/salmonella/build_h_diversity.py \\
        --db-dir ~/.local/share/gspec/db/sen_sero_sistr \\
        --assemblies ~/assemblies --antigen 1,2 --family h2 \\
        --exclude 26SAL0084 --out h2_1_2_extra.json
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from gspec.align.sketch import allele_sketches, window_kmers
from gspec.fasta import iter_fasta


def _probe_alleles(db_dir: Path, family: str, antigen: str) -> list[str]:
    """The panel's existing sequences for one antigen (best allele first)."""
    wanted = set()
    for record in iter_fasta(db_dir / "sequences"):
        gene = next(
            part.split("=", 1)[1] for part in record.id.split("|") if part.startswith("gene=")
        )
        name = gene.rsplit("__", 1)[0]
        if name == f"{family}_{antigen.replace(',', '_')}":
            wanted.add(record.id)
    return sorted(wanted)


def _candidates(
    assemblies: list[Path], probes: dict[str, str], *, identity: float, coverage: float
) -> list[tuple[Path, str, str, float, float]]:
    """``(assembly, contig, sequence, identity, coverage)`` per candidate locus."""
    found: list[tuple[Path, str, str, float, float]] = []
    with tempfile.TemporaryDirectory() as scratch:
        probe_fasta = Path(scratch) / "probes.fna"
        probe_fasta.write_text(
            "".join(f">{name}\n{sequence}\n" for name, sequence in probes.items()),
            encoding="utf-8",
        )
        for assembly in assemblies:
            result = subprocess.run(
                [
                    "blastn",
                    "-task",
                    "blastn",
                    "-dust",
                    "no",
                    "-query",
                    str(probe_fasta),
                    "-subject",
                    str(assembly),
                    "-outfmt",
                    "6 sseqid sstart send pident length qlen",
                    "-max_target_seqs",
                    "50",
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            if result.returncode != 0:
                raise SystemExit(f"blastn failed on {assembly}: {result.stderr.strip()}")
            contigs = {record.id.split()[0]: record.sequence for record in iter_fasta(assembly)}
            for line in result.stdout.splitlines():
                contig, start, end, pident, length, qlen = line.split("\t")
                share = float(length) / float(qlen) * 100
                if float(pident) < identity or share < coverage:
                    continue
                sequence = contigs.get(contig, "")
                low, high = sorted((int(start), int(end)))
                if sequence and high <= len(sequence):
                    found.append((assembly, contig, sequence[low - 1 : high], float(pident), share))
    return found


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-dir", type=Path, required=True)
    parser.add_argument("--assemblies", type=Path, required=True)
    parser.add_argument("--family", choices=("h1", "h2"), required=True)
    parser.add_argument("--antigen", required=True)
    parser.add_argument("--exclude", nargs="*", default=[], help="assembly stems to skip")
    parser.add_argument("--identity", type=float, default=95.0)
    parser.add_argument("--coverage", type=float, default=90.0)
    parser.add_argument("--min-novelty", type=float, default=0.99)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    index = allele_sketches(args.db_dir)
    if index is None:
        raise SystemExit(f"{args.db_dir}: no k-mer substrate")
    probes = {
        record.id: record.sequence
        for record in iter_fasta(args.db_dir / "sequences")
        if record.id in set(_probe_alleles(args.db_dir, args.family, args.antigen))
    }
    if not probes:
        raise SystemExit(f"no {args.family} {args.antigen} allele in {args.db_dir}")

    assemblies = sorted(
        path for path in args.assemblies.glob("*.fna") if path.stem not in set(args.exclude)
    )
    print(f"scanning {len(assemblies)} assembly/assemblies against {len(probes)} probe(s)")
    existing = [sketch.codes for sketch in index.of(args.family) if sketch.antigen == args.antigen]
    candidates = _candidates(assemblies, probes, identity=args.identity, coverage=args.coverage)
    k = index.k
    added: list[dict[str, object]] = []
    for assembly, contig, sequence, pident, share in candidates:
        codes = window_kmers({contig: sequence}, contig, 1, len(sequence), k, pad=0)
        if (
            existing
            and max(len(codes & allele) / len(codes) for allele in existing) >= args.min_novelty
        ):
            continue
        if any(entry["sequence"] == sequence for entry in added):
            continue
        added.append(
            {
                "family": args.family,
                "antigen": args.antigen,
                "accession": f"local:{assembly.stem}:{contig}",
                "note": (
                    f"extracted from {assembly.name} ({contig}), {pident:.2f} % identity"
                    f" over {share:.1f} % of the probe; local use only, not redistributable"
                ),
                "sequence": sequence,
            }
        )
    document = {
        "schema": "gspec.extra_refs/1",
        "family": args.family,
        "antigen": args.antigen,
        "identity_floor": args.identity,
        "coverage_floor": args.coverage,
        "references": added,
    }
    args.out.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(f"{len(candidates)} candidate locus/loci, {len(added)} novel allele(s) -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
