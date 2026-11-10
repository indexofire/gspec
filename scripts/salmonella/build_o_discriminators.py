#!/usr/bin/env python3
"""Build the O9/O2 discriminator pair from public-domain RefSeq annotations.

`wzx`/`wzy` cannot separate serogroup A (O:2) from D1 (O:9): the two groups share
the rfb machinery, and the actual difference sits in `rfbS` — CDP-paratose
2-epimerase, which carries tyvelose in D1 and paratose in A. Measured against the
references: the alleles are 99.29 % identical over the 840 bp gene (SeqSero2 uses
the 5' 609 bp, where they differ by a single SNP at 99.84 %).

So a discriminator exists, but it is a SNP-level one. This builder extracts it
from two NCBI RefSeq annotations (public domain):

* O:9 / D1 — ``NC_003198.1`` (*S.* Typhi CT18),
* O:2 / A  — ``NC_006511.1`` (*S.* Paratyphi A ATCC 9150).

Usage::

    python scripts/salmonella/build_o_discriminators.py \\
        --typhi ~/.cache/NC_003198.1_cds.fna \\
        --paratyphi-a ~/.cache/NC_006511.1_cds.fna \\
        --out src/gspec/data/dbs/sen_sero_sistr/o_discriminators.json \\
        --verify-against ~/genomes/          # optional: genomes to re-check

The output is a `gspec.o_discriminators/1` document: two reference sequences with
their provenance, plus the evidence floor below which the pair must NOT decide
(`min_gap`, in identity percentage points). A partial locus yields a gap of a few
tenths of a point — measured 0.3 on a 344 bp fragment — which is noise, not an
answer; the engine keeps the ambiguity in that case and says why.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from gspec.fasta import iter_fasta

#: The gene that separates O:2 from O:9, and its neighbours in the rfb cluster.
GENE = "rfbS"

#: Below this identity gap (percentage points) between the two alleles the pair
#: says nothing: at 344 bp the correct allele scored 100 % and the other 99.7 %,
#: a one-mismatch difference that arises as easily from a sequencing error as
#: from the true serotype.
MIN_GAP = 0.5


def _gene_records(path: Path) -> dict[str, str]:
    """``gene name -> sequence`` from a ``*_cds_from_genomic.fna`` export."""
    out: dict[str, str] = {}
    for record in iter_fasta(path):
        match = re.search(r"\[gene=([^\]]+)\]", record.description)
        if match is None:
            continue
        out.setdefault(match.group(1), record.sequence)
    return out


def _sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()[:16]


def _verify(sequence: str, genome: Path, label: str) -> float | None:
    """Best identity of the extracted gene against its own genome (a sanity gate)."""
    if not genome.is_file():
        return None
    with subprocess.Popen(
        [
            "blastn",
            "-task",
            "blastn",
            "-dust",
            "no",
            "-query",
            "-",
            "-subject",
            str(genome),
            "-outfmt",
            "6 pident length",
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
    ) as handle:
        assert handle.stdin is not None
        stdout, _ = handle.communicate(f">probe\n{sequence}\n")
    for line in stdout.splitlines():
        identity, _, length = line.partition("\t")
        if int(length) >= len(sequence) * 0.95:
            print(f"  {label}: {float(identity):.3f} % over {length} bp against {genome.name}")
            return float(identity)
    print(f"  {label}: no full-length hit against {genome.name}")
    return None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--typhi", type=Path, required=True, help="NC_003198.1 CDS export (O:9)")
    parser.add_argument(
        "--paratyphi-a", type=Path, required=True, help="NC_006511.1 CDS export (O:2)"
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--verify-against", type=Path, help="directory of assemblies to re-check")
    args = parser.parse_args(argv)

    pairs = {
        "D1": (
            "NC_003198.1",
            "Salmonella enterica subsp. enterica serovar Typhi str. CT18",
            args.typhi,
        ),
        "A": (
            "NC_006511.1",
            "Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC 9150",
            args.paratyphi_a,
        ),
    }
    references: list[dict[str, object]] = []
    for group, (accession, organism, path) in pairs.items():
        genes = _gene_records(path)
        if GENE not in genes:
            raise SystemExit(f"{path}: no [gene={GENE}] record")
        sequence = genes[GENE]
        genomes = (
            sorted(args.verify_against.glob(f"{accession}*.fna")) if args.verify_against else []
        )
        print(f"{group} <- {accession} {GENE} ({len(sequence)} bp)")
        identity = _verify(sequence, genomes[0], group) if genomes else None
        references.append(
            {
                "group": group,
                "gene": GENE,
                "accession": accession,
                "organism": organism,
                "source": f"NCBI RefSeq {accession} (CDS export, public domain)",
                "url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
                "length": len(sequence),
                "sha256_16": _sha256(sequence),
                "self_identity_pct": identity,
                "sequence": sequence,
            }
        )

    first, second = references
    document = {
        "schema": "gspec.o_discriminators/1",
        "name": "o9_o2_rfbS",
        "description": (
            "The A (O:2) versus D1 (O:9) discriminator: rfbS, whose paratose/"
            "tyvelose polymorphism is the only sequence difference between the two "
            "serogroups' rfb clusters. A SNP-level marker, so the pair only decides "
            "when the identity gap is decisive."
        ),
        "license": "public domain (NCBI RefSeq; extracted by gspec)",
        "extracted": "2026-11-10",
        "min_gap": MIN_GAP,
        "favours": {first["group"]: ["D1"], second["group"]: ["A"]},
        "references": references,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.out} ({args.out.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
