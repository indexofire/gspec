#!/usr/bin/env python3
"""Build differential k-mer signatures for the six S. enterica subspecies.

Subspecies is a whole-genome property — antigen genes cannot decide it. Both
SISTR (cgMLST) and SeqSero2/SalmID (clade-specific k-mers) classify against
labelled reference genomes; this builder follows SalmID's k-mer approach with
public-domain NCBI references only (no GPL content).

For each subspecies a set of DIFFERENTIAL k-mers is computed: every canonical
k-mer of its reference genome minus every k-mer of the other five references.
A query genome is then scored by how many of each subspecies' differential
k-mers it carries (recovery fraction); the highest wins. ``--max-per-set``
caps the stored signature size.

References (one per subspecies):

    I   enterica    NC_003197.2      (LT2, RefSeq)
    II  salamae     GCA_062082130.1
    IIIa arizonae   GCA_062081970.1
    IIIb diarizonae GCA_038592885.1
    IV  houtenae    GCA_061231905.1
    VI  indica      GCA_047114675.1

Usage::

    python scripts/salmonella/build_subspecies_kmers.py \
        --out src/gspec/data/dbs/sen_sero_sistr/subspecies_kmers.json

The default ``--max-per-set`` keeps the shipped sidecar small (10 000
signature k-mers per subspecies still classify the in-house 100-genome corpus
at 100 % accuracy).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from gspec.cluster.kmerdb import iter_canonical_kmer_codes

K = 31

REFERENCES = {
    "enterica": "NC_003197.2",
    "salamae": "GCA_062082130.1",
    "arizonae": "GCA_062081970.1",
    "diarizonae": "GCA_038592885.1",
    "houtenae": "GCA_061231905.1",
    "indica": "GCA_047114675.1",
}

COMP = str.maketrans("ACGT", "TGCA")


def _revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def _canonical_kmer(seq: str) -> str:
    rc = _revcomp(seq)
    return seq if seq <= rc else rc


def _kmer_hashes(sequence: str, k: int = K) -> set[int]:
    """Canonical k-mers of one sequence, as gspec's 2-bit packed codes."""
    return set(iter_canonical_kmer_codes(sequence, k))


def _download_genomes(tmp: Path) -> dict[str, str]:
    """``subspecies -> genome sequence`` for all six references."""
    genomes: dict[str, str] = {}
    # LT2 via efetch (RefSeq complete genome).
    fasta = (
        urllib.request.urlopen(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=nuccore&id={REFERENCES['enterica']}&rettype=fasta&retmode=text",
            timeout=180,
        )
        .read()
        .decode()
    )
    genomes["enterica"] = "".join(fasta.splitlines()[1:]).upper()

    for ssp, accession in REFERENCES.items():
        if ssp == "enterica":
            continue
        subprocess.run(
            [
                "datasets",
                "download",
                "genome",
                "accession",
                accession,
                "--include",
                "genome",
                "--filename",
                str(tmp / f"{ssp}.zip"),
            ],
            check=True,
            capture_output=True,
        )
        subprocess.run(
            ["unzip", "-o", "-q", str(tmp / f"{ssp}.zip"), "-d", str(tmp / ssp)],
            check=True,
        )
        fasta_path = next((tmp / ssp).rglob("*genomic.fna"))
        genomes[ssp] = "".join(
            line.strip() for line in fasta_path.read_text().splitlines() if not line.startswith(">")
        ).upper()
    return genomes


def _classify(query_kmers: set[int], signatures: dict[str, set[int]]) -> list[tuple[str, float]]:
    """Recovery fraction per subspecies, best first."""
    scored: list[tuple[str, float]] = []
    for ssp, kmer_set in signatures.items():
        if not kmer_set:
            scored.append((ssp, 0.0))
            continue
        hits = sum(1 for kmer in kmer_set if kmer in query_kmers)
        scored.append((ssp, hits / len(kmer_set)))
    return sorted(scored, key=lambda item: item[1], reverse=True)


def build(out: Path, *, max_per_set: int = 200_000) -> None:
    with TemporaryDirectory() as tmpdir:
        genomes = _download_genomes(Path(tmpdir))
        all_sets = {ssp: _kmer_hashes(seq) for ssp, seq in genomes.items()}
        signatures: dict[str, set[int]] = {}
        for ssp, own in all_sets.items():
            others = set().union(*(s for name, s in all_sets.items() if name != ssp))
            differential = own - others
            if len(differential) > max_per_set:
                differential = set(sorted(differential)[:max_per_set])
            signatures[ssp] = differential

    document = {
        "schema": "gspec.subspecies_kmers/1",
        "name": "sen_subspecies",
        "k": K,
        "description": (
            "Differential k-mer signatures for the six S. enterica subspecies,"
            " computed from public-domain NCBI reference genomes (one per"
            " subspecies). Recovery fraction = the query's share of each"
            " subspecies' signature; the highest wins."
        ),
        "license": "public domain (NCBI RefSeq/GenBank; derived by gspec)",
        "built": datetime.now(UTC).strftime("%Y-%m-%d"),
        "references": REFERENCES,
        "subspecies": {
            ssp: {
                "n": len(signatures[ssp]),
                # Packed 2-bit k-mers as hex strings, kept compact.
                "kmers_hex": [f"{kmer:0{K * 2 // 4}x}" for kmer in sorted(signatures[ssp])],
            }
            for ssp in REFERENCES
        },
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(f"{out.name}: k={K}," + " ".join(f"{ssp}={len(signatures[ssp])}" for ssp in REFERENCES))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--max-per-set", type=int, default=10_000)
    args = parser.parse_args(argv)
    build(args.out, max_per_set=args.max_per_set)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
