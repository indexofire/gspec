#!/usr/bin/env python3
"""Build the O:22/O:23 sub-factor document from public-domain NCBI records.

Serogroup G (O:13) references cannot separate factor 22 from factor 23: the two
wzx alleles differ by noise (measured on GCF_001569625.1: wzx scored 13,22 while
the dedicated O:23 genes and cgMLST point to 13,23). SeqSero2 resolves the pair
with O:22/O:23-specific genes; the O:23 side is a small sugar-biosynthesis block
— ``galE`` (UDP-glucose 4-epimerase) and a glycosyltransferase family 4 gene —
that this builder extracts from the NCBI RefSeq Cubana assembly
(``GCF_003430035.1``, O:13,23), so the sequence carries no GPL-2.0 content.

The O:22 side is not yet curated: the O:22-specific genes must come from a
verified O:13,22 reference, and the public "Diguel" assembly checked during
development (``GCA_007536395.1``) carries the O:23 block instead (its wzx gap to
13,22 is 0.32 pp, below the 0.5 pp floor, so it is ambiguous at best). Until a
curated 13,22 reference lands, the document is ONE-SIDED: presence of the O:23
genes proves O:23, their absence is read as O:22.

It also extracts the post-decision special gene ``sdf`` (Enteritidis-specific
Sdf I region, GenBank ``AF370707``). ``oafA`` (Typhimurium O5-) and
``gntR``/STM3356 (Paratyphi B d-tartrate) are natural VARIANT probes (a deletion
junction and a SNP, respectively), not standard NCBI CDS records, so they are
recorded here as provenance only and wired into ``special_checks`` when their
variant reference genomes are curated.

Usage::

    python scripts/salmonella/build_o_subfactor.py \
        --out src/gspec/data/dbs/sen_sero/o_subfactor.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from gspec.fasta import iter_fasta

#: NCBI CDS records that prove the O:23 factor (Cubana, GCF_003430035.1, O:13,23).
O23_GENES = {
    "o_disc_O23_galE": {
        "protein_id": "WP_001637675.1",
        "product": "UDP-glucose 4-epimerase GalE",
        "gene": "galE",
        "genome": "GCF_003430035.1",
    },
    "o_disc_O23_gtr": {
        "protein_id": "WP_020898978.1",
        "product": "glycosyltransferase family 4 protein",
        "gene": None,
        "genome": "GCF_003430035.1",
    },
}

#: NCBI CDS records that prove the O:22 factor (Poona, GCF_053338755.1, O:13,22).
#: The two galE/glycosyltransferase blocks are DIFFERENT genes (~62 % identical),
#: so this is a presence/absence discriminator, not an allele comparison.
O22_GENES = {
    "o_disc_O22_galE": {
        "protein_id": "WP_023214139.1",
        "product": "UDP-glucose 4-epimerase GalE",
        "gene": "galE",
        "genome": "GCF_053338755.1",
    },
    "o_disc_O22_gtr": {
        "protein_id": "WP_024154990.1",
        "product": "glycosyltransferase",
        "gene": None,
        "genome": "GCF_053338755.1",
    },
}

#: Enteritidis-specific Sdf I region (GenBank, public domain).
SDF_ACCESSION = "AF370707"

#: LT2 (NC_003197.2, public domain) carries the reference alleles for the two
#: clinically common single-locus variants: ``oafA`` (Typhimurium O5-) and the
#: ``gntR``/STM3356 d-tartrate SNP region (Paratyphi B dt+/-).
LT2_ACCESSION = "NC_003197.2"
LT2_OAF_A = "NP_461175.1"
#: The 291 bp d-tartrate SNP probe region on LT2 (1-based, forward strand);
#: SeqSero2's dt+ probe is the reverse complement of this window.
LT2_GNT_R_DT_REGION = (3525030, 3525320)
#: The single dt+ -> dt- SNP: probe position 251, G -> A.
GNT_R_DT_SNP = (251, "G", "A")

#: The enterica O:22/O:23 serovar pairs (a WHO-2007 fact), in (O:22, O:23) order.
ENTERICA_PAIRS = [
    ("Ibadan", "Mississippi"),
    ("Bahati", "Durham"),
    ("Oudwijk", "Bracknell"),
    ("Vaertan", "Ullevi"),
    ("Friedenau", "Wichita"),
    ("Diguel", "Telelkebir"),
    ("Newyork", "Okatie"),
    ("Bron", "Agbeni"),
    ("Washington", "Kintambo"),
    ("Poona", "Farmsen"),
    ("Agoueve", "Cubana"),
    ("Nimes", "Anna"),
    ("Mampong", "Picpus"),
    ("Leiden", "Fanti"),
    ("Ried", "Ajiobo"),
]

MIN_GAP = 0.5


def _sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()[:16]


def _efetch(accession: str, rettype: str = "fasta_cds_na") -> str:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={accession}&rettype={rettype}&retmode=text"
    )
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read().decode()


def _download_cds(tmp: Path, accession: str) -> dict[str, tuple[str, str]]:
    """Download one assembly's CDS and return ``protein_id -> (gene, sequence)``."""
    subprocess.run(
        [
            "datasets",
            "download",
            "genome",
            "accession",
            accession,
            "--include",
            "cds",
            "--filename",
            str(tmp / f"{accession}.zip"),
        ],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        ["unzip", "-o", "-q", str(tmp / f"{accession}.zip"), "-d", str(tmp / accession)],
        check=True,
    )
    fasta = next((tmp / accession).rglob("cds_from_genomic.fna"))
    out: dict[str, tuple[str, str]] = {}
    for record in iter_fasta(fasta):
        protein = record.description
        protein_id = (
            protein.split("[protein_id=")[-1].split("]")[0] if "[protein_id=" in protein else ""
        )
        gene = protein.split("[gene=")[-1].split("]")[0] if "[gene=" in protein else ""
        out[protein_id] = (gene, record.sequence)
    return out


def _reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return sequence.translate(table)[::-1]


def _iter_fasta_text(text: str) -> list[tuple[str, str, str]]:
    """``(id, description, sequence)`` records from an in-memory FASTA text."""
    records: list[tuple[str, str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                seqid, _, description = header.partition(" ")
                records.append((seqid, description, "".join(chunks)))
            header, chunks = line[1:].strip(), []
        else:
            chunks.append(line.strip().upper())
    if header is not None:
        seqid, _, description = header.partition(" ")
        records.append((seqid, description, "".join(chunks)))
    return records


def _lt2_oaf_a() -> dict[str, str]:
    """Full ``oafA`` from LT2 (NCBI, public domain) as a special marker."""
    text = _efetch(LT2_ACCESSION, rettype="fasta_cds_na")
    for _, description, sequence in _iter_fasta_text(text):
        if f"[protein_id={LT2_OAF_A}]" in description:
            return {
                "gene": "special_oafA",
                "protein_id": LT2_OAF_A,
                "product": "acetylation of the O-antigen (LPS)",
                "gene_name": "oafA",
                "genome": LT2_ACCESSION,
                "length": str(len(sequence)),
                "sha256_16": _sha256(sequence),
                "sequence": sequence,
            }
    raise SystemExit(f"{LT2_OAF_A} (oafA) not found in {LT2_ACCESSION} CDS")


def _lt2_gnt_r_dt_probes() -> tuple[dict[str, str], dict[str, str]]:
    """The 291 bp d-tartrate SNP probes, from the LT2 genome (public domain).

    The dt+ probe is the reverse complement of ``LT2_GNT_R_DT_REGION``; the dt-
    probe is the same sequence with the single documented SNP applied. Both are
    natural alleles of one locus, so this is an allele comparison, not a
    presence/absence test.
    """
    fasta = _efetch(LT2_ACCESSION, rettype="fasta")
    genome = "".join(fasta.splitlines()[1:]).upper()
    start, end = LT2_GNT_R_DT_REGION
    region = genome[start - 1 : end]
    plus = _reverse_complement(region)
    pos, _, _ = GNT_R_DT_SNP
    # The dt+ and dt- alleles differ by one SNP at ``pos``; flip the LT2 base
    # to the other allele for the dt- probe (auto-detecting the base keeps this
    # correct even if LT2 carries either allele at the site).
    base = plus[pos - 1]
    other = "A" if base == "G" else "G"
    minus = plus[: pos - 1] + other + plus[pos:]
    return (
        {
            "gene": "special_gntR_dtplus",
            "protein_id": "",
            "product": "STM3356 (gntR) d-tartrate SNP region, dt+ allele",
            "gene_name": "gntR",
            "genome": LT2_ACCESSION,
            "length": str(len(plus)),
            "sha256_16": _sha256(plus),
            "sequence": plus,
        },
        {
            "gene": "special_gntR_dtminus",
            "protein_id": "",
            "product": "STM3356 (gntR) d-tartrate SNP region, dt- allele",
            "gene_name": "gntR",
            "genome": LT2_ACCESSION,
            "length": str(len(minus)),
            "sha256_16": _sha256(minus),
            "sequence": minus,
        },
    )


def build(out: Path) -> None:
    references: list[dict[str, str]] = []
    o22_genes: list[str] = []
    o23_genes: list[str] = []

    with TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        cds = {
            accession: _download_cds(tmp, accession)
            for accession in {
                *(spec["genome"] for spec in O23_GENES.values()),
                *(spec["genome"] for spec in O22_GENES.values()),
            }
        }
        for genes, gene_list in ((O23_GENES, o23_genes), (O22_GENES, o22_genes)):
            for gene_id, spec in genes.items():
                gene_name, sequence = cds[spec["genome"]][spec["protein_id"]]
                gene_list.append(gene_id)
                references.append(
                    {
                        "gene": gene_id,
                        "protein_id": spec["protein_id"],
                        "product": spec["product"],
                        "gene_name": gene_name,
                        "genome": spec["genome"],
                        "length": str(len(sequence)),
                        "sha256_16": _sha256(sequence),
                        "sequence": sequence,
                    }
                )

    sdf = "".join(_efetch(SDF_ACCESSION, rettype="fasta").splitlines()[1:]).upper()
    references.append(
        {
            "gene": "special_sdf",
            "protein_id": "",
            "product": "Sdf I region (Enteritidis-specific)",
            "gene_name": "sdf-I",
            "genome": SDF_ACCESSION,
            "length": str(len(sdf)),
            "sha256_16": _sha256(sdf),
            "sequence": sdf,
        }
    )

    oaf_a = _lt2_oaf_a()
    references.append(oaf_a)
    gnt_r_plus, gnt_r_minus = _lt2_gnt_r_dt_probes()
    references.extend([gnt_r_plus, gnt_r_minus])

    document = {
        "schema": "gspec.o_subfactor/1",
        "name": "sen_sero_o_subfactor",
        "description": (
            "O:22/O:23 sub-factor discriminator (presence/absence of the"
            " O:22-specific and O:23-specific galE+glycosyltransferase blocks,"
            " which are ~62 % identical and therefore not alleles) and the"
            " Enteritidis Sdf I region, extracted from public-domain NCBI records."
        ),
        "license": "public domain (NCBI RefSeq/GenBank; extracted by gspec)",
        "extracted": datetime.now(UTC).strftime("%Y-%m-%d"),
        "min_gap": MIN_GAP,
        "o22_genes": o22_genes,
        "o23_genes": o23_genes,
        "pairs": ENTERICA_PAIRS,
        "references": references,
        "special_checks": [
            {
                "serovar": "Enteritidis",
                "markers": ["special_sdf"],
                "versus": [],
                "note": "Sdf I detected: characteristic of commonly circulating Enteritidis",
            },
            {
                "serovar": "Typhimurium",
                "markers": ["special_oafA"],
                "versus": [],
                "when": "absent",
                "note": "oafA absent/truncated: Typhimurium O5- variant (4,[5],12 antigen)",
            },
            {
                "serovar": "I 1,4,[5],12:i:-",
                "markers": ["special_oafA"],
                "versus": [],
                "when": "absent",
                "note": "oafA absent/truncated: monophasic Typhimurium O5- variant",
            },
            {
                "serovar": "Paratyphi B",
                "markers": ["special_gntR_dtminus"],
                "versus": ["special_gntR_dtplus"],
                "note": (
                    "STM3356 (gntR) d-tartrate non-fermenting allele (dt-):"
                    " Paratyphi B var. Java / typhoidal pathotype"
                ),
            },
            {
                "serovar": "Paratyphi B",
                "markers": ["special_gntR_dtplus"],
                "versus": ["special_gntR_dtminus"],
                "note": (
                    "STM3356 (gntR) d-tartrate fermenting allele (dt+):"
                    " Paratyphi B var. L(+) tartrate+"
                ),
            },
        ],
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(
        f"{out.name}: {len(o23_genes)} O:23 genes,"
        f" {len(document['pairs'])} enterica pairs,"
        f" {len(document['special_checks'])} special check(s),"
        f" {len(references)} reference sequence(s)"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    build(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
