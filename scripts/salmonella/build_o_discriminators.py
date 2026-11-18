#!/usr/bin/env python3
"""Build the O serogroup discriminator document from public-domain NCBI records.

Some serogroups cannot be separated by their wzx/wzy references at all, so each
needs a dedicated SNP-level marker. The engine consults them in order:

* **A (O:2) vs D1 (O:9)** — ``rfbS`` (CDP-paratose 2-epimerase). The two alleles
  are 99.3 % identical over 840 bp; the paratose/tyvelose polymorphism is the
  only sequence difference between the two serogroups' rfb clusters.
* **D1 (O:9,12) vs D2 (O:9,46)** — ``wbaV`` (a glycosyltransferase family 2
  gene in the rfb cluster). Neither SISTR nor SeqSero2 ships a D1 ``wzy``
  reference, and the two ``wzx`` alleles sit 0.08 pp apart (measured on a
  Gateshead genome: D1 99.46 % vs D2 99.38 %), so ``wzx``/``wzy`` cannot decide
  the pair. ``wbaV`` differs by 7 SNPs over 1002 bp and that difference is
  conserved across 14 independent Enteritidis genomes (100 % pairwise among
  them) — a clean 0.70 pp gap.

Both are SNP-level markers, so each only decides when the identity gap is
decisive; a partial locus stays ambiguous and the engine says why.

A group may carry several alleles because its own diversity can rival the
between-group gap: Typhi's ``wbaV`` is 8 SNPs from Enteritidis' while D1 vs D2 is
7, so the D1 side ships both alleles and scores by the best of them.

Every sequence is public-domain NCBI data, extracted by stable identifier
(``[gene=...]`` or ``[protein_id=...]``) so a rebuild is reproducible:

===========  =========================  =====================  ==========
group        accession                  identifier             length
===========  =========================  =====================  ==========
D1 (rfbS)    NC_003198.1  (Typhi CT18)  gene=rfbS              840
A  (rfbS)    NC_006511.1  (Paratyphi A)  gene=rfbS              840
D1 (wbaV)    GCF_000783815.2 (Enterit.)  WP_000868492.1         1002
D1 (wbaV)    NC_003198.1  (Typhi CT18)  WP_010989202.1         1002
D2 (wbaV)    GCF_009730015.1 (Gateshead) WP_076726827.1         1002
===========  =========================  =====================  ==========

Usage::

    python scripts/salmonella/build_o_discriminators.py \
        --out src/gspec/data/dbs/sen_sero_sistr/o_discriminators.json
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

#: Below this identity gap (percentage points) a pair says nothing: a partial
#: locus shows a few tenths of a point, which is sequencing noise, not an answer.
MIN_GAP = 0.5

#: The A/D1 `rfbS` pair, by gene symbol (both RefSeq complete genomes).
RFBS_PAIR = {
    "D1": ("NC_003198.1", "Salmonella enterica subsp. enterica serovar Typhi str. CT18"),
    "A": ("NC_006511.1", "Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC 9150"),
}

#: The D1/D2 `wbaV` pair, by protein_id (the gene carries no symbol in these
#: annotations). The two D1 alleles bracket the group's own diversity.
WBAV_PAIR = {
    "D1": [
        ("GCF_000783815.2", "WP_000868492.1", "RefSeq Enteritidis assembly"),
        (
            "NC_003198.1",
            "WP_010989202.1",
            "Salmonella enterica subsp. enterica serovar Typhi str. CT18",
        ),
    ],
    "D2": [
        ("GCF_009730015.1", "WP_076726827.1", "RefSeq Gateshead assembly"),
    ],
}


def _sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()[:16]


def _efetch(accession: str, rettype: str) -> str:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={accession}&rettype={rettype}&retmode=text"
    )
    with urllib.request.urlopen(url, timeout=300) as response:
        return response.read().decode()


def _parse_fasta(text: str) -> list[tuple[str, str, str]]:
    """``(id, description, sequence)`` records from in-memory FASTA text."""
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


def _cds_records(accession: str, tmp: Path) -> list[tuple[str, str, str]]:
    """One accession's CDS records, from RefSeq efetch or an assembly's export."""
    if accession.startswith("NC_"):
        return _parse_fasta(_efetch(accession, "fasta_cds_na"))
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
    found = next((tmp / accession).rglob("cds_from_genomic.fna"))
    return _parse_fasta(found.read_text(encoding="utf-8"))


def _pick(
    records: list[tuple[str, str, str]], *, gene: str | None = None, protein_id: str | None = None
) -> str:
    """The one CDS sequence carrying ``gene`` or ``protein_id``."""
    needle = f"[gene={gene}]" if gene else f"[protein_id={protein_id}]"
    for _, description, sequence in records:
        if needle in description:
            return sequence
    raise SystemExit(f"no CDS record with {needle}")


def _rfbs(tmp: Path) -> dict[str, object]:
    references: list[dict[str, object]] = []
    favours: dict[str, list[str]] = {}
    for group, (accession, organism) in RFBS_PAIR.items():
        sequence = _pick(_cds_records(accession, tmp), gene="rfbS")
        gene_id = f"o_disc_rfbS_{group}"
        favours[group] = [gene_id]
        references.append(
            {
                "gene_id": gene_id,
                "group": group,
                "gene": "rfbS",
                "accession": accession,
                "organism": organism,
                "source": f"NCBI RefSeq {accession} (CDS export, public domain)",
                "url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
                "length": len(sequence),
                "sha256_16": _sha256(sequence),
                "sequence": sequence,
            }
        )
    return {
        "name": "A/D1 rfbS",
        "min_gap_pct": MIN_GAP,
        "min_coverage": 60.0,
        "note": (
            "A (O:2) versus D1 (O:9): rfbS, whose paratose/tyvelose polymorphism is"
            " the only sequence difference between the two serogroups' rfb clusters"
        ),
        "favours": favours,
        "references": references,
    }


def _wbav(tmp: Path) -> dict[str, object]:
    references: list[dict[str, object]] = []
    favours: dict[str, list[str]] = {}
    for group, entries in WBAV_PAIR.items():
        favours[group] = []
        for accession, protein_id, organism in entries:
            sequence = _pick(_cds_records(accession, tmp), protein_id=protein_id)
            gene_id = f"o_disc_wbaV_{group}__{len(favours[group]) + 1}"
            favours[group].append(gene_id)
            references.append(
                {
                    "gene_id": gene_id,
                    "group": group,
                    "gene": "wbaV",
                    "accession": accession,
                    "protein_id": protein_id,
                    "organism": organism,
                    "source": f"NCBI RefSeq {accession} (CDS export, public domain)",
                    "url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
                    "length": len(sequence),
                    "sha256_16": _sha256(sequence),
                    "sequence": sequence,
                }
            )
    return {
        "name": "D1/D2 wbaV",
        "min_gap_pct": MIN_GAP,
        "min_coverage": 60.0,
        "note": (
            "D1 (O:9,12) versus D2 (O:9,46): wbaV, 7 SNPs apart and conserved"
            " across Enteritidis genomes; no D1 wzy reference exists and the wzx"
            " alleles differ by only 0.08 pp"
        ),
        "favours": favours,
        "references": references,
    }


#: The E1 (O:3,10) / E4 (O:1,3,19) pair, on ``wbaP`` — the O-antigen
#: initiating transferase (undecaprenyl-phosphate galactose
#: phosphotransferase), whose allele decides the first sugar of the O repeat unit
#: and therefore the serogroup. wzx/wzy cannot do it: SISTR's E1 wzx class even
#: contains O:1,3,19 references (``407|584|1,3,19|E1``), so an E4 genome matches
#: "E1" wzx at 100 %, and the two wzy alleles are 99.66 % identical.
WBAP_PAIR = {
    "E1": [
        ("GCF_002053695.1", "WP_080097323.1", "RefSeq Give assembly"),
    ],
    "E4": [
        ("GCF_002044555.1", "WP_000368702.1", "RefSeq Senftenberg assembly"),
    ],
}


def _wbap(tmp: Path) -> dict[str, object]:
    references: list[dict[str, object]] = []
    favours: dict[str, list[str]] = {}
    for group, entries in WBAP_PAIR.items():
        favours[group] = []
        for accession, protein_id, organism in entries:
            sequence = _pick(_cds_records(accession, tmp), protein_id=protein_id)
            gene_id = f"o_disc_wbaP_{group}__{len(favours[group]) + 1}"
            favours[group].append(gene_id)
            references.append(
                {
                    "gene_id": gene_id,
                    "group": group,
                    "gene": "wbaP",
                    "accession": accession,
                    "protein_id": protein_id,
                    "organism": organism,
                    "source": f"NCBI RefSeq {accession} (CDS export, public domain)",
                    "url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
                    "length": len(sequence),
                    "sha256_16": _sha256(sequence),
                    "sequence": sequence,
                }
            )
    return {
        "name": "E1/E4 wbaP",
        "min_gap_pct": MIN_GAP,
        "min_coverage": 60.0,
        "note": (
            "E1 (O:3,10) versus E4 (O:1,3,19): wbaP, the O-antigen initiating"
            " transferase, 1.22 pp apart; wzx/wzy cannot decide because SISTR's E1"
            " wzx class contains O:1,3,19 references"
        ),
        "favours": favours,
        "references": references,
    }


def build(out: Path) -> None:
    with TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        pairs = [_rfbs(tmp), _wbav(tmp), _wbap(tmp)]
    document = {
        "schema": "gspec.o_discriminators/1",
        "name": "sen_o_discriminators",
        "description": (
            "SNP-level O-serogroup discriminators for pairs whose wzx/wzy"
            " references cannot separate them: A/D1 on rfbS, D1/D2 on wbaV and"
            " E1/E4 on wbaP."
        ),
        "license": "public domain (NCBI RefSeq; extracted by gspec)",
        "extracted": datetime.now(UTC).strftime("%Y-%m-%d"),
        "min_gap": MIN_GAP,
        "pairs": pairs,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    total = sum(len(pair["references"]) for pair in pairs)
    print(f"{out.name}: {len(pairs)} pair(s), {total} reference sequence(s) -> {out}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    build(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
