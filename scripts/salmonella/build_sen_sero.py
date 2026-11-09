#!/usr/bin/env python3
"""Build the bundled ``sen_sero`` panel from annotated INSDC CDS records.

Curation, not a serotype resolver: every marker is fetched by an annotated
identifier (``gene=`` or ``locus_tag=``) from a named reference genome, written
as a gspec/v1 record, and paired with the serovar rules in
``sen_sero_serovars.tsv``. Scope is *Salmonella enterica* subsp. *enterica*
only — the other subspecies are not distinguishable by antigen genes (see
README.md).

Usage::

    python scripts/salmonella/build_sen_sero.py --dest src/gspec/data/dbs/sen_sero

The genome exports are cached under ``scripts/salmonella/.cache/`` so re-runs
are offline and reproducible; delete the cache to re-fetch from NCBI.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from gspec.dbcodec import encode_seqid

HERE = Path(__file__).resolve().parent
CACHE = HERE / ".cache"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UA = "gspec-sen-sero-curation/1 (+https://github.com/indexofire/gspec)"
RULE_IDENTITY_FLOOR = 98.0
RULE_COVERAGE_FLOOR = 80.0


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        rows = [line for line in handle if not line.lstrip().startswith("#")]
    return [row for row in csv.DictReader(rows, delimiter="\t") if row.get(next(iter(row)))]


def _cds_export(genome: str) -> Path:
    """Cached ``fasta_cds_na`` export of one genome accession."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{genome}_cds.fna"
    if path.is_file():
        return path
    url = f"{EFETCH}?db=nuccore&id={genome}&rettype=fasta_cds_na&retmode=text"
    request = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(request, timeout=120) as response:
        text = response.read().decode("utf-8")
    if not text.startswith(">"):
        raise SystemExit(f"efetch returned no FASTA for {genome}: {text[:200]!r}")
    path.write_text(text, encoding="utf-8")
    return path


def _parse_cds(path: Path) -> list[dict[str, str]]:
    """One entry per CDS: header fields plus the nucleotide sequence."""
    rows: list[dict[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                rows.append({"header": header, "sequence": "".join(chunks)})
            header, chunks = line[1:], []
        else:
            chunks.append(line.strip())
    if header is not None:
        rows.append({"header": header, "sequence": "".join(chunks)})
    return rows


def _field(header: str, key: str) -> str:
    marker = f"[{key}="
    start = header.find(marker)
    if start < 0:
        return ""
    end = header.index("]", start)
    return header[start + len(marker) : end]


def select_cds(entries: list[dict[str, str]], selector: str, *, genome: str) -> dict[str, str]:
    """The CDS whose annotation matches ``gene=<symbol>`` or ``locus_tag=<tag>``."""
    key, _, value = selector.partition("=")
    if key not in {"gene", "locus_tag"}:
        raise SystemExit(f"selector must be gene=... or locus_tag=..., got {selector!r}")
    hits = [entry for entry in entries if _field(entry["header"], key) == value]
    if not hits:
        raise SystemExit(f"{genome}: no CDS with {selector} (re-curate the reference accession)")
    if len(hits) > 1:
        ids = ", ".join(_field(hit["header"], "protein_id") or hit["header"] for hit in hits)
        raise SystemExit(f"{genome}: {selector} is ambiguous ({ids}); use a locus_tag")
    return hits[0]


#: Phase-1 flagellin sits beside the filament capping/methylase genes; the
#: phase-2 (fljB) copy sits beside the hin invertase and its fljA repressor.
_PHASE_NEIGHBOURS = {
    "phase1": {"fliD", "fliB"},
    "phase2": {"hin", "fljA"},
}


def _start(entry: dict[str, str]) -> int:
    digits = re.findall(r"\d+", _field(entry["header"], "location"))
    return int(digits[0]) if digits else 0


def check_synteny(
    entries: list[dict[str, str]], cds: dict[str, str], phase: str, *, marker: str, genome: str
) -> None:
    """The curated ORF must be the flagellin the annotation places it beside.

    Two flagellins (phase 1 and phase 2) and a second flagellin-like ORF in the
    H2 region share ``FliC/FljB family flagellin`` products, so a product string
    alone cannot identify the allele: the neighbouring gene set can. This is the
    mechanical guard against curating the wrong ORF.
    """
    expected = _PHASE_NEIGHBOURS.get(phase)
    if expected is None:
        return
    position = _start(cds)
    neighbours = {
        _field(entry["header"], "gene") for entry in entries if abs(_start(entry) - position) < 8000
    }
    if not neighbours & expected:
        raise SystemExit(
            f"{marker} ({genome}): no {'/'.join(sorted(expected))} within 8 kb of the"
            f" selected CDS; found {sorted(n for n in neighbours if n)}"
        )


def _write_fasta(path: Path, records: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            header = encode_seqid(
                record["db"], record["gene"], record["accession"], (record["role"],)
            )
            # The annotation's product rides the FASTA description: a record
            # without one reports the database description for every gene.
            product = " ".join(record["product"].split())[:120]
            handle.write(f">{header} {product}\n")
            sequence = record["sequence"]
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset : offset + 60] + "\n")


def _serotype_document(serovars: list[dict[str, str]], name: str) -> dict[str, object]:
    rules: list[dict[str, object]] = []
    for row in serovars:
        rules.append(
            {
                "phenotype": row["serovar"],
                "requires": [marker for marker in row["requires"].split(";") if marker],
                "requires_any": [],
                "excludes": [marker for marker in row["excludes"].split(";") if marker],
                "max_coverage": {},
                "min_alignment_start_pct": {},
                "max_alignment_start_pct": {},
                "coverage_floor": RULE_COVERAGE_FLOOR,
                "identity_floor": RULE_IDENTITY_FLOOR,
                "control_exempt": False,
                "notes": [
                    f"formula {row['formula']}",
                    *[note for note in row.get("notes", "").split(";") if note],
                ],
            }
        )
    return {
        "schema": "gspec.serotype/1",
        "name": name,
        "control_gene": "invA",
        "control_coverage_floor": RULE_COVERAGE_FLOOR,
        "control_identity_floor": RULE_IDENTITY_FLOOR,
        "fallback": "Salmonella subsp. enterica (serovar not in panel)",
        "ambiguous": "No prediction (markers of several serovars)",
        "unresolved": "Salmonella subsp. enterica (no complete antigen set)",
        "rules": rules,
    }


def build(dest: Path, *, name: str = "sen_sero") -> None:
    markers = _read_tsv(HERE / "sen_sero_markers.tsv")
    serovars = _read_tsv(HERE / "sen_sero_serovars.tsv")
    if not markers or not serovars:
        raise SystemExit("marker/serovar tables are empty")

    exports = {genome: _parse_cds(_cds_export(genome)) for genome in {r["genome"] for r in markers}}
    records: list[dict[str, str]] = []
    marker_rows: list[dict[str, object]] = []
    for row in markers:
        entries = exports[row["genome"]]
        cds = select_cds(entries, row["selector"], genome=row["genome"])
        check_synteny(
            entries, cds, row.get("phase", "-"), marker=row["marker_id"], genome=row["genome"]
        )
        protein_id = _field(cds["header"], "protein_id") or _field(cds["header"], "locus_tag")
        product = _field(cds["header"], "protein") or row["note"]
        records.append(
            {
                "db": name,
                "gene": row["marker_id"],
                "accession": f"{row['genome']}:{protein_id}",
                "role": row["role"],
                "product": product,
                "sequence": cds["sequence"],
            }
        )
        marker_rows.append(
            {
                "marker": row["marker_id"],
                "role": row["role"],
                "serovar": row["serovar"],
                "source": {
                    "genome": row["genome"],
                    "protein_id": protein_id,
                    "product": product,
                    "phase_check": row.get("phase", "-"),
                },
                "length": len(cds["sequence"]),
                "note": row["note"],
            }
        )

    known = {record["gene"] for record in records}
    for row in serovars:
        referenced = {
            marker for marker in (row["requires"] + ";" + row["excludes"]).split(";") if marker
        }
        missing = sorted(referenced - known)
        if missing:
            raise SystemExit(f"serovar {row['serovar']!r} names unknown markers: {missing}")

    dest.mkdir(parents=True, exist_ok=True)
    _write_fasta(dest / "sequences", records)
    (dest / "serotype.json").write_text(
        json.dumps(_serotype_document(serovars, name), indent=2) + "\n", encoding="utf-8"
    )
    (dest / "screen.json").write_text(
        json.dumps(
            {
                "schema": "gspec.screen/1",
                "contigs": {"minid": 90.0, "mincov": 60.0, "merge_fragments": True},
                "reads": {"min_breadth": 60.0, "min_identity": 90.0},
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (dest / "markers.json").write_text(
        json.dumps(
            {
                "schema": "gspec.sen_markers/1",
                "source": {
                    "tool": "NCBI E-utilities efetch (fasta_cds_na)",
                    "license": "INSDC sequence data (no copyright restrictions)",
                    "genomes": sorted(exports),
                },
                "thresholds": {
                    "rule_identity_floor": RULE_IDENTITY_FLOOR,
                    "rule_coverage_floor": RULE_COVERAGE_FLOOR,
                    "basis": "scripts/salmonella/identity_matrix.py",
                },
                "markers": marker_rows,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (dest / "bundled.json").write_text(
        json.dumps(
            {
                "name": name,
                "description": (
                    "Salmonella enterica subsp. enterica serotyping panel"
                    " (gspec-curated INSDC CDS markers)"
                ),
                "vendor": "gspec-curated (INSDC annotated CDS)",
                "dbtype": "nucl",
                "kind": "gene",
                "snapshotted": datetime.now(UTC).strftime("%Y-%m-%d"),
                "license": "INSDC sequence data (public; no redistribution restriction)",
                "note": (
                    "O-antigen/H-antigen markers extracted from annotated RefSeq CDS"
                    " (NC_003197.2 Typhimurium LT2, NC_011294.1 Enteritidis P125109,"
                    " NC_006511.1 Paratyphi A ATCC 9150, NC_003198.1 Typhi CT18);"
                    " serotype.json applies the curated serovar rules with measured"
                    " thresholds (scripts/salmonella/README.md). Subsp. enterica only:"
                    " the other subspecies cannot be separated by antigen genes."
                ),
                "guidance": {
                    "references": [
                        {
                            "citation": (
                                "Zhang S, Yin Y, Jones MB, et al. Salmonella serotype"
                                " determination utilizing high-throughput genome"
                                " sequencing data. J Clin Microbiol. 2015;53(5):1685-1692."
                            ),
                            "doi": "10.1128/JCM.00323-15",
                            "note": "SeqSero: the O/H gene approach this panel follows.",
                        },
                        {
                            "citation": (
                                "Yoshida CE, Kruczkiewicz P, Laing CR, et al. The"
                                " Salmonella In Silico Typing Resource (SISTR): an open"
                                " web-accessible tool for rapidly typing and subtyping"
                                " draft microbial genomes. PLoS One."
                                " 2016;11(8):e0160528."
                            ),
                            "doi": "10.1371/journal.pone.0160528",
                            "note": "SISTR: serovar nomenclature and validation reference.",
                        },
                    ],
                    "caveats": [
                        (
                            "Gene panel only: O:2 (group A) and O:9 (group D1) are not"
                            " separated by presence markers (rfc 98.5-99.9% and rfbS"
                            " 98.5-99.3% identical across groups); the definitive O call"
                            " needs an rfb-locus cluster stage."
                        ),
                        (
                            "Paratyphi A is provisional: its rule excludes the D1 rfbE"
                            " allele, which that serovar may carry as a pseudogene."
                        ),
                        (
                            "Vi antigen (tviB/tviC) is not screened, so Typhi reports"
                            " 9,12:d:- rather than asserting Vi."
                        ),
                        (
                            "Only subsp. enterica is typed by design; the other"
                            " subspecies need biochemical differentiation."
                        ),
                    ],
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (dest / "info.json").write_text(
        json.dumps(
            {
                "schema": "gspec.dbinfo/1",
                "name": name,
                "description": (
                    "Salmonella subsp. enterica O/H antigen panel: 4 serovars"
                    " (Typhimurium, Enteritidis, Typhi, Paratyphi A) from curated"
                    " INSDC CDS markers."
                ),
                "vendor": "gspec-curated (INSDC annotated CDS)",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"{name}: {len(records)} markers, {len(serovars)} serovar rules -> {dest}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dest", type=Path, required=True)
    parser.add_argument("--name", default="sen_sero")
    args = parser.parse_args(argv)
    build(args.dest, name=args.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
