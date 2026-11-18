#!/usr/bin/env python3
"""Build the per-accession serotype evidence table under ``~/data/sen``.

What this is, and what it is not
--------------------------------
There is **no public field that records "this genome was serotyped by
agglutination"**. INSDC records carry a submitter-asserted ``serovar``
attribute, and near-nothing about the method behind it. So a machine cannot
build a true gold-standard list from INSDC alone.

What this script does instead is collect **all 18,970 RefSeq S. enterica
assemblies** and record, per accession, every provenance signal that bears on
how much the label can be trusted:

* ``culture_collection`` — the strain is deposited in a public collection
  (NCTC / ATCC / CFSAN / LMG / CCUG / DSM). A serovar's **type or reference
  culture** is the closest thing to a gold standard that carries a genome.
* ``collected_by`` / ``submitter`` — a national reference lab (CDC, PulseNet,
  UKHSA, Institut Pasteur, RIVM …) implies routine agglutination serotyping.
* ``strain`` / ``isolate`` — a name worth chasing in the literature.
* everything else is only a submitter assertion.

The output is a table a human can filter, not a verdict. Rows are classified
into ``evidence_class`` so a downstream study can require the strong classes
only:

    reference_culture   culture_collection is populated (type/reference strain)
    reference_lab       collected_by names a national reference laboratory
    submitter_asserted  the serovar attribute alone

Usage::

    python scripts/salmonella/fetch_serotype_evidence.py --out ~/data/sen
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/28901/dataset_report"
PAGE_SIZE = 1000

#: Organisms whose submitter is a national reference laboratory: these serotype
#: routinely by agglutination, so the label has a documented process behind it.
REFERENCE_LAB_HINTS = (
    "centers for disease control",
    "pulsenet",
    "public health england",
    "ukhsa",
    "institut pasteur",
    "rivm",
    "national microbiology",
    "phac",
    "food and drug administration",
    "fsis",
    "anses",
    "iss ",
    "statens serum",
)

#: Public culture collections. A serovar's type/reference strain lives in one.
COLLECTION_HINTS = ("NCTC", "ATCC", "CFSAN", "LMG", "CCUG", "DSM", "CIP", "NBRC", "JCM")


def _attr(biosample: dict[str, object], name: str) -> str:
    """One BioSample attribute's value, ``""`` when absent."""
    for raw in biosample.get("attributes") or []:  # type: ignore[union-attr]
        if isinstance(raw, dict) and raw.get("name") == name:
            value = raw.get("value")
            return str(value).strip() if value is not None else ""
    return ""


def _evidence_class(collection: str, collector: str, submitter: str) -> str:
    if collection.strip():
        return "reference_culture"
    blob = f" {collector} {submitter} ".lower()
    if any(hint in blob for hint in REFERENCE_LAB_HINTS):
        return "reference_lab"
    return "submitter_asserted"


def fetch_all(out_dir: Path) -> Path:
    """Page the Datasets API and write one TSV row per RefSeq assembly."""
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "ncbi_salmonella_serovars.tsv"
    columns = [
        "accession",
        "serovar",
        "subspecies",
        "assembly_level",
        "evidence_class",
        "culture_collection",
        "strain",
        "isolate",
        "collected_by",
        "submitter",
        "bioproject",
        "organism_name",
    ]
    token: str | None = None
    written = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        while True:
            params = {"filters.assembly_source": "REFSEQ", "page_size": str(PAGE_SIZE)}
            if token:
                params["page_token"] = token
            url = f"{API}?{urllib.parse.urlencode(params)}"
            for attempt in range(4):
                try:
                    with urllib.request.urlopen(url, timeout=300) as response:
                        page = json.loads(response.read().decode())
                    break
                except (urllib.error.URLError, TimeoutError) as exc:  # pragma: no cover
                    if attempt == 3:
                        raise SystemExit(f"giving up on {url}: {exc}") from exc
                    time.sleep(5 * (attempt + 1))
            for report in page.get("reports", []):
                biosample = report.get("biosample") or {}
                organism = (report.get("organism") or {}).get("organism_name") or ""
                collection = _attr(biosample, "culture_collection")
                collector = _attr(biosample, "collected_by")
                submitter = (report.get("assembly_info") or {}).get("submitter") or ""
                writer.writerow(
                    [
                        report.get("accession", ""),
                        _attr(biosample, "serovar"),
                        _attr(biosample, "subspecies"),
                        (report.get("assembly_info") or {}).get("assembly_level", ""),
                        _evidence_class(collection, collector, submitter),
                        collection,
                        _attr(biosample, "strain"),
                        _attr(biosample, "isolate"),
                        collector,
                        submitter,
                        (report.get("assembly_info") or {}).get("bioproject_accession", ""),
                        organism,
                    ]
                )
                written += 1
            token = page.get("next_page_token")
            print(f"  fetched {written} rows…", file=sys.stderr)
            if not token:
                break
    print(f"{path}: {written} assemblies", file=sys.stderr)
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path.home() / "data" / "sen")
    args = parser.parse_args(argv)
    fetch_all(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
