#!/usr/bin/env python3
"""Mine published *Salmonella* papers for accession + serovar tables.

Why papers beat INSDC metadata
-----------------------------
A BioSample `serovar` attribute is a submitter assertion with no method attached.
A paper that reports its serotyping (agglutination against the Kauffmann-White
scheme) and publishes a strain table ties every accession to a **stated method**,
which is the provenance INSDC does not carry. This script harvests those tables.

How it works
------------
1. Search Europe PMC for open-access papers that are about Salmonella serotyping
   / serovar genomics and expose supplementary files.
2. Download each paper's `supplementaryFiles` archive (Europe PMC REST).
3. Parse every CSV/TSV/TXT/XLSX sheet in it and look for the two columns that
   matter: one holding accession-like tokens (`GCF_`/`GCA_`/`SAMN`/`SAMEA`/`ERS`
   / `SRR`/`ERR`/`DRR`) and one named like a serovar.
4. Emit ``accession <TAB> serovar <TAB> doi <TAB> pmcid <TAB> file <TAB> columns``
   so every row is traceable to the paper and the exact columns it came from.

XLSX is read with `zipfile` + `xml.etree` (no openpyxl in this environment).

Output: ``~/data/sen/published/``

    accessions_from_papers.tsv   one row per (accession, serovar) found
    papers.tsv                   the search result set with what was harvested
    README.md                    provenance and the caveat

Usage::

    python scripts/salmonella/harvest_published_serotypes.py --out ~/data/sen/published
    python scripts/salmonella/harvest_published_serotypes.py --serovar Hindmarsh
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
PAUSE = 0.4
USER_AGENT = "gspec-serotype-harvest/1 (+https://github.com/indexofire/gspec)"

ACCESSION_RE = re.compile(
    r"\b(GCF_\d{9,}(?:\.\d+)?|GCA_\d{9,}(?:\.\d+)?|SAMN\d{6,}|SAMEA\d{6,}|ERS\d{5,}|SRR\d{5,}|ERR\d{5,}|DRR\d{5,})\b"
)
SEROVAR_HEADER_RE = re.compile(r"sero(?:var|type|typ)", re.IGNORECASE)
#: Cells a table uses as a bucket, not a serovar. Kept out of the output: a row
#: saying "Other" is not evidence about any serovar.
_NOT_A_SEROVAR = {
    "",
    "-",
    "--",
    "n/a",
    "na",
    "none",
    "null",
    "other",
    "others",
    "unknown",
    "untypeable",
    "untyped",
    "not determined",
    "not typed",
    "not applicable",
    "pending",
    "total",
    "totals",
    "sum",
    "positive",
    "negative",
    "pos",
    "neg",
}
XLSX_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

#: A serovar value has to at least look like a name: letters, up to four words, no
#: pure codes, and none of the phrases a table uses instead of a name. Measured on
#: a 133-paper harvest, the unfiltered tail held ``/``, ``2a``, ``Yv`` and
#: "not stated explicity" — precision matters more than recall here because every
#: surviving row is meant to be checkable against the paper.
_SEROVAR_VALUE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9 ,.\[\]()'|/-]{1,40}$")


def _clean_serovar(value: str) -> str:
    """Strip the prefixes papers put in front of a serovar name."""
    text = value.strip()
    for prefix in ("Salmonella enterica ", "Salmonella ", "S. ", "S.", "serovar ", "ser. "):
        if text.lower().startswith(prefix.lower()):
            text = text[len(prefix) :].strip()
    return text


def _looks_like_serovar(value: str) -> bool:
    text = _clean_serovar(value)
    if ACCESSION_RE.search(text):
        # A table whose serovar column actually holds accessions (one transposed
        # sheet in the harvest did) must not turn into 3,690 fake serovars.
        return False
    if not text or text.lower() in _NOT_A_SEROVAR:
        return False
    if not _SEROVAR_VALUE_RE.match(text):
        return False
    if re.fullmatch(r"\d+[a-z]?", text):  # ``2a``, ``3b`` …
        return False
    if len(text.split()) > 4:  # a sentence, not a name
        return False
    lowered = text.lower()
    return not any(
        marker in lowered
        for marker in ("not stated", "not reported", "not available", "unknown", "explicity")
    )


#: Queries tried in order; each is an Europe PMC free-text expression. Open access
#: only, because the supplementary archive must be fetchable.
BASE_FILTER = "AND (OPEN_ACCESS:Y) AND HAS_SUPPL:Y"


def _get(url: str, *, binary: bool = False, attempt_limit: int = 5) -> bytes:
    for attempt in range(attempt_limit):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=180) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            # 404/403 mean "no supplement for this paper" — retrying cannot help.
            if exc.code in (403, 404):
                raise SystemExit(f"not available: {url} ({exc.code})") from exc
            if exc.code >= 500 and attempt >= 1:
                # A broken supplement endpoint stays broken; do not spend four
                # back-off sleeps on it.
                raise SystemExit(f"server error: {url} ({exc.code})") from exc
            if attempt == attempt_limit - 1:
                raise SystemExit(f"giving up on {url}: {exc}") from exc
            print(f"  retry {attempt + 1} (HTTP {exc.code})", file=sys.stderr)
            time.sleep(3 * (attempt + 1))
        except (urllib.error.URLError, TimeoutError, OSError) as exc:  # pragma: no cover
            if attempt == attempt_limit - 1:
                raise SystemExit(f"giving up on {url}: {exc}") from exc
            print(f"  retry {attempt + 1} ({type(exc).__name__})", file=sys.stderr)
            time.sleep(3 * (attempt + 1))
    raise AssertionError("unreachable")


def search(query: str, page_size: int = 25) -> list[dict[str, object]]:
    """One Europe PMC search page, core result type (carries doi/pmcid/hasSuppl)."""
    params = {
        "query": f"{query} {BASE_FILTER}",
        "format": "json",
        "pageSize": str(page_size),
        "resultType": "core",
    }
    time.sleep(PAUSE)
    payload = json.loads(_get(f"{EPMC}/search?{urllib.parse.urlencode(params)}").decode())
    hits = int(payload.get("hitCount") or 0)
    print(f"  query {query!r}: {hits} hits", file=sys.stderr)
    results = payload.get("resultList", {}).get("result", [])
    assert isinstance(results, list)
    return [r for r in results if isinstance(r, dict)]


def _xlsx_rows(blob: bytes) -> list[list[str]]:
    """Every row of the first worksheet of an .xlsx file, as strings."""
    try:
        archive = zipfile.ZipFile(io.BytesIO(blob))
    except zipfile.BadZipFile:
        return []
    shared: list[str] = []
    if "xl/sharedStrings.xml" in archive.namelist():
        root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        for item in root.iter(f"{XLSX_NS}si"):
            shared.append("".join(node.text or "" for node in item.iter(f"{XLSX_NS}t")))
    sheet = next(
        (name for name in archive.namelist() if re.fullmatch(r"xl/worksheets/sheet\d+\.xml", name)),
        None,
    )
    if sheet is None:
        return []
    rows: list[list[str]] = []
    root = ET.fromstring(archive.read(sheet))
    for row in root.iter(f"{XLSX_NS}row"):
        cells: dict[int, str] = {}
        for index, cell in enumerate(row.iter(f"{XLSX_NS}c")):
            column = cell.get("r") or ""
            letters = re.sub(r"\d", "", column)
            position = 0
            for letter in letters:
                position = position * 26 + (ord(letter.upper()) - 64)
            position = position - 1 if position else index
            kind = cell.get("t")
            if kind == "s":
                value_node = cell.find(f"{XLSX_NS}v")
                raw = value_node.text if value_node is not None else ""
                try:
                    value = shared[int(raw)]
                except (TypeError, ValueError, IndexError):
                    value = ""
            elif kind == "inlineStr":
                value = "".join(node.text or "" for node in cell.iter(f"{XLSX_NS}t"))
            else:
                value_node = cell.find(f"{XLSX_NS}v")
                value = value_node.text if value_node is not None else ""
            cells[position] = (value or "").strip()
        if cells:
            rows.append([cells.get(i, "") for i in range(max(cells) + 1)])
    return rows


def _delimited_rows(blob: bytes) -> list[list[str]]:
    """Rows of a CSV/TSV/text table, sniffing the delimiter."""
    for encoding in ("utf-8", "latin-1"):
        try:
            text = blob.decode(encoding)
            break
        except UnicodeDecodeError:
            continue
    else:  # pragma: no cover
        return []
    head = text[:4096]
    delimiter = max((",", "\t", ";"), key=head.count)
    if head.count(delimiter) < 2:
        return []
    try:
        return list(csv.reader(io.StringIO(text), delimiter=delimiter))
    except csv.Error:  # pragma: no cover
        return []


def _table_rows(name: str, blob: bytes) -> list[list[str]]:
    lowered = name.lower()
    if lowered.endswith((".xlsx", ".xlsm")):
        return _xlsx_rows(blob)
    if lowered.endswith((".csv", ".tsv", ".txt")):
        return _delimited_rows(blob)
    return []


def harvest_table(rows: list[list[str]]) -> tuple[list[tuple[str, str]], str]:
    """``([(accession, serovar)], "acc_col/serovar_col")`` for one table."""
    if not rows:
        return [], ""
    for header_index, header in enumerate(rows[:10]):
        serovar_columns = [
            index for index, cell in enumerate(header) if SEROVAR_HEADER_RE.search(cell or "")
        ]
        accession_columns = [
            index
            for index, cell in enumerate(header)
            if ACCESSION_RE.search(cell or "")
            or "accession" in (cell or "").lower()
            or "sample" in (cell or "").lower()
            or "strain" in (cell or "").lower()
        ]
        if not serovar_columns:
            continue
        found: list[tuple[str, str]] = []
        for row in rows[header_index + 1 :]:
            serovar = next(
                (
                    (row[i] or "").strip()
                    for i in serovar_columns
                    if i < len(row) and row[i].strip()
                ),
                "",
            )
            if not serovar or SEROVAR_HEADER_RE.search(serovar):
                continue
            if not _looks_like_serovar(serovar):
                continue
            serovar = _clean_serovar(serovar)
            accession = ""
            for index in accession_columns:
                if index >= len(row):
                    continue
                match = ACCESSION_RE.search(row[index] or "")
                if match:
                    accession = match.group(1)
                    break
            if accession:
                found.append((accession, serovar))
        if found:
            label = (
                f"{header[accession_columns[0]] if accession_columns else '?'}"
                f"/{header[serovar_columns[0]]}"
            )
            return found, label
    return [], ""


def harvest_paper(paper: dict[str, object], cache: Path) -> list[tuple[str, str, str]]:
    """``[(accession, serovar, "file:columns")]`` from one paper's supplement."""
    pmcid = str(paper.get("pmcid") or "")
    doi = str(paper.get("doi") or "")
    if not pmcid:
        return []
    archive_path = cache / f"{pmcid}.zip"
    if not archive_path.exists():
        try:
            time.sleep(PAUSE)
            archive_path.write_bytes(_get(f"{EPMC}/{pmcid}/supplementaryFiles"))
        except SystemExit:
            return []
    try:
        archive = zipfile.ZipFile(archive_path)
    except zipfile.BadZipFile:
        return []
    out: list[tuple[str, str, str]] = []
    for name in archive.namelist():
        if name.endswith("/"):
            continue
        try:
            blob = archive.read(name)
        except KeyError:  # pragma: no cover
            continue
        rows = _table_rows(name, blob)
        pairs, label = harvest_table(rows)
        for accession, serovar in pairs:
            out.append((accession, serovar, f"{name}:{label}"))
    if out:
        print(f"  {pmcid} ({doi}): {len(out)} rows", file=sys.stderr)
    return out


README = """# `published/` — serotype tables mined from papers

## Why this layer exists

A BioSample `serovar` attribute carries no method. A paper that says how the
strains were serotyped (agglutination against the Kauffmann-White scheme) and
publishes a strain table ties each accession to a **stated method** — the
provenance INSDC does not record. This directory holds those tables, each row
traceable to the paper and the exact columns it came from.

## Files

`accessions_from_papers.tsv`

| column | meaning |
| --- | --- |
| `accession` | `GCF_`/`GCA_` assembly, or `SAMN`/`SAMEA`/`ERS`/`SRR`/`ERR`/`DRR` sample or run |
| `serovar` | the value the paper's table gives |
| `doi` | paper DOI |
| `pmcid` | Europe PMC id (open access, so the supplement is public) |
| `source` | `supplementary-file:column-header` the row was read from |

`papers.tsv` — the search result set and how many rows each paper yielded.

## How to use it honestly

* A row means "this paper's table says accession X is serovar Y". **The paper's
  own methods section is the verification** — read it before treating a row as
  confirmed. This script cannot judge whether a paper did agglutination or just
  quoted a database.
* Accessions in papers are often **run/sample ids (`SRR`/`ERS`/`SAMN`)**, not
  assemblies: the genome may need assembling, or resolved to a `GCA_` first.
* Duplicate accessions across papers are left in place on purpose — agreement
  between two independent labs is itself evidence.

## What the first harvest found (2026-11-18)

133 open-access papers scanned, **14 with a machine-readable accession+serovar
table**, giving **4,058 (accession, serovar) rows over 106 serovars**. The
largest contributors are the O1/Paratyphi-A genomics study (PMC7530988, 1,444
rows), an Oxford-Nanopore evaluation (PMC13532270, 1,082) and the SeqSero2
validation set (PMC7147080, 94).

For the four serovars gspec cannot separate:

| serovar | rows from papers |
| --- | --- |
| Goldcoast | 52 |
| Bovismorbificans | 6 |
| Hindmarsh | **0** |
| Brikama | **0** |

So papers add real material for the O:6-bearing side but **nothing for the rare
twins** — those have to come from a laboratory, not from a publication either.

Two caveats the numbers force:

* **A table's `serovar` column is not always a serovar.** One harvested sheet had
  its columns transposed, so the serovar column held `SAMEA…` accessions — 3,690
  fake "serovars" before the guard that rejects accession-shaped values. Always
  spot-check `source` in the output.
* **Pipe-joined values (`Havana|Bron`, `Millesi|Overchurch|Tiko|Bukavu`) mean the
  paper itself reported several candidates.** They are kept verbatim: collapsing
  them would invent a decision the authors did not make, and `Yv`/`Xv`/`Ogawa`
  are *Salmonella* names the Kauffmann-White table does not list (Ogawa belongs to
  a different species group), i.e. a reminder that papers and the WHO table do not
  always agree.

## Refresh

```bash
python scripts/salmonella/harvest_published_serotypes.py --out ~/data/sen/published
python scripts/salmonella/harvest_published_serotypes.py --serovar Hindmarsh --serovar Brikama
```
"""


DEFAULT_QUERIES = (
    '"Salmonella enterica" AND serovar AND "whole genome sequencing" AND serotyping',
    '"Salmonella" AND serovar AND genomes AND "agglutination"',
    '"Salmonella" AND "Kauffmann-White" AND genome AND accession',
    '"Salmonella" AND "serovar" AND "supplementary table" AND genome AND accession',
    '"Salmonella" AND isolates AND serotyped AND "genome accession"',
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path.home() / "data" / "sen" / "published")
    parser.add_argument(
        "--serovar",
        action="append",
        default=[],
        help="also search this serovar by name (repeatable), e.g. --serovar Brikama",
    )
    parser.add_argument("--page-size", type=int, default=25)
    args = parser.parse_args(argv)

    out_dir: Path = args.out
    cache = out_dir / "supplements"
    cache.mkdir(parents=True, exist_ok=True)

    queries = list(DEFAULT_QUERIES)
    queries += [f'"Salmonella" AND "{name}" AND genome' for name in args.serovar]

    papers: dict[str, dict[str, object]] = {}
    for query in queries:
        for paper in search(query, args.page_size):
            pmcid = str(paper.get("pmcid") or "")
            if pmcid:
                papers.setdefault(pmcid, paper)

    print(f"{len(papers)} open-access papers with supplements", file=sys.stderr)
    rows: list[tuple[str, str, str, str, str]] = []
    paper_rows: list[list[str]] = []
    for index, (pmcid, paper) in enumerate(sorted(papers.items()), start=1):
        harvested = harvest_paper(paper, cache)
        for accession, serovar, source in harvested:
            rows.append((accession, serovar, str(paper.get("doi") or ""), pmcid, source))
        paper_rows.append(
            [
                pmcid,
                str(paper.get("doi") or ""),
                str(paper.get("pubYear") or ""),
                re.sub(r"<[^>]+>", "", str(paper.get("title") or "")),
                str(len(harvested)),
            ]
        )
        print(f"  [{index}/{len(papers)}] {pmcid}: {len(harvested)} rows", file=sys.stderr)

    accession_path = out_dir / "accessions_from_papers.tsv"
    with accession_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["accession", "serovar", "doi", "pmcid", "source"])
        writer.writerows(sorted(set(rows)))
    papers_path = out_dir / "papers.tsv"
    with papers_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["pmcid", "doi", "year", "title", "rows_harvested"])
        writer.writerows(sorted(paper_rows, key=lambda row: -int(row[4])))
    (out_dir / "README.md").write_text(README, encoding="utf-8")
    distinct = {row[1] for row in rows}
    print(
        f"{accession_path}: {len(set(rows))} (accession, serovar) rows"
        f" from {len(papers)} papers, {len(distinct)} distinct serovar values",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
