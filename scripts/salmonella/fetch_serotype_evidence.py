#!/usr/bin/env python3
"""Collect *Salmonella* serotype evidence into ``~/data/sen``.

Why "gold standard" needs a proxy
---------------------------------
**No INSDC field records that a genome was serotyped by agglutination.** The
`serovar` attribute is a submitter assertion. What the public record does carry
is the *provenance* around the strain, and two signals raise a label far above a
bare assertion:

1. **`culture_collection`** — the strain sits in a public collection (NCTC, ATCC,
   CFSAN, LMG, CCUG, CIP, DSM, NBRC, JCM). A serovar's **type or reference
   culture** originates from the WHO Collaborating Centre, and that is the
   closest thing to a gold standard that has a sequenced genome.
   Measured: **10,643** *S. enterica* BioSamples carry this attribute.
2. **the reference laboratory** in `collected_by` / `submitter` (CDC, PulseNet,
   UKHSA, Institut Pasteur, RIVM …) — such labs serotype by agglutination as a
   matter of routine.

Everything else is `submitter_asserted`. This script writes both layers so a
study can require the strong ones and never silently trust the rest.

Outputs (all TSV, UTF-8, no index column)
----------------------------------------
``reference_cultures.tsv``
    One row per BioSample with a `culture_collection`: biosample accession,
    assembly accession (GCF/GCA when the sample has one), serovar, collection
    and strain. The strongest public evidence class.
``ncbi_refseq_serovars.tsv``
    Every RefSeq *S. enterica* assembly with the serovar parsed out of
    `organism_name`, i.e. the broad asserted layer, useful for breadth.
``README.md``
    Provenance, column meanings and the caveat, in the same words as this
    docstring.

Usage::

    python scripts/salmonella/fetch_serotype_evidence.py --out ~/data/sen
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/28901/dataset_report"
#: NCBI allows 3 requests/second without an API key; stay just under it.
PAUSE = 0.55
BATCH = 200

SEROVAR_RE = re.compile(r"serovar\s+([^()]+?)\s*(?:\bstr\.|\b\(|$)")


def _get(url: str, *, attempt_limit: int = 6) -> bytes:
    """One E-utilities/Datasets call with retries (the API drops long reads)."""
    for attempt in range(attempt_limit):
        try:
            with urllib.request.urlopen(url, timeout=300) as response:
                return response.read()
        except (
            urllib.error.URLError,
            TimeoutError,
            json.JSONDecodeError,
            OSError,
        ) as exc:  # pragma: no cover - network flakiness
            if attempt == attempt_limit - 1:
                raise SystemExit(f"giving up on {url}: {exc}") from exc
            print(f"  retry {attempt + 1} ({type(exc).__name__})", file=sys.stderr)
            time.sleep(4 * (attempt + 1))
    raise AssertionError("unreachable")


def _json(path: str, **params: str) -> dict[str, object]:
    """One E-utilities JSON call, retrying a non-JSON body (rate limiting and
    transient errors answer with plain text, not an envelope)."""
    url = f"{EUTILS}/{path}?{urllib.parse.urlencode(params)}"
    for attempt in range(6):
        time.sleep(PAUSE)
        body = _get(url)
        try:
            payload = json.loads(body.decode())
        except (json.JSONDecodeError, UnicodeDecodeError):
            if attempt == 5:
                raise SystemExit(f"non-JSON response from {path}: {body[:200]!r}") from None
            print(f"  retry {attempt + 1} (non-JSON body)", file=sys.stderr)
            time.sleep(5 * (attempt + 1))
            continue
        assert isinstance(payload, dict)
        return payload
    raise AssertionError("unreachable")


def _search(db: str, term: str, retmax: int) -> list[str]:
    payload = _json("esearch.fcgi", db=db, term=term, retmax=str(retmax), retmode="json")
    result = payload["esearchresult"]
    assert isinstance(result, dict)
    count = int(result["count"])
    print(f"  {db}: {count} records match", file=sys.stderr)
    ids: list[str] = [str(i) for i in result["idlist"]]  # type: ignore[union-attr]
    return ids


def _summary(db: str, ids: list[str]) -> dict[str, dict[str, object]]:
    payload = _json("esummary.fcgi", db=db, id=",".join(ids), retmode="json")
    result = payload.get("result")
    assert isinstance(result, dict)
    return {
        str(uid): obj  # type: ignore[misc]
        for uid, obj in result.items()
        if uid != "uids" and isinstance(obj, dict)
    }


def _linked_uids(dbfrom: str, db: str, ids: list[str]) -> dict[str, str]:
    """``source uid -> first linked uid`` for one batched elink call."""
    time.sleep(PAUSE)
    url = f"{EUTILS}/elink.fcgi?dbfrom={dbfrom}&db={db}&id={','.join(ids)}"
    root = ET.fromstring(_get(url))  # elink answers XML only
    out: dict[str, str] = {}
    for linkset in root.iter("LinkSet"):
        source = (linkset.findtext("IdList/Id") or "").strip()
        linked = [link.findtext("Id") for link in linkset.iter("Link")]
        if source and linked and linked[0]:
            out.setdefault(source, str(linked[0]))
    return out


def _biosample_attributes(sampledata: str) -> dict[str, str]:
    """``attribute_name -> value`` out of a BioSample ``sampledata`` XML blob."""
    out: dict[str, str] = {}
    if not sampledata:
        return out
    try:
        root = ET.fromstring(sampledata)
    except ET.ParseError:  # pragma: no cover - truncated upstream payload
        return out
    for element in root.iter("Attribute"):
        name = element.get("attribute_name") or element.get("display_name") or ""
        value = (element.text or "").strip()
        if name and value:
            out.setdefault(name, value)
    return out


#: Public culture collections. A serovar's type/reference strain lives in one.
COLLECTION_TOKENS = (
    "NCTC",
    "ATCC",
    "CFSAN",
    "LMG",
    "CCUG",
    "CIP",
    "DSM",
    "NBRC",
    "JCM",
    "SGSC",
    "RIMD",
    "KACC",
)

#: The BioSample attribute is frequently filled with a placeholder rather than a
#: collection number (measured: 5,946 ``missing``, 3,891 ``not applicable``).
_JUNK = {"", "missing", "not applicable", "not provided", "not collected", "na", "n/a"}


def _real(value: str) -> str:
    return "" if value.strip().lower() in _JUNK else value.strip()


def fetch_reference_cultures(out_dir: Path) -> Path:
    """Assemblies whose strain is deposited in a public culture collection.

    Assembly-first on purpose: the deliverable is a *genome* someone can download.
    The collection token is what makes a serovar label closer to a reference
    standard — the type/reference cultures originate from the WHO Collaborating
    Centre — so a row here is the strongest public evidence class available.
    """
    path = out_dir / "reference_cultures.tsv"
    token_query = " OR ".join(f"{token}[All Fields]" for token in COLLECTION_TOKENS)
    ids = _search(
        "assembly",
        f'"Salmonella enterica"[Organism] AND ({token_query})',
        5000,
    )
    rows: list[list[str]] = []
    for start in range(0, len(ids), BATCH):
        chunk = ids[start : start + BATCH]
        assemblies = _summary("assembly", chunk)
        biosamples = {
            uid: record.get("biosampleaccn")
            for uid, record in assemblies.items()
            if record.get("biosampleaccn")
        }
        # assembly -> BioSample attributes. `esummary` takes uids only and elink is
        # not comprehensive, so the accessions are resolved with one OR-joined
        # esearch per 50 samples: a single query, no per-record round trip.
        attributes_of: dict[str, dict[str, str]] = {}
        wanted = [str(acc) for acc in biosamples.values() if acc]
        for offset in range(0, len(wanted), 50):
            part = wanted[offset : offset + 50]
            query = " OR ".join(part)
            uids = _search("biosample", query, 50)
            if not uids:
                continue
            for _uid, record in _summary("biosample", uids).items():
                attributes_of[str(record.get("accession") or "")] = _biosample_attributes(
                    str(record.get("sampledata") or "")
                )
        for record in assemblies.values():
            accession = str(record.get("assemblyaccession") or "")
            biosample = str(record.get("biosampleaccn") or "")
            attributes = attributes_of.get(biosample, {})
            collection = _real(attributes.get("culture_collection", ""))
            if not collection:
                continue  # the search token matched something else
            organism = str(record.get("organism") or "")
            serovar = _real(attributes.get("serovar", ""))
            if not serovar:
                match = SEROVAR_RE.search(organism)
                serovar = match.group(1).strip() if match else ""
            rows.append(
                [
                    accession,
                    biosample,
                    serovar,
                    collection,
                    _real(attributes.get("strain", "")),
                    _real(attributes.get("collected_by", "")),
                    str(record.get("assemblystatus") or ""),
                    organism,
                ]
            )
        print(
            f"  assemblies {min(start + BATCH, len(ids))}/{len(ids)} -> {len(rows)} kept",
            file=sys.stderr,
        )
    rows.sort(key=lambda row: (row[0], row[1]))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "assembly",
                "biosample",
                "serovar",
                "culture_collection",
                "strain",
                "collected_by",
                "assembly_level",
                "organism",
            ]
        )
        writer.writerows(rows)
    with_serovar = sum(1 for row in rows if row[2])
    print(
        f"{path}: {len(rows)} collection strains ({with_serovar} with a serovar)",
        file=sys.stderr,
    )
    return path


def fetch_refseq_serovars(out_dir: Path) -> Path:
    """Every RefSeq *S. enterica* assembly, serovar parsed from organism_name."""
    path = out_dir / "ncbi_refseq_serovars.tsv"
    columns = ["accession", "serovar", "assembly_level", "organism_name", "bioproject"]
    written = 0
    token: str | None = None
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        while True:
            params = {"filters.assembly_source": "REFSEQ", "page_size": "250"}
            if token:
                params["page_token"] = token
            time.sleep(PAUSE)
            page = json.loads(_get(f"{DATASETS}?{urllib.parse.urlencode(params)}").decode())
            for report in page.get("reports", []):
                organism = str((report.get("organism") or {}).get("organism_name") or "")
                match = SEROVAR_RE.search(organism)
                writer.writerow(
                    [
                        report.get("accession", ""),
                        match.group(1).strip() if match else "",
                        (report.get("assembly_info") or {}).get("assembly_level", ""),
                        organism,
                        (report.get("assembly_info") or {}).get("bioproject_accession", ""),
                    ]
                )
                written += 1
            token = page.get("next_page_token")
            if written % 2500 < 250:
                print(f"  refseq assemblies {written}", file=sys.stderr)
            if not token:
                break
    print(f"{path}: {written} RefSeq assemblies", file=sys.stderr)
    return path


README = """# `~/data/sen` — Salmonella serotype evidence

## What this is

Per-accession tables of *Salmonella enterica* serotype labels **with their
provenance**, built from public NCBI records. It exists because the public
record does not say which genomes were serotyped by agglutination: the `serovar`
attribute is a submitter assertion and the method is not recorded anywhere.

## Files

### `reference_cultures.tsv`

The strongest public evidence class: assemblies whose strain is deposited in a
public collection (NCTC, ATCC, CFSAN, LMG, CCUG, CIP, DSM, NBRC, JCM, SGSC, RIMD,
KACC). A serovar's **type or reference culture** originates from the WHO
Collaborating Centre, so these are the closest thing to a gold standard that has
a sequenced genome. **158 assemblies** survive the collection filter for
*S. enterica* (128 of them also carry a serovar attribute).

| column | meaning |
| --- | --- |
| `assembly` | assembly accession — the genome to download (`GCA_…` / `GCF_…`) |
| `biosample` | the BioSample behind it, for provenance |
| `serovar` | the submitter's serovar attribute (empty when unset) |
| `culture_collection` | collection + number, e.g. `ATCC:43845` |
| `strain` | strain designation |
| `collected_by` | isolating laboratory (e.g. `FDA`) |
| `assembly_level` | `Complete Genome` / `Chromosome` / `Contig` / `Scaffold` |
| `organism` | full organism string |

A note on the filter: the E-utilities text search matches the collection token
anywhere, so the candidate set is larger than the survivors; a row is kept only
when the BioSample's `culture_collection` holds a real number (not the `missing`
/ `not applicable` placeholders that 9,837 of the 10,643 matching BioSamples
carry).

### `tie_serovars_candidates.tsv`

Every RefSeq assembly whose asserted serovar is one of the four that gspec cannot
separate (`Hindmarsh`, `Bovismorbificans`, `Goldcoast`, `Brikama`) — the raw
material for a verified re-analysis, with the label's evidence class implied by
the files above. Measured: **Hindmarsh 2 (both contig-level), Bovismorbificans 63
(4 complete), Goldcoast 47 (11 complete), Brikama 0**.

### `ncbi_refseq_serovars.tsv`

The broad, *asserted* layer: every RefSeq *S. enterica* assembly (18,970 at the
time of writing) with the serovar parsed out of `organism_name`
(`Salmonella enterica subsp. enterica serovar Typhimurium str. LT2` → Typhimurium).
Use it for breadth, not for truth.

## Caveats — read before trusting a label

1. **No column says "agglutination-confirmed".** Catalogue provenance is a proxy,
   not proof. A culture-collection strain with a serovar label is *likely* to have
   been serotyped properly; nothing in INSDC guarantees it.
2. **Labels can be wrong.** Two measured examples: `GCF_004194605.3` is labelled
   *Albany* (O:8,20) yet SeqSero2 reads its O as `8` and its rfb cluster groups
   with the `6,8` genomes; the only public *Brikama* assembly
   (`GCA_021104055.1`) is a divergent contig-level genome.
3. **The collection is not evenly distributed.** Measured against SISTR's
   84,463-genome reference list: Bovismorbificans 386 genomes, Goldcoast 119,
   **Hindmarsh 1, Brikama 0**. A nearest-neighbour tool therefore cannot resolve
   rare twins, which is why its answers on such pairs are not evidence.
4. **NCBI, not gspec, owns these records.** Re-run the fetcher to refresh:
   `python scripts/salmonella/fetch_serotype_evidence.py --out ~/data/sen`.

## The verdict for the four unseparable pairs (2026-11-18)

| serovar | collection strains | RefSeq assemblies (asserted) | complete |
| --- | --- | --- | --- |
| Bovismorbificans | 0 | 63 | 4 |
| Goldcoast | 0 | 47 | 11 |
| Hindmarsh | 0 | 2 | 0 |
| **Brikama** | **0** | **0** | **0** |

**No catalogue-deposited (reference-culture) genome exists for any of the four**
— the strongest evidence class is empty for all of them. `Brikama` has no RefSeq
assembly at all and no entry in SISTR's 84,463-genome reference list. So the tie
cannot be settled from public data by provenance alone; it needs strains
serotyped by a reference laboratory, and the rare twins (Hindmarsh, Brikama)
need to be sourced from the laboratory, not from NCBI.

## The merge layer (`merged_evidence.tsv`)

`scripts/salmonella/merge_compare_serotypes.py merge` joins the four layers above
by accession (version-suffix independent) and keeps the strongest label, ranked
`verified` > `culture` > `published` > `asserted`. Two extra files appear:

* `merged_evidence.tsv` — one row per accession: winning serovar, the layer it
  came from, every source, and the full label set.
* `label_conflicts.tsv` — the rows where two sources name **different** serovars.
  Measured after the deep paper harvest: 77,249 accessions over 1,478 conflicts.
  The conflict count rises with the number of papers because one accession now
  appears in several tables, and a paper may list both a submitted and a computed
  serovar; every conflicting accession is excluded from marker analysis, because
  its own sources disagree about what it is.

Then `… compare --group NAME=serovar,serovar …` downloads the genomes, cuts the
rfb neighbourhood out of each (anchored on a panel wzx gene), clusters the CDS,
and reports whether any cluster is present in every genome of one group and none
of the other. First run (O:6 question, 11 genomes, 2,640 CDS, 2,065 clusters):
**no perfect group-splitting cluster** — matching the manual analysis, which is
the point of having the tool.

## Brikama: no genome, but the literature has the strains (2026-11-19)

Brikama is a rare serovar **whose name collides with Brikama, a town in The
Gambia** — a text search for it returns 14 open-access papers, and most are
about patients from that health centre, not the serovar. Filtering to
serovar-style mentions leaves two real leads, both recorded in
`published/brikama_leads.tsv`:

| paper | what it says |
| --- | --- |
| PMC8552638 | isolate **A030_06/19**: submitted serovar **S. Goldcoast**
(`8:r:l,w`); the authors' own in-silico call was "**S. Goldcoast or
S. Brikama**" — the exact tie gspec hits, resolved phenotypically to Goldcoast |
| PMC9440507 | **Brikama 4.6 % (n = 5)** of a Gambian NTS collection,
"overrepresented in Cluster 6" — five isolates exist in a published collection |

So the picture is:

| layer | Brikama found |
| --- | --- |
| catalogue strains (18,970 RefSeq assemblies filtered to a real `culture_collection`) | **0** |
| asserted labels (all 18,970 RefSeq assemblies) | **0** |
| machine-readable paper tables (1,199 papers mined) | **0 rows** |
| paper *text* | **2 papers, 6 isolates** — actionable, but no accessions in the text |
| SISTR's reference list (84,463 genomes / 780 serovars) | **0** |

Nothing upstream ships a Brikama genome; the strains are in the literature and in
a laboratory collection. **Chasing those isolates is the only route**, and the
PMC8552638 evidence is worth noting for the tie itself: when a laboratory had to
commit, it said **Goldcoast**.

## What would actually be gold standard

Agglutination panels from a national reference laboratory (WHO Collaborating
Centre at Institut Pasteur, CDC/PulseNet, UKHSA), or a published table that
reports the method per strain. If you have such a set, add it as a local table
here rather than inferring confirmation from catalogue status.
"""


def write_readme(out_dir: Path) -> None:
    (out_dir / "README.md").write_text(README, encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path.home() / "data" / "sen")
    parser.add_argument(
        "--skip-refseq",
        action="store_true",
        help="only rebuild reference_cultures.tsv (the RefSeq sweep takes ~20 pages)",
    )
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    fetch_reference_cultures(args.out)
    if not args.skip_refseq:
        fetch_refseq_serovars(args.out)
    write_readme(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
