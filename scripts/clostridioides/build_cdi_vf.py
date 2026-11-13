#!/usr/bin/env python3
"""Build the bundled C. difficile toxin marker database.

``cdi_vf`` is a gene database: the large glucosylating toxins TcdA/TcdB and
the binary toxin CDT subunits.  The source sequences are public INSDC
records, cut from GenBank flat files and pinned by ``protein_id``:

- strain 630 chromosome ``NC_009089.1`` — tcdA (WP_011860904.1), tcdB
  (WP_009902069.1), the PaLoc toxins;
- ``AF271719.1`` (strain CCUG 20309) — the CDT locus: cdtA (AAF81760.1),
  cdtB (AAF81761.1).

Sequences are byte-identical to the earlier locally-built panel;
annotations re-derived from INSDC.

Usage::

    python scripts/clostridioides/build_cdi_vf.py \\
        --six30 NC_009089.1.gb --cdt AF271719.1.gb \\
        --dest src/gspec/data/dbs/cdi_vf
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import UTC, datetime
from pathlib import Path

from gspec.dbcodec import encode_seqid

#: gene -> (source label, protein_id, product)
MARKERS: tuple[tuple[str, str, str, str], ...] = (
    (
        "cdtA",
        "cdt",
        "AAF81760.1",
        "cytolethal distending toxin subunit A (CDTa, ADP-ribosyltransferase)",
    ),
    ("cdtB", "cdt", "AAF81761.1", "cytolethal distending toxin subunit B (CDTb)"),
    ("tcdA", "six30", "WP_011860904.1", "glucosylating toxin TcdA (toxin A)"),
    ("tcdB", "six30", "WP_009902069.1", "glucosylating toxin TcdB (toxin B)"),
)

_BASES = "TCAG"
_AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODE = dict(
    zip(
        (a + b + c for a in _BASES for b in _BASES for c in _BASES),
        _AMINO_ACIDS,
        strict=True,
    )
)


class _Cds:
    """One CDS feature: location strand-bits plus qualifier values."""

    def __init__(self, location: str, qualifiers: dict[str, str]) -> None:
        self.complement = location.startswith("complement(")
        bounds = re.findall(r"\d+", location)
        if len(bounds) < 2:
            raise SystemExit(f"unparsable CDS location: {location}")
        self.start, self.end = int(bounds[0]), int(bounds[-1])
        self.qualifiers = qualifiers


def _parse_genbank(path: Path) -> tuple[list[_Cds], str]:
    """Return the CDS features and the ORIGIN sequence of a single-record
    GenBank flat file (qualifier values joined; continuations folded)."""
    cds: list[_Cds] = []
    sequence: list[str] = []
    in_features = False
    in_origin = False
    key = ""
    location = ""
    qualifiers: dict[str, str] = {}
    last_qualifier = ""
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if line.startswith("ORIGIN"):
            in_features = False
            in_origin = True
            continue
        if line.startswith("//"):
            break
        if in_origin:
            sequence.append(re.sub(r"[^A-Za-z]", "", line))
        elif in_features:
            if line.startswith(" " * 5) and not line.startswith(" " * 6):
                if key == "CDS" and location:
                    cds.append(_Cds(location, qualifiers))
                parts = line[5:].split(maxsplit=1)
                key = parts[0]
                location = parts[1] if len(parts) == 2 else ""
                qualifiers = {}
                last_qualifier = ""
            elif line.startswith(" " * 21) and line[21] == "/":
                name, _, value = line[22:].partition("=")
                qualifiers[name] = value.strip('"')
                last_qualifier = name
            elif line.startswith(" " * 21) and last_qualifier:
                # qualifier continuation line: drop the closing quote at wrap
                stored = qualifiers[last_qualifier].rstrip('"')
                folded = stored + " " + line[21:].strip().lstrip('"')
                qualifiers[last_qualifier] = re.sub(r"\s+", " ", folded)
    if key == "CDS" and location:
        cds.append(_Cds(location, qualifiers))
    return cds, "".join(sequence).upper()


def _cut(cds: list[_Cds], sequence: str, source: str, protein_id: str) -> _Cds:
    """Pick the one CDS whose protein_id matches; verify it translates."""
    matches = [feature for feature in cds if feature.qualifiers.get("protein_id", "") == protein_id]
    if len(matches) != 1:
        raise SystemExit(
            f"{source}: expected one CDS with protein_id={protein_id}, got {len(matches)}"
        )
    feature = matches[0]
    gene = sequence[feature.start - 1 : feature.end]
    if feature.complement:
        gene = gene.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    protein = "".join(_CODE[gene[i : i + 3]] for i in range(0, len(gene), 3))
    translation = feature.qualifiers.get("translation", "").replace(" ", "").strip('"')
    if not gene.startswith("ATG") or not protein.endswith("*") or protein[:-1] != translation:
        raise SystemExit(f"{source} {protein_id}: CDS does not translate as annotated")
    return feature


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--six30", type=Path, required=True)
    parser.add_argument("--cdt", type=Path, required=True)
    parser.add_argument("--dest", type=Path, required=True)
    args = parser.parse_args(argv)

    sources = {"six30": _parse_genbank(args.six30), "cdt": _parse_genbank(args.cdt)}
    dest: Path = args.dest
    dest.mkdir(parents=True, exist_ok=True)
    with (dest / "sequences").open("w", encoding="utf-8", newline="\n") as handle:
        for gene, source, protein_id, product in MARKERS:
            cds, sequence = sources[source]
            feature = _cut(cds, sequence, source, protein_id)
            seqid = encode_seqid("cdi_vf", gene, protein_id, ("virulence",))
            nucleotide = sequence[feature.start - 1 : feature.end]
            if feature.complement:
                nucleotide = nucleotide.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            handle.write(f">{seqid} {product}\n")
            for offset in range(0, len(nucleotide), 60):
                handle.write(nucleotide[offset : offset + 60] + "\n")
    (dest / "bundled.json").write_text(
        json.dumps(
            {
                "name": "cdi_vf",
                "description": (
                    "C. difficile toxin panel: PaLoc TcdA/TcdB + binary toxin CDT subunits"
                ),
                "vendor": "gspec-curated (INSDC)",
                "dbtype": "nucl",
                "kind": "gene",
                "snapshotted": datetime.now(UTC).strftime("%Y-%m-%d"),
                "license": "Public domain (INSDC sequences)",
                "note": (
                    "Four markers cut from public INSDC records: strain 630"
                    " chromosome NC_009089.1 (tcdA WP_011860904.1, tcdB"
                    " WP_009902069.1) and the CCUG 20309 CDT locus"
                    " AF271719.1 (cdtA AAF81760.1, cdtB AAF81761.1)."
                    " Sequences are byte-identical to the earlier"
                    " locally-built panel; annotations re-derived from"
                    " INSDC so no third-party database terms ride the"
                    " distribution. Screen-only; no typing scheme."
                ),
                "guidance": {
                    "references": [
                        {
                            "citation": (
                                "Sebaihia M, Wren BW, Mullany P, et al. The"
                                " multidrug-resistant human pathogen"
                                " Clostridium difficile has a highly mobile,"
                                " mosaic genome. Nat Genet."
                                " 2006;38(7):779-786."
                            ),
                            "url": "https://pubmed.ncbi.nlm.nih.gov/16804543/",
                            "note": "Strain 630, the genome behind the tcdA/tcdB records.",
                        },
                        {
                            "citation": (
                                "Perelle S, Gibert M, Bourlioux P, Corthier"
                                " G, Popoff MR. Production of a complete"
                                " binary toxin (actin-specific"
                                " ADP-ribosyltransferase) by Clostridium"
                                " difficile CD196. Infect Immun."
                                " 1997;65(4):1402-1407."
                            ),
                            "url": "https://pubmed.ncbi.nlm.nih.gov/9119480/",
                            "note": "The CDT locus (AF271719) behind the cdtA/cdtB records.",
                        },
                    ],
                    "caveats": [
                        (
                            "Toxin gene presence does not predict toxin"
                            " production: tcdC/regulatory changes, PaLoc"
                            " deletions and phase variation govern"
                            " expression, and RT078/BI/NAP1 strains carry"
                            " deletions in tcdC that this panel does not"
                            " assess."
                        ),
                        (
                            "Binary toxin (CDT) genes are present in a"
                            " minority of strains and can occur"
                            " independently of the PaLoc; a CDT-negative"
                            " PaLoc-positive strain still reads tcdA/tcdB"
                            " positive."
                        ),
                        (
                            "tcdA/tcdB allele variation (including the"
                            " tcdA-negative toxinotype XI and variant toxin"
                            " B genes) may lower identity below default"
                            " floors; read the identity column, not just"
                            " presence, when variant toxinotypes matter."
                        ),
                    ],
                    "usage": [
                        "gspec screen isolate.fa -d cdi_vf --output cdi.tsv",
                        "gspec screen isolate.fa -d cdi_vf --format json",
                    ],
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"wrote {dest}: {len(MARKERS)} markers")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
