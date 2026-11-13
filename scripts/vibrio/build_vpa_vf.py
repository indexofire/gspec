#!/usr/bin/env python3
"""Build the bundled V. parahaemolyticus haemolysin marker database.

``vpa_vf`` is a gene database: the species marker ``tlh`` plus the virulence
markers ``tdh`` and ``trh`` (two allele records).  The source sequences are
public INSDC records, cut from GenBank flat files:

- RIMD 2210633 chromosome 2 ``NC_004605.1`` — ``tdh`` (VPA1314) and ``tlh``
  (VPA0226; the current RefSeq re-annotation names the identical sequence a
  "type 2 secretion system lipase", protein WP_005481504.1);
- strain 10329 chromosome 2 ``CP045795.1`` — the ``trh1`` CDS (QGG36210.1);
- strain GTVSS-031 chromosome 2 ``CP176647.1`` — the ``trh2`` CDS
  (XLX46591.1).

Both trh records carry the single gene name ``trh``: the allele travels in
the accession and product, so screening output and the summary matrix key on
``trh`` alone (the record semantics of the locally-built panel this bundle
replaces — sequences byte-identical, annotations re-derived from INSDC).

Usage::

    python scripts/vibrio/build_vpa_vf.py \
        --rimd2 NC_004605.1.gb --trh1 CP045795.1.gb --trh2 CP176647.1.gb \
        --dest src/gspec/data/dbs/vpa_vf
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import UTC, datetime
from pathlib import Path

from gspec.dbcodec import encode_seqid

#: gene -> (source label, qualifier key, qualifier value, product, accession)
MARKERS: tuple[tuple[str, str, str, str, str, str], ...] = (
    (
        "tdh",
        "rimd2",
        "old_locus_tag",
        "VPA1314",
        "thermostable direct hemolysin TDH",
        "WP_005463268.1",
    ),
    (
        "trh",
        "trh1",
        "gene",
        "trh1",
        "TDH-related hemolysin TRH1 allele",
        "QGG36210.1",
    ),
    (
        "trh",
        "trh2",
        "gene",
        "trh2",
        "TDH-related hemolysin TRH2 allele",
        "XLX46591.1",
    ),
    (
        "tlh",
        "rimd2",
        "old_locus_tag",
        "VPA0226",
        "thermolabile hemolysin TLH (lecithin-dependent hemolysin)",
        "WP_005481504.1",
    ),
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


def _cut(cds: list[_Cds], sequence: str, source: str, key: str, value: str) -> _Cds:
    """Pick the one CDS whose qualifier matches; verify it translates."""
    matches = [feature for feature in cds if feature.qualifiers.get(key) == value]
    if len(matches) != 1:
        raise SystemExit(f"{source}: expected one CDS with {key}={value}, got {len(matches)}")
    feature = matches[0]
    gene = sequence[feature.start - 1 : feature.end]
    if feature.complement:
        gene = gene.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    protein = "".join(_CODE[gene[i : i + 3]] for i in range(0, len(gene), 3))
    translation = feature.qualifiers.get("translation", "").replace(" ", "").strip('"')
    if not gene.startswith("ATG") or not protein.endswith("*") or protein[:-1] != translation:
        raise SystemExit(f"{source} {key}={value}: CDS does not translate as annotated")
    return feature


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("rimd2", "trh1", "trh2"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--dest", type=Path, required=True)
    args = parser.parse_args(argv)

    sources = {
        "rimd2": _parse_genbank(args.rimd2),
        "trh1": _parse_genbank(args.trh1),
        "trh2": _parse_genbank(args.trh2),
    }
    dest: Path = args.dest
    dest.mkdir(parents=True, exist_ok=True)
    with (dest / "sequences").open("w", encoding="utf-8", newline="\n") as handle:
        for gene, source, key, value, product, accession in MARKERS:
            cds, sequence = sources[source]
            feature = _cut(cds, sequence, source, key, value)
            seqid = encode_seqid("vpa_vf", gene, accession, ("virulence",))
            nucleotide = sequence[feature.start - 1 : feature.end]
            if feature.complement:
                nucleotide = nucleotide.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            handle.write(f">{seqid} {product}\n")
            for offset in range(0, len(nucleotide), 60):
                handle.write(nucleotide[offset : offset + 60] + "\n")
    (dest / "bundled.json").write_text(
        json.dumps(
            {
                "name": "vpa_vf",
                "description": (
                    "V. parahaemolyticus haemolysin panel: tlh species marker"
                    " + tdh/trh virulence markers"
                ),
                "vendor": "gspec-curated (INSDC)",
                "dbtype": "nucl",
                "kind": "gene",
                "snapshotted": datetime.now(UTC).strftime("%Y-%m-%d"),
                "license": "Public domain (INSDC sequences)",
                "note": (
                    "tdh (VPA1314, WP_005463268.1), tlh (VPA0226,"
                    " WP_005481504.1), and the two trh allele records (TRH1"
                    " QGG36210.1 from CP045795.1, TRH2 XLX46591.1 from"
                    " CP176647.1) cut from public INSDC records: RIMD"
                    " 2210633 chromosome 2 (NC_004605.1), strain 10329"
                    " chromosome 2 (CP045795.1), strain GTVSS-031 chromosome"
                    " 2 (CP176647.1). Sequences are byte-identical to the"
                    " earlier locally-built panel; annotations re-derived"
                    " from INSDC so no third-party database terms ride the"
                    " distribution. Screen-only; no typing scheme."
                ),
                "guidance": {
                    "references": [
                        {
                            "citation": (
                                "Makino K, Oshima K, Kurokawa K, et al."
                                " Genome sequence of Vibrio parahaemolyticus:"
                                " a pathogenic mechanism distinct from that of"
                                " V. cholerae. Lancet. 2003;361(9359):743-749."
                            ),
                            "doi": "10.1016/S0140-6736(03)12659-1",
                            "note": "RIMD 2210633, the genome behind the tdh and tlh records.",
                        },
                        {
                            "citation": (
                                "Nishibuchi M, Taniguchi T, Misawa T,"
                                " Khaeomanee-iam S, Honda T, Miwatani T."
                                " Cloning and nucleotide sequence of the gene"
                                " (trh) encoding the hemolysin related to the"
                                " thermostable direct hemolysin of Vibrio"
                                " parahaemolyticus. Infect Immun."
                                " 1989;57(9):2691-2697."
                            ),
                            "url": "https://pubmed.ncbi.nlm.nih.gov/2759706/",
                            "note": "trh discovery; the two-allele framework this panel encodes.",
                        },
                        {
                            "citation": (
                                "Nishibuchi M, Kaper JB. Thermostable direct"
                                " hemolysin gene of Vibrio parahaemolyticus: a"
                                " virulence gene acquired by a marine"
                                " bacterium. Infect Immun. 1995;63(6):2093-2099."
                            ),
                            "url": "https://pubmed.ncbi.nlm.nih.gov/7768586/",
                            "note": "tdh/trh biology and strain distribution.",
                        },
                        {
                            "citation": (
                                "Nishibuchi M, Ishibashi M, Takeda Y, Kaper JB."
                                " Duplication and variation of the thermostable"
                                " direct haemolysin (tdh) gene in Vibrio"
                                " parahaemolyticus. Mol Microbiol."
                                " 1990;4(1):87-99."
                            ),
                            "url": "https://pubmed.ncbi.nlm.nih.gov/2319944/",
                            "note": "The tdh gene family: RIMD carries two copies.",
                        },
                    ],
                    "caveats": [
                        (
                            "trh is ONE gene name for two allele records: TRH1"
                            " and TRH2 are alleles of the same locus, so GENE"
                            " is 'trh' and the allele travels in the accession"
                            " and product. Typing rules and the summary matrix"
                            " therefore key on trh alone; read the accession"
                            " when the allele matters."
                        ),
                        (
                            "One contig carrying trh reports a single row even"
                            " with both allele records installed: the baseline"
                            " dedups hits sharing (contig, start, end) and"
                            " blastn's culling keeps the best-scoring subject,"
                            " so the surviving accession identifies the allele."
                        ),
                        (
                            "tdh (thermostable direct haemolysin) is the marker"
                            " of clinical isolates (Kanagawa phenomenon); trh is"
                            " associated with pandemic O3:K6 strains."
                            " Environmental isolates often carry neither. tlh is"
                            " the species marker, carried by nearly all"
                            " V. parahaemolyticus (the FDA-BAM real-time PCR"
                            " target); current RefSeq annotates the identical"
                            " sequence as a type 2 secretion system lipase"
                            " (WP_005481504.1, historically NP_799736"
                            " 'lecithin-dependent hemolysin')."
                        ),
                        (
                            "RIMD 2210633 carries two tdh copies (VPA1314 and"
                            " VPA1378) 97% identical to each other; this panel"
                            " ships VPA1314, and near-identical tdh variants in"
                            " other strains still hit it — read identity, not"
                            " just presence, when allele discrimination"
                            " matters."
                        ),
                        (
                            "Haemolysin gene presence does not by itself"
                            " establish pathogenicity: strain background and"
                            " the type III secretion systems (T3SS1,"
                            " T3SS2alpha, T3SS2beta) matter, and those are not"
                            " in this panel."
                        ),
                    ],
                    "usage": [
                        "gspec screen isolate.fa -d vpa_vf --output vpa.tsv",
                        "gspec screen isolate.fa -d vpa_vf --format json",
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
