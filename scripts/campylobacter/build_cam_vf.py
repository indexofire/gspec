#!/usr/bin/env python3
"""Build the bundled C. jejuni virulence marker database.

``cam_vf`` is a gene database: the T4SS ``virB`` panel, the CDT tripartite
toxin, the invasion antigen CiaB, and the LOS-modification markers.  The
source sequences are public INSDC records, cut from GenBank flat files and
pinned by ``protein_id``:

- NCTC 11168 chromosome ``NC_002163.1`` — ciaB, cdtA/B/C, wlaN, and the
  cgtA/neuA1 fusion (annotated ``neuA1``);
- 81-176 plasmid pVir ``AF226280.2`` (RefSeq NC_005012) — virB4, virB11;
- 81-176 plasmid pVir ``CP000550.1`` (RefSeq NC_008770) — virB8/9/10;
- ``LR134499.1`` (NCTC 13255) — ggt;
- ``JX436460.1`` (strain 108 T6SS locus) — hcp (annotated TssD);
- ``AF400048.1`` — cstII; ``AF130984.1`` (OH4384) — cgtA, cgtB.

Sequences are byte-identical to the earlier locally-built panel except
``cgtB``, whose record is strand-corrected to the coding orientation (the
old panel stored the plus-strand slice of the ``complement(5054..5959)``
CDS; blastn is strand-agnostic so screening output is unchanged).

Usage::

    python scripts/campylobacter/build_cam_vf.py \\
        --nctc11168 NC_002163.1.gb --pvir_af AF226280.2.gb \\
        --pvir_cp CP000550.1.gb --nctc13255 LR134499.1.gb \\
        --t6ss JX436460.1.gb --cstii AF400048.1.gb --los AF130984.1.gb \\
        --dest src/gspec/data/dbs/cam_vf
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import UTC, datetime
from pathlib import Path

from gspec.dbcodec import encode_seqid

#: gene -> (source label, protein_id, product).  Two cgtA records on purpose:
#: the NCTC 11168 neuA1 fusion and the OH4384 allele (panel semantics — the
#: allele travels in the accession and product).
MARKERS: tuple[tuple[str, str, str, str], ...] = (
    ("cdtA", "nctc11168", "YP_002343541.1", "cytolethal distending toxin subunit A (CDT)"),
    ("cdtB", "nctc11168", "YP_002343540.1", "cytolethal distending toxin subunit B (CDT)"),
    ("cdtC", "nctc11168", "YP_002343539.1", "cytolethal distending toxin subunit C (CDT)"),
    (
        "cgtA",
        "nctc11168",
        "YP_002344536.1",
        "bifunctional beta-1,4-N-acetylgalactosaminyltransferase/CMP-Neu5Ac"
        " synthetase (CgtA/NeuA1 fusion)",
    ),
    (
        "cgtA",
        "los",
        "AAF31769.1",
        "beta-1,4-N-acetylgalactosaminyltransferase (CgtA, OH4384 allele)",
    ),
    ("cgtB", "los", "AAF31770.1", "beta-1,3-galactosyltransferase (CgtB)"),
    ("ciaB", "nctc11168", "YP_002344312.1", "Campylobacter invasion antigen B (CiaB)"),
    ("cstII", "cstii", "AAK91725.1", "alpha-2,3-/alpha-2,8-sialyltransferase (CstII)"),
    ("ggt", "nctc13255", "VEI83882.1", "gamma-glutamyltransferase (GGT)"),
    ("hcp", "t6ss", "AGN32396.1", "type VI secretion system tube protein Hcp (TssD)"),
    ("virB10", "pvir_cp", "EAQ71778.1", "type IV secretion system protein VirB10"),
    ("virB11", "pvir_af", "AAF97748.1", "type IV secretion system protein VirB11"),
    ("virB4", "pvir_af", "AAN46946.1", "type IV secretion system protein VirB4"),
    ("virB8", "pvir_cp", "EAQ71761.1", "type IV secretion system protein VirB8"),
    ("virB9", "pvir_cp", "EAQ71741.1", "type IV secretion system protein VirB9"),
    ("wlaN", "nctc11168", "YP_002344532.1", "beta-1,3-galactosyltransferase (WlaN)"),
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
    for flag in (
        "--nctc11168",
        "--pvir-af",
        "--pvir-cp",
        "--nctc13255",
        "--t6ss",
        "--cstii",
        "--los",
    ):
        parser.add_argument(flag, type=Path, required=True)
    parser.add_argument("--dest", type=Path, required=True)
    args = parser.parse_args(argv)

    sources = {
        label: _parse_genbank(getattr(args, label.replace("-", "_")))
        for label in (
            "nctc11168",
            "pvir_af",
            "pvir_cp",
            "nctc13255",
            "t6ss",
            "cstii",
            "los",
        )
    }
    dest: Path = args.dest
    dest.mkdir(parents=True, exist_ok=True)
    with (dest / "sequences").open("w", encoding="utf-8", newline="\n") as handle:
        for gene, source, protein_id, product in MARKERS:
            cds, sequence = sources[source]
            feature = _cut(cds, sequence, source, protein_id)
            seqid = encode_seqid("cam_vf", gene, protein_id, ("virulence",))
            nucleotide = sequence[feature.start - 1 : feature.end]
            if feature.complement:
                nucleotide = nucleotide.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            handle.write(f">{seqid} {product}\n")
            for offset in range(0, len(nucleotide), 60):
                handle.write(nucleotide[offset : offset + 60] + "\n")
    (dest / "bundled.json").write_text(
        json.dumps(
            {
                "name": "cam_vf",
                "description": (
                    "C. jejuni virulence panel: T4SS virB, CDT, CiaB, GGT,"
                    " T6SS Hcp and LOS-modification markers"
                ),
                "vendor": "gspec-curated (INSDC)",
                "dbtype": "nucl",
                "kind": "gene",
                "snapshotted": datetime.now(UTC).strftime("%Y-%m-%d"),
                "license": "Public domain (INSDC sequences)",
                "note": (
                    "Sixteen markers cut from public INSDC records: NCTC"
                    " 11168 chromosome 1 NC_002163.1 (cdtA/B/C, ciaB, wlaN,"
                    " the cgtA/neuA1 fusion), 81-176 plasmid pVir via"
                    " AF226280.2 (virB4, virB11) and"
                    " CP000550.1 (virB8/9/10), NCTC 13255 via LR134499.1"
                    " (ggt), the strain 108 T6SS locus JX436460.1 (hcp,"
                    " annotated TssD), AF400048.1 (cstII) and AF130984.1"
                    " (cgtA OH4384 allele, cgtB). Sequences are"
                    " byte-identical to the earlier locally-built panel"
                    " except cgtB, strand-corrected here to the coding"
                    " orientation of complement(5054..5959); annotations"
                    " re-derived from INSDC so no third-party database"
                    " terms ride the distribution. Screen-only; no typing"
                    " scheme."
                ),
                "guidance": {
                    "references": [
                        {
                            "citation": (
                                "Parkhill J, Wren BW, Mungall K, et al. The"
                                " genome sequence of the food-borne pathogen"
                                " Campylobacter jejuni reveals hypervariable"
                                " sequences. Nature. 2000;403(6770):665-668."
                            ),
                            "url": "https://pubmed.ncbi.nlm.nih.gov/10688204/",
                            "note": "NCTC 11168: cdtABC, ciaB, wlaN and the neuA1 fusion.",
                        },
                        {
                            "citation": (
                                "Bacon DJ, Szymanski CM, Burr DH, Silver RP,"
                                " Alm RA, Guerry P. DNA sequence and"
                                " mutational analyses of the pVir plasmid of"
                                " Campylobacter jejuni 81-176. Infect Immun."
                                " 2002;70(11):6242-6250."
                            ),
                            "url": "https://pubmed.ncbi.nlm.nih.gov/12379703/",
                            "note": "The pVir T4SS virB panel (AF226280; virB8-10 on CP000550).",
                        },
                        {
                            "citation": (
                                "Pickett CL, Cottle DL, Pesci EC, Bikah G."
                                " Prevalence of cytolethal distending toxin"
                                " production in Campylobacter jejuni and"
                                " related species. Infect Immun."
                                " 1996;64(6):2070-2078."
                            ),
                            "url": "https://pubmed.ncbi.nlm.nih.gov/8675309/",
                            "note": "CDT prevalence and the cdtABC locus.",
                        },
                        {
                            "citation": (
                                "Cheng J, Hyun J, Kapoor M, et al."
                                " Multifunctionality of Campylobacter jejuni"
                                " sialyltransferase CstII: characterization of"
                                " the beta-1,4-N-acetylgalactosaminyltransfer"
                                " ase and alpha-2,8-sialyltransferase"
                                " activities. Glycobiology."
                                " 2008;18(9):686-697."
                            ),
                            "url": "https://pubmed.ncbi.nlm.nih.gov/18509108/",
                            "note": "CstII/cgtA LOS biosynthesis and ganglioside mimicry.",
                        },
                    ],
                    "caveats": [
                        (
                            "Marker presence is genetic potential, not"
                            " virulence: CDT and T4SS genes can be partial or"
                            " lost, and the LOS markers depend on phase"
                            " variable homopolymer tracts not captured by"
                            " presence screening."
                        ),
                        (
                            "cgtA is TWO records on purpose: the NCTC 11168"
                            " neuA1 fusion (1611 nt) and the OH4384 allele"
                            " (1044 nt) — key on the accession when the"
                            " distinction matters. cgtB is the coding strand"
                            " of the OH4384 locus (strand-corrected from the"
                            " earlier local panel, which stored the"
                            " plus-strand slice; blastn is strand-agnostic so"
                            " screening output is unchanged)."
                        ),
                        (
                            "virB4/virB11 come from the original pVir"
                            " submission (AF226280) and virB8/9/10 from the"
                            " 81-176 genome project's re-annotation"
                            " (CP000550): the two records differ by 5 bp"
                            " overall, so cross-hits between the two virB"
                            " source copies read near 100% identity."
                        ),
                        (
                            "cstII/cgtA/cgtB encode LOS ganglioside mimics"
                            " (GM1/GD1a) associated with Guillain-Barre"
                            " syndrome strains; absence of the gene does not"
                            " exclude a mimicking LOS made by other alleles."
                        ),
                        (
                            "This panel does not replace C. jejuni/C. coli"
                            " species identification (combine with the cpn60"
                            " or lpxA-based assays) or capsule typing (see"
                            " cje_cps)."
                        ),
                    ],
                    "usage": [
                        "gspec screen isolate.fa -d cam_vf --output cam.tsv",
                        "gspec screen isolate.fa -d cam_vf --format json",
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
