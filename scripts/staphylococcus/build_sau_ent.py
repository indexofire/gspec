#!/usr/bin/env python3
"""Build the bundled S. aureus enterotoxin gene panel (``sau_ent``).

The reference set is the ``indexofire/enterotoxin`` abricate database — 25
staphylococcal enterotoxin, enterotoxin-like and TSST-1 genes, MIT licensed
(Copyright 2023 Mark Renton). The panel adds one species control: the
S. aureus thermonuclease gene ``nuc`` from public-domain INSDC
(NC_007795.1:800034-800720), so a screen always carries a species sanity
check next to the toxin hits.

The conversion is deliberately thin — the upstream FASTA already holds the
gene sequences; only the headers are rewritten into ``gspec/v1`` and a
curated product string is attached per gene (the upstream descriptions are
mixed accession notes and NCBI record titles).

Usage::

    python scripts/staphylococcus/build_sau_ent.py \\
        --source sequences \\
        --dest src/gspec/data/dbs/sau_ent

``--source`` accepts a local path or the pinned upstream raw URL (the
default). ``--control`` defaults to the committed nuc asset next to this
script, so a rebuild is offline apart from the upstream download.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

from gspec.dbcodec import encode_seqid

#: Pinned upstream: the repository's last commit (2023-08-27). A rebuild
#: against ``main`` that differs from this is a data change and must be
#: reviewed like one.
UPSTREAM_COMMIT = "30c022eeb712"
UPSTREAM_URL = (
    "https://raw.githubusercontent.com/indexofire/enterotoxin/"
    f"{UPSTREAM_COMMIT}/enterotoxin/sequences"
)

#: gene -> product string. Fixed per gene so the output does not inherit the
#: upstream headers' accession notes and record titles.
PRODUCTS: dict[str, str] = {
    "sea": "staphylococcal enterotoxin type A (SEA)",
    "seb": "staphylococcal enterotoxin type B (SEB)",
    "sec": "staphylococcal enterotoxin type C (SEC)",
    "sed": "staphylococcal enterotoxin type D (SED)",
    "see": "staphylococcal enterotoxin type E (SEE)",
    "seg": "staphylococcal enterotoxin type G (SEG)",
    "seh": "staphylococcal enterotoxin type H (SEH)",
    "sei": "staphylococcal enterotoxin type I (SEI)",
    "sej": "staphylococcal enterotoxin type J (SEJ)",
    "sek": "staphylococcal enterotoxin type K (SEK)",
    "sel": "staphylococcal enterotoxin type L (SEL)",
    "sem": "staphylococcal enterotoxin type M (SEM)",
    "sen": "staphylococcal enterotoxin type N (SEN)",
    "sep": "staphylococcal enterotoxin type P (SEP)",
    "seq": "staphylococcal enterotoxin type Q (SEQ)",
    "ser": "staphylococcal enterotoxin type R (SER)",
    "ses": "staphylococcal enterotoxin type S (SES)",
    "set": "staphylococcal enterotoxin type T (SET)",
    "sey": "staphylococcal enterotoxin type Y (SEY)",
    "sel26": "enterotoxin-like type 26 (SEl26)",
    "sel27": "enterotoxin-like type 27 (SEl27)",
    "sel28": "enterotoxin-like type 28 (SEl28, SEl26-family variant)",
    "sel31": "enterotoxin-like type 31 (SEl31, plasmid)",
    "sel32": "enterotoxin-like type 32 (SEl32, plasmid)",
    "sel33": "enterotoxin-like type 33 (SEl33, egc recombinant)",
    "sel30": "enterotoxin-like type 30 (SEl30; identical to SEl01/SE01"
    " of Hisatsune 2017 — the same gene discovered twice)",
    "selu2": "enterotoxin-like U2 (SElU2, often called SElW)",
    "seo": "staphylococcal enterotoxin type O (SEO, egc)",
    "se02": "staphylococcal enterotoxin type 02 (SE02, validated emetic + superantigenic)",
    "selV": "enterotoxin-like type V (SElV)",
    "selX": "enterotoxin-like type X (SElX)",
    "selZ": "enterotoxin-like type Z (SElZ)",
    "tsst1": "toxic shock syndrome toxin 1 (TSST-1)",
    "seu": "staphylococcal enterotoxin U (SEU, egc)",
}

#: Curated additions cut from public-domain INSDC records (committed as
#: script assets, extracted once from the ENA records named in each file's
#: header). The first two fill gene gaps the upstream panel misses — SEU is
#: an egc-core toxin implicated in outbreaks without any classical SE
#: (Umeda et al. 2017), SElW is the egc4 recombinant (Thomas et al. 2006) —
#: the last three are intra-type variant alleles for the polymorphism
#: hotspots Lefebvre et al. (Toxins 14:249, 2022) quantified (sec: 10
#: variants, tsst1: 6). Alleles share the parent gene id and ride on the
#: accession, the cgtA/trh precedent from the other bundled panels.
EXTRAS: tuple[tuple[str, str, str, str], ...] = (
    ("seu_AY205305.fna", "seu", "AY205305:83-868", "staphylococcal enterotoxin U (SEU, egc)"),
    (
        "selu2_EF030428.fna",
        "selu2",
        "EF030428:1-771",
        "enterotoxin-like U2 (SElU2, the gene widely called SElW; egc4)",
    ),
    (
        "sec3_X51661.fna",
        "sec",
        "X51661:36-836",
        "staphylococcal enterotoxin type C-3 allele (SEC3)",
    ),
    (
        "secbov_AF217235.fna",
        "sec",
        "AF217235:c1229-2044",
        "staphylococcal enterotoxin C bovine allele (SaPIbov)",
    ),
    (
        "tst2_AF217235.fna",
        "tsst1",
        "AF217235:4242-4946",
        "toxic shock syndrome toxin 1, bovine SaPI allele",
    ),
    (
        "sei_ay158703.fna",
        "sei",
        "AY158703:82-810",
        "staphylococcal enterotoxin type I, egc3 variant allele (96.7 % to prototype)",
    ),
    (
        "seu_ay158703.fna",
        "seu",
        "AY158703:964-1749",
        "staphylococcal enterotoxin U, egc3 variant allele (97.7 % to prototype)",
    ),
    (
        "sen_ay158703.fna",
        "sen",
        "AY158703:1788-2543",
        "staphylococcal enterotoxin type N, egc3 variant allele (97.5 % to prototype)",
    ),
    (
        "sel28_MN450304.fna",
        "sel28",
        "MN450304:2472-3197",
        "enterotoxin-like type 28 (SEl28 = SEl26-family variant; TTG start)",
    ),
    (
        "se02_AP019712.fna",
        "se02",
        "AP019712:50873-51652",
        "staphylococcal enterotoxin type 02 (SE02; VALIDATED emetic"
        " + superantigenic, Tokyo 2004 outbreak; 95.4 % aa to SElZ)",
    ),
    (
        "sel31_AFTU01000025.fna",
        "sel31",
        "AFTU01000025.1:11135-11905",
        "enterotoxin-like type 31 (SEl31, plasmid cluster with SEl32; Dicks 2021 Gr8)",
    ),
    (
        "sel32_AFTU01000025.fna",
        "sel32",
        "AFTU01000025.1:10294-11025",
        "enterotoxin-like type 32 (SEl32, plasmid cluster with SEl31; Dicks 2021 Gr9)",
    ),
    (
        "sel33_FGMI01000018.fna",
        "sel33",
        "FGMI01000018.1:35991-36743",
        "enterotoxin-like type 33 (SEl33, a selw-sen egc recombinant; Dicks 2021 Gr10)",
    ),
    (
        "seo_AF285760.fna",
        "seo",
        "AF285760:370-1155",
        "staphylococcal enterotoxin type O (SEO, egc; ATT start)",
    ),
    (
        "sel30_PZ361155.fna",
        "sel30",
        "PZ361155:32-784",
        "enterotoxin-like type 30 (SEl30; identical to SEl01/SE01 of"
        " Hisatsune 2017, AP014922 — the same gene discovered twice)",
    ),
)

_CONTROL_GENE = "nuc"
_CONTROL_ACCESSION = "NC_007795.1:800034-800720"
_CONTROL_PRODUCT = "staphylococcal thermonuclease (nuc), species control"


def _read_fasta(text: str, label: str) -> list[tuple[str, str]]:
    """``(header, sequence)`` pairs from a single-line FASTA."""
    records: list[tuple[str, str]] = []
    header = ""
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:].strip(), []
        elif line.strip():
            chunks.append(line.strip())
    if header:
        records.append((header, "".join(chunks)))
    if not records:
        raise SystemExit(f"{label}: no FASTA records")
    return [(header, "".join(chunks).upper()) for header, chunks in records]


def _upstream_gene(header: str) -> tuple[str, str]:
    """``(gene, accession)`` from the upstream abricate header.

    The abricate shape is ``>DB~~~GENE~~~ACCESSION DESCRIPTION``; the accession
    is the first token of the third field (coordinates or a record title may
    follow it).
    """
    parts = header.split("~~~")
    if len(parts) != 3:
        raise SystemExit(f"unparsable abricate header: {header!r}")
    _, gene, rest = parts
    tokens = rest.split()
    accession = tokens[0] if tokens else ""
    if not gene or not accession:
        raise SystemExit(f"empty gene/accession in header: {header!r}")
    return gene, accession


def _read_text(source: str) -> str:
    if source.startswith(("http://", "https://")):
        with urllib.request.urlopen(source, timeout=60) as response:
            return response.read().decode("utf-8")
    return Path(source).read_text(encoding="utf-8")


def _validate(genes: list[str], sequences: list[str], label: str) -> None:
    if len(genes) != len(set(genes)):
        duplicates = sorted({gene for gene in genes if genes.count(gene) > 1})
        raise SystemExit(f"{label}: duplicate gene ids: {duplicates}")
    for gene, sequence in zip(genes, sequences, strict=True):
        if not sequence or set(sequence) - set("ACGTN"):
            raise SystemExit(f"{label}: {gene!r} is not a nucleotide sequence")
        if not 400 <= len(sequence) <= 4000:
            raise SystemExit(f"{label}: {gene!r} length {len(sequence)} outside the expected band")


def _write_sequence(handle, gene: str, accession: str, product: str, sequence: str) -> None:
    seqid = encode_seqid("sau_ent", gene, accession, ("virulence",))
    handle.write(f">{seqid} {product}\n")
    for offset in range(0, len(sequence), 60):
        handle.write(sequence[offset : offset + 60] + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", default=UPSTREAM_URL, help="upstream abricate sequences (path or URL)"
    )
    parser.add_argument(
        "--control",
        type=Path,
        default=Path(__file__).with_name("nuc_NC_007795.1.fna"),
        help="species-control FASTA (default: the committed nuc asset)",
    )
    parser.add_argument("--dest", type=Path, required=True, help="bundle directory to write")
    args = parser.parse_args(argv)

    upstream_text = _read_text(args.source)
    upstream = _read_fasta(upstream_text, args.source)
    records = [_upstream_gene(header) for header, _ in upstream]
    upstream_genes = [gene for gene, _ in records]
    upstream_sequences = [sequence for _, sequence in upstream]
    _validate(upstream_genes, upstream_sequences, "upstream")
    unknown = sorted(set(upstream_genes) - set(PRODUCTS))
    if unknown:
        raise SystemExit(f"no curated product string for gene(s): {unknown}")

    control = _read_fasta(args.control.read_text(encoding="utf-8"), str(args.control))
    if len(control) != 1:
        raise SystemExit(f"{args.control}: expected exactly one record")
    _validate([_CONTROL_GENE], [control[0][1]], "control")

    dest: Path = args.dest
    dest.mkdir(parents=True, exist_ok=True)
    with (dest / "sequences").open("w", encoding="utf-8", newline="\n") as handle:
        for (gene, accession), sequence in zip(records, upstream_sequences, strict=True):
            _write_sequence(handle, gene, accession, PRODUCTS[gene], sequence)
        for asset, gene, accession, product in EXTRAS:
            extra = _read_fasta((args.control.parent / asset).read_text(encoding="utf-8"), asset)
            if len(extra) != 1:
                raise SystemExit(f"{asset}: expected exactly one record")
            _validate([gene], [extra[0][1]], asset)
            _write_sequence(handle, gene, accession, product, extra[0][1])
        _write_sequence(handle, _CONTROL_GENE, _CONTROL_ACCESSION, _CONTROL_PRODUCT, control[0][1])

    (dest / "bundled.json").write_text(
        json.dumps(
            {
                "name": "sau_ent",
                "description": (
                    "S. aureus enterotoxin gene panel: SEA-SEE, SEG-SEY,"
                    " SEl26-28, SEl30-33, SElV/X/Z, SE02, SEU, SEO, SElU2 and TSST-1 + the"
                    " nuc species control, with variant alleles for SEC,"
                    " SEI, SEN, SEU and TSST-1"
                ),
                "vendor": "indexofire/enterotoxin (MIT) + gspec-curated control",
                "dbtype": "nucl",
                "kind": "gene",
                "snapshotted": datetime.now(UTC).strftime("%Y-%m-%d"),
                "license": "MIT (upstream LICENSE, Copyright 2023 Mark Renton)",
                "note": (
                    "Twenty-five enterotoxin, enterotoxin-like and TSST-1"
                    " references converted from the indexofire/enterotoxin"
                    f" abricate database (commit {UPSTREAM_COMMIT},"
                    " 2023-08-27; MIT). Sequences are the upstream records"
                    " byte-for-byte; only the headers are rewritten to"
                    " gspec/v1 and a curated product string is attached per"
                    " gene. Fourteen curated additions cut from"
                    " public-domain INSDC records (each committed as a build"
                    " asset whose FASTA header carries the citation): seu"
                    " from AY205305 (Letertre 2003; implicated in outbreaks"
                    " without classical SEs, Umeda 2017), SElU2 from"
                    " EF030428 (Thomas 2006; widely called SElW after"
                    " Collery & Smyth 2007 — the 2026 review restores"
                    " SElU2 because the 2012 'true SElW' is a different,"
                    " non-emetic protein), seo from AF285760:370-1155"
                    " (Jarraud 2001, the original egc cluster deposit; a"
                    " gap the Fanelli 2022 backtest exposed), the SEC3"
                    " allele from X51661 (Hovde 1990), the SEC-bovine and"
                    " TSST-1 bovine-SaPI alleles from AF217235 (Fitzgerald"
                    " 2001), the egc3 variant alleles for sei/seu/sen from"
                    " AY158703 (96.7-97.7 % to prototypes; Letertre),"
                    " SEl28 from MN450304:2472-3197 (egc SaPI of strain"
                    " 364P, Fanelli deposit), SEl30 from PZ361155, and"
                    " sel31/sel32/sel33 from AFTU01000025.1 and"
                    " FGMI01000018.1 (Dicks 2021, each verified 100 %"
                    " aa-identical to Table S3), and se02 from"
                    " AP019712:50873-51652 (Suzuki 2020 Tokyo 2004"
                    " outbreak — the only novel SE with validated emetic"
                    " (marmoset) and superantigenic activity; 95.4 % aa"
                    " to SElZ, listed as a distinct type by the 2026"
                    " review). SEl01/SE01 (Hisatsune 2017) is 100 %"
                    " aa-identical to sel30 — the same gene discovered"
                    " twice; sel30 already covers it."
                    " NOMENCLATURE AUDIT: the upstream"
                    " gene labels sel26/sel27 (MF370874, 744/753 bp) follow"
                    " the primary designation paper (Zhang et al., IJMM"
                    " 2018: 'designated SEl26 and SEl27') and the 2026"
                    " nomenclature-clarification review (Franjus et al.:"
                    " SEl26-1 744 nt, SEl27-1 753 nt, SEl28 = SEl26-2"
                    " 726 nt); the conflicting sel27/sel28 numbering in the"
                    " MF370874/78 description lines and MN450304's CDS"
                    " labels follows no peer-reviewed paper and is not"
                    " adopted. REMAINING GAPS (deliberate, documented in"
                    " scripts/staphylococcus/README.md): SEC2 (no complete"
                    " INSDC record; SEC1 covers at 97.4 % protein"
                    " identity), sel29p (pseudogene, absent from"
                    " contemporary strains), sel34/sel35"
                    " (bioRxiv preprint, unreviewed), and the remaining"
                    " variant alleles for seg/sek/sem/sey (Lefebvre 2022"
                    " inventory). Variant alleles share the parent gene id"
                    " and are distinguished by the accession. One species"
                    " control added, staphylococcal thermonuclease (nuc)"
                    " from public-domain INSDC"
                    f" {_CONTROL_ACCESSION}. Screen-only; no typing scheme."
                    " If any author or rights holder objects to the use or"
                    " redistribution of this database, please contact"
                    " indexofire@gmail.com and the content will be removed or"
                    " adjusted promptly upon verification."
                ),
                "guidance": {
                    "references": [
                        {
                            "citation": (
                                "Argudin MA, Mendoza MC, Rodicio MR. Food"
                                " poisoning and Staphylococcus aureus"
                                " enterotoxins. Toxins (Basel)."
                                " 2010;2(7):1751-1773."
                            ),
                            "url": "https://doi.org/10.3390/toxins2071751",
                            "note": "The SE/SEl gene family and its nomenclature.",
                        },
                        {
                            "citation": (
                                "Lefebvre D, Blanco-Valle K, Hennekinne JA,"
                                " Simon S, Fenaille F, Becher F, Nia Y."
                                " Multiplex detection of 24 staphylococcal"
                                " enterotoxins in culture supernatant using"
                                " liquid chromatography coupled to"
                                " high-resolution mass spectrometry. Toxins"
                                " (Basel). 2022;14(4):249."
                            ),
                            "url": "https://doi.org/10.3390/toxins14040249",
                            "note": (
                                "The variant inventory per gene (sec 10, seu 10,"
                                " sem 9, sen 9, sek 8, seg 6, sey 6, tsst1 6)"
                                " this panel's variant alleles mitigate."
                            ),
                        },
                    ],
                    "caveats": [
                        (
                            "Gene presence is toxigenic potential, not toxin"
                            " production: expression is strain-, niche- and"
                            " growth-dependent. Food or isolate enterotoxin"
                            " positivity still needs an immunological assay"
                            " (ELISA) or expression evidence."
                        ),
                        (
                            "The 'enterotoxin-like' types (SEl26, SEl27,"
                            " SElV, SElW, SElX, SElZ) are grouped by sequence"
                            " homology; emetic activity is not demonstrated"
                            " for all of them, so reporting them as"
                            " enterotoxins overstates the risk."
                        ),
                        (
                            "Sequence polymorphism is the family norm"
                            " (Lefebvre et al. 2022: six to ten variants per"
                            " gene for sec, seg, sek, sem, sen, seu, sey and"
                            " tsst1). This panel carries variant alleles"
                            " for sec, sei, sen, seu and tsst1; other genes"
                            " remain one reference per type, and divergence"
                            " beyond these alleles"
                            " can still fall below the default floors (80 %"
                            " identity / 80 % coverage) — key on identity"
                            " and coverage, not presence alone."
                        ),
                        (
                            "nuc is the species control: its absence argues"
                            " against S. aureus (or a failed screen), while"
                            " the toxin panel alone does not establish the"
                            " species. S. argenteus carries SElZ-like genes —"
                            " key on the accession, not the gene name alone."
                        ),
                        (
                            "Upstream references are single per-type"
                            " accessions from varied strains; alleles can"
                            " diverge below the default screen floors"
                            " (80 % identity / 80 % coverage)."
                        ),
                    ],
                    "usage": [
                        "gspec screen isolate.fa -d sau_ent --output sae.tsv",
                        "gspec screen isolate.fa -d sau_ent --format json",
                    ],
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(
        f"wrote {dest}/sequences ({len(records) + len(EXTRAS) + 1} records) and {dest}/bundled.json"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
