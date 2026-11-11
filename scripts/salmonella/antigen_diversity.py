#!/usr/bin/env python3
"""Where does the antigen panel lack reference diversity? A corpus-driven report.

The way to raise the serotype call rate is not more inference but better
references, and the corpus says where they are missing. For every genome and every
antigen family this script records the best identity reached by any reference and
reports, per antigen:

* how many genomes matched it perfectly (≥ 99 %), adequately (95-99 %), poorly
  (< 95 %), or not at all;
* the genomes whose best hit for their called antigen is below
  ``--weak-below`` (their allele is not in the panel);
* for a family with no call at all, the nearest antigen and how far away it is
  (a missing serogroup looks like this — `26SAL0087` is an O:8 genome with no O
  call).

The output is the worklist for `build_h_diversity.py` / `build_o_discriminators.py`:
"antigen X has N genomes at ≤ 96 % — extract N alleles".

The worklist is built from what each genome SHOULD carry, not from its best hit:
with a SISTR truth dump the expected O:H1:H2 tokens are known, so a genome whose
expected `h2 1,2` only reaches 96.8 % is a missing allele, while a monophasic
genome's weak `fljB` cross-reaction is correctly ignored (its expected h2 is `-`).

Usage::

    gspec screen -d sen_sero_sistr assemblies/*.fna -o table.tsv --format tsv
    sistr --no-cgmlst -f tab -o sistr.tab assemblies/*.fna
    python scripts/salmonella/antigen_diversity.py \\
        --table table.tsv --sistr-tab sistr.tab --db-dir ~/.local/share/gspec/db/sen_sero_sistr \\
        [--weak-below 98] [--out report.json]
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from gspec.tables import detect_separator, iter_rows, read_table_text

#: Identity at or above which an allele is considered represented in the panel.
STRONG = 99.0
#: Identity below which a match is too divergent to be the genome's own allele.
WEAK = 95.0


def _rows(path: Path) -> list[tuple[str, str, float, float]]:
    """``(file, gene, identity, coverage)`` from a gspec screen table."""
    text = read_table_text(path, code="INPUT_NOT_FOUND", noun="screen table")
    out: list[tuple[str, str, float, float]] = []
    for _, cells in iter_rows(text, detect_separator(text)):
        if not cells or cells[0].startswith("#") or cells[1] == "SEQUENCE":
            continue
        out.append((cells[0], cells[5], float(cells[10]), float(cells[9])))
    return out


def _antigen(gene: str) -> tuple[str, str]:
    """``h1_e_h__3`` -> ``("h1", "e,h")`` (the panel's genome ids are already decoded)."""
    name = gene.rsplit("__", 1)[0]
    family, _, antigen = name.partition("_")
    return family, antigen.replace("_", ",")


#: The formula a SISTR row carries, per family, in FIXED order O:H1:H2 (the
#: antigenic-formula column is the only place the expected tokens are stated).
_FAMILIES = ("o", "h1", "h2")


def _plain(token: str) -> str:
    """Drop WKL bracket/parenthesis notation: gene ids never carry either.

    The panel's records are keyed `h1_g__s__t` for the table's `g,[s],t`, so a
    comparison on the raw token would call every bracketed row a missing reference.
    """
    return token.replace("[", "").replace("]", "").replace("(", "").replace(")", "")


def _expected_antigens(path: Path) -> dict[str, dict[str, str]]:
    """``stem -> {family: antigen token}`` from a SISTR ``-f tab`` dump.

    Tokens carry WKL brackets (`[5]`, `(g)`); the panel's names keep them, so the
    value is used verbatim. `-` means the phase is absent for that genome, and
    families with no token are simply missing.
    """
    expected: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            stem = Path(row.get("fasta_filepath") or "").stem
            formula = (row.get("antigenic_formula") or "").strip()
            if not stem or formula.count(":") != 2:
                continue
            parts = formula.split(":")
            expected[stem] = dict(zip(_FAMILIES, (part.strip() for part in parts), strict=True))
    return expected


def _normalise_o(expected: dict[str, dict[str, str]], who_table: Path) -> dict[str, dict[str, str]]:
    """Rewrite each expected O token as the serogroup letter the panel uses.

    SISTR states the O side as an antigen list (`1,4,[5],12`) on some rows and as a
    serogroup letter on others; the panel keys on the letter, so without this every
    O row would look like a missing reference while matching at 100 %.
    """
    mapping: dict[str, str] = {}
    with who_table.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            antigen = (row.get("O_antigen") or "").strip()
            serogroup = (row.get("Serogroup") or "").strip()
            if antigen and serogroup:
                mapping.setdefault(antigen, serogroup)
    for families in expected.values():
        token = families.get("o")
        if token and token in mapping:
            families["o"] = mapping[token]
    return expected


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--sistr-tab", type=Path, required=True)
    parser.add_argument("--db-dir", type=Path, required=True)
    parser.add_argument(
        "--who-table",
        type=Path,
        help=(
            "SISTR's WHO-2007 CSV; maps an O antigen list (1,4,[5],12) to the"
            " serogroup letter the panel keys on (B). Without it every O row reads"
            " as 'absent' because the two notations differ"
        ),
    )
    parser.add_argument("--weak-below", type=float, default=98.0)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)

    best: dict[tuple[str, str], dict[str, tuple[float, float]]] = defaultdict(dict)
    for file, gene, identity, coverage in _rows(args.table):
        family, antigen = _antigen(gene)
        current = best[(file, family)].get(antigen)
        if current is None or identity > current[0]:
            best[(file, family)][antigen] = (identity, coverage)

    expected = _expected_antigens(args.sistr_tab)
    if args.who_table is not None:
        expected = _normalise_o(expected, args.who_table)
    declared = json.loads((args.db_dir / "serotype.json").read_text(encoding="utf-8"))
    antigens = {
        family["key"]: {entry["name"]: len(entry["genes"]) for entry in family["entries"]}
        for family in declared["antigens"]
    }

    per_antigen: dict[str, dict[str, dict[str, int]]] = defaultdict(lambda: defaultdict(dict))
    weak_calls: list[dict[str, object]] = []
    absent: list[dict[str, object]] = []
    missing_family: dict[str, int] = defaultdict(int)
    for (file, family), scored in sorted(best.items()):
        if family not in antigens:
            continue
        want = expected.get(Path(file).stem, {}).get(family)
        if want in (None, "-"):
            continue  # the genome is not expected to carry this phase
        if not scored:
            missing_family[family] += 1
            absent.append({"file": Path(file).name, "family": family, "expected": want})
            continue
        matches = [name for name in scored if _plain(name) == _plain(want)]
        if not matches:
            # the expected token is not even a row: the reference is missing or the
            # call would have to widen to reach it
            nearest, (nearest_identity, _) = max(scored.items(), key=lambda item: item[1][0])
            absent.append(
                {
                    "file": Path(file).name,
                    "family": family,
                    "expected": want,
                    "nearest": nearest,
                    "nearest_identity": round(nearest_identity, 2),
                }
            )
            continue
        identity, _coverage = scored[matches[0]]
        bucket = "strong" if identity >= STRONG else "adequate" if identity >= WEAK else "poor"
        per_antigen[family][want][bucket] = per_antigen[family][want].get(bucket, 0) + 1
        if identity < args.weak_below:
            weak_calls.append(
                {
                    "file": Path(file).name,
                    "family": family,
                    "expected": want,
                    "identity": round(identity, 2),
                    "alleles_in_panel": antigens[family].get(want, 0),
                }
            )

    print(f"families with no hit at all: {dict(missing_family) or 'none'}")
    print(f"expected antigens absent from the panel's rows: {len(absent)}")
    for entry in absent[:10]:
        print(f"  {entry}")
    print("\nexpected antigens matched weakly (worklist for --extra-refs / new references):")
    worklist = sorted(weak_calls, key=lambda entry: float(entry["identity"]))
    for entry in worklist[:20]:
        print(
            f"  {entry['file']:16s} {entry['family']:2s} {entry['expected']:14s}"
            f" best identity {entry['identity']:6.2f} %  ({entry['alleles_in_panel']} allele(s))"
        )
    if not worklist:
        print("  (none: every expected antigen matched a reference at or above the floor)")

    if args.out is not None:
        args.out.write_text(
            json.dumps(
                {
                    "weak_below": args.weak_below,
                    "missing_family": dict(missing_family),
                    "absent": absent,
                    "per_antigen": {family: dict(scored) for family, scored in per_antigen.items()},
                    "worklist": worklist,
                },
                indent=1,
            )
            + "\n",
            encoding="utf-8",
        )
        print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
