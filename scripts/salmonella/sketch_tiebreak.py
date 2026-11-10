#!/usr/bin/env python3
"""Can k-mer sketches settle antigen calls that alignment ranks too closely?

A prototype measurement, not shipped logic. SISTR leans on cgMLST (330 core
genes) to separate near-identical antigen evidence; gspec does not do allele
calling, but it already owns k-mer machinery. This script measures three ways of
asking "which antigen does this genome carry", per family, against the antigen
tokens in SISTR's own ``antigenic_formula`` column:

1. ``align``    — identity x coverage from a gspec screen table (what ships),
2. ``genome``   — the whole genome's k-mers against each antigen's k-mer set,
3. ``anchored`` — only the k-mers inside the aligned window, per ALLELE.

Usage::

    gspec screen -d sen_sero_sistr assemblies/*.fna -o table.tsv      # any batch
    python scripts/salmonella/sketch_tiebreak.py \\
        --db-dir ~/.local/share/gspec/db/sen_sero_sistr \\
        --sistr-tab sistr_all.tab --screen-table table.tsv [--limit N]

The measured answer (304 in-house genomes) is the reason ``anchored`` exists:
genome-wide containment is dominated by how large an antigen's reference pool
happens to be, while the anchored window asks the question that matters.
"""

from __future__ import annotations

import argparse
import csv
import struct
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from gspec.cluster.kmerdb import (
    count_fasta_kmers,
    count_sequence_kmers,
    load_kmer_set,
    read_kmerdb,
    set_hit_fraction,
)
from gspec.fasta import iter_fasta
from gspec.tables import detect_separator, iter_rows, read_table_text

FAMILIES = ("o", "h1", "h2")


def _sistr_antigens(path: Path) -> dict[str, tuple[str, str, str]]:
    """``stem -> (o, h1, h2)`` from SISTR's antigenic-formula column."""
    truth: dict[str, tuple[str, str, str]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            formula = (row.get("antigenic_formula") or "").strip()
            stem = Path(row.get("fasta_filepath") or "").stem
            if not stem:
                continue
            parts = formula.split(":")
            truth[stem] = (parts[0], parts[1], parts[2]) if len(parts) == 3 else ("", "", "")
    return truth


def _antigen_of(gene: str) -> str:
    """``h1_e_h__3`` or a gspec header -> ``h1_e,h`` (family_preserved)."""
    if "|" in gene:
        gene = next(part.split("=", 1)[1] for part in gene.split("|") if part.startswith("gene="))
    antigen = gene.rsplit("__", 1)[0]
    family, _, name = antigen.partition("_")
    return f"{family}_{name.replace('_', ',')}"


def _screen_rows(path: Path) -> dict[str, list[list[str]]]:
    """``file -> rows`` from a (batch) gspec screen table."""
    text = read_table_text(path, code="INPUT_NOT_FOUND", noun="screen table")
    rows: dict[str, list[list[str]]] = defaultdict(list)
    for _, cells in iter_rows(text, detect_separator(text)):
        if not cells or cells[0].startswith("#") or cells[1] == "SEQUENCE":
            continue
        rows[cells[0]].append(cells)
    return rows


def _allele_sets(db_dir: Path) -> dict[str, list[tuple[str, set[int]]]]:
    """``family -> [(antigen, k-mers)]`` one entry per reference ALLELE."""
    antigen_of: dict[str, str] = {}
    for record in iter_fasta(db_dir / "sequences"):
        gene = next(
            part.split("=", 1)[1] for part in record.id.split("|") if part.startswith("gene=")
        )
        antigen_of[record.id] = _antigen_of(gene)
    out: dict[str, list[tuple[str, set[int]]]] = defaultdict(list)
    for path in (db_dir / "kmers").glob("*.full.u64"):
        raw = path.read_bytes()
        gene = path.name[: -len(".full.u64")]
        family, _, name = antigen_of.get(gene, "__").partition("_")
        if family:
            out[family].append((name, set(struct.unpack(f"<{len(raw) // 8}Q", raw))))
    return out


def _window_kmers(
    contigs: dict[str, str], sequence_field: str, start: int, end: int, k: int
) -> set[int]:
    query: set[int] = set()
    for contig in sequence_field.split(","):
        sequence = contigs.get(contig, "")
        if sequence:
            query |= count_sequence_kmers(sequence[max(0, start - 60) : end + 60], k)
    return query


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-dir", type=Path, required=True)
    parser.add_argument("--sistr-tab", type=Path, required=True)
    parser.add_argument("--screen-table", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args(argv)

    document = read_kmerdb(args.db_dir)
    if document is None:
        raise SystemExit(f"{args.db_dir}: no kmer.json")
    k = document.k
    differential = {name: load_kmer_set(args.db_dir, name) for name in document.sets}
    alleles = _allele_sets(args.db_dir)
    truth = _sistr_antigens(args.sistr_tab)
    rows_by_file = _screen_rows(args.screen_table)

    tally: Counter[str] = Counter()
    disagreements: list[str] = []
    for index, (file_path, rows) in enumerate(sorted(rows_by_file.items())):
        if args.limit and index >= args.limit:
            break
        genome_path = Path(file_path)
        expected = truth.get(genome_path.stem)
        if expected is None:
            continue
        contigs = {record.id.split()[0]: record.sequence for record in iter_fasta(genome_path)}
        genome = count_fasta_kmers(genome_path, k)
        align: dict[str, tuple[float, str]] = {}
        best_windows: dict[str, list[tuple[float, set[int]]]] = defaultdict(list)
        for cells in rows:
            antigen = _antigen_of(cells[5])
            family, _, name = antigen.partition("_")
            if family not in FAMILIES:
                continue
            score = float(cells[9]) / 100.0 * float(cells[10]) / 100.0
            if score > align.get(family, (0.0, ""))[0]:
                align[family] = (score, name)
            best_windows[family].append(
                (score, _window_kmers(contigs, cells[1], int(cells[2]), int(cells[3]), k))
            )
        # Only the best few windows per family: the table holds hundreds of rows
        # per genome, and unioning them all would sweep up the whole chromosome
        # (and every paralogue) instead of the locus in question.
        windows: dict[str, set[int]] = {
            family: set().union(*(w for _, w in sorted(entries, reverse=True)[:5]))
            for family, entries in best_windows.items()
        }
        genome_best: dict[str, tuple[float, str]] = {}
        anchored_best: dict[str, tuple[float, str]] = {}
        for family in FAMILIES:
            for set_name in document.sets:
                if set_name.split("_", 1)[0] != family:
                    continue
                fraction = set_hit_fraction(genome, differential[set_name])
                if fraction > genome_best.get(family, (0.0, ""))[0]:
                    genome_best[family] = (fraction, set_name.split("_", 1)[1])
            for name, codes in alleles.get(family, ()):
                if not codes or not windows.get(family):
                    continue
                overlap = len(windows[family] & codes) / len(codes)
                if overlap > anchored_best.get(family, (0.0, ""))[0]:
                    anchored_best[family] = (overlap, name)
        for position, family in enumerate(FAMILIES):
            want = expected[position]
            if not want or want == "-":
                continue
            if family == "o":
                # SISTR's O column mixes serogroup letters (B) with O antigen
                # numbers (1,4,[5],12) for the same group, so an O token
                # comparison measures notation, not agreement. Score the phases.
                continue
            tally["rows"] += 1
            for label, table in (
                ("align", align),
                ("genome", genome_best),
                ("anchored", anchored_best),
            ):
                hit = table.get(family, (0.0, ""))[1] == want
                tally[label] += hit
            if anchored_best.get(family, (0.0, ""))[1] != want and len(disagreements) < 8:
                disagreements.append(
                    f"{genome_path.stem} {family}: sistr={want}"
                    f" align={align.get(family, (0, '-'))[1]}"
                    f" anchored={anchored_best.get(family, (0, '-'))[1]}"
                )

    rows_total = tally["rows"] or 1
    print(f"antigen tokens compared: {rows_total}")
    for label, title in (
        ("align", "alignment (identity x coverage, ships today)"),
        ("genome", "whole-genome k-mers vs each antigen's set"),
        ("anchored", "anchored window k-mers, per allele"),
    ):
        print(f"  {title:46s} {100 * tally[label] / rows_total:5.1f}%")
    if disagreements:
        print("\nexamples where the anchored variant missed:")
        for line in disagreements:
            print("  " + line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
