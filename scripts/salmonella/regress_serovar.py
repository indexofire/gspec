#!/usr/bin/env python3
"""Measure a Salmonella database's serovar accuracy against a truth table.

SISTR ships ``data/genomes-to-serovar.txt`` (one ``assembly_id<TAB>serovar``
line per genome), which makes a real regression set possible: screen every
assembly, designate it, and diff the call against the expected serovar.

Usage::

    python scripts/salmonella/regress_serovar.py \\
        --db sen_sero_sistr --assemblies /path/to/fastas \\
        --expected ~/.mamba/envs/sistr/…/data/genomes-to-serovar.txt \\
        --datadir ~/.local/share/gspec/db

The assembly directory is matched by file stem (``GCA_000006945.2.fna`` ->
``GCA_000006945.2``); unmatched expected rows are reported as "missing", not
counted as errors. Output is a per-serovar confusion table plus the overall
accuracy, so a call can be judged against its alternatives rather than a bare
percentage.

This harness measures; it does not ship data. Point ``--db`` at a locally
imported database (`import_upstream.py … --dest <datadir>/NAME`), never at a
bundled one you did not curate.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

FASTA_SUFFIXES = (".fna", ".fa", ".fasta", ".fas", ".fna.gz", ".fa.gz", ".fasta.gz")


def _expected_map(path: Path) -> dict[str, str]:
    """``assembly_id -> serovar`` (the SISTR truth table's two columns)."""
    expected: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2 or not parts[0] or parts[0].startswith("#"):
            continue
        expected[parts[0].strip()] = parts[1].strip()
    return expected


def _assemblies(root: Path, expected: dict[str, str]) -> list[tuple[str, Path]]:
    """The files whose stem is in the truth table, in truth-table order."""
    found: dict[str, Path] = {}
    for path in sorted(root.rglob("*")):
        if not path.is_file() or not path.name.lower().endswith(FASTA_SUFFIXES):
            continue
        stem = path.name
        for suffix in FASTA_SUFFIXES:
            if stem.lower().endswith(suffix):
                stem = stem[: -len(suffix)]
                break
        found.setdefault(stem, path)
    return [(name, found[name]) for name in expected if name in found]


def _designate(
    database: str, datadir: Path | None, files: list[Path]
) -> dict[str, tuple[str, tuple[str, ...]]]:
    """``path str -> phenotype`` for one batch (screen to a table, then type)."""
    from gspec.engines import OutputFormat
    from gspec.screening.contigs import run_screen
    from gspec.typing.ops import run_typing

    table = run_screen(
        files,
        database,
        datadir,
        80.0,
        80.0,
        1,
        1,
        None,
        True,  # quiet
        False,  # noheader
        False,  # nopath: the table must carry the path we key the result by
        False,
        OutputFormat.tsv,
    )
    rendered = run_typing([], datadir, OutputFormat.json, True, stdin_text=table)
    document = json.loads(rendered)
    calls: dict[str, tuple[str, tuple[str, ...]]] = {}
    for entry in document.get("files", []):
        phenotypes = entry.get("phenotypes") or {}
        call = phenotypes.get("serotype") or {}
        candidates = tuple(
            str(candidate.get("phenotype"))
            for candidate in (call.get("ambiguous") or [])
            if isinstance(candidate, dict) and candidate.get("phenotype")
        )
        calls[str(entry.get("file"))] = (str(call.get("phenotype") or ""), candidates)
    return calls


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, help="database name to screen against")
    parser.add_argument("--assemblies", type=Path, required=True)
    parser.add_argument("--expected", type=Path, required=True)
    parser.add_argument("--datadir", type=Path)
    parser.add_argument(
        "--min-accuracy",
        type=float,
        default=0.0,
        help="exit non-zero when accuracy falls below this fraction",
    )
    args = parser.parse_args(argv)

    expected = _expected_map(args.expected)
    pairs = _assemblies(args.assemblies, expected)
    if not pairs:
        raise SystemExit(f"no assembly in {args.assemblies} matches a stem in {args.expected}")
    calls = _designate(args.db, args.datadir, [path for _, path in pairs])

    hits = 0
    listed = 0
    confusion: Counter[tuple[str, str]] = Counter()
    missing: list[str] = []
    for name, path in pairs:
        phenotype, candidates = calls.get(str(path), ("", ()))
        truth = expected[name]
        if phenotype == truth:
            hits += 1
        elif truth in candidates:
            # Ambiguous but right: the upstream tools pipe-join these names, so a
            # truth inside the candidate set is information actually delivered.
            listed += 1
        elif not phenotype and not candidates:
            missing.append(name)
        else:
            confusion[(truth, phenotype or "|".join(candidates))] += 1

    total = len(pairs)
    print(
        f"assemblies: {total}  exact: {hits}  ambiguous-but-listed: {listed}"
        f"  wrong: {total - hits - listed - len(missing)}  no-call: {len(missing)}"
        f"  | exact accuracy: {100 * hits / total:.1f}%"
        f"  | contains-truth: {100 * (hits + listed) / total:.1f}%"
    )
    if confusion:
        print("\nconfusions (expected -> called):")
        for (truth, call), count in confusion.most_common(20):
            print(f"  {count:5d}  {truth} -> {call}")
    if missing:
        print(f"\nno-call examples: {', '.join(missing[:10])}")
    accuracy = (hits + listed) / total
    if accuracy < args.min_accuracy:
        print(f"FAIL: accuracy {accuracy:.3f} < --min-accuracy {args.min_accuracy}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
