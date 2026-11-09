#!/usr/bin/env python3
"""Pairwise identity of the ``sen_sero`` markers — the threshold evidence.

The panel's rule floors come from measured separations, not from copied
parameters: a marker pair that is >97% identical cannot separate serovars at a
97% threshold, and that is how `rfc`/`rfbS` were ruled out as O-group markers.

Usage::

    python scripts/salmonella/identity_matrix.py [--markers TSV] [--cache DIR]

Requires the CDS export cache left by ``build_sen_sero.py`` (delete the cache
and run that script first, or pass ``--cache``). Prints a symmetric matrix of
percent identity plus the closest non-self pair per marker, which is the number
the rule floors must clear.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _read_markers(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        rows = [line for line in handle if not line.lstrip().startswith("#")]
    return [row for row in csv.DictReader(rows, delimiter="\t") if row.get("marker_id")]


def _cds_index(cache: Path, genome: str) -> list[dict[str, str]]:
    """One entry per CDS: ``protein_id``/``locus_tag``/``gene`` plus the sequence."""
    path = cache / f"{genome}_cds.fna"
    if not path.is_file():
        raise SystemExit(f"missing {path}: run build_sen_sero.py once to populate the cache")
    entries: list[dict[str, str]] = []
    key: str | None = None
    chunks: list[str] = []

    def field(line: str, name: str) -> str:
        match = re.search(rf"\[{name}=([^\]]*)\]", line)
        return match.group(1) if match else ""

    def flush() -> None:
        if key is None:
            return
        entries.append(
            {
                "protein_id": field(key, "protein_id"),
                "locus_tag": field(key, "locus_tag"),
                "gene": field(key, "gene"),
                "sequence": "".join(chunks),
            }
        )

    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            flush()
            key, chunks = line, []
        else:
            chunks.append(line.strip())
    flush()
    return entries


def _resolve(entries: list[dict[str, str]], selector: str, *, marker: str) -> dict[str, str]:
    field_name, _, value = selector.partition("=")
    hits = [entry for entry in entries if entry.get(field_name) == value]
    if len(hits) != 1:
        raise SystemExit(f"{marker}: {selector} resolved to {len(hits)} CDS")
    return hits[0]


def _from_cache(markers: list[dict[str, str]], cache: Path) -> dict[str, str]:
    """Marker sequences resolved through the cached CDS exports (build-time view)."""
    resolved: dict[str, str] = {}
    for row in markers:
        entries = _cds_index(cache, row["genome"])
        resolved[row["marker_id"]] = _resolve(entries, row["selector"], marker=row["marker_id"])[
            "sequence"
        ]
    return resolved


def _from_bundle(bundle: Path) -> dict[str, str]:
    """Marker sequences read from the built bundle (offline, what ships)."""
    sequences: dict[str, str] = {}
    gene: str | None = None
    chunks: list[str] = []
    for line in (bundle / "sequences").read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if gene is not None:
                sequences[gene] = "".join(chunks)
            match = re.search(r"\|gene=([^|]*)", line)
            gene, chunks = (match.group(1) if match else None), []
        else:
            chunks.append(line.strip())
    if gene is not None:
        sequences[gene] = "".join(chunks)
    return sequences


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--markers", type=Path, default=HERE / "sen_sero_markers.tsv")
    parser.add_argument("--cache", type=Path, default=HERE / ".cache")
    parser.add_argument(
        "--bundle",
        type=Path,
        default=HERE.parents[1] / "src/gspec/data/dbs/sen_sero",
        help="Measure the built bundle instead of fetching (works offline)",
    )
    args = parser.parse_args(argv)

    resolved = (
        _from_bundle(args.bundle)
        if args.bundle.is_dir() and not args.cache.is_dir()
        else _from_cache(_read_markers(args.markers), args.cache)
    )

    names = list(resolved)
    width = max(len(name) for name in names) + 2
    print("".rjust(width) + "".join(f"{name[:11]:>12s}" for name in names))
    closest: dict[str, tuple[float, str]] = {}
    for left in names:
        row = f"{left:<{width}s}"
        for right in names:
            if left == right:
                row += f"{100.0:12.1f}"
                continue
            a, b = resolved[left], resolved[right]
            if min(len(a), len(b)) < 0.6 * max(len(a), len(b)):
                row += f"{'--':>12s}"
                continue
            identity = 100.0 * difflib.SequenceMatcher(None, a, b, autojunk=False).ratio()
            row += f"{identity:12.1f}"
            if left not in closest or identity > closest[left][0]:
                closest[left] = (identity, right)
        print(row)

    print("\nclosest non-self pair (the rule floor must sit above it):")
    for name in names:
        if name in closest:
            identity, other = closest[name]
            print(f"  {name:<10s} {identity:6.1f}%  vs {other}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
