#!/usr/bin/env python3
"""Build the whole ``~/data/sen`` serotype-evidence dataset, in order.

This is the one entry point. Each step is also runnable on its own (its own
``--help`` explains the options), but the steps feed each other, so the order
matters:

1. **`fetch_serotype_evidence.py`** — the two public layers plus the curated
   tie-pair and literature-lead tables:
   ``reference_cultures.tsv`` (strains deposited in a public collection),
   ``ncbi_refseq_serovars.tsv`` (every RefSeq assembly, submitter labels),
   ``tie_serovars_candidates.tsv`` (assemblies of the serovars the panel cannot
   separate — the serovar list comes from `naming_priority.json`),
   ``brikama_leads.tsv`` (rendered from the curated `brikama_leads.json`).
2. **`harvest_published_serotypes.py`** — mines open-access papers for
   accession+serovar tables into ``published/``. Incremental: the per-query
   search index and the per-paper row caches mean a re-run only touches papers
   it has not seen.
3. **`merge_compare_serotypes.py merge`** — joins the layers by accession,
   keeps the strongest label and writes ``merged_evidence.tsv`` +
   ``label_conflicts.tsv``.

Then, optionally, ``compare`` answers a marker question for a serovar pair and
``--prune`` drops the cached supplements that produced no table.

Every artifact under the output directory is produced by one of these three
scripts — nothing is hand-edited into the dataset.

Usage::

    # the full pipeline (network, ~20 min when the caches are warm)
    python scripts/salmonella/build_sen_dataset.py

    # add a serovar to the paper search, and refresh everything downstream
    python scripts/salmonella/build_sen_dataset.py --serovar Brikama

    # ask whether a marker splits two serovar groups (needs at least two --group)
    python scripts/salmonella/build_sen_dataset.py \\
        --group "factor6=Goldcoast,Bovismorbificans" \\
        --group "no_factor6=Albany,Hindmarsh"
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _run(script: str, *arguments: str) -> None:
    command = [sys.executable, str(HERE / script), *arguments]
    print(f"\n$ {' '.join(command)}", file=sys.stderr)
    result = subprocess.run(command, check=False)
    if result.returncode != 0:
        raise SystemExit(f"{script} failed with exit code {result.returncode}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path.home() / "data" / "sen")
    parser.add_argument(
        "--skip-refseq",
        action="store_true",
        help="reuse the existing ncbi_refseq_serovars.tsv instead of re-sweeping 19k assemblies",
    )
    parser.add_argument(
        "--skip-papers",
        action="store_true",
        help="do not re-run the paper harvest (keeps whatever published/ holds)",
    )
    parser.add_argument(
        "--serovar",
        action="append",
        default=[],
        help="add a serovar keyword to the paper search (repeatable)",
    )
    parser.add_argument("--max-papers", type=int, default=400)
    parser.add_argument("--max-serovar-papers", type=int, default=2000)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--min-layer", default="published")
    parser.add_argument("--per-serovar", type=int, default=4)
    parser.add_argument(
        "--group",
        action="append",
        default=[],
        metavar="NAME=serovar,serovar",
        help="also run the marker comparison for these groups (needs two or more)",
    )
    parser.add_argument(
        "--prune",
        action="store_true",
        help="afterwards, delete cached supplements whose paper produced no rows",
    )
    args = parser.parse_args(argv)

    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)

    fetch_arguments = ["--out", str(out)]
    if args.skip_refseq:
        fetch_arguments.append("--skip-refseq")
    _run("fetch_serotype_evidence.py", *fetch_arguments)

    if not args.skip_papers:
        harvest_arguments = [
            "--out",
            str(out / "published"),
            "--max-papers",
            str(args.max_papers),
            "--max-serovar-papers",
            str(args.max_serovar_papers),
            "--jobs",
            str(args.jobs),
        ]
        for serovar in args.serovar:
            harvest_arguments += ["--serovar", serovar]
        _run("harvest_published_serotypes.py", *harvest_arguments)

    _run("merge_compare_serotypes.py", "merge", "--root", str(out))

    if len(args.group) >= 2:
        compare_arguments = ["compare", "--root", str(out), "--min-layer", args.min_layer]
        compare_arguments += ["--per-serovar", str(args.per_serovar)]
        for group in args.group:
            compare_arguments += ["--group", group]
        _run("merge_compare_serotypes.py", *compare_arguments)

    if args.prune:
        _run(
            "harvest_published_serotypes.py",
            "--out",
            str(out / "published"),
            "--prune-supplements",
        )

    print("\ndataset complete:", file=sys.stderr)
    for path in sorted(out.rglob("*")):
        if path.is_file() and path.suffix in {".tsv", ".md", ".json"}:
            print(f"  {path.relative_to(out)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
