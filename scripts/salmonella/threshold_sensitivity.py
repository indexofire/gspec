#!/usr/bin/env python3
"""How close does a batch sit to each decision threshold? (Are the knobs safe?)

The panel's numbers were set against real data, so the honest question is whether a
*different* batch invalidates them. This script answers it per batch by reporting the
observed distribution against each threshold, instead of trusting the calibration.

The parameters fall into three kinds, and only the middle kind is at risk:

**Definitional** (a fact, not a fit) — they encode what a word means:
`SAME_LOCUS` (0.5: same physical gene), `PHASE_END_BASES` (150 bp ends),
`_PHASE_ALTERNATIVE` (0.98), `MIN_SKETCH`/`SKETCH_MARGIN`, the k-mer k (27 reads,
21 sketches), the union/`_distinct_loci` overlap rule. A new batch cannot invalidate
them; only a change in what we mean would.

**Reporting** — `margin` (0.01) decides whether a family is reported as one name or a
candidate list, and widening never changes *which* row ranks first. It is the most
batch-sensitive number here (measured: 77 % of multi-antigen families have a top-two
gap inside 3x margin) but by construction it cannot turn a right call into a wrong
one, only a decided call into a listed one.

**Acceptance** — `identity_floor` (80), `coverage_floor` (60),
`partial_coverage_floor` (45), `partial_identity_floor` (95), `phase_floor` (0.75).
These decide whether a locus counts at all, so they CAN change a call. They are safe
when the observed values leave an empty band around them: measured on the in-house
corpus, phase scores are either ≤ 0.691 (cross-reaction) or ≥ 0.936 (present), so
`phase_floor` sits in a 0.245-wide desert and would survive a large shift in assembly
quality. Coverage and identity floors are denser at their boundaries (5-8 % of hits
within ±2-5 pp), so they are the first to re-check on a new batch.

Usage::

    gspec screen -d sen_sero_sistr assemblies/*.fna -o table.tsv --format tsv
    python scripts/salmonella/threshold_sensitivity.py \\
        --table table.tsv --db-dir ~/.local/share/gspec/db/sen_sero_sistr [--out report.json]
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from gspec.tables import detect_separator, iter_rows, read_table_text

#: Observations within this distance of a threshold are called "on the edge".
SCORE_BAND = 0.05
IDENTITY_BAND = 2.0
COVERAGE_BAND = 5.0


def _best(path: Path) -> dict[tuple[str, str], dict[str, tuple[float, float]]]:
    """``(genome, family) -> antigen -> (identity, coverage)`` best hit per antigen."""
    text = read_table_text(path, code="INPUT_NOT_FOUND", noun="screen table")
    best: dict[tuple[str, str], dict[str, tuple[float, float]]] = defaultdict(dict)
    for _, cells in iter_rows(text, detect_separator(text)):
        if not cells or cells[0].startswith("#") or cells[1] == "SEQUENCE":
            continue
        root = cells[5].rsplit("__", 1)[0]
        family = root.partition("_")[0]
        antigen = root.partition("_")[2].replace("_", ",")
        identity, coverage = float(cells[10]), float(cells[9])
        current = best[(Path(cells[0]).stem, family)].get(antigen)
        if current is None or identity * coverage > current[0] * current[1]:
            best[(Path(cells[0]).stem, family)][antigen] = (identity, coverage)
    return best


def antigen_score(identity: float, coverage: float) -> float:
    return identity * coverage / 10_000.0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--db-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)

    document = json.loads((args.db_dir / "serotype.json").read_text(encoding="utf-8"))
    best = _best(args.table)
    genomes = len({stem for stem, _ in best})
    print(f"batch: {genomes} genome(s), {sum(len(v) for v in best.values())} antigen/row(s)\n")

    report: dict[str, object] = {"genomes": genomes}
    phases = [
        max(antigen_score(i, c) for i, c in scored.values())
        for (_, family), scored in best.items()
        if family in ("h1", "h2")
    ]
    floor = float(document.get("phase_floor", 0.75))
    below = [value for value in phases if value < floor]
    above = [value for value in phases if value >= floor]
    edge = [value for value in phases if abs(value - floor) <= SCORE_BAND]
    print(f"phase_floor = {floor}: {len(phases)} phase family/families")
    print(f"  below threshold: {len(below)} (max {max(below):.3f})" if below else "  below: none")
    print(f"  above threshold: {len(above)} (min {min(above):.3f})" if above else "  above: none")
    gap = (min(above) - max(below)) if below and above else None
    print(f"  empty band around the threshold: {gap:.3f}" if gap is not None else "  no band")
    print(
        f"  observations within ±{SCORE_BAND}: {len(edge)}"
        f" ({'re-check on this batch' if edge else 'threshold is not load-bearing here'})"
    )
    report["phase_floor"] = {
        "value": floor,
        "below_max": max(below) if below else None,
        "above_min": min(above) if above else None,
        "edge": len(edge),
    }

    margin = float(document.get("margin", 0.01))
    gaps: list[float] = []
    for scored in best.values():
        if len(scored) < 2:
            continue
        values = sorted((antigen_score(i, c) for i, c in scored.values()), reverse=True)
        gaps.append(values[0] - values[1])
    near = [value for value in gaps if value <= 3 * margin]
    print(f"\nmargin = {margin} (reporting only: widening adds candidates, never reorders)")
    print(f"  families with >1 antigen: {len(gaps)}")
    print(
        f"  top-two gap within 3x margin: {len(near)}"
        f" ({100 * len(near) / max(1, len(gaps)):.1f} %) — expected, this is a"
    )
    print("  panel of near-identical references; it changes the call FORM, not the row")
    report["margin"] = {"value": margin, "families": len(gaps), "near": len(near)}

    identity_floor = float(document["identity_floor"])
    coverage_floor = float(document["coverage_floor"])
    hits = [pair for scored in best.values() for pair in scored.values()]
    ident_edge = [pair for pair in hits if abs(pair[0] - identity_floor) <= IDENTITY_BAND]
    cov_edge = [pair for pair in hits if abs(pair[1] - coverage_floor) <= COVERAGE_BAND]
    print(
        f"\nidentity_floor = {identity_floor}: {len(ident_edge)}/{len(hits)} hit(s)"
        f" within ±{IDENTITY_BAND} pp ({100 * len(ident_edge) / max(1, len(hits)):.1f} %)"
    )
    print(
        f"coverage_floor = {coverage_floor}: {len(cov_edge)}/{len(hits)} hit(s)"
        f" within ±{COVERAGE_BAND} pp ({100 * len(cov_edge) / max(1, len(hits)):.1f} %)"
    )
    print("  (acceptance thresholds: denser boundaries mean re-check these first)")
    report["identity_floor"] = {"value": identity_floor, "edge": len(ident_edge)}
    report["coverage_floor"] = {"value": coverage_floor, "edge": len(cov_edge)}

    partial_coverage = float(document.get("partial_coverage_floor", 45.0))
    partial_identity = float(document.get("partial_identity_floor", 95.0))
    partial = [pair for pair in hits if pair[0] >= partial_identity and pair[1] < coverage_floor]
    accepted = [pair for pair in partial if pair[1] >= partial_coverage]
    print(
        f"\npartial rule (identity ≥ {partial_identity}, coverage ≥ {partial_coverage}):"
        f" {len(accepted)} of {len(partial)} truncated-near-identical hit(s) accepted"
    )
    if accepted:
        coverages = sorted(pair[1] for pair in accepted)
        print(
            f"  their coverage: {coverages[0]:.1f}-{coverages[-1]:.1f} %"
            f" (floor is {partial_coverage})"
        )
    report["partial"] = {"accepted": len(accepted), "candidates": len(partial)}

    if args.out is not None:
        args.out.write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
        print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
