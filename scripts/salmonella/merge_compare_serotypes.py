#!/usr/bin/env python3
"""Merge serotype evidence and test a serovar pair for a discriminating marker.

Two jobs, because they always run together:

``merge``
    Join every evidence layer under ``~/data/sen`` by accession, keep the most
    trustworthy label, and surface the rows where two sources disagree. The
    layers, strongest first:

    ================  ======  ==================================================
    layer             weight  meaning
    ================  ======  ==================================================
    ``verified``          4  a laboratory's own agglutination result
                             (``verified_manual.tsv`` — you supply this)
    ``culture``           3  strain deposited in a public collection: the
                             serovar's reference culture lives there
    ``published``         2  a paper's supplementary accession+serovar table
    ``asserted``          1  the INSDC submitter's ``serovar`` attribute
    ================  ======  ==================================================

``compare``
    For two groups of serovars, download the genomes each layer can supply, cut
    the rfb neighbourhood out of every one (anchored on a panel wzx gene),
    cluster the CDS by identity, and ask a single question: **is there a cluster
    whose presence — or whose allele — tracks the group?** A marker only counts
    when it is present in *every* genome of one group and in *no* genome of the
    other, and its within-group identity is at least as high as the panel's own
    working markers (≈100 %).

Usage::

    # 1. merge whatever is on disk (``verified_manual.tsv`` is optional)
    python scripts/salmonella/merge_compare_serotypes.py merge

    # 2. ask whether O:6 tracks the two pairs gspec cannot separate
    python scripts/salmonella/merge_compare_serotypes.py compare \\
        --group "factor6=Goldcoast,Bovismorbificans" \\
        --group "no_factor6=Albany,Hindmarsh,Brikama"

The ``compare`` step needs ``datasets`` and ``blastn`` on PATH (both already
required by gspec's own database build).
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import subprocess
import sys
import zipfile
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

DEFAULT_ROOT = Path.home() / "data" / "sen"

#: Evidence layers and their weight. A higher weight wins a label conflict.
LAYER_WEIGHT = {"verified": 4, "culture": 3, "published": 2, "asserted": 1}

#: One row per accession: where the label came from and how much it is worth.
EVIDENCE_COLUMNS = [
    "accession",
    "base",
    "serovar",
    "layer",
    "sources",
    "n_labels",
    "conflict",
    "assembly_level",
]

BASE_RE = re.compile(
    r"^((?:GCF|GCA)_\d{9,}|(?:SAMN|SAMEA)\d{6,}|(?:ERS|SRR|ERR|DRR)\d{5,})(?:\.\d+)?$"
)


# --------------------------------------------------------------------------- #
# merge
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Evidence:
    """One label for one accession, with the layer that asserted it."""

    accession: str
    serovar: str
    layer: str
    source: str
    assembly_level: str = ""


def _base(accession: str) -> str:
    """``GCF_000783815.2`` -> ``GCF_000783815`` so versions still join."""
    match = BASE_RE.match(accession.strip())
    return match.group(1) if match else accession.strip()


def _clean(serovar: str) -> str:
    text = " ".join(serovar.split())
    for prefix in ("Salmonella enterica subsp. enterica ", "Salmonella enterica ", "Salmonella "):
        if text.lower().startswith(prefix.lower()):
            text = text[len(prefix) :]
    return text.strip()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def collect_evidence(root: Path) -> list[Evidence]:
    """Every (accession, serovar, layer) the four sources can supply."""
    out: list[Evidence] = []

    for row in _read_tsv(root / "verified_manual.tsv"):
        accession = (row.get("accession") or row.get("assembly") or "").strip()
        serovar = _clean(row.get("serovar") or "")
        if accession and serovar:
            out.append(Evidence(accession, serovar, "verified", "verified_manual.tsv"))

    for row in _read_tsv(root / "reference_cultures.tsv"):
        accession = (row.get("assembly") or row.get("biosample") or "").strip()
        serovar = _clean(row.get("serovar") or "")
        if accession and serovar:
            out.append(
                Evidence(
                    accession,
                    serovar,
                    "culture",
                    f"reference_cultures.tsv:{row.get('culture_collection', '')}",
                    row.get("assembly_level", ""),
                )
            )

    for row in _read_tsv(root / "published/accessions_from_papers.tsv"):
        accession = (row.get("accession") or "").strip()
        serovar = _clean(row.get("serovar") or "")
        if accession and serovar:
            out.append(Evidence(accession, serovar, "published", f"doi:{row.get('doi', '')}"))

    for row in _read_tsv(root / "ncbi_refseq_serovars.tsv"):
        accession = (row.get("accession") or "").strip()
        serovar = _clean(row.get("serovar") or "")
        if accession and serovar:
            out.append(
                Evidence(
                    accession,
                    serovar,
                    "asserted",
                    "ncbi_refseq_serovars.tsv",
                    row.get("assembly_level", ""),
                )
            )
    return out


def merge_evidence(evidence: Iterable[Evidence]) -> list[dict[str, object]]:
    """One row per accession, best label kept, conflicts marked.

    A conflict is not an error to hide: two sources naming different serovars for
    one accession is exactly the signal that the label cannot be trusted, so the
    row keeps every candidate label in ``conflict`` and the caller can drop it.
    """
    by_base: dict[str, list[Evidence]] = defaultdict(list)
    for item in evidence:
        by_base[_base(item.accession)].append(item)

    rows: list[dict[str, object]] = []
    for base, items in sorted(by_base.items()):
        labels: dict[str, int] = defaultdict(int)
        for item in items:
            labels[item.serovar] = max(labels[item.serovar], LAYER_WEIGHT[item.layer])
        # Best label wins; ties break alphabetically so the output is stable.
        winner = sorted(labels.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        best = max(items, key=lambda item: LAYER_WEIGHT[item.layer])
        levels = [item.assembly_level for item in items if item.assembly_level]
        rows.append(
            {
                "accession": best.accession,
                "base": base,
                "serovar": winner,
                "layer": best.layer,
                "sources": ";".join(sorted({item.source for item in items})),
                "n_labels": len(items),
                "conflict": ";".join(sorted(labels)),
                "assembly_level": max(levels, key=len) if levels else "",
            }
        )
    return rows


def write_merge(root: Path) -> Path:
    rows = merge_evidence(collect_evidence(root))
    path = root / "merged_evidence.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=EVIDENCE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    conflicts = [row for row in rows if len(str(row["conflict"]).split(";")) > 1]
    by_layer: dict[str, int] = defaultdict(int)
    for row in rows:
        by_layer[str(row["layer"])] += 1
    print(
        f"{path}: {len(rows)} accessions"
        f" ({', '.join(f'{k}={v}' for k, v in sorted(by_layer.items()))});"
        f" {len(conflicts)} with conflicting labels",
        file=sys.stderr,
    )
    if conflicts:
        conflict_path = root / "label_conflicts.tsv"
        with conflict_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=EVIDENCE_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(conflicts)
        print(f"{conflict_path}: {len(conflicts)} rows to review", file=sys.stderr)
    return path


# --------------------------------------------------------------------------- #
# compare
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Genome:
    """One genome to analyse, with the label the evidence supports."""

    accession: str
    serovar: str
    group: str
    layer: str


def select_genomes(
    rows: Iterable[dict[str, object]],
    groups: dict[str, list[str]],
    *,
    minimum_layer: str,
    per_serovar: int = 4,
) -> list[Genome]:
    """Up to ``per_serovar`` accessions per wanted serovar, best layer first.

    A marker has to hold across a *group*, not one genome, so several genomes per
    serovar are the point. Conflicting labels are skipped outright: a genome whose
    own sources disagree cannot support a marker decision. Only assemblies are
    useful (``datasets`` fetches GCF/GCA), so `SRR`/`SAMN` rows are dropped.
    """
    wanted = {serovar.lower(): group for group, serovars in groups.items() for serovar in serovars}
    floor = LAYER_WEIGHT[minimum_layer]
    candidates: dict[str, list[tuple[int, int, int, str, str, str]]] = defaultdict(list)
    for row in rows:
        serovar = str(row["serovar"]).lower()
        group = wanted.get(serovar)
        if group is None:
            continue
        if len(str(row["conflict"]).split(";")) > 1:
            continue
        layer = str(row["layer"])
        if LAYER_WEIGHT[layer] < floor:
            continue
        accession = str(row["accession"])
        if not accession.startswith(("GCF_", "GCA_")):
            continue
        # Prefer RefSeq ("complete" second), then the strongest evidence layer.
        complete = 1 if "complete" in str(row["assembly_level"]).lower() else 0
        refseq = 1 if accession.startswith("GCF_") else 0
        # The candidate carries its own serovar and layer: reading them back from
        # a loop variable later silently substituted the last row parsed.
        candidates[serovar].append(
            (-LAYER_WEIGHT[layer], -refseq, -complete, accession, layer, str(row["serovar"]))
        )

    chosen: list[Genome] = []
    for serovar, entries in sorted(candidates.items()):
        seen: set[str] = set()
        for _weight, _refseq, _complete, accession, layer, label in sorted(entries)[:per_serovar]:
            if accession in seen:
                continue
            seen.add(accession)
            chosen.append(Genome(accession, label, wanted[serovar], layer))
    return sorted(chosen, key=lambda genome: (genome.group, genome.serovar, genome.accession))


def _run(command: list[str], **kwargs: Any) -> str:
    result = subprocess.run(command, capture_output=True, text=True, check=False, **kwargs)
    if result.returncode != 0:
        raise SystemExit(f"{' '.join(command)} failed:\n{result.stderr[-2000:]}")
    return result.stdout


def download(genomes: list[Genome], work: Path) -> dict[str, Path]:
    """Fetch genome + CDS for each accession; returns ``base -> genome fasta``."""
    work.mkdir(parents=True, exist_ok=True)
    data = work / "ncbi_dataset" / "data"
    if not data.is_dir():
        accessions = "\n".join(genome.accession for genome in genomes) + "\n"
        accession_file = work / "accessions.txt"
        accession_file.write_text(accessions, encoding="utf-8")
        _run(
            [
                "datasets",
                "download",
                "genome",
                "accession",
                "--inputfile",
                str(accession_file),
                "--include",
                "genome,cds",
                "--filename",
                str(work / "set.zip"),
            ]
        )
        with zipfile.ZipFile(work / "set.zip") as archive:
            archive.extractall(work)
    found: dict[str, Path] = {}
    for directory in sorted(data.iterdir()):
        genome = next((p for p in directory.glob("*genomic.fna") if "cds_from" not in p.name), None)
        if genome is not None:
            found[directory.name] = genome
    return found


def rfb_window(genome: Path, anchor_fasta: Path, flank: int) -> tuple[str, int, int]:
    """``(sequence, low, high)`` of the rfb neighbourhood around the anchor hit."""
    database = genome.with_suffix(".anchor")
    _run(
        [
            "makeblastdb",
            "-in",
            str(genome),
            "-dbtype",
            "nucl",
            "-out",
            str(database),
            "-logfile",
            "/dev/null",
        ]
    )
    out = _run(
        [
            "blastn",
            "-query",
            str(anchor_fasta),
            "-db",
            str(database),
            "-outfmt",
            "6 sseqid sstart send pident",
            "-max_target_seqs",
            "1",
            "-dust",
            "no",
        ]
    ).strip()
    if not out:
        raise SystemExit(f"{genome}: no anchor hit")
    contig, start, end, _identity = out.split("\n")[0].split("\t")
    low, high = sorted((int(start), int(end)))
    chunks: list[str] = []
    inside = False
    for line in genome.read_text().splitlines():
        if line.startswith(">"):
            inside = line[1:].split()[0] == contig
        elif inside:
            chunks.append(line.strip().upper())
    sequence = "".join(chunks)
    from_coordinate = max(0, low - flank)
    to_coordinate = high + flank
    return sequence[from_coordinate:to_coordinate], from_coordinate, to_coordinate


def cds_in_window(cds_path: Path, low: int, high: int) -> list[dict[str, str]]:
    """CDS records of one assembly whose coordinates fall inside the window."""
    records: list[dict[str, str]] = []
    for block in cds_path.read_text().split(">")[1:]:
        first = block.split("\n", 1)[0]
        match = re.search(r"\[location=(?:complement\()?<?(\d+)\.\.>?(\d+)", first)
        if not match:
            continue
        a, b = int(match.group(1)), int(match.group(2))
        start, end = min(a, b), max(a, b)
        if end < low or start > high:
            continue
        sequence = "".join(block.split("\n")[1:]).upper()
        product = (re.search(r"\[protein=([^\]]+)\]", first) or [None, ""])[1]
        gene = (re.search(r"\[gene=([^\]]+)\]", first) or [None, ""])[1]
        records.append(
            {
                "start": str(start),
                "end": str(end),
                "product": product,
                "gene": gene,
                "sequence": sequence,
            }
        )
    return sorted(records, key=lambda record: int(record["start"]))


def cluster_cds(fasta: Path, work: Path) -> dict[str, set[str]]:
    """Single-linkage clusters of CDS at >=95 % identity over >=80 % length."""
    database = work / "cds"
    _run(
        [
            "makeblastdb",
            "-in",
            str(fasta),
            "-dbtype",
            "nucl",
            "-out",
            str(database),
            "-logfile",
            "/dev/null",
        ]
    )
    out = _run(
        [
            "blastn",
            "-query",
            str(fasta),
            "-db",
            str(database),
            "-outfmt",
            "6 qseqid sseqid pident length qlen slen",
            "-dust",
            "no",
            "-evalue",
            "1e-20",
            "-max_target_seqs",
            "5000",
        ]
    )
    identifiers = [
        line[1:].split()[0] for line in fasta.read_text().splitlines() if line.startswith(">")
    ]
    parent = {identifier: identifier for identifier in identifiers}

    def find(item: str) -> str:
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    for line in out.splitlines():
        if not line:
            continue
        query, subject, identity, length, query_length, subject_length = line.split("\t")
        if (
            float(identity) >= 95.0
            and int(length) / max(int(query_length), int(subject_length)) >= 0.80
        ):
            left, right = find(query), find(subject)
            if left != right:
                parent[right] = left
    clusters: dict[str, set[str]] = defaultdict(set)
    for identifier in identifiers:
        clusters[find(identifier)].add(identifier)
    return clusters


def _identity(first: str, second: str) -> float:
    length = min(len(first), len(second))
    if length < 200:
        return 0.0
    mismatches = sum(1 for i in range(length) if first[i] != second[i])
    return 100.0 * (length - mismatches) / length


def analyse(
    genomes: list[Genome],
    fastas: dict[str, Path],
    work: Path,
    *,
    anchor_fasta: Path,
    flank: int,
) -> list[dict[str, object]]:
    """Cluster the rfb CDS and score each cluster against the groups.

    A cluster is a candidate marker only if it is present in *every* genome of
    one group and in *no* genome of the other, **and** its sequences are
    conserved inside the group — the panel's working markers are ~100 % identical
    within group, so anything much looser is strain variation, not a serogroup
    signal.
    """
    groups = sorted({genome.group for genome in genomes})
    totals: dict[str, int] = defaultdict(int)
    for genome in genomes:
        totals[genome.group] += 1

    fasta_path = work / "rfb_cds.fna"
    annotations: dict[str, dict[str, str]] = {}
    group_of: dict[str, str] = {}
    written = 0
    with fasta_path.open("w", encoding="utf-8") as handle:
        for genome in genomes:
            directory = work / "ncbi_dataset" / "data" / genome.accession
            genome_fasta = fastas.get(genome.accession)
            cds_path = directory / "cds_from_genomic.fna"
            if genome_fasta is None or not cds_path.is_file():
                print(f"  skip {genome.accession} (no genome or CDS)", file=sys.stderr)
                continue
            _sequence, low, high = rfb_window(genome_fasta, anchor_fasta, flank)
            for record in cds_in_window(cds_path, low, high):
                identifier = f"{genome.accession}|{record['start']}"
                annotations[identifier] = record
                group_of[identifier] = genome.group
                handle.write(f">{identifier}\n{record['sequence']}\n")
                written += 1
    print(f"  {written} CDS in {len(annotations)} records", file=sys.stderr)
    clusters = cluster_cds(fasta_path, work)
    genome_group = {genome.accession: genome.group for genome in genomes}

    rows: list[dict[str, object]] = []
    for members in clusters.values():
        present: dict[str, int] = defaultdict(int)
        per_genome: dict[str, str] = {}
        for identifier in members:
            accession = identifier.split("|")[0]
            present[group_of[identifier]] += 1
            per_genome.setdefault(accession, annotations[identifier]["sequence"])
        representative = sorted(members, key=lambda item: int(item.split("|")[1]))[0]
        annotation = annotations[representative]
        within: list[float] = []
        between: list[float] = []
        accessions = sorted(per_genome)
        for index, first in enumerate(accessions):
            for second in accessions[index + 1 :]:
                score = _identity(per_genome[first], per_genome[second])
                same = genome_group[first] == genome_group[second]
                (within if same else between).append(score)
        rows.append(
            {
                "n": len(members),
                "present": ";".join(
                    f"{group}={present.get(group, 0)}/{totals[group]}" for group in groups
                ),
                "gene": annotation.get("gene", ""),
                "product": annotation.get("product", ""),
                "representative": representative,
                "within_pct": round(sum(within) / len(within), 2) if within else 0.0,
                "between_pct": round(sum(between) / len(between), 2) if between else 0.0,
            }
        )
    return rows


def verdict(rows: list[dict[str, object]], groups: list[str]) -> list[str]:
    """The clusters whose presence splits the groups perfectly."""
    lines: list[str] = []
    for row in rows:
        presence = dict(item.split("=") for item in str(row["present"]).split(";"))
        counts = {group: presence.get(group, "0/0") for group in groups}
        full = [
            group for group, value in counts.items() if value.split("/")[0] == value.split("/")[1]
        ]
        empty = [group for group, value in counts.items() if value.split("/")[0] == "0"]
        if len(full) == 1 and len(empty) == len(groups) - 1:
            lines.append(
                f"  {row['gene'] or '-'} / {str(row['product'])[:48]:48}"
                f" {counts}  <== only in {full[0]}"
            )
    return lines


def write_report(
    path: Path, genomes: list[Genome], groups: dict[str, list[str]], rows: list[dict[str, object]]
) -> None:
    group_names = sorted(groups)
    lines = ["# rfb marker comparison", ""]
    lines.append("## Genomes analysed")
    lines.append("")
    lines.append("| serovar | group | accession | evidence |")
    lines.append("| --- | --- | --- | --- |")
    for genome in genomes:
        lines.append(f"| {genome.serovar} | {genome.group} | {genome.accession} | {genome.layer} |")
    lines.append("")
    lines.append("## Clusters that split the groups perfectly")
    lines.append("")
    found = verdict(rows, group_names)
    lines.extend(
        found or ["  (none — no cluster is present in one group and absent from the other)"]
    )
    lines.append("")
    lines.append(
        "A cluster is only a candidate if it is `N/N` for one group and `0/N` for the"
        " other; `within %` is the mean identity between genomes of the same group and"
        " `between %` across groups. The panel's working markers are ~100 % within group,"
        " so a within-group figure well below that is strain variation, not a serogroup"
        " signal."
    )
    lines.append("")
    lines.append("## All clusters")
    lines.append("")
    lines.append("| n | presence | within % | between % | gene | product |")
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for row in sorted(rows, key=lambda item: -int(item["n"])):
        lines.append(
            f"| {row['n']} | {row['present']} | {row['within_pct']} | {row['between_pct']}"
            f" | {row['gene'] or '-'} | {str(row['product'])[:52]} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("merge", "compare"))
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument(
        "--group",
        action="append",
        default=[],
        metavar="NAME=serovar,serovar",
        help="compare mode: a group name and its serovars (repeat for each group)",
    )
    parser.add_argument(
        "--min-layer",
        default="published",
        choices=sorted(LAYER_WEIGHT),
        help="weakest evidence layer a compare genome may come from (default published)",
    )
    parser.add_argument(
        "--anchor-gene",
        default="o_C2-C3",
        help="panel gene id used to locate the rfb cluster (default o_C2-C3)",
    )
    parser.add_argument("--flank", type=int, default=14000)
    parser.add_argument(
        "--per-serovar",
        type=int,
        default=4,
        help="how many genomes per serovar to analyse (default 4)",
    )
    parser.add_argument("--datadir", type=Path, default=Path.home() / ".local/share/gspec/db")
    args = parser.parse_args(argv)

    args.root.mkdir(parents=True, exist_ok=True)
    merged = write_merge(args.root)
    if args.command == "merge":
        return 0

    if len(args.group) < 2:
        raise SystemExit("compare needs at least two --group NAME=serovar,serovar arguments")
    groups: dict[str, list[str]] = {}
    for spec in args.group:
        name, _, serovars = spec.partition("=")
        groups[name.strip()] = [item.strip() for item in serovars.split(",") if item.strip()]

    rows = list(csv.DictReader(io.StringIO(merged.read_text(encoding="utf-8")), delimiter="\t"))
    genomes = select_genomes(
        rows, groups, minimum_layer=args.min_layer, per_serovar=args.per_serovar
    )
    if not genomes:
        raise SystemExit(
            "no genome qualified — lower --min-layer or add the serovars to verified_manual.tsv"
        )
    work = args.root / "analysis"
    print(f"{len(genomes)} genomes selected; downloading…", file=sys.stderr)
    fastas = download(genomes, work)

    anchor_source = args.datadir / "sen_sero_sistr" / "sequences"
    anchor_fasta = work / "anchor.fna"
    if not anchor_fasta.is_file():
        found = False
        chunks: list[str] = []
        inside = False
        for line in anchor_source.read_text().splitlines():
            if line.startswith(">"):
                inside = f"gene={args.anchor_gene}|" in line
                if inside:
                    found = True
            elif inside:
                chunks.append(line.strip())
        if not found:
            raise SystemExit(f"anchor gene {args.anchor_gene!r} not found in {anchor_source}")
        anchor_fasta.write_text(f">{args.anchor_gene}\n{''.join(chunks)}\n", encoding="utf-8")

    clusters = analyse(genomes, fastas, work, anchor_fasta=anchor_fasta, flank=args.flank)
    report = args.root / "analysis_report.md"
    write_report(report, genomes, groups, clusters)
    found = verdict(clusters, sorted(groups))
    print(f"{report}: {len(clusters)} clusters", file=sys.stderr)
    print(
        "\n".join(found) if found else "no perfect group-splitting cluster found",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
