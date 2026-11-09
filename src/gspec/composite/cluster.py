"""The `db` composite stage on a CLUSTER database: roster loci, then designation."""

import json
from pathlib import Path

from gspec.composite.calls import Call
from gspec.composite.models import CompositeStage
from gspec.composite.stages import cache_path, inputs_digest
from gspec.db import core as db
from gspec.engines import OutputFormat
from gspec.screening.contigs import run_screen


def cluster_cache_key(db_name: str, stage: CompositeStage, paths: list[Path]) -> str:
    """The cluster stage's cache key: database, both gates, and the inputs."""
    cluster_cov = stage.min_cluster_cov if stage.min_cluster_cov is not None else 96.0
    gene_cov = stage.min_gene_cov if stage.min_gene_cov is not None else 90.0
    return f"{db_name}-C{cluster_cov:g}-G90-g{gene_cov:g}-{inputs_digest(paths)}"


def cluster_calls(
    database: db.Database,
    paths: list[Path],
    stage: CompositeStage,
    scheme: str,
    *,
    datadir: Path | None,
    threads: int,
    jobs: int,
    cache_dir: Path | None = None,
) -> dict[str, Call]:
    """Per-file designation of a cluster database (integrated typing)."""
    key = cluster_cache_key(database.name, stage, paths)
    text: str | None = None
    if cache_dir is not None:
        cached = cache_path(cache_dir, key, "json")
        if cached.is_file():
            text = cached.read_text(encoding="utf-8")
    if text is None:
        text = run_screen(
            paths,
            database.name,
            datadir,
            80.0,
            80.0,  # the cluster engine's default gene-engine pair (guard default)
            threads,
            jobs,
            None,
            True,
            False,
            True,
            False,
            OutputFormat.json,
            database=database,
            min_gene_cov=stage.min_gene_cov if stage.min_gene_cov is not None else 90.0,
            min_gene_id=90.0,
            min_cluster_cov=(stage.min_cluster_cov if stage.min_cluster_cov is not None else 96.0),
        )
        if cache_dir is not None:
            cache_path(cache_dir, key, "json").write_text(text, encoding="utf-8")
    document = json.loads(text)
    calls: dict[str, Call] = {}
    for entry in document.get("files", []):
        phenotypes = entry.get("phenotypes")
        if phenotypes is None:
            # gspec.cluster/1 (single-scheme document): the best-locus call.
            best = entry.get("best") or {}
            detail = best.get("phenotype_detail") or {}
            calls[entry["file"]] = Call(
                phenotype=best.get("phenotype"),
                score=detail.get("score", 0.0) or 0.0,
                confidence=detail.get("confidence", "low"),
                scheme=scheme,
            )
            continue
        call = phenotypes.get(scheme)
        if call is None:
            continue
        calls[entry["file"]] = Call(
            phenotype=call.get("phenotype"),
            score=call.get("score", 0.0) or 0.0,
            confidence=call.get("confidence", "low"),
            scheme=scheme,
        )
    return calls
