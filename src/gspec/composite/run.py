"""Running a ``gspec.composite/1`` document: screen, type, fuse.

The pipeline is a second-stage orchestration over the existing contracts:
each stage screens its database through :func:`gspec.screening.run_screen`
(the same path ``gspec screen`` uses, so per-database ``screen.json``
defaults — including the ``gspec.screen/2`` coverage model — apply), then
each file is designated through the database's own typing document. Nothing
about thresholds or phenotype logic lives here: the pipeline only decides
*which stage's accepted call wins* and keeps every stage's evidence.

Grouping: stages that name the same database share one screen run, so a
locus database carrying several tier schemes (``tier_a``/``tier_b``/
``tier_c``) is aligned once and typed once — three stages for the price of
one. Cluster stages honour per-stage ``min_cluster_cov``/``min_gene_cov``
overrides (one locus database, several gates) and therefore screen per
distinct gate pair.
"""

from pathlib import Path

from gspec.cluster.core import require_typing
from gspec.composite.calls import Call, FileCall
from gspec.composite.cluster import cluster_calls
from gspec.composite.gene import gene_calls, screen_gene_table
from gspec.composite.kmer import kmer_calls
from gspec.composite.models import CompositeDocument, CompositeStage
from gspec.composite.pcr import base_type, load_primers, primer_products, primer_refine
from gspec.composite.stages import gene_params, stage_database, stage_scheme
from gspec.db import core as db
from gspec.screening.defaults import screen_defaults
from gspec.typing.models import TypingDocument


def apply_conflicts(
    document: CompositeDocument,
    winner: tuple[CompositeStage, Call],
    calls: dict[str, Call | None],
    stages: dict[str, CompositeStage],
) -> tuple[CompositeStage, Call]:
    """Swap a winning call for a declared override stage's call.

    Published assays carry documented cross-reactivity exceptions ("the HS41
    marker also amplifies in HS55 strains: trust the locus there"). The pipeline
    applies them between the cascade and any PCR refinement, and only when both
    competing calls were accepted on their own gates.
    """
    stage, call = winner
    base = base_type(call.phenotype or "")
    for rule in document.conflicts:
        when_stage = rule.when.stage
        if when_stage is not None and when_stage != stage.name:
            continue
        if base_type(rule.when.phenotype) != base:
            continue
        target = base_type(rule.override.phenotype)
        candidates = (
            [rule.override.stage]
            if rule.override.stage is not None
            else list(stages)  # document order (dicts insert in stage order)
        )
        for name in candidates:
            if name is None:
                continue
            override_stage = stages.get(name)
            override_call = calls.get(name)
            if override_stage is None or not accept(override_call, override_stage):
                continue
            assert override_call is not None  # accept guarantees it
            if base_type(override_call.phenotype or "") != target:
                continue
            return override_stage, override_call
    return winner


def accept(call: Call | None, stage: CompositeStage) -> bool:
    """Whether a stage's call is usable: designated, and clearing the stage gates.

    A stage that fell through to its own fallback (``confidence == "low"``, e.g.
    a marker scheme answering ``marker-negative``) designates nothing: accepting
    it would stop the pipeline on a layer that explicitly found no evidence, so
    fallbacks always fall through to the next stage.
    """
    if call is None or call.phenotype is None:
        return False
    if call.confidence == "low":
        return False
    if call.score < stage.min_score:
        return False
    return not (stage.require_confidence == "high" and call.confidence != "high")


def run_composite(
    paths: list[Path],
    document: CompositeDocument,
    *,
    datadir: Path | None,
    threads: int = 1,
    jobs: int = 1,
    quiet: bool = True,
    cache_dir: Path | None = None,
) -> list[FileCall]:
    """Screen every stage's database, designate every file, then fuse in stage
    order (first accepted call wins; every stage's call stays as evidence).

    ``kind="pcr"`` stages run after the database cascade: they refine the
    winning call (override, or bracket a complex type) from the primer products
    that amplified, and their evidence lands next to the database calls.
    """
    db_stages = tuple(stage for stage in document.stages if stage.kind == "db")
    kmer_stages = tuple(stage for stage in document.stages if stage.kind == "kmer")
    evidence_stages = tuple(stage for stage in document.stages if stage.kind in ("db", "kmer"))
    primer_stages = tuple(stage for stage in document.stages if stage.kind == "pcr")
    databases: dict[str, db.Database] = {}
    typing_docs: dict[str, TypingDocument] = {}
    schemes: dict[str, str] = {}
    evidence: dict[str, dict[str, Call]] = {stage.name: {} for stage in document.stages}

    # One screen per (db, cluster gate pair); gene stages of the same db share
    # a single table, cluster stages re-screen only when the gates differ.
    for stage in db_stages:
        database = databases.get(stage.db)
        if database is None:
            database = stage_database(document, stage.db, datadir)
            databases[stage.db] = database
            typing_docs[stage.db] = require_typing(database)
        schemes[stage.name] = stage_scheme(typing_docs[stage.db], stage)

    gene_tables: dict[str, str] = {}
    for stage in db_stages:
        database = databases[stage.db]
        scheme = schemes[stage.name]
        if database.kind == "cluster":
            evidence[stage.name] = cluster_calls(
                database,
                paths,
                stage,
                scheme,
                datadir=datadir,
                threads=threads,
                jobs=jobs,
                cache_dir=cache_dir,
            )
            continue
        table = gene_tables.get(stage.db)
        if table is None:
            defaults = screen_defaults(database)
            table = screen_gene_table(
                database,
                paths,
                gene_params(database, defaults, threads),
                datadir=datadir,
                jobs=jobs,
                quiet=quiet,
                cache_dir=cache_dir,
            )
            gene_tables[stage.db] = table
        evidence[stage.name] = gene_calls(database, table, typing_docs[stage.db], scheme)

    for stage in kmer_stages:
        database = stage_database(document, stage.db, datadir)
        evidence[stage.name] = kmer_calls(database, stage, paths, jobs=jobs)

    staged: dict[str, dict[str, Call | None]] = {}
    for path in paths:
        label = str(path)
        staged[label] = {stage.name: evidence[stage.name].get(label) for stage in evidence_stages}

    stage_by_name = {item.name: item for item in evidence_stages}
    cascade: dict[str, tuple[CompositeStage, Call] | None] = {}
    for label in (str(path) for path in paths):
        call_map = staged[label]
        winner: tuple[CompositeStage, Call] | None = None
        for stage in evidence_stages:
            call = call_map[stage.name]
            if accept(call, stage):
                assert call is not None
                winner = apply_conflicts(document, (stage, call), call_map, stage_by_name)
                break
        cascade[label] = winner

    # PCR products per stage: {file: {pair: [sizes]}}, scanning only the pairs
    # this file's call could need.
    products_by_stage: dict[str, dict[str, dict[str, list[int]]]] = {}
    for stage in primer_stages:
        primers = load_primers(stage, document)
        phenotypes = {
            label: (winner[1].phenotype if winner is not None else None)
            for label, winner in cascade.items()
        }
        products_by_stage[stage.name] = primer_products(
            stage, paths, primers, jobs, phenotypes=phenotypes
        )

    def primer_evidence(stage: CompositeStage, label: str) -> tuple[tuple[str, str], ...]:
        products = products_by_stage.get(stage.name, {}).get(label, {})
        return tuple((name, f"{sizes[0]}bp") for name, sizes in sorted(products.items()))

    results: list[FileCall] = []
    for label in (str(path) for path in paths):
        call_map = staged[label]
        winner = cascade[label]
        if winner is None:
            phenotype, tier, stage_name, score = document.fallback, "-", "-", 0.0
        else:
            stage, call = winner
            phenotype = call.phenotype or document.fallback
            tier, stage_name, score = stage.tier, stage.name, call.score

        # PCR refinement: the assay is authoritative for the call it names.
        for pcr_stage in primer_stages:
            refined = primer_refine(
                pcr_stage,
                (phenotype, stage_name, score),
                products_by_stage.get(pcr_stage.name, {}).get(label, {}),
            )
            if refined is None:
                continue
            phenotype, stage_name, score = refined
            tier = pcr_stage.tier
            break

        evidence_items: tuple[tuple[str, Call | str | None], ...] = tuple(
            (name, call_map[name]) for name in call_map
        )
        for pcr_stage in primer_stages:
            products = primer_evidence(pcr_stage, label)
            evidence_items += (
                (pcr_stage.name, "+".join(f"{name}({size})" for name, size in products) or None),
            )
        results.append(
            FileCall(
                file=label,
                phenotype=phenotype,
                tier=tier,
                stage=stage_name,
                score=score,
                evidence=evidence_items,
            )
        )
    return results
