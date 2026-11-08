"""K-mer stages: an alignment-free evidence source inside a pipeline.

A ``kind: "kmer"`` stage scores a genome by differential k-mer sets (the
``gspec.kmerdb/1`` substrate cluster databases already use for SNP-level pairs)
and turns the best set's hit fraction into a call. It is the pipeline's third
method family next to alignment (``db``) and PCR (``pcr``), so a scheme can mix
blastn, minimap2, k-mers and in-silico PCR in one answer.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from gspec.composite_models import read_composite
from gspec.kmerdb import KmerDbDocument, KmerSetSpec, materialize_kmer_sets

SEQUENCE_A = "ACGT" * 60 + "TTGG" * 30
SEQUENCE_B = "GGCC" * 60 + "AATT" * 30
SEQUENCE_C = "TTAA" * 120


def _kmer_db(tmp_path: Path) -> Path:
    """A minimal installed database carrying a k-mer sidecar."""
    db_dir = tmp_path / "datadir" / "kmer_db"
    db_dir.mkdir(parents=True)
    (db_dir / "sequences").write_text(
        f">gspec|db=kmer_db|gene=A|acc=a|func= A\n{SEQUENCE_A}\n"
        f">gspec|db=kmer_db|gene=B|acc=b|func= B\n{SEQUENCE_B}\n",
        encoding="utf-8",
    )
    (db_dir / "gspec-manifest.json").write_text(
        json.dumps(
            {
                "schema": "gspec.manifest/1",
                "name": "kmer_db",
                "source_urls": ["local"],
                "fetched_at": "2026-01-01T00:00:00Z",
                "sha256": "0" * 64,
                "n_records": 2,
                "dbtype": "nucl",
                "kind": "gene",
            }
        ),
        encoding="utf-8",
    )
    references = tmp_path / "references.fna"
    references.write_text(
        f">A\n{SEQUENCE_A}\n>B\n{SEQUENCE_B}\n>C\n{SEQUENCE_C}\n", encoding="utf-8"
    )
    document = KmerDbDocument.model_validate(
        {
            "schema": "gspec.kmerdb/1",
            "k": 21,
            "sets": {
                "A.diff": KmerSetSpec(of="A", against=("B", "C")).model_dump(),
                "B.diff": KmerSetSpec(of="B", against=("A", "C")).model_dump(),
            },
        }
    )
    (db_dir / "kmer.json").write_text(
        json.dumps(document.model_dump(by_alias=True), indent=1), encoding="utf-8"
    )
    materialize_kmer_sets(db_dir, document, references)
    return db_dir.parent


def _genome(tmp_path: Path, sequence: str, name: str = "asm.fna") -> Path:
    path = tmp_path / name
    path.write_text(f">ctg1\n{sequence}\n", encoding="utf-8")
    return path


def _pipeline(tmp_path: Path, *, min_fraction: float = 0.5, sets: list[dict] | None = None) -> Path:
    document = {
        "schema": "gspec.composite/1",
        "name": "kmer_demo",
        "fallback": "untypeable",
        "stages": [
            {
                "name": "snps",
                "kind": "kmer",
                "db": "kmer_db",
                "tier": "K",
                "sets": sets
                or [
                    {"set": "A.diff", "phenotype": "TypeA", "min_fraction": min_fraction},
                    {"set": "B.diff", "phenotype": "TypeB", "min_fraction": min_fraction},
                ],
            }
        ],
    }
    path = tmp_path / "composite.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    return path


def test_kmer_stage_calls_from_the_matching_set(tmp_path: Path) -> None:
    from gspec.composite_ops import run_composite

    datadir = _kmer_db(tmp_path)
    document = read_composite(_pipeline(tmp_path))
    result = run_composite([_genome(tmp_path, SEQUENCE_A)], document, datadir=datadir)
    assert result[0].phenotype == "TypeA"
    assert result[0].stage == "snps"
    assert result[0].score > 0.5


def test_kmer_stage_separates_two_references(tmp_path: Path) -> None:
    from gspec.composite_ops import run_composite

    datadir = _kmer_db(tmp_path)
    document = read_composite(_pipeline(tmp_path))
    genomes = [
        _genome(tmp_path, SEQUENCE_A, "a.fna"),
        _genome(tmp_path, SEQUENCE_B, "b.fna"),
    ]
    results = run_composite(genomes, document, datadir=datadir)
    assert [result.phenotype for result in results] == ["TypeA", "TypeB"]


def test_min_fraction_rejects_a_weak_match(tmp_path: Path) -> None:
    from gspec.composite_ops import run_composite

    datadir = _kmer_db(tmp_path)
    document = read_composite(_pipeline(tmp_path, min_fraction=0.99))
    # C is in no set, so nothing crosses a 0.99 floor.
    result = run_composite([_genome(tmp_path, SEQUENCE_C)], document, datadir=datadir)[0]
    assert result.phenotype == "untypeable"


def test_unknown_set_is_rejected_at_load(tmp_path: Path) -> None:
    from gspec.composite_ops import run_composite

    datadir = _kmer_db(tmp_path)
    document = read_composite(_pipeline(tmp_path, sets=[{"set": "nope.diff", "phenotype": "X"}]))
    with pytest.raises(Exception) as exc_info:
        run_composite([_genome(tmp_path, SEQUENCE_A)], document, datadir=datadir)
    assert getattr(exc_info.value, "code", "") == "COMPOSITE_MALFORMED"


def test_database_without_a_kmer_sidecar_is_a_typed_error(tmp_path: Path) -> None:
    from gspec.composite_ops import run_composite

    datadir = _kmer_db(tmp_path)
    (datadir / "kmer_db" / "kmer.json").unlink()
    document = read_composite(_pipeline(tmp_path))
    with pytest.raises(Exception) as exc_info:
        run_composite([_genome(tmp_path, SEQUENCE_A)], document, datadir=datadir)
    assert getattr(exc_info.value, "code", "") == "KMERDB_MALFORMED"


def test_kmer_stage_needs_sets(tmp_path: Path) -> None:
    from gspec.composite_models import CompositeDocument

    with pytest.raises(ValueError):
        CompositeDocument.model_validate(
            {
                "schema": "gspec.composite/1",
                "name": "x",
                "stages": [{"name": "s", "kind": "kmer", "db": "kmer_db", "sets": []}],
            }
        )


def test_cluster_stage_cache_key_covers_each_gate(tmp_path: Path) -> None:
    """A cluster stage's cache key must change with its gates, not only its db."""
    from gspec.composite_models import CompositeStage
    from gspec.composite_ops import _cluster_cache_key

    genomes = [tmp_path / "a.fna"]
    genomes[0].write_text(">c\nACGT\n", encoding="utf-8")
    stage = CompositeStage.model_validate({"name": "c", "db": "locus_db"})
    strict = CompositeStage.model_validate(
        {"name": "c", "db": "locus_db", "min_cluster_cov": 90.0, "min_gene_cov": 90.0}
    )
    lenient = CompositeStage.model_validate(
        {"name": "c", "db": "locus_db", "min_cluster_cov": 60.0, "min_gene_cov": 78.0}
    )
    keys = {
        _cluster_cache_key("locus_db", stage, genomes),
        _cluster_cache_key("locus_db", strict, genomes),
        _cluster_cache_key("locus_db", lenient, genomes),
        _cluster_cache_key("other_db", stage, genomes),
    }
    assert len(keys) == 4
