"""Anchored allele sketches: the alignment-free tiebreak for antigen panels.

The measured facts these tests pin (see `scripts/salmonella/README.md`):

- a window's k-mers identify its allele far better than a BLAST HSP does, so a
  family whose alignment scores sit inside the margin can still be decided;
- a *whole genome* against an antigen's k-mer set cannot (it saturates), which is
  why the API takes a window;
- a tie the sketch cannot separate must stay a tie — the candidates are the
  honest answer, never a guess.
"""

from __future__ import annotations

import json
from pathlib import Path

from gspec.align.sketch import (
    SKETCH_MARGIN,
    AlleleSketch,
    SketchIndex,
    antigen_scores,
    window_kmers,
)
from gspec.cluster.core import GeneCall
from gspec.typing.serotype_engine import evaluate_serotype
from gspec.typing.serotype_models import (
    AntigenEntry,
    AntigenFamily,
    SerotypeDocumentV2,
    SerotypeFormula,
)

K = 21
GENE_A = "ACGTTGCAAGGCTTACGATCGATCGTAGCTAGGCTAACGTTAGGCATCGATCGTAGCTAGGCTAAC"
GENE_B = "ACGTTGCAAGGCTTACGATCGATCGTAGCTAGGCTAACGTTAGGCATCGATCGTAGCAAGGCTAAC"


def _codes(sequence: str) -> frozenset[int]:
    return window_kmers({"contig1": sequence}, "contig1", 1, len(sequence), K, pad=0)


def test_a_window_scores_the_allele_it_came_from() -> None:
    """Given a window that is allele A, When scored, Then A wins outright."""
    sketches = [
        AlleleSketch("h1", "a", _codes(GENE_A)),
        AlleleSketch("h1", "b", _codes(GENE_B)),
    ]
    scores = antigen_scores(_codes(GENE_A), sketches)
    assert scores["a"] == 1.0
    assert scores["b"] < scores["a"]


def test_an_antigen_scores_by_its_best_allele() -> None:
    """Two alleles of one antigen: the better one represents the antigen."""
    sketches = [
        AlleleSketch("h1", "a", _codes(GENE_A)),
        AlleleSketch("h1", "a", _codes(GENE_B[:40])),
        AlleleSketch("h1", "b", _codes(GENE_B)),
    ]
    scores = antigen_scores(_codes(GENE_A), sketches)
    assert scores["a"] == 1.0
    assert scores["a"] > scores["b"]


def test_a_window_with_no_shared_k_mers_scores_nothing() -> None:
    """An unrelated window names no antigen at all (no zero-score entries)."""
    unrelated = "TTTTAAAA" * 20
    assert antigen_scores(_codes(unrelated), [AlleleSketch("h1", "a", _codes(GENE_A))]) == {}


def test_the_index_groups_sketches_by_family() -> None:
    index = SketchIndex(
        k=K,
        families=(
            ("h1", (AlleleSketch("h1", "a", _codes(GENE_A)),)),
            ("o", (AlleleSketch("o", "B", _codes(GENE_B)),)),
        ),
    )
    assert [sketch.antigen for sketch in index.of("h1")] == ["a"]
    assert index.of("h2") == ()


def _call(identity: float, sketch: float | None = None) -> GeneCall:
    return GeneCall(
        gene_id="g",
        start=0,
        end=0,
        strand="+",
        coverage_pct=100.0,
        identity_pct=identity,
        verdict="present",
        sketch_pct=sketch,
    )


def _document() -> SerotypeDocumentV2:
    """Two H1 antigens inside the margin, each with its own serovar row."""
    return SerotypeDocumentV2(
        name="panel",
        control_gene=None,
        identity_floor=80.0,
        coverage_floor=60.0,
        margin=0.01,
        antigens=(
            AntigenFamily(key="o", entries=(AntigenEntry(name="C2", genes=("o_C2",)),)),
            AntigenFamily(
                key="h1",
                entries=(
                    AntigenEntry(name="e,h", genes=("h1_eh",)),
                    AntigenEntry(name="e,n,x", genes=("h1_enx",)),
                ),
            ),
            AntigenFamily(key="h2", entries=(AntigenEntry(name="1,2,7", genes=("h2_127",)),)),
        ),
        formulas=(
            SerotypeFormula(o="C2", h1="e,h", h2="1,2,7", serovar="Newport"),
            SerotypeFormula(o="C2", h1="e,n,x", h2="1,2,7", serovar="Other"),
        ),
    )


def test_the_sketch_breaks_a_margin_tie() -> None:
    """Given two antigens 0.005 apart in alignment but far apart in sketch,
    When the family is ranked, Then the sketch winner is the call — not a
    widened candidate list."""
    calls = {
        "o_C2": _call(100.0, 100.0),
        "h1_eh": _call(99.5, 100.0),
        "h1_enx": _call(100.0, 60.0),  # alignment winner, sketch loser
        "h2_127": _call(100.0, 100.0),
    }
    call = evaluate_serotype(calls, _document())
    assert call.phenotype == "Newport"
    assert "formula C2:e,h:1,2,7" in call.notes


def test_a_tie_the_sketch_cannot_separate_stays_a_tie() -> None:
    """Given both antigens at the same sketch score, When the family is ranked,
    Then the candidates stay visible rather than one being guessed."""
    calls = {
        "o_C2": _call(100.0, 100.0),
        "h1_eh": _call(99.5, 100.0),
        "h1_enx": _call(100.0, 100.0 - 100 * SKETCH_MARGIN / 2),
        "h2_127": _call(100.0, 100.0),
    }
    call = evaluate_serotype(calls, _document())
    assert call.phenotype is None
    assert call.confidence == "ambiguous"
    assert {entry.phenotype for entry in call.ambiguous} == {"Newport", "Other"}


def test_without_any_sketch_the_tie_widens_as_before() -> None:
    """A database that ships no sketches keeps the previous behaviour exactly."""
    calls = {
        "o_C2": _call(100.0),
        "h1_eh": _call(99.5),
        "h1_enx": _call(100.0),
        "h2_127": _call(100.0),
    }
    call = evaluate_serotype(calls, _document())
    assert call.phenotype is None
    assert {entry.phenotype for entry in call.ambiguous} == {"Newport", "Other"}


def test_typing_reads_the_sketch_from_a_screen_document(tmp_path: Path) -> None:
    """Given a screen JSON report, When typing parses it, Then each gene call
    carries its anchored sketch score; the TSV twin cannot (frozen columns)."""
    from gspec.typing.input import read_typing_input

    document = {
        "schema": "gspec.report/1",
        "tool": {"name": "gspec", "version": "0.0.0"},
        "created_at": "2026-01-01T00:00:00Z",
        "params": {"db": "panel", "minid": 80.0, "mincov": 60.0, "threads": 1},
        "files": [
            {
                "file": "genome.fna",
                "hits": [
                    {
                        "sequence": "ctg1",
                        "start": 1,
                        "end": 1500,
                        "strand": "+",
                        "gene": "h1_e_h",
                        "coverage": "1-1500/1500",
                        "coverage_map": "==========",
                        "gaps": "0/0",
                        "coverage_pct": 100.0,
                        "identity_pct": 99.5,
                        "database": "panel",
                        "accession": "upstream:e,h",
                        "product": "antigen eh",
                        "resistance": "antigen",
                        "sketch_pct": 100.0,
                    }
                ],
            }
        ],
    }
    path = tmp_path / "screen.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    parsed = read_typing_input([path])
    (calls,) = parsed.files
    assert calls.calls["h1_e_h"].sketch_pct == 100.0
    assert calls.calls["h1_e_h"].subject_len == 1500
