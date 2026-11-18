"""gspec.subspecies_kmers/1: differential k-mer signatures for the six
S. enterica subspecies.

Subspecies is a whole-genome property the antigen panel cannot decide. The
signatures here follow SalmID's clade-specific k-mer idea but are derived from
public-domain NCBI references (one per subspecies). Each set is that subspecies'
canonical k-mers minus every other subspecies', so the sets are disjoint by
construction and a query genome is classified by the highest recovery fraction.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, cast

DOC_PATH = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "gspec"
    / "data"
    / "dbs"
    / "sen_sero_sistr"
    / "subspecies_kmers.json"
)

SUBSPECIES = ("enterica", "salamae", "arizonae", "diarizonae", "houtenae", "indica")


def _load() -> dict[str, Any]:
    return cast(dict[str, Any], json.loads(DOC_PATH.read_text(encoding="utf-8")))


def _signatures(document: dict[str, Any]) -> dict[str, frozenset[int]]:
    """``subspecies -> canonical k-mer codes`` (the document stores hex strings)."""
    payload = cast(dict[str, Any], document["subspecies"])
    out: dict[str, frozenset[int]] = {}
    for name in SUBSPECIES:
        entry = cast(dict[str, Any], payload[name])
        hex_codes = cast(list[str], entry["kmers_hex"])
        out[name] = frozenset(int(code, 16) for code in hex_codes)
    return out


def test_document_covers_all_six_subspecies() -> None:
    document = _load()
    assert document["schema"] == "gspec.subspecies_kmers/1"
    assert document["k"] == 31
    subspecies = cast(dict[str, Any], document["subspecies"])
    assert set(subspecies) == set(SUBSPECIES)
    for name in SUBSPECIES:
        entry = cast(dict[str, Any], subspecies[name])
        assert cast(int, entry["n"]) > 0


def test_differential_sets_are_disjoint() -> None:
    """A k-mer may belong to at most one subspecies (it is differential)."""
    signatures = _signatures(_load())
    memberships: dict[int, int] = {}
    for codes in signatures.values():
        for code in codes:
            memberships[code] = memberships.get(code, 0) + 1
    shared = {code: count for code, count in memberships.items() if count > 1}
    assert not shared, f"{len(shared)} k-mers shared across subspecies"


def test_classify_recovers_the_right_subspecies() -> None:
    """A genome carrying only one subspecies' k-mers scores that subspecies top."""
    signatures = _signatures(_load())
    query = set(list(signatures["enterica"])[:1000])
    scored = sorted(
        ((name, len(query & codes) / len(codes)) for name, codes in signatures.items()),
        key=lambda item: item[1],
        reverse=True,
    )
    assert scored[0][0] == "enterica"
    assert scored[0][1] > 0.0
    # The other five must score zero (the query has no member of their sets).
    assert all(score == 0.0 for _, score in scored[1:])


def test_subspecies_smoke(tmp_path: Path) -> None:
    """A short random assembly still classifies without error (every recovery
    is ~0, but the command path must run end to end)."""
    query = tmp_path / "query.fna"
    query.write_text(
        ">q\n" + "A" * 500 + "C" * 500 + "G" * 500 + "T" * 500 + "\n",
        encoding="utf-8",
    )
    from gspec.cluster.kmerdb import count_fasta_kmers
    from gspec.subspecies import classify_kmers

    document = _load()
    signatures = _signatures(document)
    scored = classify_kmers(count_fasta_kmers(query, cast(int, document["k"])), signatures)
    assert len(scored) == 6
    assert all(0.0 <= item.recovery <= 1.0 for item in scored)
