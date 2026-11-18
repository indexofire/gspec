"""Salmonella subspecies identification from differential k-mer signatures.

Subspecies is a whole-genome property: the antigen panel cannot decide it.
SISTR uses cgMLST and SeqSero2/SalmID uses clade-specific k-mers; gspec follows
the latter with differential k-mer sets derived from public-domain NCBI
references (``scripts/salmonella/build_subspecies_kmers.py``). A query genome is
classified by the highest recovery fraction — the share of one subspecies'
differential k-mers it carries.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from gspec.cluster.kmerdb import count_fasta_kmers
from gspec.errors import DatabaseError

SUBSPECIES_FILENAME = "subspecies_kmers.json"

#: The six S. enterica subspecies, with their Kauffmann-White roman numerals.
SUBSPECIES = (
    ("enterica", "I"),
    ("salamae", "II"),
    ("arizonae", "IIIa"),
    ("diarizonae", "IIIb"),
    ("houtenae", "IV"),
    ("indica", "VI"),
)


@dataclass(frozen=True)
class SubspeciesScore:
    """One subspecies' recovery fraction for a query genome."""

    name: str
    numeral: str
    recovery: float


def load_subspecies_signatures(path: Path) -> tuple[int, dict[str, frozenset[int]]]:
    """``(k, {subspecies: canonical k-mer codes})`` from a signature document."""
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != "gspec.subspecies_kmers/1":
        raise DatabaseError(
            f"{path}: unsupported subspecies signature schema {document.get('schema')!r}",
            code="SUBSPECIES_MALFORMED",
            context={"file": str(path)},
        )
    k = int(document["k"])
    signatures = {
        name: frozenset(int(code, 16) for code in payload["kmers_hex"])
        for name, payload in document["subspecies"].items()
    }
    for name, _ in SUBSPECIES:
        if name not in signatures:
            raise DatabaseError(
                f"{path}: missing subspecies {name!r}",
                code="SUBSPECIES_MALFORMED",
                context={"file": str(path)},
            )
    return k, signatures


def classify_kmers(query: set[int], signatures: dict[str, frozenset[int]]) -> list[SubspeciesScore]:
    """Recovery fraction per subspecies, best first."""
    numeral = dict(SUBSPECIES)
    scored = [
        SubspeciesScore(
            name=name,
            numeral=numeral[name],
            recovery=sum(1 for code in codes if code in query) / len(codes),
        )
        for name, codes in signatures.items()
    ]
    return sorted(scored, key=lambda item: item.recovery, reverse=True)


def classify_fasta(
    fasta: Path, signatures: dict[str, frozenset[int]], k: int
) -> list[SubspeciesScore]:
    """Classify one assembly FASTA against the signature sets."""
    return classify_kmers(count_fasta_kmers(fasta, k), signatures)


__all__ = [
    "SUBSPECIES",
    "SUBSPECIES_FILENAME",
    "SubspeciesScore",
    "classify_fasta",
    "classify_kmers",
    "load_subspecies_signatures",
]
