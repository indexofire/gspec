"""Result containers for one composite run: a stage's call and the fused file call."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Call:
    """One stage's designation of one file."""

    phenotype: str | None
    score: float
    confidence: str
    scheme: str


@dataclass(frozen=True)
class FileCall:
    """The fused result: the winning stage, its call, and every stage's raw call."""

    file: str
    phenotype: str

    tier: str
    stage: str
    score: float
    evidence: tuple[tuple[str, Call | str | None], ...]
