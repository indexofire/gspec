"""Composite result rendering: the gspec.composite_result/1 document and its table."""

import json
from datetime import UTC, datetime
from typing import Any

from gspec.composite.calls import Call, FileCall
from gspec.composite.models import CompositeDocument
from gspec.formats.json import utc_timestamp


def evidence_cell(value: Call | str | None) -> str:
    """One stage's evidence as a table cell: a call, a PCR product list, or ``-``."""
    if isinstance(value, str):
        return value
    if value is None or not value.phenotype:
        return "-"
    return f"{value.phenotype}({value.score:.4f})"


def render_composite_tsv(
    results: list[FileCall], document: CompositeDocument, *, noheader: bool = False
) -> str:
    """One row per file: the fused call plus every stage's raw call as evidence."""
    lines: list[str] = []
    if not noheader:
        lines.append("FILE\tPIPELINE\tPHENOTYPE\tTIER\tSTAGE\tSCORE\tEVIDENCE")
    for result in results:
        evidence = "|".join(f"{name}={evidence_cell(value)}" for name, value in result.evidence)
        lines.append(
            "\t".join(
                [
                    result.file,
                    document.name,
                    result.phenotype,
                    result.tier,
                    result.stage,
                    f"{result.score:.4f}",
                    evidence,
                ]
            )
        )
    return "\n".join(lines) + "\n"


def evidence_document(value: Call | str | None) -> Any:
    """One stage's evidence in the JSON result: call object, PCR string, or null."""
    if isinstance(value, str):
        return {"products": value.split("+") if value else []}
    if value is None or not value.phenotype:
        return None
    return {
        "phenotype": value.phenotype,
        "score": round(value.score, 4),
        "confidence": value.confidence,
        "scheme": value.scheme,
    }


def render_composite_json(
    results: list[FileCall], document: CompositeDocument, *, now: datetime | None = None
) -> str:
    """``gspec.composite_result/1`` — the fused call plus per-stage evidence."""
    stamp = utc_timestamp(now or datetime.now(UTC))
    payload = {
        "schema": COMPOSITE_RESULT_SCHEMA,
        "pipeline": document.name,
        "created_at": stamp,
        "stages": [
            {
                "name": stage.name,
                "kind": stage.kind,
                "db": stage.db,
                "scheme": stage.scheme,
                "tier": stage.tier,
                "min_score": stage.min_score,
                "require_confidence": stage.require_confidence,
                "primers": stage.primers,
            }
            for stage in document.stages
        ],
        "files": [
            {
                "file": result.file,
                "phenotype": result.phenotype,
                "tier": result.tier,
                "stage": result.stage,
                "score": round(result.score, 4),
                "evidence": {name: evidence_document(value) for name, value in result.evidence},
            }
            for result in results
        ],
    }
    return json.dumps(payload, indent=2) + "\n"


COMPOSITE_RESULT_SCHEMA = "gspec.composite_result/1"
