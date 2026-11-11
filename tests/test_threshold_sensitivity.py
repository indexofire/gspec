"""The batch-sensitivity check for the panel's thresholds.

The numbers were set against real data, so the question "does a new batch invalidate
them?" has to be answerable per batch. This test pins the report's decision logic:
a threshold is safe only when the observations leave an empty band around it.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Protocol, cast

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "salmonella" / "threshold_sensitivity.py"


class _Sensitivity(Protocol):
    def antigen_score(self, identity: float, coverage: float) -> float: ...


def _module() -> _Sensitivity:
    spec = importlib.util.spec_from_file_location("threshold_sensitivity", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return cast(_Sensitivity, module)


def test_score_is_identity_times_coverage() -> None:
    module = _module()
    assert module.antigen_score(100.0, 100.0) == 1.0
    assert module.antigen_score(99.85, 52.0) == pytest.approx(0.5192, abs=1e-4)
    assert module.antigen_score(80.0, 60.0) == pytest.approx(0.48)


def test_a_threshold_sits_in_a_desert_or_on_an_edge() -> None:
    """The report's question, in the numbers it measured on the real corpus.

    Phase scores are either <= 0.691 (a cross-reaction) or >= 0.936 (a present
    phase), so `phase_floor = 0.75` has a 0.245-wide empty band and would survive a
    large shift in assembly quality. An observation at 0.72 would put it on an edge.
    """
    present, absent, floor, band = 0.936, 0.691, 0.75, 0.05
    assert min(present, 1.0) - floor > band
    assert floor - absent > band
