"""The decision tool for "should we add a similarity channel?".

The script answers one question from a corpus: does a Jaccard threshold put every
same-serovar pair above every different-serovar pair? Its measured answer on 125
Salmonella assemblies was "no" (the most similar different-serovar pair, 0.997, is
closer than the most similar same-serovar pair, 0.979), which is why gspec does not
carry a reference-genome sketch channel.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Protocol, cast

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "salmonella" / "sketch_separability.py"


class _Separability(Protocol):
    """The script surface these tests exercise (it lives outside the package)."""

    def jaccard(self, first: frozenset[int], second: frozenset[int]) -> float: ...

    def sketch(self, path: Path, *, size: int, k: int) -> frozenset[int]: ...

    def separable(self, within: list[float], between: list[float]) -> bool: ...


def _module() -> _Separability:
    spec = importlib.util.spec_from_file_location("sketch_separability", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return cast(_Separability, module)


def test_jaccard_bounds() -> None:
    module = _module()
    assert module.jaccard(frozenset({1, 2, 3}), frozenset({1, 2, 3})) == 1.0
    assert module.jaccard(frozenset({1, 2}), frozenset({3, 4})) == 0.0
    assert module.jaccard(frozenset(), frozenset({1})) == 0.0


def test_sketch_is_bounded_and_keeps_the_smallest_codes() -> None:
    import tempfile

    module = _module()
    sequence = "ACGT" * 500
    with tempfile.TemporaryDirectory() as scratch:
        path = Path(scratch) / "tiny.fna"
        path.write_text(f">contig\n{sequence}\n", encoding="utf-8")
        codes = module.sketch(path, size=50, k=21)
    assert len(codes) <= 50
    assert codes == frozenset(sorted(codes)[:50])


def test_separable_needs_a_clean_gap() -> None:
    module = _module()
    assert module.separable([0.9, 0.95], [0.5, 0.7]) is True
    # one different-serovar pair above the closest same-serovar pair: no threshold
    assert module.separable([0.9, 0.95], [0.5, 0.97]) is False
    assert module.separable([], [0.5]) is False
