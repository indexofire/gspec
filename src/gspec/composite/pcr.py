"""The `pcr` composite stage: primer-set loading, the per-file in-silico PCR scan
(gseek.ispcr) and the product-size rules that refine a base call.
"""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

from gspec.composite.models import CompositeDocument, CompositeStage
from gspec.errors import DatabaseError, DependencyError


def primer_path(stage: CompositeStage, document: CompositeDocument) -> Path:
    """A stage's primer file, resolved next to the pipeline document."""
    path = Path(stage.primers or "")
    if path.is_absolute() or document.source_dir is None:
        return path
    return document.source_dir / path


def load_primers(stage: CompositeStage, document: CompositeDocument) -> tuple[Any, ...]:
    """Load a PCR stage's primer pairs through gseek (lazy, typed error)."""
    try:
        import gseek
    except ImportError as exc:  # pragma: no cover - dependency gate
        raise DependencyError(
            "pcr stages need the gseek package with ispcr support (install gseek>=0.10.5)",
            code="MISSING_DEPENDENCY",
            context={"stage": stage.name},
        ) from exc
    if not hasattr(gseek, "read_primers"):
        raise DependencyError(
            "the installed gseek has no in-silico PCR support (need gseek>=0.10.5)",
            code="MISSING_DEPENDENCY",
            context={"stage": stage.name, "gseek": getattr(gseek, "__version__", "?")},
        )
    pairs = tuple(gseek.read_primers(primer_path(stage, document)))
    if stage.report_all_products:
        return pairs
    needed = {product.name for rule in stage.rules for product in rule.products}
    known = {pair.name for pair in pairs}
    missing = sorted(needed - known)
    if missing:
        raise DatabaseError(
            f"pcr stage {stage.name!r} names product(s) the primer set lacks: {', '.join(missing)}",
            code="COMPOSITE_MALFORMED",
            context={"stage": stage.name, "products": ", ".join(missing)},
        )
    return tuple(pair for pair in pairs if pair.name in needed)


def primer_worker(
    payload: tuple[str, tuple[Any, ...], int, int],
) -> tuple[str, dict[str, list[int]]]:
    """One file's in-silico PCR, in a worker process.

    Module-level (not a closure) so the payload pickles: primer matching is a
    CPU-bound byte scan that does not release the GIL, so threads bought
    nothing here and processes are what makes ``--jobs`` pay.
    """
    import gseek

    path, primers, max_mismatches, max_edits = payload
    if not primers:
        # Nothing this file's call could need: skip the scan entirely.
        return path, {}
    result = gseek.ispcr(Path(path), primers, max_mismatches=max_mismatches, max_edits=max_edits)
    products: dict[str, list[int]] = {}
    for amplicon in result.amplicons:
        products.setdefault(amplicon.pair, []).append(amplicon.size)
    return path, products


def needed_products(stage: CompositeStage, phenotype: str | None) -> set[str]:
    """The products this call could possibly need.

    A rule with ``when`` can only fire on a matching base type, so a genome
    called HS2 does not need the HS4A/HS4B/HS5_31 scans. Unconditional rules
    (the assay's own overrides) are always scanned.
    """
    base = base_type(phenotype) if phenotype else None
    needed: set[str] = set()
    for rule in stage.rules:
        if rule.when and (base is None or base not in {base_type(item) for item in rule.when}):
            continue
        needed.update(product.name for product in rule.products)
    return needed


def primer_products(
    stage: CompositeStage,
    paths: list[Path],
    primers: tuple[Any, ...],
    jobs: int,
    *,
    phenotypes: dict[str, str | None] | None = None,
) -> dict[str, dict[str, list[int]]]:
    """Positive products per file: ``{file: {pair: [sizes]}}``.

    ``jobs > 1`` fans the scan out over worker PROCESSES — the primer search is
    a pure-Python byte scan holding the GIL, so a thread pool would serialise
    it (measured 1.01x); processes give the expected ~Nx.
    """
    payloads = []
    for path in paths:
        subset = primers
        if not stage.report_all_products and phenotypes is not None:
            needed = needed_products(stage, phenotypes.get(str(path)))
            subset = tuple(pair for pair in primers if pair.name in needed)
        payloads.append((str(path), subset, stage.max_mismatches, stage.max_edits))
    if jobs <= 1 or len(payloads) <= 1:
        return dict(primer_worker(payload) for payload in payloads)
    try:
        with ProcessPoolExecutor(max_workers=min(jobs, len(payloads))) as executor:
            return dict(executor.map(primer_worker, payloads))
    except (OSError, RuntimeError):  # pragma: no cover - restricted sandboxes
        # A platform that forbids extra processes must not break the pipeline:
        # fall back to the sequential scan (same results, slower).
        return dict(primer_worker(payload) for payload in payloads)


def base_type(phenotype: str) -> str:
    """The call without its parenthetical bracket (``HS4c (… complex)`` → ``HS4c``)."""
    return phenotype.split(" (", 1)[0].strip()


def primer_refine(
    stage: CompositeStage,
    current: tuple[str, str, float],
    products: dict[str, list[int]],
) -> tuple[str, str, float] | None:
    """Apply the stage's rules to one file's call; None when nothing applies."""
    phenotype = current[0]
    base = base_type(phenotype)
    for rule in stage.rules:
        if rule.when and base not in {base_type(when) for when in rule.when}:
            continue
        if not all(
            any(
                product.min_size <= size <= product.max_size
                for size in products.get(product.name, ())
            )
            for product in rule.products
        ):
            continue
        size = next(
            size
            for size in products[rule.products[0].name]
            if rule.products[0].min_size <= size <= rule.products[0].max_size
        )
        return (rule.set, stage.name, float(size))
    return None
