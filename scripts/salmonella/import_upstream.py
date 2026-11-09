#!/usr/bin/env python3
"""Import an upstream Salmonella antigen database into a gspec database.

Two upstream layouts are understood, and both are treated as *inputs to a
conversion*, never as content this repository ships:

- **SISTR** (Apache-2.0): ``data/antigens/{wzx,wzy,fliC,fljB}.fasta`` plus
  ``data/Salmonella-serotype_serogroup_antigen_table-WHO_2007.csv``. Record ids
  end in ``|<antigen name>`` — the convention SISTR's own
  ``get_antigen_name()`` splits on.
- **SeqSero2** (GPL-2.0): ``seqsero2_db/H_and_O_and_specific_genes.fasta``,
  whose ids are ``<gene>_<antigen>_<strain>_<epitope>__<len>`` (fliC/fljB) or
  ``O-<group>_<gene>_...__<len>`` (O antigen) plus ``tyr-O-2``/``tyr-O-9`` style
  special genes.

Licence reality (read before importing): SISTR's database is Apache-2.0 and may
be redistributed with attribution; SeqSero2's is GPL-2.0, so its content must
NOT travel with the MIT wheel — importing it produces a *local* database in a
datadir, which is the same relationship a user has with the upstream tool. The
generated ``info.json``/manifest carries the upstream licence, the source URL,
the conversion date and the objection clause.

Usage::

    # SISTR (dir or unpacked tarball root)
    python scripts/salmonella/import_upstream.py sistr \
        --source /path/to/SISTR_V_1.1.3_db --datadir ~/.local/share/gspec/db --name sal_sistr

    # SeqSero2 FASTA
    python scripts/salmonella/import_upstream.py seqsero2 \
        --source /path/to/H_and_O_and_specific_genes.fasta \\
        --formulas /path/to/serotype_table.csv \\
        --datadir ~/.local/share/gspec/db --name sal_seqsero

The formula table is required either way: SeqSero2 resolves serovar names in
code, gspec resolves them from data, so ``--formulas`` takes a CSV whose header
names a serovar column plus O/H1/H2 columns (case-insensitive; ``--formulas``
accepts SISTR's WHO-2007 CSV unchanged).
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections.abc import Iterator
from datetime import UTC, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from gspec.dbcodec import encode_seqid

FALLBACK_PHASE = "-"
#: SeqSero2 marks genes that only sharpen a call and must not become antigens
#: (e.g. the "not in 3,10" discriminator, or the O:1,3,19 short probe).
_SEQSERO2_SKIP = re.compile(r"not_in|partial$|_full$")


def _read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header, chunks = line[1:].strip(), []
        else:
            chunks.append(line.strip().upper())
    if header is not None:
        yield header, "".join(chunks)


def _gene_id(name: str) -> str:
    """A safe gspec gene id (no tabs/newlines, alias-free)."""
    cleaned = re.sub(r"[^A-Za-z0-9_.:-]+", "_", name.strip())
    return cleaned.strip("_") or "unnamed"


def _sistr_antigen(record_id: str) -> str:
    return record_id.rsplit("|", 1)[-1].strip()


def _seqsero2_antigen(record_id: str) -> tuple[str, str] | None:
    """(family, antigen) for one SeqSero2 record id, or None to skip it."""
    base = record_id.split("__")[0]
    if _SEQSERO2_SKIP.search(base):
        return None
    # SeqSero2 ids are ``<gene>_<antigen>_<strain>_<epitope>__<len>``; the
    # antigen is the one token right after the gene name (``fliC_g,m_...``).
    tokens = base.split("_")
    if tokens[0] == "fliC":
        return "h1", tokens[1]
    if tokens[0] == "fljB":
        return "h2", tokens[1]
    if base.startswith(("tyr-O-", "tyr_O-")):
        return "o", base.split("-", 2)[-1]
    reduced = re.sub(r"^O[-:]", "", base)
    reduced = re.sub(r"_(wzx|wzy|wzt|wba[A-Za-z]*|gene\d+).*$", "", reduced)
    if not reduced or not re.match(r"^[0-9]", reduced):
        return None
    return "o", reduced.replace("_", "")


def _detect_family(path: Path) -> str | None:
    stem = path.stem.lower()
    if stem.startswith("fljb"):
        return "h2"
    if stem.startswith("flic"):
        return "h1"
    if stem.startswith(("wzx", "wzy")):
        return "o"
    return None


def _formula_rows(path: Path) -> list[dict[str, str]]:
    """CSV rows normalized to ``{serovar, o, h1, h2}`` (header names guessed)."""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise SystemExit(f"{path}: empty formula table")
        keys = {name.strip().lower(): name for name in reader.fieldnames}

        def pick(*candidates: str) -> str | None:
            for candidate in candidates:
                for lower, original in keys.items():
                    if candidate in lower:
                        return original
            return None

        serovar_col = pick("serovar")
        o_col = pick("serogroup", "o_antigen", "o antigen", "o group")
        h1_col = pick("h1", "phase 1", "antigen_h1")
        h2_col = pick("h2", "phase 2", "antigen_h2")
        if serovar_col is None or o_col is None or h1_col is None or h2_col is None:
            raise SystemExit(
                f"{path}: could not find serovar/O/H1/H2 columns in {reader.fieldnames}"
            )
        rows = []
        for row in reader:
            serovar = (row.get(serovar_col) or "").strip()
            if not serovar:
                continue
            rows.append(
                {
                    "serovar": serovar,
                    "o": (row.get(o_col) or FALLBACK_PHASE).strip() or FALLBACK_PHASE,
                    "h1": (row.get(h1_col) or FALLBACK_PHASE).strip() or FALLBACK_PHASE,
                    "h2": (row.get(h2_col) or FALLBACK_PHASE).strip() or FALLBACK_PHASE,
                }
            )
        return rows


def _write_database(
    *,
    name: str,
    source_label: str,
    licence: str,
    source_url: str,
    markers: dict[str, list[tuple[str, str]]],
    formulas: list[dict[str, str]],
    dest: Path,
) -> None:
    """Write ``sequences`` + ``serotype.json`` + ``info.json`` for one database."""
    dest.mkdir(parents=True, exist_ok=True)
    families = [
        {
            "key": key,
            "entries": [
                {"name": antigen, "genes": [gene]}
                for gene, antigen in sorted(markers[key], key=lambda item: item[1])
            ],
        }
        for key in ("o", "h1", "h2")
        if markers.get(key)
    ]
    document = {
        "schema": "gspec.serotype/2",
        "name": name,
        "control_gene": None,
        "identity_floor": 80.0,
        "coverage_floor": 60.0,
        "margin": 0.01,
        "antigens": families,
        "formulas": formulas,
        "fallback": "Not Salmonella",
        "unresolved": "Salmonella subsp. enterica (no antigen set resolved)",
        "partial": "Salmonella {formula} (serovar not in the formula table)",
    }
    (dest / "serotype.json").write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    (dest / "info.json").write_text(
        json.dumps(
            {
                "schema": "gspec.dbinfo/1",
                "name": name,
                "description": (
                    f"{source_label}: {sum(len(v) for v in markers.values())} antigen"
                    f" references, {len(formulas)} serovar formulas"
                ),
                "vendor": f"converted from {source_label}",
                "note": (
                    f"Converted by scripts/salmonella/import_upstream.py on"
                    f" {datetime.now(UTC).strftime('%Y-%m-%d')} from {source_url}."
                    f" Upstream licence: {licence}. This is a converted copy, not"
                    " gspec-curated content: verify the upstream terms before"
                    " redistributing it. If any author or rights holder objects to"
                    " the use or redistribution of this database, please contact"
                    " indexofire@gmail.com and the content will be removed or"
                    " adjusted promptly upon verification."
                ),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"{name}: {sum(len(v) for v in markers.values())} antigen refs, {len(formulas)} formulas")


def import_sistr(source: Path, dest: Path, *, name: str, formulas: Path | None) -> None:
    antigen_dir = source / "antigens" if (source / "antigens").is_dir() else source
    files = sorted(antigen_dir.glob("*.fasta"))
    if not files:
        raise SystemExit(f"{source}: no antigens/*.fasta found")
    table = formulas or (source / "Salmonella-serotype_serogroup_antigen_table-WHO_2007.csv")
    if not table.is_file():
        raise SystemExit(f"{table}: formula table not found (pass --formulas)")
    pairs: list[tuple[str, str, str]] = []
    markers: dict[str, list[tuple[str, str]]] = {"o": [], "h1": [], "h2": []}
    for path in files:
        family = _detect_family(path)
        if family is None:
            continue
        for record_id, sequence in _read_fasta(path):
            antigen = _sistr_antigen(record_id)
            gene = _gene_id(f"{family}_{antigen}")
            markers[family].append((gene, antigen))
            pairs.append((gene, antigen, sequence))
    _write_fasta(dest, name, pairs)
    _write_database(
        name=name,
        source_label="SISTR antigen database",
        licence="Apache-2.0",
        source_url="https://zenodo.org/records/15078586",
        markers=markers,
        formulas=_formula_rows(table),
        dest=dest,
    )


def import_seqsero2(source: Path, dest: Path, *, name: str, formulas: Path) -> None:
    if not source.is_file():
        raise SystemExit(f"{source}: not a file")
    markers: dict[str, list[tuple[str, str]]] = {"o": [], "h1": [], "h2": []}
    pairs: list[tuple[str, str, str]] = []
    for record_id, sequence in _read_fasta(source):
        parsed = _seqsero2_antigen(record_id)
        if parsed is None:
            continue
        family, antigen = parsed
        gene = _gene_id(f"{family}_{antigen}")
        markers[family].append((gene, antigen))
        pairs.append((gene, antigen, sequence))
    _write_fasta(dest, name, pairs)
    _write_database(
        name=name,
        source_label="SeqSero2 antigen database",
        licence="GPL-2.0 (code and database; do not redistribute with the MIT wheel)",
        source_url="https://github.com/denglab/SeqSero2/blob/main/seqsero2_db",
        markers=markers,
        formulas=_formula_rows(formulas),
        dest=dest,
    )


def _write_fasta(dest: Path, name: str, pairs: list[tuple[str, str, str]]) -> None:
    """gspec/v1 FASTA: id carries db/gene/accession/function, description the antigen."""
    dest.mkdir(parents=True, exist_ok=True)
    with (dest / "sequences").open("w", encoding="utf-8", newline="\n") as handle:
        for gene, antigen, sequence in pairs:
            header = encode_seqid(name, gene, f"upstream:{antigen}", ("antigen",))
            handle.write(f">{header} antigen {antigen}\n")
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset : offset + 60] + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("upstream", choices=("sistr", "seqsero2"))
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--dest", type=Path, required=True, help="database directory to create")
    parser.add_argument("--name", required=True)
    parser.add_argument("--formulas", type=Path, help="serovar formula CSV (O, H1, H2, serovar)")
    args = parser.parse_args(argv)
    args.dest.mkdir(parents=True, exist_ok=True)
    if args.upstream == "sistr":
        import_sistr(args.source, args.dest, name=args.name, formulas=args.formulas)
    else:
        if args.formulas is None:
            raise SystemExit("seqsero2 import needs --formulas (SeqSero2 names serovars in code)")
        import_seqsero2(args.source, args.dest, name=args.name, formulas=args.formulas)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
