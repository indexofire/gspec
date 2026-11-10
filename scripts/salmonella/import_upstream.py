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
from typing import cast

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from gspec.dbcodec import encode_seqid

HERE = Path(__file__).resolve().parent

FALLBACK_PHASE = "-"
#: Reference alleles of one antigen that are at least this identical are the same
#: evidence: keeping them all makes a reads table split its reads across several
#: references (each then misses the breadth floor) without adding resolution.
# Exact duplicates only. Collapsing at 99 % used to keep one representative per
# antigen, which silently deleted the allele a genome actually carries: the top
# hit then fell to a *sibling* antigen one SNP away (Enteritidis' fliC `g,m` lost
# to `g,q` at 99.93 %), and with `top_result()`-style scoring that sibling wins on
# bitscore noise. An exact allele must survive so the call can be perfect and the
# sibling cannot outrank it. The reads path is unaffected: the k-mer substrate
# already reduces each antigen to one representative.
REFERENCE_IDENTITY_FLOOR = 100.0
#: SISTR's ``{SEROGROUP,H1,H2}_*_SIMILARITY_GROUPS``, ported verbatim from
#: ``sistr/src/serovar_prediction/constants.py`` (SISTR is Apache-2.0; the
#: attribution in ``bundled.json`` names it). They are the upstream answer to
#: WKL bracket tokens: two tokens in a class are the same antigen, so a genome
#: whose allele matches one of them still resolves the row.
SIMILARITY_GROUPS = HERE / "sistr_similarity_groups.json"
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


def _reference_ids(
    markers: dict[str, dict[str, list[str]]], sequences: dict[str, str]
) -> dict[str, dict[str, list[str]]]:
    """Assign deterministic gene ids: ``<family>_<antigen>``, then ``__2``, ``__3``.

    Numbering is per antigen in input order, so a rebuild produces the same ids
    (a global counter shifted them whenever an unrelated antigen changed). Two
    antigen names can sanitise to the same token; the later name (by sort order)
    gets a ``__g<n>`` prefix tag so both stay unique and stable.
    """
    bases: dict[str, list[str]] = {}
    for family, antigens in markers.items():
        for antigen in antigens:
            bases.setdefault(_gene_id(f"{family}_{antigen}"), []).append(antigen)
    tagged: dict[tuple[str, str], str] = {}
    for base, antigens in bases.items():
        for index, antigen in enumerate(sorted(antigens)):
            family = base.split("_", 1)[0]
            tagged[(family, antigen)] = base if index == 0 else f"{base}__g{index + 1}"
    out: dict[str, dict[str, list[str]]] = {}
    for family, antigens in markers.items():
        out[family] = {}
        for antigen, old_genes in antigens.items():
            prefix = tagged[(family, antigen)]
            ids = [prefix] + [f"{prefix}__{n}" for n in range(2, len(old_genes) + 1)]
            for new, old in zip(ids, old_genes, strict=True):
                sequences[new] = sequences.pop(old)
            out[family][antigen] = ids
    return out


def _gene_id(name: str) -> str:
    """A safe gspec gene id (no tabs/newlines, alias-free)."""
    cleaned = re.sub(r"[^A-Za-z0-9_.:-]+", "_", name.strip())
    return cleaned.strip("_") or "unnamed"


def _norm(text: str) -> str:
    """Collapse whitespace (upstream tokens carry trailing/double spaces)."""
    return " ".join(text.split())


def _sistr_antigen(record_id: str) -> str:
    return _norm(record_id.rsplit("|", 1)[-1])


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


def _sistr_formulas(path: Path, *, subspecies: str | None) -> list[dict[str, str]]:
    """SISTR's WHO-2007 table as formula rows.

    Two SISTR semantics are preserved: ``subspecies`` narrows the rows (the
    antigen genes cannot decide the subspecies, so ``enterica`` is the default
    scope), and ``can_h2_be_missing`` adds a phase-off alias row so a genome
    caught in phase 2-off still resolves to the same serovar.
    """
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows: list[dict[str, str]] = []
        for row in reader:
            row_spp = (row.get("subspecies") or "").strip()
            if subspecies is not None and row_spp and row_spp != subspecies:
                # A table without the column (older exports, hand-made fixtures)
                # keeps every row rather than silently dropping all of them.
                continue
            serovar = (row.get("Serovar") or "").strip()
            if not serovar:
                continue
            o = (row.get("Serogroup") or FALLBACK_PHASE).strip() or FALLBACK_PHASE
            h1 = (row.get("H1") or FALLBACK_PHASE).strip() or FALLBACK_PHASE
            h2 = (row.get("H2") or FALLBACK_PHASE).strip() or FALLBACK_PHASE
            rows.append({"o": o, "h1": h1, "h2": h2, "serovar": serovar})
            if (
                row.get("can_h2_be_missing") or ""
            ).strip().upper() == "TRUE" and h2 != FALLBACK_PHASE:
                rows.append({"o": o, "h1": h1, "h2": FALLBACK_PHASE, "serovar": serovar})
    return rows


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


def _write_bundle_metadata(
    dest: Path,
    *,
    name: str,
    source_label: str,
    licence: str,
    source_url: str,
    attribution: str,
    dropped_serovars: set[str] | None = None,
    dbtype: str = "nucl",
) -> None:
    """``bundled.json`` + ``screen.json``: make ``dest`` a wheel-eligible bundle.

    Only valid for an upstream dataset whose licence permits redistribution
    (AGENTS.md §6 gate 1/2); the attribution is kept verbatim because the
    Apache-2.0 NOTICE requirement travels with a conversion of the content.
    """
    (dest / "bundled.json").write_text(
        json.dumps(
            {
                "name": name,
                "description": f"{source_label}: O/H antigen panel with its formula table",
                "vendor": f"converted from {source_label}",
                "dbtype": dbtype,
                "kind": "gene",
                "snapshotted": datetime.now(UTC).strftime("%Y-%m-%d"),
                "license": licence,
                "attribution": attribution,
                "note": (
                    f"Converted by gspec's scripts/salmonella/import_upstream.py from"
                    f" {source_url} (retrieved "
                    f"{datetime.now(UTC).strftime('%Y-%m-%d')}). The antigen references"
                    " and the O:H1:H2 formula table are the upstream project's data,"
                    " converted to gspec/v1 records and gspec.serotype/2; gspec claims"
                    " no authorship of the content and does not imply upstream"
                    " endorsement of this conversion."
                ),
                "guidance": {
                    "references": [
                        {
                            "citation": (
                                "Yoshida CE, Kruczkiewicz P, Laing CR, Lingohr EJ,"
                                " Gannon VPJ, Nash JHE, Taboada EN. The Salmonella In"
                                " Silico Typing Resource (SISTR): an open"
                                " web-accessible tool for rapidly typing and"
                                " subtyping draft microbial genomes. PLoS One."
                                " 2016;11(8):e0160528."
                            ),
                            "doi": "10.1371/journal.pone.0160528",
                            "note": (
                                "SISTR: the antigen references and the"
                                " O:H1:H2 formula table this database was"
                                " converted from."
                            ),
                        },
                        {
                            "citation": (
                                "Grimont PAD, Weill FX. Antigenic formulae of the"
                                " Salmonella serovars (Kauffmann-White-Le Minor"
                                " scheme), WHO Collaborating Centre for Reference and"
                                " Research on Salmonella, Institut Pasteur."
                            ),
                            "note": (
                                "The antigenic-formula scheme the table encodes;"
                                " bracket and parenthesis notation in the tokens is"
                                " the scheme's."
                            ),
                        },
                    ],
                    "caveats": [
                        "Converted upstream data: see `attribution` for the upstream"
                        " project, version and URL.",
                        "Antigen determination is a ranking step: a margin below the"
                        " measured gap between near-identical alleles is required, so"
                        " re-calibrate after an upstream refresh.",
                        *(
                            [
                                f"{len(dropped_serovars)} upstream serovar rows were"
                                " dropped at conversion because no reference allele"
                                " can reach their formula: "
                                + ", ".join(sorted(dropped_serovars)[:10])
                                + ("…" if len(dropped_serovars) > 10 else "")
                            ]
                            if dropped_serovars
                            else []
                        ),
                    ],
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (dest / "screen.json").write_text(
        json.dumps(
            {
                # screen/2 is required for the evidence-model knobs: under the
                # /1 schema `dedup`/`coverage` are not part of the document and
                # were silently ignored, so the panel still deduped per query and
                # dropped sibling antigen references (measured: H2 vanished on
                # real genomes and group-B calls landed on monophasic rows).
                "schema": "gspec.screen/2",
                # An antigen panel is a set of near-identical references: the
                # default per-query dedup would keep only the best-matching
                # sibling and hide the rest, so the family candidates collapse to
                # one name (measured on real genomes: H2 vanished and group-B
                # calls landed on monophasic rows). Subject dedup keeps every
                # reference's own HSPs, which is what the ranking step needs.
                #
                # ``coverage: union`` is not optional either. The parity path
                # passes ``-culling_limit 1``, which keeps one HSP per *query
                # interval*; on a panel where dozens of references align to one
                # locus that silently reduced 1552 HSPs to 1 (measured on
                # Enteritidis: only ``h1_g_q`` survived and ``g,m`` -- the allele
                # the genome actually carries at 100 % identity -- was dropped).
                # Union coverage drops culling and aggregates per subject.
                "contigs": {
                    "minid": 80.0,
                    # Screen LOW and let the serotype document decide presence:
                    # antigen genes sit in repeat-rich regions, so short-read
                    # assemblies break them across contigs. A real D1 genome held
                    # wzx at 99.85 % identity over 52 % union coverage (split
                    # across contig00291/contig00348) and a 60 % screen floor
                    # silently erased the serogroup. `serotype.json` carries the
                    # strict floor plus a near-identity partial rule
                    # (partial_coverage_floor/partial_identity_floor).
                    "mincov": 45.0,
                    "coverage": "union",
                    "merge_fragments": False,
                    "dedup": "subject",
                    # Anchored allele sketches: the alignment-free evidence the
                    # serotype ranker uses to break near-ties. Measured on 60
                    # genomes: 78.7 % -> 94.7 % agreement with SISTR's H tokens.
                    "allele_sketch": True,
                },
                "reads": {"min_breadth": 60.0, "min_identity": 80.0},
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def _merge_cross_labels(
    groups: dict[str, list[list[str]]], cross_labels: Path
) -> dict[str, list[list[str]]]:
    """Fold audited cross-labelled classes into the similarity groups.

    A pair of labels that cannot be told apart by k-mer containment (`g,q` versus
    `g,m,q`, `i` versus `g,t`) must land in one widening class, otherwise the
    top-scoring name's rows are the only ones reachable and a genome carrying the
    sibling name loses its serovar rows entirely.
    """
    audited = json.loads(cross_labels.read_text(encoding="utf-8"))
    added = 0
    for family, classes in audited["groups"].items():
        existing = [list(group) for group in groups.get(family, [])]
        for component in classes:
            members = [name for name in component if name]
            if len(members) < 2:
                continue
            overlapped = [group for group in existing if set(group) & set(members)]
            merged = set(members)
            for group in overlapped:
                merged.update(group)
                existing.remove(group)
            existing.append(sorted(merged))
            added += 1
        groups[family] = existing
    print(f"  cross-labelled classes merged: {added}")
    return groups


def _add_o_discriminators(
    document: Path,
    markers: dict[str, dict[str, list[str]]],
    pairs: list[tuple[str, str, str]],
) -> dict[str, object]:
    """Add the curated O9/O2 marker pair and return its ``serotype/2`` block.

    The markers are NOT antigen entries: they must never widden the O family's
    candidate list or be reported as an antigen. They live as standalone records
    whose gene ids the document points at, and the engine only consults them when
    the O tie actually spans their names.
    """
    spec = json.loads(document.read_text(encoding="utf-8"))
    favours: dict[str, str] = {}
    for reference in spec["references"]:
        gene = f"o_disc_{reference['group']}"
        pairs.append((gene, f"discriminator:{reference['group']}", reference["sequence"]))
        for antigen in spec["favours"][reference["group"]]:
            favours[antigen] = gene
    return {
        "favours": favours,
        "min_gap_pct": spec["min_gap"],
        "min_coverage": 60.0,
        "note": spec["description"],
    }


def _write_database(
    *,
    name: str,
    source_label: str,
    licence: str,
    source_url: str,
    markers: dict[str, dict[str, list[str]]],
    formulas: list[dict[str, str]],
    dest: Path,
    similarity_groups: dict[str, list[list[str]]] | None = None,
    o_discriminator: dict[str, object] | None = None,
) -> None:
    """Write ``sequences`` + ``serotype.json`` + ``info.json`` for one database."""
    dest.mkdir(parents=True, exist_ok=True)
    families = [
        {
            "key": key,
            "entries": [
                {"name": antigen, "genes": genes} for antigen, genes in sorted(markers[key].items())
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
        **({"similarity_groups": similarity_groups} if similarity_groups else {}),
        **({"o_discriminator": o_discriminator} if o_discriminator else {}),
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
                    f"{source_label}:"
                    f" {sum(len(genes) for fam in markers.values() for genes in fam.values())}"
                    f" antigen references for"
                    f" {sum(len(fam) for fam in markers.values())} antigens,"
                    f" {len(formulas)} serovar formulas"
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
    print(
        f"{name}:"
        f" {sum(len(genes) for fam in markers.values() for genes in fam.values())} refs,"
        f" {sum(len(fam) for fam in markers.values())} antigens,"
        f" {len(formulas)} formulas"
    )


def _collapse_references(
    markers: dict[str, dict[str, list[str]]],
    sequences: dict[str, str],
    *,
    floor: float = REFERENCE_IDENTITY_FLOOR,
) -> tuple[dict[str, dict[str, list[str]]], dict[str, str]]:
    """Drop near-identical references of the same antigen (keep the longest).

    The floor defaults to 100 %, i.e. only exact duplicates of one antigen are
    dropped. Anything looser deletes the allele a genome carries and hands the
    call to a sibling antigen that differs by a SNP (see the floor constant); the
    *reads* path is not a reason to collapse, because the k-mer substrate already
    keeps one representative per antigen.
    """
    import difflib

    # A 15-mer containment prefilter (0.99 ** 15 ~= 0.86) rejects the clearly
    # different pairs before the quadratic alignment, which took the whole panel
    # from minutes to seconds; the survivors are confirmed by alignment.
    prefilter_k = 15
    prefilter_floor = (floor / 100.0) ** prefilter_k

    def fingerprints(sequence: str) -> set[str]:
        return {
            sequence[i : i + prefilter_k]
            for i in range(len(sequence) - prefilter_k + 1)
            if set(sequence[i : i + prefilter_k]) <= set("ACGT")
        }

    kept_markers: dict[str, dict[str, list[str]]] = {key: {} for key in markers}
    kept_sequences: dict[str, str] = {}
    dropped = 0
    for family, antigens in markers.items():
        for antigen, genes in antigens.items():
            keep: list[str] = []
            prints: dict[str, set[str]] = {}
            for gene in sorted(genes, key=lambda name: len(sequences[name]), reverse=True):
                sequence = sequences[gene]
                print_set = fingerprints(sequence)
                near = [
                    other
                    for other in keep
                    if prints[other]
                    and len(print_set & prints[other]) / len(print_set | prints[other])
                    >= prefilter_floor
                ]
                if any(
                    difflib.SequenceMatcher(
                        None, sequence, sequences[other], autojunk=False
                    ).ratio()
                    * 100.0
                    >= floor
                    for other in near
                ):
                    dropped += 1
                    continue
                keep.append(gene)
                prints[gene] = print_set
                kept_sequences[gene] = sequence
            kept_markers[family][antigen] = keep
    print(
        f"  references: {sum(len(v) for fam in markers.values() for v in fam.values())}"
        f" -> {len(kept_sequences)} (collapsed {dropped} at >={floor:.0f}% identity)"
    )
    return kept_markers, kept_sequences


def _renumber(
    markers: dict[str, dict[str, list[str]]], sequences: dict[str, str]
) -> tuple[dict[str, dict[str, list[str]]], list[tuple[str, str, str]]]:
    """Compact gene ids after collapsing: ``<family>_<antigen>``, ``__2``, ``__3``.

    The pre-collapse numbering left gaps (``h2_1_2__65``) because it counted the
    upstream records that were then dropped; renumbering the survivors keeps the
    ids short and stable for a given input.
    """
    markers_out: dict[str, dict[str, list[str]]] = {}
    pairs: list[tuple[str, str, str]] = []
    for family, antigens in markers.items():
        markers_out[family] = {}
        used: set[str] = set()
        bases = {_gene_id(f"{family}_{antigen}"): antigen for antigen in antigens}
        collisions = {
            base: sorted(name for name, other in bases.items() if other == antigen)
            for base, antigen in bases.items()
            if sum(1 for other in bases.values() if other == antigen) > 1
        }
        for antigen, genes in antigens.items():
            base = _gene_id(f"{family}_{antigen}")
            if base in collisions:
                position = collisions[base].index(antigen)
                base = base if position == 0 else f"{base}__g{position + 1}"
            ids = [base] + [f"{base}__{n}" for n in range(2, len(genes) + 1)]
            # Last line of defence: an antigen name can sanitise onto another
            # name's id even after the collision tag above.
            for position, gene_id in enumerate(ids):
                unique = gene_id
                suffix = 2
                while unique in used:
                    unique = f"{gene_id}__u{suffix}"
                    suffix += 1
                used.add(unique)
                ids[position] = unique
            markers_out[family][antigen] = ids
            for gene, new in zip(genes, ids, strict=True):
                sequences[new] = sequences[gene]
                pairs.append((new, antigen, sequences[new]))
    return markers_out, pairs


def _drop_unreachable(
    families: dict[str, dict[str, list[str]]],
    formulas: list[dict[str, str]],
    *,
    groups: dict[str, list[list[str]]] | None = None,
) -> tuple[list[dict[str, str]], set[str]]:
    """Keep only rows an antigen reference can actually reach.

    A row whose token has no reference sequence (directly or through a
    similarity class) can never fire, so shipping it is dead weight that also
    inflates the apparent coverage. The dropped serovars are returned so the
    bundle can record them instead of losing them silently.
    """
    known = {key: set(family) | {FALLBACK_PHASE} for key, family in families.items()}

    def reachable(family: str, token: str) -> bool:
        if token in known[family]:
            return True
        for group in (groups or {}).get(family, []):
            if token in group and any(member in known[family] for member in group):
                return True
        return False

    kept: list[dict[str, str]] = []
    dropped: set[str] = set()
    for row in formulas:
        if reachable("o", row["o"]) and reachable("h1", row["h1"]) and reachable("h2", row["h2"]):
            kept.append(row)
        else:
            dropped.add(row["serovar"])
    return kept, dropped


def _report_reachability(
    families: dict[str, dict[str, list[str]]],
    formulas: list[dict[str, str]],
    *,
    groups: dict[str, list[list[str]]] | None = None,
) -> None:
    """How many formula rows can actually be reached from the declared antigens.

    A row whose O/H1/H2 token has no reference sequence is unreachable: the
    genome can never resolve that token, so the row is dead weight. Printing the
    number (and the worst offenders) keeps the gap visible instead of shipping a
    table that silently never fires.
    """
    known = {key: set(family) | {FALLBACK_PHASE} for key, family in families.items()}

    def reachable_token(family: str, token: str) -> bool:
        if token in known[family]:
            return True
        for group in (groups or {}).get(family, []):
            if token in group and any(member in known[family] for member in group):
                return True
        return False

    dead = [
        row
        for row in formulas
        if not (
            reachable_token("o", row["o"])
            and reachable_token("h1", row["h1"])
            and reachable_token("h2", row["h2"])
        )
    ]
    total = len(formulas)
    if total == 0:
        print("  formula rows: 0 (the table produced no rows — check the input columns)")
        return
    reachable = total - len(dead)
    print(f"  formula rows: {total} reachable: {reachable} ({100 * reachable / total:.1f}%)")
    if dead:
        sample = sorted({(row["o"], row["h1"], row["h2"], row["serovar"]) for row in dead})[:5]
        print(f"  unreachable rows reference antigens with no reference sequence, e.g. {sample}")


def import_sistr(
    source: Path,
    dest: Path,
    *,
    name: str,
    formulas: Path | None,
    bundle: bool = False,
    subspecies: str | None = "enterica",
    extra_h_from: Path | None = None,
    kmer_k: int = 0,
    collapse_floor: float = REFERENCE_IDENTITY_FLOOR,
    o_discriminators: Path | None = None,
    cross_labels: Path | None = None,
    extra_refs: Path | None = None,
) -> None:
    antigen_dir = source / "antigens" if (source / "antigens").is_dir() else source
    files = sorted(antigen_dir.glob("*.fasta"))
    if not files:
        raise SystemExit(f"{source}: no antigens/*.fasta found")
    table = formulas or (source / "Salmonella-serotype_serogroup_antigen_table-WHO_2007.csv")
    if not table.is_file():
        raise SystemExit(f"{table}: formula table not found (pass --formulas)")
    pairs: list[tuple[str, str, str]] = []
    markers: dict[str, dict[str, list[str]]] = {"o": {}, "h1": {}, "h2": {}}
    for path in files:
        family = _detect_family(path)
        if family is None:
            continue
        for record_id, sequence in _read_fasta(path):
            antigen = _sistr_antigen(record_id)
            if not antigen:
                continue
            # Every upstream allele reference is kept: an antigen covered by
            # several sequences scores by the best of them, exactly as the
            # upstream top-hit step does.
            genes = markers[family].setdefault(antigen, [])
            gene = (
                _gene_id(f"{family}_{antigen}")
                if not genes
                else (f"{_gene_id(f'{family}_{antigen}')}__{len(genes) + 1}")
            )
            genes.append(gene)
            pairs.append((gene, antigen, sequence))
    if extra_h_from is not None:
        # A hybrid: SISTR's O groups and formula table, plus a second project's
        # H alleles. Both projects name the H antigens in WKL notation, so the
        # extra references only widen each antigen's evidence; the O family must
        # stay SISTR's (its serogroup letters are what the table keys on).
        if bundle:
            raise SystemExit(
                "a hybrid using SeqSero2 H alleles carries GPL-2.0 content and"
                " cannot ride the MIT wheel (AGENTS.md §6 gate 2); drop --bundle"
            )
        added = 0
        for record_id, sequence in _read_fasta(extra_h_from):
            parsed = _seqsero2_antigen(record_id)
            if parsed is None:
                continue
            family, antigen = parsed
            if family not in {"h1", "h2"}:
                continue
            genes = markers[family].setdefault(antigen, [])
            gene = (
                _gene_id(f"{family}_{antigen}")
                if not genes
                else (f"{_gene_id(f'{family}_{antigen}')}__{len(genes) + 1}")
            )
            genes.append(gene)
            pairs.append((gene, antigen, sequence))
            added += 1
        print(f"  extra H references from {extra_h_from.name}: {added}")
    if extra_refs is not None:
        # Locally extracted alleles (build_h_diversity.py): they widen one
        # antigen's reference diversity when upstream lacks the variant a genome
        # carries. Never bundled — an extracted allele is neither gspec-curated
        # nor a redistributable dataset, so `--bundle` is refused with them.
        if bundle:
            raise SystemExit(
                "extracted alleles (--extra-refs) are local-only: they are neither"
                " gspec-curated nor a redistributable upstream dataset (AGENTS.md §6)"
            )
        document = json.loads(extra_refs.read_text(encoding="utf-8"))
        for reference in document["references"]:
            family = str(reference["family"])
            antigen = str(reference["antigen"])
            genes = markers[family].setdefault(antigen, [])
            gene = (
                _gene_id(f"{family}_{antigen}")
                if not genes
                else f"{_gene_id(f'{family}_{antigen}')}__{len(genes) + 1}"
            )
            genes.append(gene)
            pairs.append((gene, antigen, str(reference["sequence"])))
        print(f"  extra references from {extra_refs.name}: {len(document['references'])}")
    sequences = {gene: sequence for gene, _, sequence in pairs}
    markers, kept = _collapse_references(markers, sequences, floor=collapse_floor)
    markers, pairs = _renumber(markers, kept)
    # After the renumber: the discriminators are standalone records, not antigen
    # alleles, so the collapse and the allele renumbering must not see them.
    discriminator = (
        _add_o_discriminators(o_discriminators, markers, pairs)
        if o_discriminators is not None
        else None
    )
    if discriminator is not None:
        print(
            f"  O discriminator: {len(cast(dict[str, str], discriminator['favours']))}"
            f" group(s), min gap {discriminator['min_gap_pct']} pp"
        )
    rows = _sistr_formulas(table, subspecies=subspecies)
    groups = json.loads(SIMILARITY_GROUPS.read_text(encoding="utf-8"))
    if cross_labels is not None:
        groups = _merge_cross_labels(groups, cross_labels)
    rows, dropped_serovars = _drop_unreachable(markers, rows, groups=groups)
    _report_reachability(markers, rows, groups=groups)
    if dropped_serovars:
        print(
            f"  dropped {len(dropped_serovars)} unreachable serovar row(s), e.g."
            f" {sorted(dropped_serovars)[:5]}"
        )
    _write_fasta(dest, name, pairs)
    _write_database(
        name=name,
        source_label="SISTR antigen database",
        licence="Apache-2.0",
        source_url="https://zenodo.org/records/15078586",
        markers=markers,
        formulas=rows,
        dest=dest,
        similarity_groups=groups,
        o_discriminator=discriminator,
    )
    if kmer_k:
        _write_kmerdb(
            dest,
            name=name,
            markers=markers,
            sequences={gene: sequence for gene, _, sequence in pairs},
            k=kmer_k,
        )
    if bundle:
        _write_bundle_metadata(
            dest,
            name=name,
            source_label="SISTR antigen database",
            licence="Apache-2.0",
            source_url="https://zenodo.org/records/15078586",
            dropped_serovars=dropped_serovars,
            attribution=(
                "SISTR (Salmonella In Silico Typing Resource), PHAC-NML"
                " (https://github.com/phac-nml/sistr_cmd), database"
                " SISTR_V_1.1.3_db (Zenodo record 15078586, Apache-2.0)."
                " Converted to gspec format."
            ),
        )


def import_seqsero2(
    source: Path,
    dest: Path,
    *,
    name: str,
    formulas: Path,
    bundle: bool = False,
    collapse_floor: float = REFERENCE_IDENTITY_FLOOR,
) -> None:
    if not source.is_file():
        raise SystemExit(f"{source}: not a file")
    markers: dict[str, dict[str, list[str]]] = {"o": {}, "h1": {}, "h2": {}}
    pairs: list[tuple[str, str, str]] = []
    for record_id, sequence in _read_fasta(source):
        parsed = _seqsero2_antigen(record_id)
        if parsed is None:
            continue
        family, antigen = parsed
        genes = markers[family].setdefault(antigen, [])
        gene = (
            _gene_id(f"{family}_{antigen}")
            if not genes
            else (f"{_gene_id(f'{family}_{antigen}')}__{len(genes) + 1}")
        )
        genes.append(gene)
        pairs.append((gene, antigen, sequence))
    sequences = {gene: sequence for gene, _, sequence in pairs}
    markers, kept = _collapse_references(markers, sequences, floor=collapse_floor)
    markers, pairs = _renumber(markers, kept)
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
    if bundle:
        raise SystemExit(
            "SeqSero2 is GPL-2.0: its content cannot ride the MIT wheel"
            " (AGENTS.md §6 gate 2). Import into a datadir instead."
        )


def _write_kmerdb(
    dest: Path,
    *,
    name: str,
    markers: dict[str, dict[str, list[str]]],
    sequences: dict[str, str],
    k: int = 27,
) -> None:
    """One differential k-mer set per antigen (`gspec.kmerdb/1`).

    The set is the antigen's representative reference minus every other antigen of
    the same family, which is what makes reads callable by containment: reads are
    counted once and scored per antigen, so near-identical references of one
    antigen no longer split the evidence the way alignment does.
    """
    sets: dict[str, dict[str, object]] = {}
    for family, antigens in markers.items():
        reps = {
            antigen: max(genes, key=lambda gene: len(sequences[gene]))
            for antigen, genes in antigens.items()
        }
        for antigen, gene in reps.items():
            # One representative per other antigen keeps the sidecar small; the
            # subtraction is therefore approximate (documented in the caveats).
            against = [
                encode_seqid(name, member, f"upstream:{other}", ("antigen",))
                for other, member in reps.items()
                if other != antigen
            ]
            # Set names carry the antigen verbatim (``h2_1,2``), so a caller can
            # report the antigen as the table spells it; ids stay sanitized.
            sets[f"{family}_{antigen}"] = {
                "of": encode_seqid(name, gene, f"upstream:{antigen}", ("antigen",)),
                "against": against,
                "n": 0,
            }
    (dest / "kmer.json").write_text(
        json.dumps({"schema": "gspec.kmerdb/1", "k": k, "sets": sets}, indent=1) + "\n",
        encoding="utf-8",
    )
    print(f"  kmer.json: {len(sets)} differential sets (k={k})")


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
    parser.add_argument(
        "--extra-h-from",
        type=Path,
        help=(
            "SISTR only: add another project's fliC/fljB references (e.g. "
            "SeqSero2's H_and_O_and_specific_genes.fasta). The O family stays "
            "SISTR's because the formula table keys on its serogroup letters; a "
            "hybrid carries GPL-2.0 content, so --bundle is refused."
        ),
    )
    parser.add_argument(
        "--subspecies",
        default="enterica",
        help=(
            "SISTR only: restrict the formula table to one subspecies"
            " (antigen genes cannot decide the subspecies, so enterica is the"
            " default; pass 'any' to keep every row)"
        ),
    )
    parser.add_argument(
        "--kmerdb",
        type=int,
        default=0,
        metavar="K",
        help=(
            "also write kmer.json with one differential k-mer set per antigen at"
            " this k (<=31); the sets are materialized by `db build --kmerdb`"
        ),
    )
    parser.add_argument(
        "--extra-refs",
        type=Path,
        help=(
            "gspec.extra_refs/1 document from build_h_diversity.py: additional"
            " alleles for one antigen, extracted from local assemblies (local use"
            " only; --bundle is refused)"
        ),
    )
    parser.add_argument(
        "--cross-labels",
        type=Path,
        help=(
            "SISTR only: gspec.cross_labels/1 document from audit_cross_labels.py;"
            " its classes are merged into the document's similarity_groups"
        ),
    )
    parser.add_argument(
        "--o-discriminators",
        type=Path,
        help=(
            "SISTR only: gspec.o_discriminators/1 document (built by"
            " build_o_discriminators.py) whose marker pair separates serogroups the"
            " antigen references cannot, e.g. A (O:2) from D1 (O:9)"
        ),
    )
    parser.add_argument(
        "--collapse-identity",
        type=float,
        default=REFERENCE_IDENTITY_FLOOR,
        metavar="PCT",
        help=(
            "drop same-antigen references at or above this identity (default"
            " 100, i.e. exact duplicates only; a looser floor deletes the allele"
            " a genome carries and lets a sibling antigen win)"
        ),
    )
    parser.add_argument(
        "--bundle",
        action="store_true",
        help=(
            "also write bundled.json/screen.json so the result can ride the wheel"
            " (permissive upstream licences only; SeqSero2 is refused)"
        ),
    )
    args = parser.parse_args(argv)
    args.dest.mkdir(parents=True, exist_ok=True)
    if args.upstream == "sistr":
        import_sistr(
            args.source,
            args.dest,
            name=args.name,
            formulas=args.formulas,
            bundle=args.bundle,
            subspecies=None if args.subspecies == "any" else args.subspecies,
            extra_h_from=args.extra_h_from,
            kmer_k=args.kmerdb,
            collapse_floor=args.collapse_identity,
            o_discriminators=args.o_discriminators,
            cross_labels=args.cross_labels,
            extra_refs=args.extra_refs,
        )
    else:
        if args.formulas is None:
            raise SystemExit("seqsero2 import needs --formulas (SeqSero2 names serovars in code)")
        import_seqsero2(
            args.source,
            args.dest,
            name=args.name,
            formulas=args.formulas,
            bundle=args.bundle,
            collapse_floor=args.collapse_identity,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
