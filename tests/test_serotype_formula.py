"""gspec.serotype/2: ranked antigen families + the formula table.

This is the shape a full Salmonella panel needs (tens of O groups, hundreds of
H alleles, thousands of O:H1:H2 rows). The tests pin the two steps that make it
work: the per-family ranking (near-identical alleles resolved by score and a
margin) and the formula lookup (monophasic rows use "-", so an unresolved phase
must look up as "-").

The importer tests run the upstream *parsers* on tiny inline fixtures: the real
SISTR/SeqSero2 databases are never committed (Apache-2.0 attribution and
GPL-2.0 separation, see the script docstring).
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from gspec.cluster.core import GeneCall
from gspec.typing.serotype_engine import evaluate_serotype
from gspec.typing.serotype_models import (
    SerotypeDocumentV2,
    SerotypeFormula,
    read_serotype_document,
)

IMPORTER = Path(__file__).resolve().parents[1] / "scripts" / "salmonella" / "import_upstream.py"


def _import(*args: str) -> None:
    """Run the conversion CLI exactly as a curator would."""
    result = subprocess.run(
        [sys.executable, str(IMPORTER), *args], capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr


def _call(gene: str, *, identity: float = 100.0, coverage: float = 100.0) -> GeneCall:
    return GeneCall(
        gene_id=gene,
        start=1,
        end=600,
        strand="+",
        coverage_pct=coverage,
        identity_pct=identity,
        verdict="present",
    )


def _document(**overrides: object) -> SerotypeDocumentV2:
    payload: dict[str, object] = {
        "schema": "gspec.serotype/2",
        "name": "sal_test",
        "control_gene": "invA",
        "identity_floor": 90.0,
        "coverage_floor": 60.0,
        "margin": 0.01,
        "antigens": [
            {
                "key": "o",
                "entries": [
                    {"name": "4", "genes": ["o_4"]},
                    {"name": "9", "genes": ["o_9"]},
                    {"name": "2", "genes": ["o_2"]},
                ],
            },
            {
                "key": "h1",
                "entries": [
                    {"name": "i", "genes": ["h1_i"]},
                    {"name": "g,m", "genes": ["h1_gm"]},
                    {"name": "a", "genes": ["h1_a"]},
                ],
            },
            {"key": "h2", "entries": [{"name": "1,2", "genes": ["h2_12"]}]},
        ],
        "formulas": [
            {"o": "4", "h1": "i", "h2": "1,2", "serovar": "Typhimurium"},
            {"o": "9", "h1": "g,m", "h2": "-", "serovar": "Enteritidis"},
            {"o": "2", "h1": "a", "h2": "-", "serovar": "Paratyphi A"},
        ],
    }
    payload.update(overrides)
    return SerotypeDocumentV2.model_validate(payload)


def test_complete_formula_is_named() -> None:
    calls = {
        "invA": _call("invA"),
        "o_4": _call("o_4"),
        "h1_i": _call("h1_i"),
        "h2_12": _call("h2_12"),
    }
    call = evaluate_serotype(calls, _document())
    assert (call.phenotype, call.confidence) == ("Typhimurium", "high")
    assert "formula 4:i:1,2" in call.notes


def test_monophasic_serovar_looks_up_the_dash_row() -> None:
    """Enteritidis is tabulated with h2 "-", and a genome without the phase
    resolves to "-": the two must be the same token."""
    calls = {"invA": _call("invA"), "o_9": _call("o_9"), "h1_gm": _call("h1_gm")}
    call = evaluate_serotype(calls, _document())
    assert call.phenotype == "Enteritidis"
    assert "formula 9:g,m:-" in call.notes


def test_near_identical_alleles_are_ranked_not_merged() -> None:
    """O:2 and O:9 references sit at ~99 % identity: the higher score wins and
    a narrow margin keeps them apart (the SeqSero2 tyr comparison, in data)."""
    calls = {
        "invA": _call("invA"),
        "o_2": _call("o_2", identity=100.0),
        "o_9": _call("o_9", identity=98.9),
        "h1_a": _call("h1_a"),
    }
    call = evaluate_serotype(calls, _document())
    assert call.phenotype == "Paratyphi A"


def test_margin_tie_makes_the_family_ambiguous() -> None:
    """Two alleles inside the margin leave the family undecided, so the call
    falls back to the partial/unresolved path instead of picking one."""
    calls = {
        "invA": _call("invA"),
        "o_2": _call("o_2", identity=100.0),
        "o_9": _call("o_9", identity=99.9),  # closer than margin 0.01
        "h1_a": _call("h1_a"),
    }
    call = evaluate_serotype(calls, _document(margin=0.01))
    # The O family is undecided, so the lookup widens to every serogroup; only
    # one row pairs H:a, so the table still resolves it (SISTR's widening order).
    assert call.phenotype == "Paratyphi A"


def test_unknown_formula_reports_the_partial_template() -> None:
    calls = {"invA": _call("invA"), "o_4": _call("o_4"), "h1_a": _call("h1_a")}
    call = evaluate_serotype(calls, _document())
    assert call.phenotype == "Salmonella 4:a:- (serovar not in the formula table)"


def test_missing_control_gene_falls_back() -> None:
    calls = {"o_4": _call("o_4"), "h1_i": _call("h1_i"), "h2_12": _call("h2_12")}
    call = evaluate_serotype(calls, _document())
    assert (call.phenotype, call.confidence) == ("Not Salmonella", "low")


def _widened_document() -> SerotypeDocumentV2:
    """Two rows whose H1 tokens only meet through a similarity class."""
    return _document(
        antigens=[
            {
                "key": "o",
                "entries": [{"name": "B", "genes": ["o_B"]}],
            },
            {
                "key": "h1",
                "entries": [
                    {"name": "z13,z28", "genes": ["h1_z13_z28"]},
                    {"name": "l,z13,[z28]", "genes": []}
                    if False
                    else {"name": "l,z13,[z28]", "genes": ["h1_never"]},
                ],
            },
            {"key": "h2", "entries": [{"name": "1,2", "genes": ["h2_12"]}]},
        ],
        formulas=[
            {"o": "B", "h1": "l,z13,[z28]", "h2": "1,2", "serovar": "Widened"},
            {"o": "B", "h1": "z13,z28", "h2": "1,2", "serovar": "Exact"},
        ],
        similarity_groups={"h1": (("z13,z28", "l,z13,[z28]"),)},
    )


def test_the_exact_formula_row_is_the_call() -> None:
    """The table carries the observed formula: that row IS the call, and the
    similarity class cannot import sibling rows the antigens never named.

    Measured on a Muenchen genome (`C2-C3:d:1,2`): widening the ``1,2`` token
    to its class also reached the `1,5` rows, so a four-candidate ambiguous
    call buried a formula the table carries verbatim."""
    calls = {
        "invA": _call("invA"),
        "o_B": _call("o_B"),
        "h1_z13_z28": _call("h1_z13_z28"),
        "h2_12": _call("h2_12"),
    }
    call = evaluate_serotype(calls, _widened_document())
    assert call.phenotype == "Exact"
    assert call.confidence == "high"
    assert not call.ambiguous


def test_widening_resolves_a_bracket_token() -> None:
    """Only the bracket token has a reference allele; its class mate resolves it."""
    document = _widened_document()
    document = document.model_copy(
        update={"formulas": tuple(row for row in document.formulas if row.serovar != "Exact")}
    )
    calls = {
        "invA": _call("invA"),
        "o_B": _call("o_B"),
        "h1_z13_z28": _call("h1_z13_z28"),
        "h2_12": _call("h2_12"),
    }
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "Widened"


def test_widening_to_several_serovars_is_ambiguous() -> None:
    """Widening can match more than one row; gspec reports the candidates instead
    of pipe-joining them (SISTR's behaviour) or picking one."""
    document = _widened_document().model_copy(
        update={
            "formulas": (
                # No row carries the allele's own token, so widening decides; two
                # rows share the widened key, which is genuine ambiguity.
                SerotypeFormula(o="B", h1="l,z13,[z28]", h2="1,2", serovar="Widened"),
                SerotypeFormula(o="B", h1="l,z13,[z28]", h2="1,2", serovar="Other"),
            )
        }
    )
    calls = {
        "invA": _call("invA"),
        "o_B": _call("o_B"),
        "h1_z13_z28": _call("h1_z13_z28"),
        "h2_12": _call("h2_12"),
    }
    call = evaluate_serotype(calls, document)
    assert call.phenotype is None
    assert call.confidence == "ambiguous"
    assert {entry.phenotype for entry in call.ambiguous} == {"Widened", "Other"}


def test_two_serovars_sharing_one_formula_keep_the_primary() -> None:
    """The table can list several serovars for one formula (57 keys do): the
    first-declared name is the call — the WHO table lists the common serovar
    first (`C2-C3:d:1,2` is Muenchen, then Virginia) — and the twin rides in
    RUNNER_UP with a note, instead of refusing to name the formula at all."""
    document = _widened_document().model_copy(
        update={
            "formulas": (
                SerotypeFormula(o="B", h1="z13,z28", h2="1,2", serovar="Exact"),
                SerotypeFormula(o="B", h1="z13,z28", h2="1,2", serovar="Twin"),
            )
        }
    )
    calls = {
        "invA": _call("invA"),
        "o_B": _call("o_B"),
        "h1_z13_z28": _call("h1_z13_z28"),
        "h2_12": _call("h2_12"),
    }
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "Exact"
    assert call.runner_up is not None and call.runner_up.phenotype == "Twin"
    assert any("formula shared with Twin" in note for note in call.notes)


def _g_subfactor_document() -> SerotypeDocumentV2:
    """The Diguel (O:22) / Telelkebir (O:23) pair: one G:d:e,n,z15 formula."""
    return _document(
        antigens=[
            {"key": "o", "entries": [{"name": "G", "genes": ["o_G"]}]},
            {"key": "h1", "entries": [{"name": "d", "genes": ["h1_d"]}]},
            {"key": "h2", "entries": [{"name": "e,n,z15", "genes": ["h2_en_z15"]}]},
        ],
        formulas=[
            {"o": "G", "h1": "d", "h2": "e,n,z15", "serovar": "Diguel"},
            {"o": "G", "h1": "d", "h2": "e,n,z15", "serovar": "Telekebir"},
        ],
        o_subfactor={
            "o22_genes": ["o22_g1", "o22_g2", "o22_g3"],
            "o23_genes": ["o23_g1", "o23_g2", "o23_g3"],
            "pairs": [["Diguel", "Telekebir"]],
            "min_gap": 0.5,
        },
    )


def test_o22_o23_subfactor_breaks_a_shared_formula() -> None:
    """Diguel and Telelkebir share G:d:e,n,z15; only the dedicated O:22/O:23
    genes separate them, and they decide the tie the antigen references cannot."""
    calls = {
        "invA": _call("invA"),
        "o_G": _call("o_G"),
        "h1_d": _call("h1_d"),
        "h2_en_z15": _call("h2_en_z15"),
        "o23_g1": _call("o23_g1"),
        "o23_g2": _call("o23_g2"),
        "o23_g3": _call("o23_g3"),
    }
    call = evaluate_serotype(calls, _g_subfactor_document())
    assert call.phenotype == "Telekebir"
    assert call.runner_up is None
    assert any("O:22/O:23 sub-factor" in note for note in call.notes)


def test_o22_o23_subfactor_stays_ambiguous_when_the_gap_is_small() -> None:
    """One weak gene on each side is noise, not a verdict: the call stays with
    the formula's ambiguity instead of a lucky pick."""
    calls = {
        "invA": _call("invA"),
        "o_G": _call("o_G"),
        "h1_d": _call("h1_d"),
        "h2_en_z15": _call("h2_en_z15"),
        "o22_g1": _call("o22_g1", coverage=60.0),
        "o23_g1": _call("o23_g1", coverage=60.0),
    }
    call = evaluate_serotype(calls, _g_subfactor_document())
    assert call.phenotype == "Diguel"
    assert call.runner_up is not None and call.runner_up.phenotype == "Telekebir"
    assert any("not separated" in note for note in call.notes)


def test_special_check_annotates_the_decided_serovar() -> None:
    """oafA O5- is not an antigen: it annotates a Typhimurium call it can prove."""
    document = _document(
        special_checks=[
            {
                "serovar": "Typhimurium",
                "markers": ["oafA_o5minus"],
                "versus": ["oafA_full"],
                "note": "oafA deletion: Typhimurium O5- variant",
            }
        ]
    )
    calls = {
        "invA": _call("invA"),
        "o_4": _call("o_4"),
        "h1_i": _call("h1_i"),
        "h2_12": _call("h2_12"),
        "oafA_o5minus": _call("oafA_o5minus"),
    }
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "Typhimurium"
    assert "oafA deletion: Typhimurium O5- variant" in call.notes


def test_special_check_when_absent_annotates_a_deletion() -> None:
    """``when: absent`` fires when the marker gene is missing entirely — the
    oafA-deletion O5- signal — and stays silent when the gene is present."""
    document = _document(
        special_checks=[
            {
                "serovar": "Typhimurium",
                "markers": ["oafA"],
                "versus": [],
                "when": "absent",
                "note": "oafA absent: O5- variant",
            }
        ]
    )
    base = {
        "invA": _call("invA"),
        "o_4": _call("o_4"),
        "h1_i": _call("h1_i"),
        "h2_12": _call("h2_12"),
    }
    call = evaluate_serotype(dict(base), document)
    assert call.phenotype == "Typhimurium"
    assert "oafA absent: O5- variant" in call.notes

    call_with_gene = evaluate_serotype({**base, "oafA": _call("oafA")}, document)
    assert "oafA absent: O5- variant" not in call_with_gene.notes


def test_non_enterica_subspecies_renders_numeric_formula() -> None:
    """A non-enterica subspecies is reported as ``NUMERAL O:H1:H2``, not as a
    serovar name from the enterica-only table."""
    document = _document(
        o_numeric={"indica": {"C1": "6,7"}, "enterica": {"X": "47"}},
        antigens=[
            {"key": "o", "entries": [{"name": "C1", "genes": ["o_C1"]}]},
            {"key": "h1", "entries": [{"name": "z41", "genes": ["h1_z41"]}]},
            {"key": "h2", "entries": [{"name": "1,7", "genes": ["h2_1_7"]}]},
        ],
        formulas=[{"o": "C1", "h1": "z41", "h2": "1,7", "serovar": "Koumra"}],
    )
    calls = {
        "invA": _call("invA"),
        "o_C1": _call("o_C1"),
        "h1_z41": _call("h1_z41"),
        "h2_1_7": _call("h2_1_7"),
        "subspecies_indica": _call("subspecies_indica"),
    }
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "VI 6,7:z41:1,7"
    assert any("subspecies indica" in note for note in call.notes)


def test_special_check_versus_compares_two_alleles() -> None:
    """The gntR dt SNP is an allele comparison: whichever probe scores higher
    names the variant, and the other is suppressed."""
    document = _document(
        special_checks=[
            {
                "serovar": "Paratyphi B",
                "markers": ["gntR_dtminus"],
                "versus": ["gntR_dtplus"],
                "note": "dt- allele (Paratyphi B var. Java)",
            },
            {
                "serovar": "Paratyphi B",
                "markers": ["gntR_dtplus"],
                "versus": ["gntR_dtminus"],
                "note": "dt+ allele (Paratyphi B var. L(+) tartrate+)",
            },
        ],
        formulas=[
            {"o": "4", "h1": "b", "h2": "1,2", "serovar": "Paratyphi B"},
        ],
        antigens=[
            {"key": "o", "entries": [{"name": "4", "genes": ["o_4"]}]},
            {"key": "h1", "entries": [{"name": "b", "genes": ["h1_b"]}]},
            {"key": "h2", "entries": [{"name": "1,2", "genes": ["h2_12"]}]},
        ],
    )
    base = {
        "invA": _call("invA"),
        "o_4": _call("o_4"),
        "h1_b": _call("h1_b"),
        "h2_12": _call("h2_12"),
    }
    call = evaluate_serotype({**base, "gntR_dtminus": _call("gntR_dtminus")}, document)
    assert call.phenotype == "Paratyphi B"
    assert any("dt- allele" in n for n in call.notes)
    assert not any("dt+ allele" in n for n in call.notes)

    call2 = evaluate_serotype({**base, "gntR_dtplus": _call("gntR_dtplus")}, document)
    assert any("dt+ allele" in n for n in call2.notes)
    assert not any("dt- allele" in n for n in call2.notes)


def test_a_near_tie_still_names_the_formula_with_a_visible_margin_note() -> None:
    """The margin keeps a sibling antigen alive (`g,m` at 0.9993 against `g,q`
    at 0.9987 — one SNP of difference on a real Enteritidis genome): the exact
    formula row is still the call, the sibling rides along as a note, and the
    winning allele's own identity/coverage is in the per-family note so the
    call is auditable without opening the screen table."""
    document = _document(
        antigens=[
            {"key": "o", "entries": [{"name": "9", "genes": ["o_9"]}]},
            {
                "key": "h1",
                "entries": [
                    {"name": "g,m", "genes": ["h1_gm"]},
                    {"name": "g,q", "genes": ["h1_gq"]},
                ],
            },
            {"key": "h2", "entries": [{"name": "1,2", "genes": ["h2_12"]}]},
        ],
    )
    calls = {
        "invA": _call("invA"),
        "o_9": _call("o_9", identity=99.85),
        "h1_gm": _call("h1_gm", identity=99.93),
        "h1_gq": _call("h1_gq", identity=99.87),
    }
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "Enteritidis"
    assert any(note.startswith("H1 near-tie within 0.01: g,q") for note in call.notes)
    assert any(note.startswith("H1: g,m via h1_gm (99.93% id, 100.00% cov)") for note in call.notes)


def test_loader_dispatches_on_the_schema_tag(tmp_path: Path) -> None:
    path = tmp_path / "serotype.json"
    path.write_text(json.dumps(_document().model_dump(by_alias=True)), encoding="utf-8")
    document = read_serotype_document(path)
    assert isinstance(document, SerotypeDocumentV2)

    path.write_text(json.dumps({"schema": "gspec.serotype/9", "name": "x"}), encoding="utf-8")
    with pytest.raises(Exception, match="unsupported serotype schema tag"):
        read_serotype_document(path)


def test_imported_bundle_is_discoverable_and_attributed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``--bundle`` produces a wheel-eligible bundle: `gspec.bundled` discovers it,
    its metadata parses, and the Apache-2.0 attribution is carried verbatim
    (AGENTS.md §6 gate 1) while SeqSero2 (GPL-2.0) is refused outright."""
    from gspec.db.bundled import bundled_databases

    sistr = tmp_path / "sistr"
    (sistr / "antigens").mkdir(parents=True)
    (sistr / "antigens" / "wzy.fasta").write_text(
        ">ref|gene|O:4\n" + "ACGT" * 200 + "\n", encoding="utf-8"
    )
    (sistr / "antigens" / "fliC.fasta").write_text(
        ">ref|gene|i\n" + "TTTT" * 200 + "\n", encoding="utf-8"
    )
    table = tmp_path / "table.csv"
    table.write_text("Serovar,Serogroup,H1,H2\nTyphimurium,O:4,i,-\n", encoding="utf-8")

    root = tmp_path / "bundles"
    dest = root / "sen_sero"
    _import(
        "sistr",
        "--source",
        str(sistr),
        "--dest",
        str(dest),
        "--name",
        "sen_sero",
        "--formulas",
        str(table),
        "--bundle",
    )
    monkeypatch.setattr("gspec.db.bundled.BUNDLED_ROOT", root)
    found = {database.name: database for database in bundled_databases()}
    assert "sen_sero" in found
    metadata = found["sen_sero"].metadata
    assert metadata.license == "Apache-2.0"
    assert "PHAC-NML" in metadata.attribution
    assert "15078586" in metadata.attribution
    assert (dest / "screen.json").is_file()

    result = subprocess.run(
        [
            sys.executable,
            str(IMPORTER),
            "seqsero2",
            "--source",
            str(sistr / "antigens" / "fliC.fasta"),
            "--dest",
            str(tmp_path / "db_seqsero"),
            "--name",
            "sal_seqsero",
            "--formulas",
            str(table),
            "--bundle",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "GPL-2.0" in result.stderr


def test_importer_parses_both_upstream_layouts(tmp_path: Path) -> None:
    """The parsers handle SISTR's ``...|<antigen>`` ids and SeqSero2's
    ``fliC_<antigen>_...``/``O-<group>_wzy__<len>`` ids, and build a /2 document."""
    sistr = tmp_path / "sistr"
    (sistr / "antigens").mkdir(parents=True)
    (sistr / "antigens" / "wzy.fasta").write_text(
        ">ref|gene|O:4\n" + "ACGT" * 200 + "\n>ref|gene|O:9\n" + "ACGT" * 200 + "\n",
        encoding="utf-8",
    )
    (sistr / "antigens" / "fliC.fasta").write_text(
        ">ref|gene|i\n" + "TTTT" * 200 + "\n", encoding="utf-8"
    )
    (sistr / "antigens" / "fljB.fasta").write_text(
        ">ref|gene|1,2\n" + "GGGG" * 200 + "\n", encoding="utf-8"
    )
    table = sistr / "table.csv"
    table.write_text(
        "Serovar,Serogroup,H1,H2\nTyphimurium,O:4,i,1,2\n",
        encoding="utf-8",
    )
    dest = tmp_path / "db_sistr"
    _import(
        "sistr",
        "--source",
        str(sistr),
        "--dest",
        str(dest),
        "--name",
        "sal_sistr",
        "--formulas",
        str(table),
    )
    document = json.loads((dest / "serotype.json").read_text(encoding="utf-8"))
    assert document["schema"] == "gspec.serotype/2"
    assert [entry["name"] for entry in document["antigens"][0]["entries"]] == ["O:4", "O:9"]
    assert (dest / "sequences").read_text(encoding="utf-8").count(">") == 4

    seqsero = tmp_path / "H_and_O_and_specific_genes.fasta"
    seqsero.write_text(
        ">fliC_i_Typhimuriumstr_AUQT01000018_r,i__1476\n" + "AAAA" * 200 + "\n"
        ">fljB_1,2_Saintpaulstr_AYDS01000021_1__1511\n" + "CCCC" * 200 + "\n"
        ">O-9,46_wzy__1191\n" + "GGGG" * 200 + "\n"
        ">tyr-O-2__608\n" + "TTTT" * 200 + "\n"
        ">O-1,3,19_not_in_3,10__130\n" + "ACGT" * 30 + "\n",
        encoding="utf-8",
    )
    dest2 = tmp_path / "db_seqsero"
    _import(
        "seqsero2",
        "--source",
        str(seqsero),
        "--dest",
        str(dest2),
        "--name",
        "sal_seqsero",
        "--formulas",
        str(table),
    )
    document2 = json.loads((dest2 / "serotype.json").read_text(encoding="utf-8"))
    families = {
        family["key"]: [entry["name"] for entry in family["entries"]]
        for family in document2["antigens"]
    }
    assert families["h1"] == ["i"]
    assert families["h2"] == ["1,2"]
    assert families["o"] == ["2", "9,46"]
    assert document2["formulas"][0]["serovar"] == "Typhimurium"


def test_a_disagreeing_hit_set_is_not_a_perfect_call() -> None:
    """SISTR's rule: a family is perfect only when a 100 %-id/100 %-cov hit
    exists AND every hit names the same antigen. A cross-labelled reference
    (`1,2` vs `1,2,7`) makes the names disagree, so the call is downgraded."""
    document = _document(
        antigens=[
            {"key": "o", "entries": [{"name": "B", "genes": ["o_B"]}]},
            {
                "key": "h1",
                "entries": [
                    {"name": "i", "genes": ["h1_i"]},
                    {"name": "r,[i]", "genes": ["h1_r_i"]},
                ],
            },
            {"key": "h2", "entries": [{"name": "1,2", "genes": ["h2_12"]}]},
        ],
        formulas=[{"o": "B", "h1": "i", "h2": "1,2", "serovar": "Typhimurium"}],
    )
    # Both H1 alleles hit at 100 %/100 %: the names disagree.
    calls = {
        "invA": _call("invA"),
        "o_B": _call("o_B"),
        "h1_i": _call("h1_i"),
        "h1_r_i": _call("h1_r_i"),
        "h2_12": _call("h2_12"),
    }
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "Typhimurium"
    # Named, but two H1 names hold identical evidence: resolved, not certain.
    assert call.confidence == "medium"
    assert any("not a perfect match" in note for note in call.notes)


def test_a_single_perfect_allele_keeps_high_confidence() -> None:
    calls = {
        "invA": _call("invA"),
        "o_B": _call("o_B"),
        "h1_i": _call("h1_i"),
        "h2_12": _call("h2_12"),
    }
    call = evaluate_serotype(
        calls,
        _document(
            antigens=[
                {"key": "o", "entries": [{"name": "B", "genes": ["o_B"]}]},
                {"key": "h1", "entries": [{"name": "i", "genes": ["h1_i"]}]},
                {"key": "h2", "entries": [{"name": "1,2", "genes": ["h2_12"]}]},
            ],
            formulas=[{"o": "B", "h1": "i", "h2": "1,2", "serovar": "Typhimurium"}],
        ),
    )
    assert call.confidence == "high"


def test_a_near_tie_inside_the_margin_widens_every_name() -> None:
    """Given two antigens within the document's margin, When the lookup runs,
    Then both names' rows are reported.

    Committing to the winner discards the right name: on real data Newport's
    `e,h` scored 0.9103 against `e,n,x`'s 0.9152, and the winner's rows have no
    serovar, so the serovar the upstream tool reports was unreachable.
    """
    document = _document(
        antigens=[
            {"key": "o", "entries": [{"name": "C2", "genes": ["o_C2"]}]},
            {
                "key": "h1",
                "entries": [
                    {"name": "e,h", "genes": ["h1_eh"]},
                    {"name": "e,n,x", "genes": ["h1_enx"]},
                ],
            },
            {"key": "h2", "entries": [{"name": "1,2,7", "genes": ["h2_127"]}]},
        ],
        formulas=[
            {"o": "C2", "h1": "e,h", "h2": "1,2,7", "serovar": "Newport"},
            # `e,n,x` with `1,2,7` is not a tabulated serovar.
        ],
    )
    calls = {
        "invA": _call("invA"),
        "o_C2": _call("o_C2"),
        "h1_eh": _call("h1_eh", identity=99.5),
        "h1_enx": _call("h1_enx", identity=100.0),  # the winner, by 0.005
        "h2_127": _call("h2_127"),
    }
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "Newport"
    # The observed formula is not in the table (the name came from the wide
    # pass), which is `medium` evidence: named, but not read off a row.
    assert call.confidence == "medium"
    assert any("widening" in note for note in call.notes)


def test_a_high_identity_partial_hit_counts_as_assembled_breakage() -> None:
    """A locus the assembly failed to close is evidence, not absence.

    Real case: a D1 genome carried `wzx` at 99.85 % identity over 52 % union
    coverage because the rfb cluster was split across two contigs; the coverage
    floor discarded it and the serogroup vanished (SeqSero2 called O:9).
    """
    document = _document()
    calls = {
        "invA": _call("invA"),
        "o_9": _call("o_9", coverage=52.0, identity=99.85),
        "h1_gm": _call("h1_gm"),
    }
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "Enteritidis"
    # Truncated evidence names the serovar but can never reach `high`.
    assert call.confidence == "medium"


def test_a_divergent_partial_hit_is_still_rejected() -> None:
    """The weaker rule needs near-identity: a 52 %-covered 85 % hit is noise.

    The serogroup then stays undecided — the formula records ``-`` and the
    lookup widens, which is the honest answer — instead of being pinned to 9.
    """
    calls = {
        "invA": _call("invA"),
        "o_9": _call("o_9", coverage=52.0, identity=85.0),
        "h1_gm": _call("h1_gm"),
    }
    call = evaluate_serotype(calls, _document())
    assert "formula -:g,m:-" in call.notes


def test_the_discriminator_settles_a_serogroup_the_references_cannot() -> None:
    """Given O:2 and O:9 at the same score, When rfbS separates them, Then the
    call follows the marker (SeqSero2's `tyr` comparison, curated from RefSeq)."""
    document = _document(
        antigens=[
            {
                "key": "o",
                "entries": [
                    {"name": "A", "genes": ["o_A"]},
                    {"name": "D1", "genes": ["o_D1"]},
                ],
            },
            {"key": "h1", "entries": [{"name": "g,m", "genes": ["h1_gm"]}]},
            {"key": "h2", "entries": [{"name": "1,2", "genes": ["h2_12"]}]},
        ],
        formulas=[
            {"o": "A", "h1": "g,m", "h2": "-", "serovar": "Paratyphi A"},
            {"o": "D1", "h1": "g,m", "h2": "-", "serovar": "Enteritidis"},
        ],
        o_discriminator={
            "favours": {"A": "o_disc_A", "D1": "o_disc_D1"},
            "min_gap_pct": 0.5,
            "min_coverage": 60.0,
        },
    )
    calls = {
        "invA": _call("invA"),
        "o_A": _call("o_A", identity=99.4),
        "o_D1": _call("o_D1", identity=99.3),
        "o_disc_A": _call("o_disc_A", identity=99.2),
        "o_disc_D1": _call("o_disc_D1", identity=100.0),
        "h1_gm": _call("h1_gm"),
    }
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "Enteritidis"
    assert any("D1 over A" in note for note in call.notes)


def test_a_snp_level_gap_without_a_marker_stays_ambiguous() -> None:
    """A partial locus shows a few tenths of a point: that is noise, so the
    marker abstains and the ambiguity stays visible."""
    calls = {
        "invA": _call("invA"),
        "o_A": _call("o_A", identity=99.4),
        "o_D1": _call("o_D1", identity=99.3),
        "o_disc_A": _call("o_disc_A", identity=99.9, coverage=30.0),
        "o_disc_D1": _call("o_disc_D1", identity=100.0, coverage=30.0),
        "h1_gm": _call("h1_gm"),
    }
    document = _document(
        antigens=[
            {
                "key": "o",
                "entries": [
                    {"name": "A", "genes": ["o_A"]},
                    {"name": "D1", "genes": ["o_D1"]},
                ],
            },
            {"key": "h1", "entries": [{"name": "g,m", "genes": ["h1_gm"]}]},
            {"key": "h2", "entries": [{"name": "1,2", "genes": ["h2_12"]}]},
        ],
        formulas=[
            {"o": "A", "h1": "g,m", "h2": "-", "serovar": "Paratyphi A"},
            {"o": "D1", "h1": "g,m", "h2": "-", "serovar": "Enteritidis"},
        ],
        o_discriminator={
            "favours": {"A": "o_disc_A", "D1": "o_disc_D1"},
            "min_gap_pct": 0.5,
            "min_coverage": 60.0,
        },
    )
    call = evaluate_serotype(calls, document)
    assert any("not separated" in note for note in call.notes)


def test_a_monophasic_call_does_not_borrow_the_other_phase() -> None:
    """Given h1 `i` and no h2 at all, When the table is searched, Then only the
    `i:-` row matches.

    Substituting "-" for a phase that WAS resolved matched every `B:-:-` serovar:
    on the in-house corpus `I 1,4,[5],12:i:-` (a fifth of the genomes) arrived with
    11 unrelated candidates and never as a single call.
    """
    document = _document(
        antigens=[
            {"key": "o", "entries": [{"name": "B", "genes": ["o_B"]}]},
            {"key": "h1", "entries": [{"name": "i", "genes": ["h1_i"]}]},
            {"key": "h2", "entries": [{"name": "1,2", "genes": ["h2_12"]}]},
        ],
        formulas=[
            {"o": "B", "h1": "i", "h2": "-", "serovar": "I 1,4,[5],12:i:-"},
            {"o": "B", "h1": "i", "h2": "1,2", "serovar": "Typhimurium"},
        ],
    )
    calls = {"invA": _call("invA"), "o_B": _call("o_B"), "h1_i": _call("h1_i")}
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "I 1,4,[5],12:i:-"
    # An absent H2 is a monophasic signal, not weak evidence: the row is read off
    # the table with 100 %-id/full-coverage O and H1, so this is a `high` call.
    assert call.confidence == "high"


def test_a_weak_phase_hit_is_absent_not_named() -> None:
    """Below the naming bar a phase is absent: monophasic genomes still carry
    weak `fljB` cross-reactions (a 0.69 `l,z13` hit invented a second phase)."""
    document = _document(
        antigens=[
            {"key": "o", "entries": [{"name": "B", "genes": ["o_B"]}]},
            {"key": "h1", "entries": [{"name": "i", "genes": ["h1_i"]}]},
            {"key": "h2", "entries": [{"name": "1,2", "genes": ["h2_12"]}]},
        ],
        formulas=[
            {"o": "B", "h1": "i", "h2": "-", "serovar": "I 1,4,[5],12:i:-"},
            {"o": "B", "h1": "i", "h2": "1,2", "serovar": "Typhimurium"},
        ],
        phase_floor=0.75,
    )
    calls = {
        "invA": _call("invA"),
        "o_B": _call("o_B"),
        "h1_i": _call("h1_i"),
        "h2_12": _call("h2_12", identity=92.0, coverage=75.0),  # 0.69
    }
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "I 1,4,[5],12:i:-"


def test_confidence_high_needs_a_table_row_whole_evidence_and_no_conflict() -> None:
    """The ladder's top tier: exact row, full-length, one name per family.

    Calibrated on the 100-genome corpus, where this tier holds 22 genomes and
    every one of them agrees with both reference tools — the tier means what a
    lab needs it to mean, which the old binary `perfect` flag did not (it
    required a 100 %-identity hit *and* no other antigen scoring at all, so
    0/100 genomes ever reached `high`).
    """
    calls = {
        "invA": _call("invA"),
        "o_4": _call("o_4", identity=99.5),
        "h1_i": _call("h1_i", identity=99.2),
        "h2_12": _call("h2_12"),
    }
    call = evaluate_serotype(calls, _document())
    assert (call.phenotype, call.confidence) == ("Typhimurium", "high")


def test_confidence_drops_to_medium_when_the_formula_row_is_shared() -> None:
    """A one-to-many formula row names a serovar by declaration order, so the
    twin keeps the call out of the top tier even with perfect evidence."""
    document = _widened_document().model_copy(
        update={
            "formulas": (
                SerotypeFormula(o="B", h1="z13,z28", h2="1,2", serovar="Exact"),
                SerotypeFormula(o="B", h1="z13,z28", h2="1,2", serovar="Twin"),
            )
        }
    )
    calls = {
        "invA": _call("invA"),
        "o_B": _call("o_B"),
        "h1_z13_z28": _call("h1_z13_z28"),
        "h2_12": _call("h2_12"),
    }
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "Exact"
    assert call.confidence == "medium"


def test_confidence_drops_to_low_when_the_weakest_antigen_is_soft() -> None:
    """Below `MEDIUM_SCORE_FLOOR` the name is still produced, but as `low`: this
    is the tier that held the weak-H2 `P:k:k` genome and the 0.74-scoring one."""
    calls = {
        "invA": _call("invA"),
        "o_4": _call("o_4"),
        "h1_i": _call("h1_i", identity=99.0),
        # 95 % id over 85 % coverage = 0.8075, under MEDIUM_SCORE_FLOOR (0.90)
        # but still above the document's 90 %/60 % presence floors.
        "h2_12": _call("h2_12", identity=95.0, coverage=85.0),
    }
    call = evaluate_serotype(calls, _document())
    assert call.phenotype == "Typhimurium"
    assert call.confidence == "low"


def test_an_identical_antigen_in_both_phases_drops_the_weaker_phase() -> None:
    """H1 and H2 cannot be the same antigen: no WHO-2007 row has H1 == H2.

    Measured case: a genome arrived as `P:k:k` because `h2 k` hit at 87.15 %
    identity while `h1 k` was perfect — a `fljB` cross-reaction on a monophasic
    genome. The weaker phase is dropped, so the formula stays a real one.
    """
    document = _document(
        antigens=[
            {"key": "o", "entries": [{"name": "P", "genes": ["o_P"]}]},
            {"key": "h1", "entries": [{"name": "k", "genes": ["h1_k"]}]},
            {"key": "h2", "entries": [{"name": "k", "genes": ["h2_k"]}]},
        ],
        formulas=[
            {"o": "P", "h1": "k", "h2": "-", "serovar": "Inverness"},
            {"o": "P", "h1": "k", "h2": "k", "serovar": "NotASerovar"},
        ],
    )
    calls = {
        "invA": _call("invA"),
        "o_P": _call("o_P"),
        "h1_k": _call("h1_k"),
        # above the document's 90 % identity floor, so it IS a phase — and
        # it duplicates H1, which is what the rule must reject.
        "h2_k": _call("h2_k", identity=95.0),
    }
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "Inverness"
    assert any("repeats H1" in note for note in call.notes)
    assert "formula P:k:-" in call.notes


def test_the_weaker_phase_is_the_one_dropped() -> None:
    """The rule is asymmetric on purpose: the strong phase must survive."""
    document = _document(
        antigens=[
            {"key": "o", "entries": [{"name": "P", "genes": ["o_P"]}]},
            {"key": "h1", "entries": [{"name": "k", "genes": ["h1_k"]}]},
            {"key": "h2", "entries": [{"name": "k", "genes": ["h2_k"]}]},
        ],
        formulas=[{"o": "P", "h1": "-", "h2": "k", "serovar": "OtherPhaseWins"}],
    )
    calls = {
        "invA": _call("invA"),
        "o_P": _call("o_P"),
        "h1_k": _call("h1_k", identity=95.0),
        "h2_k": _call("h2_k"),
    }
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "OtherPhaseWins"
    assert any("repeats H2" in note for note in call.notes)


def test_an_unambiguous_serogroup_reports_its_numeric_o() -> None:
    """`P` is O:38 and nothing else, so the letter can be translated for a lab
    comparing against SeqSero2's numeric notation."""
    document = _document(
        o_numeric={"enterica": {"P": "38"}},
        antigens=[
            {"key": "o", "entries": [{"name": "P", "genes": ["o_P"]}]},
            {"key": "h1", "entries": [{"name": "k", "genes": ["h1_k"]}]},
            {"key": "h2", "entries": [{"name": "junk", "genes": ["h2_junk"]}]},
        ],
        formulas=[{"o": "P", "h1": "k", "h2": "-", "serovar": "Inverness"}],
    )
    calls = {"invA": _call("invA"), "o_P": _call("o_P"), "h1_k": _call("h1_k")}
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "Inverness"
    assert "O group P = O:38" in call.notes


def test_a_multi_factor_serogroup_gets_no_numeric_claim() -> None:
    """`C2-C3` covers 6,8 / 8,20 / 8 / 6,8,20, so there is no single number to
    print and the note must stay absent."""
    document = _document(
        o_numeric={"enterica": {"P": "38"}},
        formulas=[{"o": "4", "h1": "i", "h2": "1,2", "serovar": "Typhimurium"}],
    )
    calls = {
        "invA": _call("invA"),
        "o_4": _call("o_4"),
        "h1_i": _call("h1_i"),
        "h2_12": _call("h2_12"),
    }
    call = evaluate_serotype(calls, document)
    assert not any(note.startswith("O group") for note in call.notes)


def test_the_partial_formula_prefers_the_numeric_o() -> None:
    """A formula the table does not carry renders with the numeric O, which is
    what makes `Salmonella 38:k:-` comparable to SeqSero2's `I 38:k:-`."""
    document = _document(
        o_numeric={"enterica": {"P": "38"}},
        antigens=[
            {"key": "o", "entries": [{"name": "P", "genes": ["o_P"]}]},
            {"key": "h1", "entries": [{"name": "k", "genes": ["h1_k"]}]},
            {"key": "h2", "entries": [{"name": "junk", "genes": ["h2_junk"]}]},
        ],
    )
    calls = {"invA": _call("invA"), "o_P": _call("o_P"), "h1_k": _call("h1_k")}
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "Salmonella 38:k:- (serovar not in the formula table)"


def test_merge_evidence_keeps_the_strongest_layer_and_flags_conflicts() -> None:
    """The evidence merge decides which of several labels for one accession wins.

    A lab's own agglutination outranks a culture-collection strain, which
    outranks a paper table, which outranks the INSDC submitter attribute. Two
    different labels on one accession are kept in ``conflict`` and the row is
    flagged: a genome whose own sources disagree cannot support a marker call.
    """
    # Loaded by path: other tests change the process working directory, so the
    # repo root is not reliably on sys.path by the time this runs.
    import importlib.util

    script = Path(__file__).resolve().parents[1] / "scripts" / "salmonella"
    script = script / "merge_compare_serotypes.py"
    spec = importlib.util.spec_from_file_location("merge_compare_serotypes", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    # `@dataclass` resolves annotations through sys.modules, so register first.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    Evidence, merge_evidence = module.Evidence, module.merge_evidence

    rows = merge_evidence(
        [
            Evidence("GCF_000000001.1", "Enteritidis", "asserted", "ncbi"),
            Evidence("GCF_000000001.2", "Enteritidis", "published", "doi:10.1/x"),
            Evidence("GCF_000000001.2", "Kentucky", "asserted", "ncbi"),
            Evidence("GCF_000000002.1", "Typhimurium", "culture", "NCTC:1"),
            Evidence("GCF_000000002.1", "Typhimurium", "verified", "lab"),
        ]
    )
    by_base = {str(row["base"]): row for row in rows}
    # Version suffixes collapse onto one base accession.
    assert set(by_base) == {"GCF_000000001", "GCF_000000002"}
    # The published label outranks the asserted one, and the clash is visible.
    assert by_base["GCF_000000001"]["serovar"] == "Enteritidis"
    assert by_base["GCF_000000001"]["conflict"] == "Enteritidis;Kentucky"
    assert by_base["GCF_000000001"]["n_labels"] == 3
    # The verified layer wins without a conflict.
    assert by_base["GCF_000000002"]["serovar"] == "Typhimurium"
    assert by_base["GCF_000000002"]["layer"] == "verified"
    assert by_base["GCF_000000002"]["conflict"] == "Typhimurium"
