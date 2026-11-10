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


def test_widening_applies_before_the_lookup() -> None:
    """The name is widened to its class *before* the lookup (SISTR's `isin` over
    the class), so both rows of the class are reachable: the call is the candidate
    set, not whichever row happened to match by string."""
    calls = {
        "invA": _call("invA"),
        "o_B": _call("o_B"),
        "h1_z13_z28": _call("h1_z13_z28"),
        "h2_12": _call("h2_12"),
    }
    call = evaluate_serotype(calls, _widened_document())
    assert call.confidence == "ambiguous"
    assert {entry.phenotype for entry in call.ambiguous} == {"Exact", "Widened"}


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


def test_two_serovars_sharing_one_formula_are_ambiguous() -> None:
    """The table can list several serovars for one formula; an exact hit then
    reports them all rather than letting dict order pick one."""
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
    assert call.confidence == "ambiguous"
    assert {entry.phenotype for entry in call.ambiguous} == {"Exact", "Twin"}


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
    dest = root / "sen_sero_sistr"
    _import(
        "sistr",
        "--source",
        str(sistr),
        "--dest",
        str(dest),
        "--name",
        "sen_sero_sistr",
        "--formulas",
        str(table),
        "--bundle",
    )
    monkeypatch.setattr("gspec.db.bundled.BUNDLED_ROOT", root)
    found = {database.name: database for database in bundled_databases()}
    assert "sen_sero_sistr" in found
    metadata = found["sen_sero_sistr"].metadata
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
    assert call.confidence == "low"  # not perfect, but still named
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
    assert call.confidence == "low"  # the observed formula is not in the table
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
    assert call.confidence == "low"  # truncated evidence is never perfect


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
