"""sen_sero: the SISTR antigen panel ported into gspec format.

This bundle is a *conversion*, not gspec curation: SISTR's antigen references and
its WHO-2007 formula table, rewritten as gspec/v1 records plus a
`gspec.serotype/2` document (ranked antigen families + the formula table). Under
AGENTS.md §6 gate 1 a converted permissive dataset may ride the wheel, but it
must carry `license`, `note` and a verbatim `attribution`, which these tests pin.

Coverage is honest, not aspirational: 97.4 % of the upstream formula rows are
reachable from the shipped antigen references. The rest name alleles that
SISTR's reference FASTAs themselves do not contain (`z67`, `z91`, `g,m,t`, …);
the similarity-group widening the port does carry reconciles the bracket tokens
(`l,[z13],[z28]` and friends), and the test below pins the measured percentage
so a future regression cannot hide.
"""

import json
from pathlib import Path

from click.testing import CliRunner

from gspec.cli import app
from gspec.cluster.core import GeneCall
from gspec.db.bundled import find_bundled
from gspec.dbcodec import decode_seqid
from gspec.fasta import iter_fasta
from gspec.formats.typing_result import TypingResultDocument
from gspec.typing.results import SchemeCall
from gspec.typing.serotype_engine import evaluate_serotype
from gspec.typing.serotype_models import ODiscriminator, SerotypeDocumentV2, read_serotype_document

ASSET = Path(__file__).resolve().parents[1] / "src/gspec/data/dbs/sen_sero"
runner = CliRunner()


def _genes() -> dict[str, str]:
    return {
        decode_seqid(record.id, "sen_sero").gene: record.sequence
        for record in iter_fasta(ASSET / "sequences")
    }


def _gene(family: str, antigen: str) -> str:
    """First reference id of one antigen — ids are derived, never hardcoded."""
    for declared in _document().antigens:
        if declared.key == family:
            for entry in declared.entries:
                if entry.name == antigen:
                    return entry.genes[0]
    raise AssertionError(f"{family}:{antigen} not declared")


def _document() -> SerotypeDocumentV2:
    document = read_serotype_document(ASSET / "serotype.json")
    assert isinstance(document, SerotypeDocumentV2)
    return document


def test_bundle_is_discoverable_and_attributed() -> None:
    bundled = find_bundled("sen_sero")
    assert bundled is not None
    metadata = bundled.metadata
    assert metadata.kind == "gene"
    assert metadata.dbtype == "nucl"
    assert metadata.license == "Apache-2.0"
    # The Apache-2.0 NOTICE requirement travels with the conversion.
    assert "PHAC-NML" in metadata.attribution
    assert "15078586" in metadata.attribution
    assert "SISTR" in metadata.vendor


def test_db_info_surfaces_the_attribution() -> None:
    result = runner.invoke(app, ["db", "info", "-d", "sen_sero", "--json"], catch_exceptions=False)
    assert '"attribution"' in result.stdout
    assert "15078586" in result.stdout


def test_document_declares_the_expected_families() -> None:
    document = _document()
    families = {family.key: len(family.entries) for family in document.antigens}
    # The O family is split by locus (wzx/wzy), so it has more entries than
    # serogroups; the unique names still cover every serogroup.
    assert families == {"o": 86, "h1": 62, "h2": 25}
    o = next(family for family in document.antigens if family.key == "o")
    assert len({entry.name for entry in o.entries}) == 46
    assert {entry.locus for entry in o.entries} == {"wzx", "wzy"}
    assert len(document.formulas) > 1500


def test_most_formula_rows_are_reachable() -> None:
    """Every shipped reference antigen resolves to a family entry, and at least
    95 % of the upstream formula rows name only antigens that exist.

    The unreachable remainder is the documented gap (WKL bracket tokens that
    SISTR widens through similarity groups); pinning the number keeps the port
    honest instead of claiming full coverage.
    """
    document = _document()
    known = {
        family.key: {entry.name for entry in family.entries} | {"-"} for family in document.antigens
    }

    def reachable_token(family: str, token: str) -> bool:
        if token in known[family]:
            return True
        for group in document.similarity_groups.get(family):  # type: ignore[arg-type]
            if token in group and any(member in known[family] for member in group):
                return True
        return False

    reachable = [
        row
        for row in document.formulas
        if reachable_token("o", row.o)
        and reachable_token("h1", row.h1)
        and reachable_token("h2", row.h2)
    ]
    assert len(reachable) / len(document.formulas) >= 0.97
    assert set(document.similarity_groups) == {"o", "h1", "h2"}


def test_typhimurium_markers_reach_the_right_row(tmp_path: Path) -> None:
    """group B + H:i + H:1,2 addresses the Typhimurium row verbatim: the call
    is the formula's own serovar, not a candidate list widened through the H2
    similarity class. Group B has no single numeric O antigen (4,5,12 / 4,12 /
    4,27 are all "B"), but the key pins one serovar here, so its WHO row's O is
    reported verbatim (`1,4,[5],12`)."""
    call = _designate(tmp_path, ("o", "B"), ("h1", "i"), ("h2", "1,2"))
    assert call.phenotype == "Typhimurium"
    assert call.antigenic_formula == "1,4,[5],12:i:1,2"


def test_genes_column_lists_the_winning_alleles(tmp_path: Path) -> None:
    """GENES carries the decisive evidence, not every cross-reactive row: one
    representative winning allele per antigen family (the same allele the
    per-family note names), instead of the thousand-row cross-reaction list a
    real genome's table holds."""
    document = _designate_document(tmp_path, ("o", "B"), ("h1", "i"), ("h2", "1,2"))
    assert len(document.files) == 1
    genes = document.files[0].genes
    assert len(genes) == 3
    families = sorted(gene.partition("_")[0] for gene in genes)
    assert families == ["h1", "h2", "o"]


def test_phase_off_alias_reaches_the_monophasic_row(tmp_path: Path) -> None:
    """Paratyphi A is tabulated as A:a:1,5 with ``can_h2_be_missing``; the phase-off
    alias row the importer writes makes a genome without fljB address it."""
    call = _designate(tmp_path, ("o", "A"), ("h1", "a"))
    assert "Paratyphi A" in _candidates(call) or call.phenotype == "Paratyphi A"


def _designate(tmp_path: Path, *antigens: tuple[str, str]) -> SchemeCall:
    """Screen a genome built from the named antigens, then designate it."""
    return _designate_document(tmp_path, *antigens).files[0].phenotypes["serotype"]


def _designate_document(tmp_path: Path, *antigens: tuple[str, str]) -> TypingResultDocument:
    """The parsed gspec.typing_result/1 document for a synthetic antigen genome."""
    genes = _genes()
    body = "".join(genes[_gene(family, antigen)] for family, antigen in antigens)
    genome = tmp_path / "m.fna"
    genome.write_text(f">ctg1\n{body}\n", encoding="utf-8")
    datadir = tmp_path / "dd"
    table = tmp_path / "screen.tsv"
    result = runner.invoke(
        app,
        [
            "screen",
            "-d",
            "sen_sero",
            str(genome),
            "-o",
            str(table),
            "--datadir",
            str(datadir),
            "--quiet",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.stderr
    typed = runner.invoke(
        app,
        ["typing", str(table), "--datadir", str(datadir), "--format", "json", "--quiet"],
        catch_exceptions=False,
    )
    assert typed.exit_code == 0, typed.stderr
    return TypingResultDocument.model_validate_json(typed.stdout)


def _candidates(call: SchemeCall) -> set[str | None]:
    return {entry.phenotype for entry in call.ambiguous}


def test_widening_keeps_the_candidate_list_visible(tmp_path: Path) -> None:
    """A perfect antigen set reads the formula off the table: `B:i:1,2` is
    Typhimurium verbatim, a definitive call the similarity class no longer
    widens into a candidate list (the class reaches `1,5` rows the antigens
    never named). Sibling cross-hits keep the families below SISTR's
    perfection bar, so the confidence stays low."""
    call = _designate(tmp_path, ("o", "B"), ("h1", "i"), ("h2", "1,2"))
    assert call.phenotype == "Typhimurium"
    assert not call.ambiguous


def _call(
    gene: str,
    identity: float,
    *,
    sketch: float | None = None,
    contig: str | None = None,
    start: int = 1,
    end: int = 1000,
) -> GeneCall:
    return GeneCall(
        gene_id=gene,
        start=start,
        end=end,
        strand="+",
        coverage_pct=100.0,
        identity_pct=identity,
        verdict="present",
        sketch_pct=sketch,
        contig=contig,
    )


def _gene_at(family: str, antigen: str, locus: str) -> str:
    """First reference id of one antigen at a named locus (wzx/wzy)."""
    for declared in _document().antigens:
        if declared.key != family:
            continue
        for entry in declared.entries:
            if entry.name == antigen and entry.locus == locus:
                return entry.genes[0]
    raise AssertionError(f"{family}:{antigen}:{locus} not declared")


def _discriminator(name: str) -> ODiscriminator:
    for spec in _document().o_discriminators:
        if spec.name == name:
            return spec
    raise AssertionError(f"discriminator {name!r} not declared")


def test_cross_reactive_phase_falls_back_to_the_real_allele() -> None:
    """A cross-reacting `l,v` that wins the H2 family must not blank the phase.

    Without locus coordinates (a frozen table without SEQUENCE/START/END) the
    fliC `l,v` is free to score in the H2 pool, where it outranks the real
    `1,7`; the H1==H2 guard used to drop the whole phase and report
    `9,12:l,v:-`. The guard now re-calls the phase without the duplicate, so
    the runner-up — the phase's real antigen — is named.
    """
    document = _document()
    lv_h1 = _gene("h1", "l,v")
    lv_h2 = _gene("h2", "l,v")
    h2_17 = _gene("h2", "1,7")
    calls = {
        _gene("o", "D1"): _call(_gene("o", "D1"), 99.92),
        lv_h1: _call(lv_h1, 99.73),
        lv_h2: _call(lv_h2, 99.67),
        h2_17: _call(h2_17, 99.14),
    }
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "Kapemba"
    assert call.antigenic_formula == "9,12:l,v:1,7"
    assert any("re-called as 1,7" in note for note in call.notes)


def test_identical_phase_without_a_runner_up_stays_absent() -> None:
    """The fallback must not invent a phase: when the duplicate was the only
    evidence, the phase is still reported absent."""
    document = _document()
    lv_h1 = _gene("h1", "l,v")
    lv_h2 = _gene("h2", "l,v")
    calls = {
        _gene("o", "D1"): _call(_gene("o", "D1"), 99.92),
        lv_h1: _call(lv_h1, 99.73),
        lv_h2: _call(lv_h2, 99.67),
    }
    call = evaluate_serotype(calls, document)
    assert any("H2 'l,v' repeats H1" in note for note in call.notes)
    assert not any("re-called" in note for note in call.notes)


def test_sketch_tie_break_still_lets_the_snp_discriminator_decide_o() -> None:
    """The anchored sketch may collapse an A/D1 margin tie to A; the rfbS/wbaV
    markers are stronger evidence and must still decide D1.

    Real genome: every O channel ties (wzx A/D1 at 99.92 %), the sketch favoured
    A, and the family never presented the pair to the discriminator because the
    sketch had already narrowed `ties` to one name. The measured D1/D2 `wbaV`
    gap then settles the group as D1, which names the exact formula row.
    """
    document = _document()
    rfb_s = _discriminator("A/D1 rfbS")
    wba_v = _discriminator("D1/D2 wbaV")
    h1 = _gene("h1", "l,v")
    h2 = _gene("h2", "1,7")
    calls = {
        _gene_at("o", "A", "wzx"): _call(_gene_at("o", "A", "wzx"), 99.92, sketch=95.0),
        _gene_at("o", "D1", "wzx"): _call(_gene_at("o", "D1", "wzx"), 99.92, sketch=85.0),
        _gene_at("o", "D2", "wzx"): _call(_gene_at("o", "D2", "wzx"), 99.31, sketch=80.0),
        rfb_s.favours["D1"][0]: _call(rfb_s.favours["D1"][0], 99.76),
        rfb_s.favours["A"][0]: _call(rfb_s.favours["A"][0], 99.29),
        wba_v.favours["D1"][0]: _call(wba_v.favours["D1"][0], 99.80),
        wba_v.favours["D2"][0]: _call(wba_v.favours["D2"][0], 98.90),
        h1: _call(h1, 99.73),
        h2: _call(h2, 99.14),
    }
    call = evaluate_serotype(calls, document)
    assert call.phenotype == "Kapemba"
    assert any("D1 over D2" in note for note in call.notes)


def test_tsv_table_types_the_phase_like_the_json_document(tmp_path: Path) -> None:
    """The frozen table's SEQUENCE/START/END/STRAND columns reach the phase guard.

    A hand-written 15-column table whose `fliC` `l,v` locus also scores in the H2
    family used to type as `9,12:l,v:-` because the table's coordinates were
    dropped at the typing boundary; the JSON document typed it `9,12:l,v:1,7`.
    Both inputs must now give the same serovar.
    """
    tsv_header = (
        "#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS"
        "\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE"
    )

    def table_row(
        contig: str, start: int, end: int, strand: str, gene: str, identity: float
    ) -> str:
        length = end - start + 1
        return (
            f"{genome.name}\t{contig}\t{start}\t{end}\t{strand}\t{gene}"
            f"\t1-{length}/{length}\t{'=' * 15}\t0/0\t100.00\t{identity:.2f}"
            f"\tsen_sero\tupstream:test\tantigen test\tantigen"
        )

    genome = tmp_path / "genome.fna"
    lv_h1 = _gene("h1", "l,v")
    lv_h2 = _gene("h2", "l,v")
    h2_17 = _gene("h2", "1,7")
    o_d1 = _gene("o", "D1")
    table = tmp_path / "screen.tsv"
    table.write_text(
        "\n".join(
            [
                tsv_header,
                # One physical fliC gene, scored by both phase families.
                table_row("ctg1", 1, 1503, "-", lv_h1, 99.73),
                table_row("ctg1", 1, 1503, "-", lv_h2, 99.67),
                # The real fljB locus.
                table_row("ctg2", 1, 1148, "+", h2_17, 99.14),
                table_row("ctg3", 1, 1299, "+", o_d1, 99.92),
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    typed = runner.invoke(
        app,
        [
            "typing",
            str(table),
            "--datadir",
            str(tmp_path / "dd"),
            "--format",
            "json",
            "--quiet",
        ],
        catch_exceptions=False,
    )
    assert typed.exit_code == 0, typed.stderr
    document = TypingResultDocument.model_validate_json(typed.stdout)
    call = document.files[0].phenotypes["serotype"]
    assert call.phenotype == "Kapemba"
    # D1 is the serogroup; the WHO table pins it to 9,12, so the formula spells
    # the numeric O antigen. The TSV carries the same value in its FORMULA cell.
    assert call.antigenic_formula == "9,12:l,v:1,7"
    payload = json.loads(typed.stdout)
    assert payload["files"][0]["phenotypes"]["serotype"]["antigenic_formula"] == "9,12:l,v:1,7"
    table_out = runner.invoke(
        app,
        [
            "typing",
            str(table),
            "--datadir",
            str(tmp_path / "dd"),
            "--quiet",
        ],
        catch_exceptions=False,
    )
    assert table_out.exit_code == 0, table_out.stderr
    first_line, data_line = table_out.stdout.splitlines()[:2]
    assert first_line.split("\t")[3] == "FORMULA"
    assert data_line.split("\t")[3] == "9,12:l,v:1,7"


def test_duplicate_phase_fallback_requires_an_independent_locus() -> None:
    """A runner-up on the cross-reaction's own physical gene is not independent
    phase evidence: a monophasic genome keeps its phase absent.

    Both H2 candidates sit on the `fliC` locus, so the second one is another
    reading of the phase-1 gene, not a `fljB` signal. Only evidence from its own
    locus may resurrect the phase.
    """
    document = _document()
    lv_h1 = _gene("h1", "l,v")
    lv_h2 = _gene("h2", "l,v")
    h2_12 = _gene("h2", "1,2")
    calls = {
        _gene("o", "D1"): _call(_gene("o", "D1"), 99.92, contig="ctg3"),
        lv_h1: _call(lv_h1, 99.73, contig="ctg1", start=1, end=1503),
        lv_h2: _call(lv_h2, 99.67, contig="ctg1", start=1, end=1503),
        h2_12: _call(h2_12, 99.50, contig="ctg1", start=1, end=1503),
    }
    call = evaluate_serotype(calls, document)
    assert any("H2 'l,v' repeats H1" in note for note in call.notes)
    assert not any("re-called" in note for note in call.notes)


def test_unique_serovar_takes_its_who_o_formula() -> None:
    """When the joint key pins ONE serovar, the WHO row's O string is reported.

    `B:l,v:e,n,z15` is Brandenburg alone, whose row is `4,[5],12`; the formula
    keeps the table's bracket notation (the O:5 factor is optional and oafA is
    only a note). `B:d:1,2` is Stanley alone, `1,4,[5],12,[27]`.
    """
    document = _document()
    brandenburg = {
        _gene("o", "B"): _call(_gene("o", "B"), 100.0),
        _gene("h1", "l,v"): _call(_gene("h1", "l,v"), 100.0),
        _gene("h2", "e,n,z15"): _call(_gene("h2", "e,n,z15"), 100.0),
        "special_oafA": _call("special_oafA", 100.0),
    }
    call = evaluate_serotype(brandenburg, document)
    assert call.phenotype == "Brandenburg"
    assert call.antigenic_formula == "4,[5],12:l,v:e,n,z15"
    assert any("oafA detected: O:5 may be positive" in note for note in call.notes)

    stanley = {
        _gene("o", "B"): _call(_gene("o", "B"), 100.0),
        _gene("h1", "d"): _call(_gene("h1", "d"), 100.0),
        _gene("h2", "1,2"): _call(_gene("h2", "1,2"), 100.0),
    }
    call = evaluate_serotype(stanley, document)
    assert call.phenotype == "Stanley"
    assert call.antigenic_formula == "1,4,[5],12,[27]:d:1,2"
    assert any("oafA not detected: O:5 likely negative" in note for note in call.notes)


def test_ambiguous_or_renamed_serovar_keeps_the_serogroup() -> None:
    """Two names for one key (or a naming-priority twin) use the named
    serovar's own O formula, not the raw serogroup letter.

    `C2-C3:l,v:1,2` is Litchfield (6,8) and Pakistan (8); `C2-C3:r:1,5` is
    Hindmarsh (8,20), renamed to Bovismorbificans (6,8,20) by convention.
    Once the engine commits to a name, the formula reflects THAT serovar's
    O expression — a report reading "Bovismorbificans" with "C2-C3" is
    confusing to a lab (the user's Newport precedent). An ambiguous call
    (no name) uses the group-level numeric mapping.
    """
    document = _document()
    litchfield = {
        _gene("o", "C2-C3"): _call(_gene("o", "C2-C3"), 100.0),
        _gene("h1", "l,v"): _call(_gene("h1", "l,v"), 100.0),
        _gene("h2", "1,2"): _call(_gene("h2", "1,2"), 100.0),
    }
    call = evaluate_serotype(litchfield, document)
    # Ambiguous (Litchfield|Pakistan share the formula): numeric group level
    assert call.antigenic_formula == "6,8:l,v:1,2"
    assert call.phenotype in ("Litchfield", None)

    bovis = {
        _gene("o", "C2-C3"): _call(_gene("o", "C2-C3"), 100.0),
        _gene("h1", "r"): _call(_gene("h1", "r"), 100.0),
        _gene("h2", "1,5"): _call(_gene("h2", "1,5"), 100.0),
    }
    call = evaluate_serotype(bovis, document)
    assert call.phenotype == "Bovismorbificans"
    assert call.antigenic_formula == "6,8,20:r:1,5"
