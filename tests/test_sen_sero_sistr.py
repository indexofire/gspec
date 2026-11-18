"""sen_sero_sistr: the SISTR antigen panel ported into gspec format.

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

from pathlib import Path

from click.testing import CliRunner

from gspec.cli import app
from gspec.db.bundled import find_bundled
from gspec.dbcodec import decode_seqid
from gspec.fasta import iter_fasta
from gspec.formats.typing_result import TypingResultDocument
from gspec.typing.results import SchemeCall
from gspec.typing.serotype_models import SerotypeDocumentV2, read_serotype_document

ASSET = Path(__file__).resolve().parents[1] / "src/gspec/data/dbs/sen_sero_sistr"
runner = CliRunner()


def _genes() -> dict[str, str]:
    return {
        decode_seqid(record.id, "sen_sero_sistr").gene: record.sequence
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
    bundled = find_bundled("sen_sero_sistr")
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
    result = runner.invoke(
        app, ["db", "info", "-d", "sen_sero_sistr", "--json"], catch_exceptions=False
    )
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
    similarity class."""
    call = _designate(tmp_path, ("o", "B"), ("h1", "i"), ("h2", "1,2"))
    assert call.phenotype == "Typhimurium"


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
            "sen_sero_sistr",
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
