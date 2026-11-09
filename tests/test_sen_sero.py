"""sen_sero: the bundled Salmonella subsp. enterica O/H antigen panel.

The panel is a gene database, so the two-stage contract applies: `screen`
reports the markers, `typing` designates. These tests pin the committed
asset's shape, the rule↔marker consistency, and the end-to-end call on
synthetic genomes assembled from the bundle's own markers (zero network).

The negative cases matter as much as the positive ones: a monophasic serovar's
rule excludes `fljB_12`, so adding that marker must stop the call; and the
provisional Paratyphi A rule fires only when the D1 `rfbE` allele is absent —
the README explains why that is a limitation of a presence panel, not a feature.
"""

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from gspec.cli import app
from gspec.db.bundled import find_bundled
from gspec.dbcodec import decode_seqid
from gspec.fasta import iter_fasta

ASSET = Path(__file__).resolve().parents[1] / "src/gspec/data/dbs/sen_sero"
TABLE = Path(__file__).resolve().parents[1] / "scripts/salmonella/sen_sero_markers.tsv"
runner = CliRunner()


def _marker_sequences() -> dict[str, str]:
    return {
        decode_seqid(record.id, "sen_sero").gene: record.sequence
        for record in iter_fasta(ASSET / "sequences")
    }


def _genome(tmp_path: Path, *markers: str) -> Path:
    sequences = _marker_sequences()
    path = tmp_path / f"{'_'.join(markers)}.fna"
    path.write_text(
        f">ctg1\n{''.join(sequences[marker] for marker in markers)}\n", encoding="utf-8"
    )
    return path


@pytest.fixture(scope="module")
def datadir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """One materialized datadir for the module (the bundle builds on first use)."""
    root = tmp_path_factory.mktemp("sen_sero_dd")
    result = runner.invoke(
        app,
        ["db", "setup", "--datadir", str(root), "--quiet"],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.stderr
    return root


def _designate(tmp_path: Path, datadir: Path, *markers: str) -> dict[str, str]:
    """Screen a synthetic genome and designate: the typing row as a dict."""
    genome = _genome(tmp_path, *markers)
    table = tmp_path / "screen.tsv"
    screened = runner.invoke(
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
    assert screened.exit_code == 0, screened.stderr
    typed = runner.invoke(
        app,
        ["typing", str(table), "--datadir", str(datadir), "--format", "tsv", "--quiet"],
        catch_exceptions=False,
    )
    assert typed.exit_code == 0, typed.stderr
    lines = typed.stdout.splitlines()
    header = lines[0].split("\t")
    return dict(zip(header, lines[1].split("\t"), strict=True))


def test_panel_markers_are_mutually_separable() -> None:
    """Curation gate: every marker pair must be < 95% identical.

    A marker that is ~99% identical to another cannot separate serovars at a
    98% rule floor; adding one would quietly turn calls ambiguous. The measured
    closest pair is fliC_r vs fliC_i at ~88%, so the margin is real. Pairs that
    only differ by a few residues (fliC_gm vs fliC_gp: 3 nt in 1518) are the
    documented reason serovars such as Dublin stay out of this panel.
    """
    import difflib

    sequences = _marker_sequences()
    names = list(sequences)
    worst = (0.0, "", "")
    for index, left in enumerate(names):
        for right in names[index + 1 :]:
            a, b = sequences[left], sequences[right]
            if min(len(a), len(b)) < 0.6 * max(len(a), len(b)):
                continue  # different gene families; blastn will not confuse them
            identity = 100.0 * difflib.SequenceMatcher(None, a, b, autojunk=False).ratio()
            if identity > worst[0]:
                worst = (identity, left, right)
    assert worst[0] < 95.0, f"{worst[1]} and {worst[2]} are {worst[0]:.1f}% identical"


def test_bundle_is_discoverable_and_complete() -> None:
    bundled = find_bundled("sen_sero")
    assert bundled is not None
    assert bundled.metadata.kind == "gene"
    assert bundled.metadata.dbtype == "nucl"
    assert "subsp. enterica" in bundled.metadata.description
    for sidecar in ("sequences", "serotype.json", "screen.json", "bundled.json"):
        assert (ASSET / sidecar).is_file()


def test_serotype_rules_only_name_panel_markers() -> None:
    document = json.loads((ASSET / "serotype.json").read_text(encoding="utf-8"))
    known = set(_marker_sequences())
    assert document["control_gene"] in known
    assert [rule["phenotype"] for rule in document["rules"]] == [
        "Salmonella Typhimurium",
        "Salmonella Enteritidis",
        "Salmonella Typhi",
        "Salmonella Heidelberg",
        "Salmonella Paratyphi A",
    ]
    for rule in document["rules"]:
        referenced = set(rule["requires"]) | set(rule["requires_any"]) | set(rule["excludes"])
        assert referenced <= known, (rule["phenotype"], sorted(referenced - known))
        assert not set(rule["requires"]) & set(rule["excludes"])


def test_committed_markers_match_the_curated_table() -> None:
    rows = [
        line.split("\t")[0]
        for line in TABLE.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ][1:]  # drop the header row
    assert rows == list(_marker_sequences())


def test_typhimurium_markers_designate_typhimurium(tmp_path: Path, datadir: Path) -> None:
    row = _designate(tmp_path, datadir, "invA", "rfbJ_B", "fliC_i", "fljB_12")
    assert row["PHENOTYPE"] == "Salmonella Typhimurium"
    assert row["CONFIDENCE"] == "high"
    assert "4,5,12:i:1,2" in row["NOTES"]


def test_enteritidis_markers_designate_enteritidis(tmp_path: Path, datadir: Path) -> None:
    row = _designate(tmp_path, datadir, "invA", "rfbE_D1", "fliC_gm")
    assert row["PHENOTYPE"] == "Salmonella Enteritidis"
    assert "9,12:g,m:-" in row["NOTES"]


def test_typhi_markers_designate_typhi(tmp_path: Path, datadir: Path) -> None:
    """Typhi shares group D1 with Enteritidis; only the H1 allele separates them."""
    row = _designate(tmp_path, datadir, "invA", "rfbE_D1", "fliC_d")
    assert row["PHENOTYPE"] == "Salmonella Typhi"
    assert "9,12,[Vi]:d:-" in row["NOTES"]


def test_diphasic_marker_breaks_the_monophasic_call(tmp_path: Path, datadir: Path) -> None:
    """Given Enteritidis markers plus fljB_12, Then Enteritidis is not called:
    its rule excludes the phase-2 allele (a diphasic D1 strain is a different
    serovar, and this panel says so instead of guessing)."""
    row = _designate(tmp_path, datadir, "invA", "rfbE_D1", "fliC_gm", "fljB_12")
    assert row["PHENOTYPE"] != "Salmonella Enteritidis"
    assert "not in panel" in row["PHENOTYPE"] or "no complete" in row["PHENOTYPE"]


def test_heidelberg_markers_designate_heidelberg(tmp_path: Path, datadir: Path) -> None:
    """Group B shares rfbJ with Typhimurium; H:r + H:1,2 is what names Heidelberg."""
    row = _designate(tmp_path, datadir, "invA", "rfbJ_B", "fliC_r", "fljB_12")
    assert row["PHENOTYPE"] == "Salmonella Heidelberg"
    assert "4,[5],12:r:1,2" in row["NOTES"]


def test_heidelberg_without_phase_two_is_not_heidelberg(tmp_path: Path, datadir: Path) -> None:
    row = _designate(tmp_path, datadir, "invA", "rfbJ_B", "fliC_r")
    assert row["PHENOTYPE"] != "Salmonella Heidelberg"


def test_paratyphi_a_is_provisional_and_says_so(tmp_path: Path, datadir: Path) -> None:
    """The O:2 rule fires only without the D1 rfbE allele, and carries the caveat."""
    row = _designate(tmp_path, datadir, "invA", "rfbS_A", "fliC_a")
    assert row["PHENOTYPE"] == "Salmonella Paratyphi A"
    assert "rfb cluster stage" in row["NOTES"]


def test_missing_control_gene_suppresses_every_call(tmp_path: Path, datadir: Path) -> None:
    """Without invA the profile reports unresolved: the control gate runs first."""
    row = _designate(tmp_path, datadir, "rfbJ_B", "fliC_i", "fljB_12")
    assert row["PHENOTYPE"] != "Salmonella Typhimurium"
