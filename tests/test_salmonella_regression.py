"""Regression harness + the hybrid (SISTR O + SeqSero2 H) import.

Two things are pinned here:

* the accuracy harness (`scripts/salmonella/regress_serovar.py`) reports the
  right number on a set of synthetic assemblies whose serovars are known;
* `import_upstream.py sistr --extra-h-from <SeqSero2 fasta>` really does widen the
  H panels — and refuses `--bundle`, because the result carries GPL-2.0 content.

The hybrid test needs the local upstream files, so it skips when they are absent
(CI does not have them); the harness test runs everywhere on data generated from
the committed bundle.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from gspec.dbcodec import decode_seqid
from gspec.fasta import iter_fasta
from gspec.typing.serotype_models import SerotypeDocumentV2, read_serotype_document

ASSET = Path(__file__).resolve().parents[1] / "src/gspec/data/dbs/sen_sero_sistr"
REGRESS = Path(__file__).resolve().parents[1] / "scripts/salmonella/regress_serovar.py"
IMPORTER = Path(__file__).resolve().parents[1] / "scripts/salmonella/import_upstream.py"
SISTR_DATA = Path.home() / ".mamba/envs/sistr/lib/python3.9/site-packages/sistr/data"
SEQSERO2_FASTA = Path.home() / ".mamba/envs/seqsero/seqsero2_db/H_and_O_and_specific_genes.fasta"


def _genes() -> dict[str, str]:
    return {
        decode_seqid(record.id, "sen_sero_sistr").gene: record.sequence
        for record in iter_fasta(ASSET / "sequences")
    }


def _gene(document: SerotypeDocumentV2, family: str, antigen: str) -> str:
    for declared in document.antigens:
        if declared.key == family:
            for entry in declared.entries:
                if entry.name == antigen:
                    return entry.genes[0]
    raise AssertionError(f"{family}:{antigen} not declared")


@pytest.mark.parametrize(
    "name,antigens,expected",
    [
        ("typhimurium", (("o", "B"), ("h1", "i"), ("h2", "1,2")), "Typhimurium"),
        ("heidelberg", (("o", "B"), ("h1", "r"), ("h2", "1,2")), "Heidelberg"),
        ("enteritidis", (("o", "D1"), ("h1", "g,m")), "Enteritidis"),
    ],
)
def test_harness_reports_full_accuracy_on_known_assemblies(
    tmp_path: Path, name: str, antigens: tuple[tuple[str, str], ...], expected: str
) -> None:
    """Given assemblies built from the panel's own antigens for known serovars,
    When the harness scores them, Then it reports 100 % and no confusions."""
    document = read_serotype_document(ASSET / "serotype.json")
    assert isinstance(document, SerotypeDocumentV2)
    genes = _genes()
    assemblies = tmp_path / "assemblies"
    assemblies.mkdir()
    body = "".join(genes[_gene(document, family, antigen)] for family, antigen in antigens)
    (assemblies / f"{name}.fna").write_text(f">ctg1\n{body}\n", encoding="utf-8")
    truth = tmp_path / "expected.txt"
    truth.write_text(f"{name}\t{expected}\n", encoding="utf-8")

    result = subprocess.run(
        [
            sys.executable,
            str(REGRESS),
            "--db",
            "sen_sero_sistr",
            "--assemblies",
            str(assemblies),
            "--expected",
            str(truth),
            "--min-accuracy",
            "1.0",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr + result.stdout
    # Exact or ambiguous-but-listed: the harness reports both, and the truth must
    # be among what was delivered.
    assert "contains-truth: 100.0%" in result.stdout


def test_harness_reports_a_confusion_instead_of_hiding_it(tmp_path: Path) -> None:
    """A wrong truth row is reported as (expected -> called), and the exit code
    follows --min-accuracy."""
    document = read_serotype_document(ASSET / "serotype.json")
    assert isinstance(document, SerotypeDocumentV2)
    genes = _genes()
    assemblies = tmp_path / "assemblies"
    assemblies.mkdir()
    body = "".join(
        genes[_gene(document, family, antigen)]
        for family, antigen in (("o", "B"), ("h1", "i"), ("h2", "1,2"))
    )
    (assemblies / "strain.fna").write_text(f">ctg1\n{body}\n", encoding="utf-8")
    truth = tmp_path / "expected.txt"
    truth.write_text("strain\tDublin\n", encoding="utf-8")

    result = subprocess.run(
        [
            sys.executable,
            str(REGRESS),
            "--db",
            "sen_sero_sistr",
            "--assemblies",
            str(assemblies),
            "--expected",
            str(truth),
            "--min-accuracy",
            "0.5",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 1
    assert "contains-truth: 0.0%" in result.stdout
    assert "Dublin ->" in result.stdout


@pytest.mark.skipif(
    not (SISTR_DATA / "antigens").is_dir() or not SEQSERO2_FASTA.is_file(),
    reason="local SISTR + SeqSero2 reference files absent",
)
def test_hybrid_import_widens_the_h_panels_and_refuses_bundling(tmp_path: Path) -> None:
    """SISTR's O groups + formula table with SeqSero2's H alleles: the H panels
    grow, reachability stays, and `--bundle` is refused (GPL-2.0 content)."""
    table = SISTR_DATA / "Salmonella-serotype_serogroup_antigen_table-WHO_2007.csv"
    dest = tmp_path / "sal_hybrid"
    result = subprocess.run(
        [
            sys.executable,
            str(IMPORTER),
            "sistr",
            "--source",
            str(SISTR_DATA),
            "--formulas",
            str(table),
            "--extra-h-from",
            str(SEQSERO2_FASTA),
            "--dest",
            str(dest),
            "--name",
            "sal_hybrid",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    document = json.loads((dest / "serotype.json").read_text(encoding="utf-8"))
    families = {family["key"]: len(family["entries"]) for family in document["antigens"]}
    base = json.loads((ASSET / "serotype.json").read_text(encoding="utf-8"))
    base_families = {family["key"]: len(family["entries"]) for family in base["antigens"]}
    assert families["h1"] >= base_families["h1"]
    assert families["h2"] >= base_families["h2"]
    assert families["o"] == base_families["o"]

    refused = subprocess.run(
        [
            sys.executable,
            str(IMPORTER),
            "sistr",
            "--source",
            str(SISTR_DATA),
            "--formulas",
            str(table),
            "--extra-h-from",
            str(SEQSERO2_FASTA),
            "--dest",
            str(tmp_path / "bundled_hybrid"),
            "--name",
            "sal_hybrid",
            "--bundle",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert refused.returncode != 0
    assert "GPL-2.0" in refused.stderr
