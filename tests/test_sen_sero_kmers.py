"""k-mer antigen calling from reads (AGENTS: the reads path, alignment-free).

The panel ships one differential k-mer set per antigen (`kmer.json`, k=27).
Counting reads once and scoring each set answers "which O group / H1 / H2" with a
much wider margin than mapping, because an antigen's near-identical references no
longer split the evidence (the alignment path measured `h1_i` at 80 % breadth
against a 75 % sibling; here the same reads give 1.00 against a 0.11 runner-up).
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

from gspec.dbcodec import decode_seqid
from gspec.fasta import iter_fasta

HERE = Path(__file__).resolve().parents[1]
ASSET = HERE / "src/gspec/data/dbs/sen_sero_sistr"
SCRIPT = HERE / "scripts/salmonella/kmers_call_reads.py"


def _module() -> object:
    spec = importlib.util.spec_from_file_location("kmers_call_reads", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["kmers_call_reads"] = module
    spec.loader.exec_module(module)
    return module


def _genes() -> dict[str, str]:
    return {
        decode_seqid(record.id, "sen_sero_sistr").gene: record.sequence
        for record in iter_fasta(ASSET / "sequences")
    }


def _antigen_gene(family: str, antigen: str) -> str:
    document = json.loads((ASSET / "serotype.json").read_text(encoding="utf-8"))
    for declared in document["antigens"]:
        if declared["key"] == family:
            for entry in declared["entries"]:
                if entry["name"] == antigen:
                    return str(entry["genes"][0])
    raise AssertionError(f"{family}:{antigen} not declared")


def _reads(tmp_path: Path, *antigens: tuple[str, str]) -> Path:
    genes = _genes()
    body = "".join(genes[_antigen_gene(family, antigen)] for family, antigen in antigens)
    reads: list[str] = []
    offset = 0
    while offset + 150 <= len(body):
        window = body[offset : offset + 150]
        reads += [window, str(window.translate(str.maketrans("ACGT", "TGCA"))[::-1])]
        offset += 25
    path = tmp_path / "reads.fq"
    path.write_text(
        "".join(f"@r{n}\n{read}\n+\n{'I' * len(read)}\n" for n, read in enumerate(reads)),
        encoding="utf-8",
    )
    return path


def test_kmer_calling_separates_the_three_families(tmp_path: Path) -> None:
    """Given reads from the panel's own Typhimurium antigens, When they are
    scored by k-mer containment, Then each family's antigen wins by a wide
    margin — the reads are counted once, so nothing is split."""
    module = _module()
    fastq = _reads(tmp_path, ("o", "B"), ("h1", "i"), ("h2", "1,2"))
    calls = module.call(ASSET, [fastq])  # type: ignore[attr-defined]
    assert calls == {"o": "B", "h1": "i", "h2": "1,2"}


def test_kmer_calling_resolves_a_near_identical_allele_pair(tmp_path: Path) -> None:
    """O:2 vs O:9 is the alignment problem child (98.5-99.2 % identical). With
    k-mers the correct group still wins outright on reads from group A."""
    module = _module()
    fastq = _reads(tmp_path, ("o", "A"), ("h1", "a"))
    calls = module.call(ASSET, [fastq])  # type: ignore[attr-defined]
    assert calls["o"] == "A"
    assert calls["h1"] == "a"


def test_kmer_manifest_declares_every_antigen(tmp_path: Path) -> None:
    """Every antigen has a set, and its name carries the antigen verbatim."""
    document = json.loads((ASSET / "kmer.json").read_text(encoding="utf-8"))
    assert document["schema"] == "gspec.kmerdb/1"
    assert document["k"] == 27
    antigens = json.loads((ASSET / "serotype.json").read_text(encoding="utf-8"))["antigens"]
    declared = sum(len(family["entries"]) for family in antigens)
    assert len(document["sets"]) == declared
    assert "h2_1,2" in document["sets"]
    entry = document["sets"]["h2_1,2"]
    assert entry["of"].startswith("gspec|db=sen_sero_sistr|gene=h2_1_2")
    # One representative per *other antigen of the same family*: h2 declares 25.
    h2_antigens = next(len(family["entries"]) for family in antigens if family["key"] == "h2")
    assert len(entry["against"]) == h2_antigens - 1


def test_kmers_are_absent_from_a_wrong_serovar(tmp_path: Path) -> None:
    """Reads from group B leave the group A set (nearly) empty: the fraction is
    what discriminates, not merely the presence of any k-mer."""
    from gspec.cluster.kmerdb import (
        count_fastq_kmers,
        load_kmer_set,
        materialize_kmer_sets,
        read_kmerdb,
        set_hit_fraction,
    )

    database = tmp_path / "db"
    database.mkdir()
    (database / "sequences").write_text(
        (ASSET / "sequences").read_text(encoding="utf-8"), encoding="utf-8"
    )
    (database / "kmer.json").write_text(
        (ASSET / "kmer.json").read_text(encoding="utf-8"), encoding="utf-8"
    )
    document = read_kmerdb(database)
    assert document is not None
    materialize_kmer_sets(database, document, database / "sequences")

    fastq = _reads(tmp_path, ("o", "B"))
    query = count_fastq_kmers(fastq, k=27)
    group_b = set_hit_fraction(query, load_kmer_set(database, "o_B"))
    group_a = set_hit_fraction(query, load_kmer_set(database, "o_A"))
    assert group_b > 0.5
    assert group_a < 0.2


def test_cli_reads_route_through_the_kmer_substrate(tmp_path: Path) -> None:
    """`screen --r1` on a gene database that ships a k-mer substrate answers from
    containment (no minimap2), and `typing` designates from those rows."""
    from click.testing import CliRunner

    from gspec.cli import app
    from gspec.db.bundled import materialize_bundled

    datadir = tmp_path / "dd"
    datadir.mkdir()
    materialize_bundled("sen_sero_sistr", datadir, quiet=True)
    assert (datadir / "sen_sero_sistr" / "kmers").is_dir()

    fastq = _reads(tmp_path, ("o", "B"), ("h1", "i"), ("h2", "1,2"))
    table = tmp_path / "screen.tsv"
    runner = CliRunner()
    screened = runner.invoke(
        app,
        [
            "screen",
            "-d",
            "sen_sero_sistr",
            "--r1",
            str(fastq),
            "--read-type",
            "map-ont",
            "-D",
            str(datadir),
            "-o",
            str(table),
            "--quiet",
        ],
        catch_exceptions=False,
    )
    assert screened.exit_code == 0, screened.stderr
    rows = [line.split("\t") for line in table.read_text().splitlines()[1:]]
    genes = {row[1]: row for row in rows}
    assert {"o_B", "h1_i", "h2_1_2"} <= set(genes)
    assert all(genes[gene][5] == "yes" for gene in ("o_B", "h1_i", "h2_1_2"))

    typed = runner.invoke(
        app,
        ["typing", str(table), "--datadir", str(datadir), "--format", "tsv", "--quiet"],
        catch_exceptions=False,
    )
    assert typed.exit_code == 0, typed.stderr
    assert typed.stdout.splitlines()[1].split("\t")[2] == "Typhimurium"
