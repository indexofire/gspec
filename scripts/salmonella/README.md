# sen_sero curation — Salmonella enterica subsp. enterica serotyping panel

`sen_sero` is a **gspec-curated gene panel** for *Salmonella enterica* subsp.
*enterica* molecular serotyping, built from annotated INSDC CDS records.

## Scope (deliberate)

Only **subsp. enterica** is targeted. The antigenic formula of the other five
subspecies (salamae, arizonae, diarizonae, houtenae, indica) can be read the
same way, but *which* subspecies a genome belongs to cannot be decided from
antigen genes — it needs biochemical tests. A panel that claimed to type them
would be wrong about the subspecies it printed, so nothing here does.

## Full port route: importing an upstream antigen database

`import_upstream.py` converts a full upstream database into a gspec database with
the `gspec.serotype/2` shape (ranked antigen families + the O:H1:H2 formula
table). Two layouts are understood:

| upstream | files | licence | may the wheel carry it? |
|---|---|---|---|
| SISTR | `data/antigens/{wzx,wzy,fliC,fljB}.fasta` + the WHO-2007 CSV | **Apache-2.0** | **yes** with attribution — §6 gate 1 now admits a converted permissive dataset, and `--bundle` records it |
| SeqSero2 | `seqsero2_db/H_and_O_and_specific_genes.fasta` (+ its own formula logic) | **GPL-2.0** | **no** — importing it locally is the only permitted use |

```bash
python scripts/salmonella/import_upstream.py sistr \
    --source /path/to/SISTR_V_1.1.3_db --dest ~/.local/share/gspec/db/sal_sistr \
    --name sal_sistr            # reads the WHO-2007 CSV from the same directory

python scripts/salmonella/import_upstream.py seqsero2 \
    --source H_and_O_and_specific_genes.fasta --formulas who2007.csv \
    --dest ~/.local/share/gspec/db/sal_seqsero --name sal_seqsero
```

What the conversion produces, and why the shape is `/2`:

- each antigen family becomes a **ranked list** of allele references — the
  "top hit wins" step SeqSero2 and SISTR both implement, expressed as data;
- the formula table becomes `formulas[]` (`O:H1:H2 -> serovar`, `"-"` for a
  phase the formula does not carry), which is the only viable representation for
  thousands of rows (one rule per serovar would not scale);
- `info.json` records the upstream licence, URL, conversion date and the
  objection clause, so `gspec db info` shows what a user is actually running.

The importer is a *conversion*, not curation. `--bundle` additionally writes
`bundled.json` + `screen.json` so the result can ride the wheel, and records the
upstream attribution verbatim (AGENTS.md §6 gate 1 requires `license`, `note`
and `attribution` for a non-curated bundle). `seqsero2 --bundle` is refused: its
content is GPL-2.0 and must stay local.

Building the full bundled panel is one command once the upstream database is on
disk:

```bash
python scripts/salmonella/import_upstream.py sistr \
    --source /path/to/SISTR_V_1.1.3_db \
    --dest src/gspec/data/dbs/sen_sero_sistr --name sen_sero_sistr --bundle
python -m pytest tests/test_bundled.py tests/test_serotype_formula.py
```

(The SISTR database is distributed as the Zenodo tarball
`SISTR_V_1.1.3_db.tar.gz`; zenodo.org answered **403** from the environment this
importer was developed in, so the committed bundle is generated where the file is
reachable — the conversion and the bundle metadata are exercised by
`tests/test_serotype_formula.py` on an inline fixture.)

## How the panel is built

```
sen_sero_markers.tsv     curated markers: id, role, genome, selector, note
sen_sero_serovars.tsv    serovar rules: serovar, formula, requires/requires_any/excludes
build_sen_sero.py        fetch (cached) -> sequences + markers.json + serotype.json
                         + screen.json + bundled.json + info.json  in --dest
identity_matrix.py       pairwise identity of the panel, the evidence behind the
                         thresholds (run it whenever a marker is added)
```

```bash
python scripts/salmonella/build_sen_sero.py --dest src/gspec/data/dbs/sen_sero
python scripts/salmonella/identity_matrix.py     # prints the calibration matrix
```

The builder caches each reference genome's `fasta_cds_na` export under
`scripts/salmonella/.cache/` (gitignored) so re-runs are offline and
deterministic; delete the cache to re-fetch.

Every marker is selected by an **annotated identifier** (`gene=` symbol or
`locus_tag=`), never by hand-copied coordinates, and the built record carries
`<genome>:<protein_id>` as its accession so any call can be traced back.

## Threshold calibration (measured, not copied)

`identity_matrix.py` over the current panel reports, among others:

| pair | identity | consequence |
|---|---|---|
| `fliC_i` vs `fliC_gm` / `fliC_a` / `fliC_d` | 41–76 % | H1 phase-1 alleles separate cleanly; the rule floor (98 %) is far above |
| `fljB_12` vs the Paratyphi A H2-region flagellin | 97.0 % | H2 presence needs a floor above 97 %, otherwise a monophasic genome "carries" H2 |
| `rfc` group B vs D1 vs A | 98.5–99.9 % | **`rfc` cannot resolve O groups** — it is not an O marker here |
| `rfbS` group A vs D1 | 98.5–99.3 % | why the D1 marker here is `rfbE` (a different gene) rather than `rfbS`; it is also why O:2 vs O:9 is not a presence question |
| `tyr_A` (group A) vs `rfbE_D1` (group D1) | 98.5–99.2 % | **near-identical by design**: separated by the rule-level `outranks` ranking, never by a presence floor — the mechanism SeqSero2 and SISTR both use |
| `rfbJ` (group B, CDP-abequose synthase) vs everything else | 13–36 % | group B (O:4) *is* separable at the gene level |
| `fliC_i` (H:i) vs `fliC_r` (H:r) | 87.8–90.4 % | closest H1 pair in the panel; still far below the 98 % floor |
| `fliC_b` (H:b) / `fliC_eh` (H:e,h) vs the rest | ≤77.6 % / ≤76.9 % | the two group-B additions separate on their own H1 allele |
| `fljB` H:1,2 of Typhimurium vs Heidelberg | 99.9 % | one shared H2 marker covers both group-B serovars |
| `fliC_gm` (H:g,m) vs `fliC_gp` (H:g,p) | **99.80 %** (3 nt in 1518) | **inseparable by template identity at any sane floor** — the reason Dublin/Gallinarum/Pullorum are out of scope |
| `invA` vs everything | ≤25 % | control gene is unambiguous |

The panel therefore ships **seven serovars**: Typhimurium, Heidelberg,
Paratyphi B and Saintpaul (group B, `rfbJ` + one H1 allele each + the shared
`fljB_12`), Enteritidis and Typhi (group D1, `rfbE`), and Paratyphi A (group A,
`tyr_A`, ranked against the D1 copy). The three generic "H not determined" rows
answer the partial cases (group B alone, group A alone, group D1 alone), and a
monophasic variant of a diphasic serovar lands there instead of being misnamed
(Paratyphi B var. Java is the pinned test).

`tests/test_sen_sero.py::test_panel_markers_are_mutually_separable` turns the
matrix into a gate: every marker pair must stay below 95 % identity, so a future
curation pass cannot silently add a marker that makes calls ambiguous.

Within the shipped panel the closest non-self pair is `fliC_i` vs `fliC_a` at
76.4 %, so the 98 % floor has ample margin; it is set that high because a
genome outside the panel (Paratyphi A's H2-region flagellin) sits at 97.0 %.
Rule floors are therefore `identity_floor: 98`, `coverage_floor: 80`; the DB
defaults (`screen.json`) are `minid 90`, `mincov 60`. These come from the matrix
and from blastn behaviour on the synthetic fixtures in
`tests/test_sen_sero.py` — **re-calibrate against real genomes before trusting a
call** (see caveats below).

## How the reference tools separate O:2 from O:9 (and what that means here)

Neither SeqSero2 nor SISTR decides O:2 vs O:9 by gene presence — both use
**allele-level best-hit ranking** against near-identical references:

- **SeqSero2** (GPL-2.0, read for method only) keeps two references,
  `tyr-O-2__608` and `tyr-O-9__609`, in
  `seqsero2_db/H_and_O_and_specific_genes.fasta`. After the O:9,46 markers
  establish the group family, `decide_O_type_and_get_special_genes()` compares
  the *scores* of the best hit to each tyr allele and takes the larger
  (`if O2 > O9: O_choice = "O-2"`), with an explicit tie-break note favouring
  O:9. It is a strict `>` comparison, so the effective margin is zero.
- **SISTR** (Apache-2.0, same) BLASTs against `data/antigens/{wzx,wzy,fliC,fljB}.fasta`,
  whose record IDs end in the antigen name (`...|O:2`, `...|O:9`), takes the
  **top hit** per antigen (`get_antigen_name(qseqid)`), and looks the triplet up
  in the WHO 2007 serogroup/antigen table. Again: ranking, then a formula table.

`sen_sero` applies that step inside the naming profile: `gspec.serotype/1` rules
gained an `outranks` gate — a ranked allele must be present **and** score
strictly higher (identity × coverage, the `best_hit` metric) than each named
competitor — and the panel declares it in **both directions**:

```json
{"phenotype": "Salmonella Typhi", "requires": ["rfbE_D1", "fliC_d"],
 "excludes": ["fljB_12"], "outranks": {"rfbE_D1": ["tyr_A"]}}
{"phenotype": "Salmonella Paratyphi A", "requires": ["tyr_A", "fliC_a"],
 "outranks": {"tyr_A": ["rfbE_D1"]}}
```

Both generic group rows carry the same gate, so a tie (a synthetic genome
carrying both alleles at 100 %) reports the unresolved string instead of
inheriting SeqSero2's O:9 tie-break. `tests/test_serotype_engine.py` pins the
five gate outcomes (higher wins, lower blocked, tie unresolved, absent
competitor ignores, competitor below floor ignores);
`tests/test_sen_sero.py::test_the_o2_o9_pair_is_decided_by_ranking` pins the
wiring, and `test_both_tyr_alleles_tied_leave_the_group_undecided` pins the tie.

What this panel therefore does **not** need is a different alignment of the
`rfb` cluster: measured over the whole locus, group A and D1 are ~99 % identical
and differ by copy number/arrangement (A carries five duplicated
glycosyltransferase ORFs between `rfbE` and `rfbN`; D1 carries one), so cluster
identity, k-mer sets and gene presence all fail for the same reason. Ranking two
allele references is the correct — and cheapest — mechanism.

## Known limitations of this first slice

1. **O groups beyond B/D1/A are not in the panel.** Group C1/C2, E1 and the
   rest need their own allele references curated the same way (wzx/wzy allele
   ranking, as SISTR does). The `rfbS_A` marker stays in `sequences` as evidence
   of the group-A backbone but gates no rule any more.
2. **Paratyphi A is provisional.** Its rule requires `rfbS` + `fliC_a` and
   excludes the D1 `rfbE`; that exclusion fires only when the D1 allele is
   absent or divergent, which is not guaranteed — the test suite asserts the
   unresolved outcome rather than pretending it works.
3. **Vi antigen is not screened** (`tviB`/`tviC`), so Typhi and Paratyphi C
   report their O:H formula, not "Vi positive".
4. **Phase variation**: only `fljB` presence/absence is used (H2 on/off). The
   `hin` invertase state is not read, so a phase-OFF culture of a diphasic
   serovar is reported by its H1+H2 genes, not by the phase it expressed.
5. **Alleles that differ by a few residues are not separable.** Enteritidis
   (H:g,m) and Dublin (H:g,p) differ in 3 of 1518 nt of `fliC`; Gallinarum and
   Pullorum are non-motile variants of the same group. Distinguishing them needs
   a different signal (epitope-level rules or SISTR-style subtyping), so they are
   deliberately absent instead of being guessed at.
6. The builder verifies **synteny** for every flagellin marker (phase 1 beside
   `fliD`/`fliB`, phase 2 beside `hin`/`fljA`) — the two ORFs share a product
   string, so the neighbourhood, not the annotation text, identifies the allele.
7. Serovars not in `sen_sero_serovars.tsv` return the panel's `fallback`
   string with the markers that were detected, never a guessed name.

## Provenance

Sequences are extracted from RefSeq/GenBank annotated CDS (`NC_003197.2`,
`NC_011294.1`, `NC_006511.1`, `NC_003198.1`) via NCBI E-utilities. INSDC
sequence data carries no copyright restrictions; the selection, naming and
threshold calibration are gspec's own, which is what the bundled-database
curation gate requires. `bundled.json` carries the accessions, the retrieval
date and the caveats so `gspec db info -d sen_sero` shows them to a user.
