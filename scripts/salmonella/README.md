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
    --dest src/gspec/data/dbs/sen_sero --name sen_sero --bundle
python -m pytest tests/test_bundled.py tests/test_serotype_formula.py
```

The port was run against a local SISTR install
(`…/site-packages/sistr/data`, SISTR db 1.1.3) and ships as `sen_sero`:

| | count |
|---|---|
| O serogroups (wzx + wzy) | 46 |
| H1 alleles (fliC) | 62 |
| H2 alleles (fljB) | 25 |
| antigen reference sequences | 1,549 |
| subsp. enterica formula rows (incl. 52 phase-off aliases) | 1,577 |
| reference sequences after collapsing at ≥99 % identity | 621 |
| **reachable formula rows** | **1,536 (97.4 %)** |

Two conversion steps keep the panel usable:

- **References are collapsed at ≥99 % identity** (1,549 → 621). Contig screening
  is unaffected — any of them would hit — but the *reads* path aggregates per
  gene, so ~10 near-identical `fliC` references per antigen split the reads and
  each reference falls under the breadth floor (measured: two `o_B` references at
  100 %/98 % breadth, two `h1_i` references at 80 %/75 %). Gene ids are
  renumbered after collapsing (`<family>_<antigen>`, then `__2`, `__3`) so a
  rebuild does not shift them.
- **SISTR's similarity groups are ported** (`sistr_similarity_groups.json`,
  Apache-2.0, acknowledged in `bundled.json`). The `gspec.serotype/2` lookup is
  layered: the observed formula is matched **verbatim first** (several serovars
  sharing the exact row keep the first-declared name as the call, the twins in
  RUNNER_UP), and a token only widens to its similarity class when the exact
  formula has no row — several widened matches stay an ambiguous call with the
  candidates, not pipe-joined. That took reachability from 95.4 % to 97.4 %,
  and the layering keeps a widened sibling from burying a formula the table
  carries verbatim (a hand-confirmed Muenchen arrived as four candidates
  before it: the `1,2` class also reaches `1,5` rows).

The remaining 41 rows name alleles that SISTR's reference FASTAs themselves do
not contain (`z67`, `z91`, `g,m,t`, …). `gspec db info -d sen_sero` shows
the Apache-2.0 attribution; `tests/test_sen_sero.py` pins the numbers and
exercises Typhimurium, the phase-off alias (Paratyphi A) and the partial
formula end to end.

(The SISTR database is distributed as the Zenodo tarball
`SISTR_V_1.1.3_db.tar.gz`; zenodo.org answered **403** from the environment this
importer was developed in, so the committed bundle is generated where the file is
reachable — the conversion and the bundle metadata are exercised by
`tests/test_serotype_formula.py` on an inline fixture.)

## Hybrid and accuracy routes

**Hybrid (`--extra-h-from`).** SISTR's O serogroups and WHO-2007 formula table
key on serogroup *letters*; SeqSero2 names the O antigens by *number* (`O-9,46`,
`O-4`), so its O family cannot address the table (measured overlap: 0/35
serogroups). Its H alleles, however, use the same WKL notation (H1 overlap
44/55, H2 13/17), and 1,467 of 1,525 enterica rows have both phases covered by
them. `import_upstream.py sistr --extra-h-from <SeqSero2 fasta>` therefore keeps
SISTR's O family and adds SeqSero2's `fliC`/`fljB` references: measured 46 O + 66
H1 + 25 H2 antigens, 1,851 upstream references collapsed to 637, reachability
97.5 %. `--bundle` is refused for a hybrid (GPL-2.0 content).

**Accuracy (`regress_serovar.py`).** Point it at a directory of assemblies and
SISTR's `data/genomes-to-serovar.txt`:

```bash
python scripts/salmonella/regress_serovar.py --db sal_hybrid \
    --assemblies /path/to/fastas \
    --expected ~/.mamba/envs/sistr/lib/python3.9/site-packages/sistr/data/genomes-to-serovar.txt \
    --datadir ~/.local/share/gspec/db --min-accuracy 0.9
```

It prints accuracy, a per-serovar confusion table and the no-call list, and exits
non-zero below `--min-accuracy`. The runner resolves the database's own
`screen.json` (minid/mincov/dedup/coverage) — passing the CLI defaults explicitly
overrides the panel's evidence model and silently changes the result.

Measured on eight RefSeq complete genomes (NC_003197.2 Typhimurium, NC_011294.1
Enteritidis, NC_003198.1 Typhi, NC_010102.1 Paratyphi B, NC_011080.1 Newport,
NC_011083.1 Heidelberg, NC_011205.1 Dublin, NC_011149.1 Agona): **2/8 = 25 %**.
The failures are informative, not random:

- group-B genomes call the *monophasic* formula (`I 1,4,[5],12:i:-`,
  `I 1,4,[5],12:r:-`) — the H2 antigen is not detected on real data even though
  `fljB` is present, so the `"-"` row wins;
- `Enteritidis -> Moscow` and `Dublin -> Naestved` are H1 near-allele confusions
  (`g,m` vs `g,q`);
- Paratyphi B and Newport return no call.

Diagnosis, now precise:

1. **Upstream cross-labels H2 alleles.** SISTR's `fljB.fasta` holds the same
   sequence under `1,2` and `1,2,7` (max identity between the two antigens in the
   shipped panel: 99.54 %), and real LT2 fljB matches the `1,2,7` reference at
   100 % — so the H2 family resolves to the wrong antenna name and the monophasic
   row wins.
2. **The lookup commits to one name per family.** SISTR's robustness comes from
   doing the opposite twice: it accepts *several* predicted antigens per family
   (`isin`), and it replaces an **unknown** family with the table's whole column
   (`sg = list(df['Serogroup'].unique())`). On these assemblies the O group is
   frequently absent — SISTR's own answers are `-:i:1,2`, `-:r:l,w` — which gspec
   currently turns into a partial call.
3. `screen.json` had to be `gspec.screen/2`: under `/1` the `dedup`/`coverage`
   keys are not part of the document and were silently ignored, so the panel still
   deduped per query and dropped sibling references.

**Landed:** the O group now widens when it is undecided — an unknown serogroup
looks up against every serogroup the panel declares, which is SISTR's
`sg = list(df['Serogroup'].unique())`. A genome without an O call that carries
H:i + H:1,2 therefore designates Typhimurium instead of a partial formula. The H
phases are deliberately **not** widened: their absence is a real monophasic
signal and the `"-"` alias rows already encode it (Paratyphi B var. Java stays a
group-B partial).

### How the two tools handle cross-antigen >= 99 % references

Researched (source read, no data copied):

- **SISTR** — `blast_wrapper.BlastReader.top_result()`: it collects every antigen
  name among the hits for that gene (`{id.split("|")[-1] for id in df['qseqid']}`)
  and only calls a **perfect match** when a 100 %-identity/100 %-coverage hit
  exists *and* that name set has exactly one member. Otherwise it takes the
  highest-bitscore row and flags the call truncated. The near-identical siblings
  are then reconciled in `get_serovar()`, which replaces the predicted name with
  its curated `*_SIMILARITY_GROUPS` class **before** the table lookup and matches
  with `isin(...)` — i.e. *many names means "not perfect"*, and *the name — not a
  candidate set — is widened*. Several serovars can match, and SISTR pipe-joins
  them.
- **SeqSero2** — one best hit per gene with explicit coverage/identity thresholds,
  and the O:2/O:9 pair is settled by comparing the two `tyr` hit scores (strict
  `>`, O:9 on a tie), as documented above.

The faithful implementation is three-fold, and all three are landed:

1. **Per-family agreement.** `_FamilyCall.perfect`/`.ambiguous` record SISTR's
   `top_result()` rule: a call is perfect only when a hit scores full marks *and*
   every hit names the same antigen. A serovar found from a non-perfect family set
   is still named, but with `confidence: low` and a `(not a perfect match)` note.
2. **Widening before the lookup.** Every resolved name is expanded to its
   `similarity_groups` class *before* the table match (SISTR's `isin`), so
   `1,2` also reaches the `1,2,7` / `1,5` rows. A missing O group still widens to
   every serogroup; a missing H phase stays `"-"` (the monophasic alias). Several
   matches produce an `ambiguous` call with the candidate list — the structured
   form of SISTR's pipe-joined `Typhimurium|Lagos`. **Later revised (2026-11-16):
   the widening now runs only when the observed formula has no exact row** — see
   the layered-lookup bullet above; an exact row is a definitive call (the
   first-declared name), and the near-ties ride as notes.
3. **The evidence model of the panel, not the parity default.** Two settings had
   been silently discarding the answer:

   - **`coverage: "union"`** (screen.json). The parity path passes
     `-culling_limit 1`, which keeps one HSP per *query interval*. Around a
     single ~1.5 kb antigen locus dozens of references align, so 1552 HSPs became
     1: Enteritidis screening reported only `h1_g,q` at 99.93 % and dropped
     `h1_g,m` — the allele its genome actually carries — entirely. Union coverage
     drops the culling flag and aggregates per subject.
   - **collapse only exact duplicates** (`import_upstream.py`, floor 100 %). The
     old 99 % same-antigen collapse kept one representative per antigen and
     deleted the carried allele; the top hit then fell to a sibling antigen one
     SNP away and won on bitscore noise. 1549 upstream references are now 1435
     instead of 621.

## Three defects a three-way comparison exposed

Running the whole `26SAL` corpus against **both** references at once (not just
SISTR) found three real bugs that the SISTR-only regression had hidden:

1. **`_formula_keys` substituted `"-"` for a phase that WAS resolved.** A
   monophasic genome (`h1 i`, no `h2`) therefore also matched every `B:-:-` row:
   `I 1,4,[5],12:i:-` — a fifth of the corpus — arrived as 11 unrelated candidates
   and never as a single call. The exact key already covers a missing phase,
   because the caller passes `"-"` for it.
2. **No naming bar for a phase.** A monophasic genome still carries weak `fljB`
   cross-reactions; a 0.69 `l,z13` hit invented a second phase. `serotype/2`
   gained `phase_floor` (0.75): below it a phase is *absent*, not weak. Presence
   floors stay lower — this is the bar for NAMING a phase.
3. **The cross-label merge was harmful where it mattered.** `audit_cross_labels.py`
   had merged `i` with `g,t`, which pulled it into an 18-member class; SISTR's own
   H1 groups contain no `i` at all. The audit stays as a report, `--cross-labels`
   is optional and **not** part of the shipped bundle.

The harness was also silently measuring a different engine: `regress_serovar.py`
never passed `allele_sketch` from the panel's screen defaults, so every corpus run
ran without the sketch/phase evidence. That is why the earlier ablation showed "no
difference" — both arms were the same engine.

## Three-way comparison: gspec, SISTR, SeqSero2 (125 `26SAL` genomes)

Absolute numbers (SeqSero2 calls `-m k -t 4`, SISTR `--no-cgmlst -f tab`):

| metric | value |
| --- | --- |
| SISTR ∩ SeqSero2 | 92.0 % |
| SISTR ∩ gspec | **97.6 %** |
| SeqSero2 ∩ gspec | 92.8 % |
| all three agree | 92.0 % |
| gspec single nomination | 26 % (SISTR: 30 %) |
| gspec misses the two-reference consensus | **5 genomes** |

Normalised for naming convention (subspecies prefix, `1,4,[5],12` versus
`4,[5],12`) — without that normalisation SISTR∩SeqSero2 reads 79 %, which is a
formatting artefact, not a disagreement.

The five real misses: one H2 reference-diversity case (`26SAL0084`, `1,7` vs
`1,2`), one "no call at all" (`26SAL0087`, O:8), and three where gspec's candidate
set is disjoint from both references (`26SAL0172`, `0501`, `0511`). Everything else
is either exact agreement or a superset of SISTR's set.

## Accuracy, measured

`regress_serovar.py` runs the same pipeline the CLI runs (it resolves the panel's
own `screen.json`; a hard-coded `run_screen` call screens with the parity defaults
and would measure the collapse above, not the panel). It also reads
`sistr --no-cgmlst -f tab` output directly, because SISTR reports a formula's whole
serovar set in one `|`-joined cell.

```bash
sistr --no-cgmlst -t 4 -f tab -o sistr_all.tab assemblies/*.fna
python scripts/salmonella/regress_serovar.py --db sen_sero \
    --assemblies assemblies --sistr-tab sistr_all.tab --datadir <datadir>
```

Corpora, truth from `sistr_cmd` 1.1.3:

| corpus | exact | contains truth | before this work |
| --- | --- | --- | --- |
| 304 in-house Salmonella assemblies (25/26SAL) | 26 % | **99.7 %** | 25 % (8 genomes) |
| 8 RefSeq complete genomes | 0 | **7/8** | 2/8 |
| 8 `26SAL` assemblies with hand-checked truth | 3/8 | **7/8** | 2/8 |

"contains truth" counts an `ambiguous` call whose candidate list holds the right
serovar, which is the information SISTR's `Typhimurium|Lagos` strings deliver —
and SISTR itself names a *single* serovar for only about a quarter of this corpus
(the rest are pipe-joined sets), so `exact` understates both tools equally.

Two genomes out of 304 disagree, and in both gspec names a serovar where SISTR
falls back to the formula: `C1:e,h:-` -> Braenderup and `F:f,g:-` -> Rissen. Both
look like SISTR failing to name (Rissen *is* `6,7:f,g:-`); they are recorded here
rather than counted as wins.

## Scripts

| script | what it does |
| --- | --- |
| `import_upstream.py` | build a panel from SISTR/SeqSero2 (`--bundle`, `--kmerdb`, `--o-discriminators`, `--cross-labels`, `--extra-refs`) |
| `build_o_discriminators.py` | extract the O9/O2 `rfbS` pair from RefSeq annotations |
| `audit_cross_labels.py` | find references whose k-mer sets overlap across antigen labels |
| `build_h_diversity.py` | extract novel H alleles from local assemblies |
| `regress_serovar.py` | measure accuracy against `sistr --no-cgmlst -f tab` |
| `sketch_tiebreak.py` | measure the anchored-sketch signal itself |

## Current state after the fixes (spot check)

The five-miss list above predates the formula-key, phase-floor, group-typo and
phase-rework commits, so it is stale. Re-measured on the first 20 `26SAL` genomes
with the shipped panel:

| metric | value |
| --- | --- |
| single-name calls | 8/20 (40 %) |
| candidate set contains SISTR's truth | 11/20 |
| no call where SISTR also has none (`-:-:-`) | 1/20 |
| genuinely wrong | **0/20** |

The full 125-genome re-measurement needs a longer window than the interactive
session allows (each genome costs ~18 s with the anchored sketch enabled); run it
with `regress_serovar.py --assemblies ... --sistr-tab ...` on a machine that can
keep the job alive.

## 2026-11-16: exact-formula-first designation (re-measured on the same 20)

The widening doctrine above refused to name a formula the table carries
verbatim, which the in-house batch made acute: a hand-confirmed Muenchen
(`C2-C3:d:1,2`) arrived as four candidates (the `1,2` class also reaches the
`1,5` rows), and ten D1:`g,m` genomes as thirteen candidates each. The lookup is
now layered, in order:

1. **Exact pass.** The observed `O:H1:H2` keys a row the table carries: that
   row is the call (first-declared name the phenotype, twins in RUNNER_UP with
   a `formula shared with ...` note).
2. **Equal-tie pass.** A family with two antigens at *identical* evidence
   (alignment score and sketch tied; class mates excluded — see the cross-label
   below) contributes all its tied names, and the rows those names key are the
   candidates: a tight set, never the margin siblings' rows or the class
   cross-product.
3. **Widening pass.** The pre-2026-11-16 behaviour (similarity classes,
   undecided-O, monophasic fallback), only for formulas no row carries.

Margin near-ties (strictly below the winner) no longer block the call: they
ride as notes (`H1 near-tie within 0.01: g,q`), which is what SISTR and
SeqSero2 effectively do — both commit to the top antigen name instead of
reporting a margin sibling's serovar.

**The `g,m`/`g,q` cross-label, measured and registered.** The ten D1:`g,m`
genomes still paired `{Enteritidis, Moscow}` after the layering: references
`h1_g_m` and `h1_g_q` are 99.868 % identical (1 mismatch, 1 gap over 1519 bp —
the same allele under two labels), so every channel tied. The audit had missed
the pair because k=21 containment lands at 0.972 (two differences destroy ~42
k-mers), under its 0.98 floor; it is hand-registered in `cross_labels.json`
with the measurement, folded into the shipped panel's `similarity_groups`, and
`import_upstream.py --cross-labels` reproduces the merge on the next rebuild.

Re-measured on the same 20 genomes (single-name 2/20 before the layering):

| metric | before | after |
| --- | --- | --- |
| single-name calls | 2/20 | **16/20** |
| ambiguous with >10 candidates | 15/20 | 0/20 (the 3 remaining ambiguous are the documented rfb assembly gaps: `0501`, `0511`, `0520`) |
| GENES ids on a ranked panel | 575–1,271 | 2–3 (one winning allele per family) |

The hand-confirmed Muenchen calls Muenchen (Virginia in RUNNER_UP); the
`26SAL0516` batch of D1:`g,m:-` genomes calls Enteritidis.

## Do the thresholds survive a different batch?

The panel's numbers were set against real data, so this is a fair challenge. The
parameters split three ways, and `threshold_sensitivity.py` answers the question per
batch instead of trusting the calibration.

**Definitional** (a fact, not a fit): `SAME_LOCUS` 0.5, `PHASE_END_BASES` 150,
`_PHASE_ALTERNATIVE` 0.98, `MIN_SKETCH`/`SKETCH_MARGIN`, the k-mer k, and the union/
`_distinct_loci` overlap rule. A new batch cannot invalidate these; only a change in
what the words mean would.

**Reporting**: `margin` 0.01 decides whether a family is reported as one name or a
candidate list. Widening never reorders rows, so it cannot turn a right call into a
wrong one — only a decided call into a listed one. It is the *most* batch-sensitive
number here (77 % of multi-antigen families have a top-two gap inside 3x margin) and
it is also the safest to move.

**Acceptance** (these can change a call): `identity_floor` 80, `coverage_floor` 60,
`partial_coverage_floor` 45, `partial_identity_floor` 95, `phase_floor` 0.75. They
are safe while the observations leave an empty band around them. Measured on the
in-house corpus:

| threshold | evidence | verdict |
| --- | --- | --- |
| `phase_floor` 0.75 | phase scores are ≤ 0.691 (cross-reaction) or ≥ 0.936 (present): a **0.245-wide empty band**, 0/248 observations within ±0.05 | survives a large shift in assembly quality |
| `identity_floor` 80 | 5.4 % of hits within ±2 pp | re-check first on a new batch |
| `coverage_floor` 60 | 7.5 % of hits within ±5 pp | re-check first on a new batch |
| partial rule (95 / 45) | 57 accepted hits land at 51.4–55.5 % coverage, i.e. above the floor | safe, but the margin is data-shape dependent |

So the answer is: no, a new batch does not invalidate them wholesale — but the
acceptance thresholds must be re-measured, which is exactly what the tool prints. Run
it on every new batch and re-calibrate whichever threshold reports observations inside
its band.

## Case-by-case diagnosis of the remaining misses

Five genomes disagreed with both references in the 26SAL comparison. Each was
traced to its raw hits rather than guessed at, and the split matters: **most are
assembly problems, not database content.**

| genome | SISTR / SeqSero2 | what the hits say | class |
| --- | --- | --- | --- |
| `26SAL0084` | `D1:l,v:l,v` / Mendoza | `wzx` D1 99.85 % but only **52 %** of the gene (two contigs); `fljB` matches `1,7` at 99.61 % and `1,2` at 96.83 % | assembly break + reference diversity |
| `26SAL0087` | `C2-C3:r:-` / `8:r:-` | O `C2-C3` 99.58 % full length, H1 `r,[i]`/`r`, no H2 — the table simply has **no serovar** for `C2-C3:r:-` | **not a miss** (SeqSero2 agrees) |
| `26SAL0172` | `C2-C3:l,v:l,v` / Litchfield | O `C2-C3` 100 %, H1 `l,v` 100 %, H2 `1,2` 99.94 % on its own contig — the phase guard now keeps them apart | **fixed by the phase rework** |
| `26SAL0501` | Indiana / `-:z:1,7` | O group B reference hits at 99.65 % but over only **287 of 1293 bp (22 %)**; H1 `z` and H2 `1,7` are perfect | assembly gap in `rfb` |
| `26SAL0511` | `Typhimurium|Lagos` / `-:i:1,2` | O group B reaches only 80.6 % over 227 bp — the *marker* is absent from the assembly; H1 `i` and H2 `1,2` are perfect | assembly gap in `rfb` |

> `26SAL0511` is re-diagnosed in the 2026-11-16 open-question section below:
> per-gene LB probes show nearly the whole `rfb` cluster present at 99.6–100 %
> group B, with only `rfbU`/`rfbX` (the panel's markers) in the gap — the
> cluster is not absent, the marker is.

Three of five are missing or broken `rfb` sequence: no reference can fix that, and
merging fragments cannot either (the bases are not in the assembly). The response is
assembly QC — flag "O family absent while H is perfect" and ask for long reads or a
re-assembly — which is exactly the signal `antigen_diversity.py` surfaces as
"family with no hit" and `--mincov`/`partial_coverage_floor` already encode as a
distinction between *absent* and *truncated*.

## Two decision tools (run these before building anything)

**`sketch_separability.py`** answers "would a genomic-similarity channel separate
serovars?" from a corpus and a truth table, and exits 1 when it would not. Measured
on 125 in-house assemblies (k=21, bottom-1000 MinHash):

| | min | median | max |
| --- | --- | --- | --- |
| within a serovar (135 pairs) | 0.630 | 0.864 | 0.979 |
| between serovars (7615 pairs) | 0.016 | 0.519 | **0.997** |

The most similar different-serovar pair is closer than the most similar
same-serovar pair, and 5.7 % of the between pairs sit above the closest within
pair: no threshold exists. `k` and sketch size are not the limit either — bottom-k
error is `sqrt(J(1-J)/size)` (1.9 % at s=1000, J=0.9) and `J ~ p**k/(2-p**k)` makes
k=21 resolve a 0.5 % ANI gap by 8 Jaccard points. The overlap is biological
(accessory genome, phage, and clonal variants such as Typhimurium versus its
monophasic `1,4,[5],12:i:-` form), so a reference-genome sketch is not a path to
serovar designation. It remains useful for batch QC (duplicates, contamination) and
for flagging "antigen says X but the genome is far from every X reference".

**`antigen_diversity.py`** turns a screen table plus a SISTR truth dump into a
worklist: for every genome and family it looks up the identity of the antigen the
genome *should* carry (not its best hit, so a monophasic genome's weak `fljB`
cross-reaction is correctly ignored) and reports what is weak or absent.

Two things came out of running it on the 26SAL corpus:

* it found an **upstream data typo**: SISTR's H1 similarity group spells its last
  member `'g,s,tg,t'`, so `g,[s],t` never reached the rows keyed `g,s,t`.
  `import_upstream.py` now corrects it (`SIMILARITY_GROUP_FIXES`, printed at build
  time) and Senftenberg is reachable again — one genome recovered, documented
  because it deviates from SISTR's constants on purpose.
* it flagged `h2 e,n,z15` at 86.9 % for three Braenderup/Chester-lineage genomes.
  Twenty extra alleles were extracted from the other 301 assemblies
  (`build_h_diversity.py`) and the panel grew 1435 → 1455 — **and nothing changed**:
  the class widening already reached the right row, because `e,n,x` and `e,n,z15`
  share a SISTR class. A weak identity is a candidate, not a defect; check whether
  widening already covers it before adding data.

## Sketch tiebreak: shipped

SISTR separates near-identical antigen evidence with cgMLST (330 core genes).
gspec does not do allele calling and will not grow an MLST engine, but it already
owns k-mer machinery (`cluster/kmerdb.py`), so a sketch can decide a family whose
alignment scores sit inside the margin. Measured with `sketch_tiebreak.py` over 60
in-house genomes (94 H1/H2 antigen tokens, truth = SISTR's `antigenic_formula`):

| method | agrees with SISTR |
| --- | --- |
| alignment, identity x coverage | 78.7 % |
| whole-genome k-mers vs each antigen's set | 53.2 % |
| **anchored window k-mers, per allele** | **94.7 %** |

The middle row is the trap: comparing a 4.8 Mb genome against a 1.5 kb antigen's
k-mer set saturates, so the winner is decided by how large that antigen's reference
pool happens to be — it is *worse* than alignment (it even prefers `r,[i]` over
`i`). What works is anchoring: take the genomic window the alignment already found,
count its k-mers, and score each reference **allele** by
`|window ∩ allele| / |allele|`; the best allele names the antigen. That ignores
gapped/extended HSPs and asks only "whose k-mers are here".

It ships as an evidence field, not as a second engine:

- the panel's `screen.json` sets `contigs.allele_sketch: true`; every hit gains a
  `sketch_pct` in the `gspec.report/1` document (`--format json`);
- `gspec typing` accepts that JSON and, when a family's top two alignment scores
  are inside `margin`, keeps only the antigens within `SKETCH_MARGIN` (5 points) of
  the best sketch. A tie the sketch cannot separate stays a tie, and a database
  without k-mer sketches behaves exactly as before.

**What it changes in the aggregate.** On the 304-genome corpus the serovar metric
does *not* move (99.3 % before and after): the widening strategy already put the
truth in the candidate list, so a corpus where every tie spans interchangeable
serovars cannot show the gain. The gain is antigen-level fidelity (the table
above) and fewer wrong-antigen calls on genomes whose competing antigens lead to
*disjoint* serovar sets — e.g. `fljB` `1,6` vs `1,7`, which now resolves
(`24SAL0292` -> `London|Winterthur`, matching SISTR's own set).

## Phase cross-talk: `fliC` and `fljB` share `l,v`

SISTR's antigen pools are labelled per gene family, but the same antigen name can
appear in both: `l,v` is a phase-1 antigen in `fliC.fasta` (30 records) and a
phase-2 antigen in `fljB.fasta` (34 records), and the alleles are ~99 % identical.
A genome's single `fliC` gene therefore scores in *both* families, and the
phase-2 score can outrank the genome's real `fljB`.

Measured on `26SAL0084` (Mendoza, `9:l,v:1,2`): one contig hit `h1 l,v` at 97.14 %
and `h2 l,v` at 99.53 %, while the actual `fljB` on another contig matched `1,2`
at 96.83 %. The phase-2 pool's copy won, so gspec reported `l,v` for H2 — and so
does SISTR, which has no guard either.

`_separate_phases` now keeps one physical gene out of both families. Hits are
grouped into loci (same contig, overlapping interval) and a locus claimed by both
phases is handed to the family that owns a *competitive* locus elsewhere — within
2 % of the contested score, because alleles of one antigen sit inside ~1 %. Both
conditions are deliberately narrow, and each was calibrated against a regression:
a 10 % band let a stray 0.86 hit stand in for Braenderup's real `e,h` gene (0.913)
and the call fell to `z39`, and an earlier version that ignored the antigen name
deleted `e,h` even though the locus's second name was `e,n,z15` — cross-reactivity
between flagellin alleles, not a duplicated reference. With the guard,
`26SAL0084` reports `l,v` for H1 and `1,7` for H2 (the residual `1,2` vs `1,7` gap
is the cross-labelled-reference problem, 2.8 points apart), and the 304-genome
corpus stays at 99.3 %.

## Phase affinity: which flagellin is in this locus

The phase guard (see "Phase cross-talk" below) used to infer the phase from score
ordering. The direct evidence is in the ends: `fliC` and `fljB` share ~99 % of their
middle — that is why one gene scores in both pools — but their 5'/3' head and tail
are conserved *within* a phase, which is exactly what SeqSero2's head/tail markers
exploit.

Screen now anchors a `phase_pct` on every H hit: the fraction of the locus' end
k-mers (150 bp each side) found in **any** allele of that family. Measured on a real
phase conflict (one locus, both pools): `fliC` 86.1 % against `fljB` 74.2 %, a
12-point gap that names phase 1 outright. A same-antigen comparison gave only 5
points — the family-wide comparison is the one that works, because the ends are
shared across antigens of a phase.

When a contested locus is reassigned, the losing family keeps the antigen **name**
as a low-confidence candidate instead of erasing it, and only while its surviving
locus is weak (< 0.9). That matters: Braenderup's `e,h` references fit its `fljB`
locus at 91 % and its own `fliC` at 86 %, so deleting the demoted name dropped the
correct serovar from the candidate list entirely.

## Reference diversity: extracting new alleles (`--extra-refs`)

A pool can be complete enough for rows to be reachable and still miss the variant a
genome carries: one Mendoza genome's `fljB` matched `1,7` at 99.61 % and `1,2` at
96.83 %, so the call landed on `1,7` while SeqSero2 (with its own reference set)
said `1,2`. `build_h_diversity.py` extracts candidate alleles from assemblies by
homology (no annotation), keeps only alleles not already represented (k-mer
containment < 0.99 of an existing one), and writes `gspec.extra_refs/1`;
`import_upstream.py --extra-refs` appends them. Measured on the 304-genome corpus:
26 novel `1,2` alleles, panel 1435 -> 1461. `--bundle` is refused — an extracted
allele is neither gspec-curated nor a redistributable dataset.

That specific genome still reports `1,7`: its `fljB` is 99.61 % identical to a
genuine `1,7` reference (verified 100 % against the upstream `1,7` record that
produced it), so no `1,2` allele from this corpus outranks it. With SeqSero2 saying
`1,2` and SISTR saying `l,v`, that genome's H2 is genuinely contested — recorded
here rather than tuned away.

## Curated O9/O2 discriminator (public-domain, own extraction)

`wzx`/`wzy` cannot separate serogroup A (O:2) from D1 (O:9): the two share the rfb
machinery and differ in `rfbS`, the paratose/tyvelose epimerase. Measured on the
references: the alleles are 99.29 % identical over 840 bp, and SeqSero2's `tyr-O-9`
marker turns out to be exactly the 5' 609 bp of `rfbS`, where the two differ by a
single SNP (99.84 %).

`build_o_discriminators.py` extracts that pair from two NCBI RefSeq annotations —
`NC_003198.1` (Typhi CT18, O:9/D1) and `NC_006511.1` (Paratyphi A ATCC 9150,
O:2/A), both public domain — records the provenance, and writes
`scripts/salmonella/o_discriminators.json`. `import_upstream.py --o-discriminators`
adds the two sequences as standalone records (never as antigen entries, so they
cannot widen a candidate list) and writes an `o_discriminator` block into
`serotype.json`.

The engine only consults them when the O tie actually spans their names, and only
when the gap is decisive: `min_gap_pct` is 0.5, because a partial locus shows a few
tenths of a point (measured 0.3 on a 344 bp fragment) and that is noise. On
`26SAL0084` the gap was 0.60 pp and the marker chose D1 over A — the same serogroup
SeqSero2 reported.

## Systematic cross-labelled references

`audit_cross_labels.py` finds upstream references that carry two antigen labels by
comparing each antigen's representative k-mer set (containment ≥ 0.98 of the
shorter one). On the SISTR panel it reports five classes:

```
h1: g,m == g,m,q == g,q     g,p == g,p,s == g,p,u     g,t == i
h2: z35 == z52
```

These are folded into the document's `similarity_groups`
(`import_upstream.py --cross-labels`), so both names' rows stay reachable and the
answer is an honest ambiguity instead of whichever label won on bitscore. The audit
also flagged a record whose antigen name is empty (`h1__f_g_t`), which the merge
skips. The `g,m == g,q` pair is hand-registered (2026-11-16): `h1_g_m` and
`h1_g_q` are 99.868 % identical (1 mismatch, 1 gap over 1519 bp — one allele
under two labels), but the audit's representative comparison and its 0.98 k=21
floor both miss a pair with two differences, and ten D1:`g,m` genomes paired
`{Enteritidis, Moscow}` until the merge landed.

The `1,2` versus `1,7` case is **not** cross-labelling: k-mer containment between
them is 0.80, i.e. genuinely different antigens. A genome whose `fljB` matches
`1,7` at 99.61 % and `1,2` at 96.83 % is therefore a *reference diversity* gap — the
`1,2` pool does not carry that allele — and the fix is more `1,2` references, not a
merger.

## Case study: a serogroup lost to an assembly break (26SAL0084)

`26SAL0084` came out as `-:l,v:l,v` while SeqSero2 reported `9:l,v:1,2`
(Mendoza). Both failures were checked against the references directly, and
**neither needs extra database sequence**:

- **O group absent.** Aligning SISTR's `wzx` D1 references to the genome gives
  100 %/99.6 % identity — but in two blocks of 409 bp and 266 bp, on
  `contig00348` and `contig00291`. The `rfb` cluster is repeat-rich and this
  short-read assembly broke it in two, so the gene's union coverage is 52 % and
  the panel's 60 % screen floor erased it. The reference is fine; the *assembly*
  is fragmented. Fixes landed: the panel screens at `mincov 45` and lets the
  serotype document decide presence, and `gspec.serotype/2` gained a
  near-identity partial rule (`partial_coverage_floor` 45 /
  `partial_identity_floor` 95) that accepts a truncated hit **and ranks it by the
  identity it shows**, because truncation is an artefact rather than evidence of
  absence. A partial hit can never be perfect.
- **H2 near-miss, not a gap.** The genome matches the panel's `h2_1_2` alleles at
  98.6 % over the full gene, so the right allele is already in the pool. It loses
  to `l,v` at 99.67 % and `1,7` at 99.61 % — cross-labelled upstream references
  one SNP apart, 1.07 points above the correct allele, i.e. just outside the
  document's 0.01 `margin`. SISTR gets this genome wrong the same way (it also
  reported `l,v`), which is what a margin calibration can and cannot fix.

### Why SeqSero2 called O9 where gspec and SISTR did not

SeqSero2 does not simply align a full-length `wzx`. Its O decision uses

1. **purpose-built, sometimes partial markers** (`O-9,46_wbaV` at 1002 bp,
   `O-9,46_wzy_partial` at 216 bp, `tyr-O-9` at 609 bp) instead of full genes, so a
   fragmented locus still matches at high coverage — measured on this genome:
   `O-9,46_wbaV` 77 %, `tyr-O-9` 56 %;
2. a **relative** coverage gate (`> 15 % of the highest O hit`) rather than an
   absolute floor; and
3. a dedicated **functional test** to split O9 from O2, comparing `tyr-O-9` with
   `tyr-O-2` — the two are allelic, so the wzx/wzy references cannot separate them.

Items 1–2 are why the truncated locus still counted; item 3 is why the A/D1
ambiguity resolves upstream and cannot here (SISTR ships no `tyr` reference, and
SeqSero2's database is GPL-2.0). Extracting a `tyr`/`prt` pair from public-domain
genomes would be a legitimate curated addition if this class of error matters.

The lesson recorded here: before adding reference sequence, check whether the
genome already matches one at high identity — an assembly break and a
cross-labelled sibling both look like "no hit" from the outside.

## Open question (2026-11-16): O missing at the panel markers — assembly gaps, not database divergence

**Status: recorded for the next discussion; deliberately not implemented.**
Pending wet-lab validation (gap PCR / long reads) from the curator.

Three in-house genomes whose O family came out absent were dissected with
per-gene probes taken from the LT2 reference cluster (`NC_003197.2`,
2,157,000–2,182,000): `gnd → rfbP rfbK rfbM rfbN rfbU rfbV rfbX rfbJ rfbH rfbG
rfbF rfbI rfbC rfbA rfbD rfbB → galF`. A probe present at ≥90 % identity over
≥80 % coverage means the gene is assembled; the panel's O markers are `rfbX`
(the wzx/flippase — the shipped `o_B` reference aligns to LT2 `rfbX` exactly,
2,168,623–2,169,915) and its allele variants.

| genome | cluster in the assembly | gap | consequence |
| --- | --- | --- | --- |
| `26SAL0501` | O marker present as a **287/1293 bp (22.2 %) remnant** at 99.65–100 %, and **only** group B hits at ≥70 % id | none (truncation) | group B is inferable from the remnant — the case SeqSero2's relative gate (≥15 % of the best O hit) converts into a call |
| `26SAL0511` | **almost the whole cluster, 99.6–100 % identical to LT2 group B**: `gnd/rfbP/rfbK/rfbM/rfbN` on `contig00012`, `rfbV` 100 % on `contig00075`, `rfbJ…rfbB` on `contig00049` | **rfbU + rfbX (~2.4 kb)** — exactly the genes the panel uses as O markers | the isolate is almost certainly O:4/B, but the marker sits in the assembly gap |
| `26SAL0520` | only the cluster ends: `rfbP` (`contig00137`), `rfbK/rfbM` (`contig00121`), `rfbB/rfbD/rfbA` next to `galF` (`contig00004`) | **rfbN…rfbC (~10 kb)**, including every group-specific gene | no tool can infer O from this assembly |

**Not database divergence.** For `26SAL0520` there is **zero** O homology even
at 50 % identity / e ≤ 1e-3 (a present-but-divergent locus lights up related
groups: the D1-positive control hits six). For `26SAL0511` the surviving genes
are 99.9–100 % identical to group B, so the sequence is the panel's group —
only the marker gene is missing from the assembly.

**On adding `rfbU`/`rfbV` as fallback markers — declined for now.** A
two-genome check (group B vs Enteritidis D1: `rfbU`/`rfbV`/`rfbX` do not hit
the D1 cluster at ≥80 % id, while `rfbH` does at 99.4 %) is a hypothesis, not
evidence of group specificity. Accepting the markers requires a formal
comparative-genomics study: all 46 O groups as a negative set with multiple
strains per group, within-group diversity for the candidate genes, a
calibrated identity/coverage gate, and an end-to-end re-measurement on the
`26SAL` corpus so the change cannot introduce cross-group false calls. Until
that exists, the shipped behaviour (O absent → formula + H evidence + candidate
list, flagged as assembly QC) is the honest one.

**How to settle "true deletion vs assembly/sequencing artefact"** (cheapest
first): (1) map the original reads back to the LT2 cluster — reads spanning
the gap with normal depth prove the DNA carries the region and the assembler
lost it; (2) PCR across the gap — `rfbN`↔`rfbJ` for `26SAL0511` (~2.4 kb,
standard polymerase), `rfbM`↔`rfbC` for `26SAL0520` (~10 kb, long-range);
(3) long-read sequencing (ONT/HiFi) as the final arbiter, especially for the
10 kb gap; (4) compare a positive slide-agglutination result against the
molecular finding. If a deletion is confirmed **and** agglutination is
positive, that is a result worth recording on its own: it would say the
serological determinant does not live in `rfbN…rfbI`, which is also the reason
a single-gene remnant must not be promoted into a marker without the study
above.

## The union evidence model and paralogues

`coverage: "union"` sums a reference's HSP blocks, which is what a divergent or
split locus needs. But a subject hit perfectly at one site *and* weakly at a
second (a paralogous remnant, common in these assemblies) used to average into a
mediocre identity — a 98.9 % fljB block merged with an 80.6 % remnant reported
90.6 % and the antigen was lost. Blocks whose reference spans overlap by more than
half of the shorter one are now treated as alternate placements of one locus and
only the better is kept (`align/blast._distinct_loci`), so identity comes from one
locus while disjoint blocks (a locus split across contigs) still add up.

Additionally, a family name is no longer committed to when a sibling antigen is
within the document's `margin` (0.01): every name inside the margin is widened and
looked up (`_FamilyCall.ties`), and a serovar whose *observed* formula is absent
from the table carries a "matched after widening" note plus `confidence: low`.

## Screening cost: cache the sketch per window (2026-11-16)

The anchored-sketch annotation is a pure function of a hit's window and its
antigen family, but a locus produces one hit row per near-identical reference
allele — measured on one 4.8 Mb genome: 1,070 hits over 33 distinct windows,
and the per-hit recomputation spent 62 s of a 68 s screen (91 % of the run;
profile: `align/blast._annotate_sketch`). The annotation now caches per
`(window, family)`: the same genome screens in 9 s, 5.9 s of which is blastn
itself, with byte-identical sketches and phases, and a 20-genome batch finishes
in 41 s at four-way parallelism. Remaining headroom if it ever matters: the
`allele_sketches` index load (1,443 small `.full.u64` reads, ~0.5 s per process,
`lru_cache`d across files) and an inverted k-mer index for `antigen_scores`
(~1.7 s per genome today).

## Bundle size

`sen_sero` ships 2.9 MiB (sequences + `serotype.json` + `kmer.json` +
`bundled.json`); the packed k-mer sets are computed by `db setup`. They used to be
committed (8.2 MiB, the largest single item in the wheel) because
`scripts/salmonella/kmers_call_reads.py` materialized them into the source bundle
— it now uses a scratch directory. Bundled root: 20.6 MiB of the 30 MiB budget.

## Reads by k-mer (no alignment)

`import_upstream.py … --kmerdb 27` writes `kmer.json` (`gspec.kmerdb/1`): one
**differential k-mer set per antigen** — the antigen's representative reference
minus one representative of every other antigen in the same family (133 sets for
the SISTR panel). Reads are counted once and each set is scored by containment
(`set_hit_fraction`), which sidesteps the alignment path's weakness: with several
near-identical references per antigen, mapping splits the reads (`h1_i` 80 % next
to its sibling's 75 %), while k-mers give the true antigen 1.00 against a 0.11
runner-up on the same reads.

```bash
python scripts/salmonella/kmers_call_reads.py \
    --db-dir src/gspec/data/dbs/sen_sero --fastq reads.fq
#   o: B (1.00 … runner-up 0.00)   h1: i (1.00 … 0.11)   h2: 1,2 (1.00 … 0.07)
```

Caveats, measured: 128/133 antigens have a non-empty differential set — the five
that do not (`f,g,[s]`, `f,g,s`, `g,m,s`, `g,q`, `z4,z23,z32`) are wholly
contained in another antigen of their family, so k-mers cannot separate them
either. The subtraction uses one representative per other antigen, so the sets
are slightly permissive.

**The CLI uses it.** `bundled.py` now ships `kmer.json` inside the bundle and
materializes the sets on first use, and `run_screen_reads` routes a database with
a k-mer substrate away from minimap2: `gspec screen -d sen_sero --r1
reads.fq` returns one row per representative antigen (`BREADTH%` = the set's
containment fraction, `IDENTITY%` = the mash-style ANI estimate), which
`gspec typing` designates exactly like a mapped table. On the synthetic reads that
is `o_B`/`h1_i`/`h2_1_2` at 100.00 each → `Typhimurium`, where the alignment path
was splitting H1 evidence 80 %/75 %.

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

## 2026-11-17: the naming-tie problem (open, needs verified serotype data)

Measured on the in-house 100-genome *Salmonella enterica* corpus with the
shipped `sen_sero` panel and the three SNP-level O discriminators
(`rfbS` A/D1, `wbaV` D1/D2, `wbaP` E1/E4): **92/100 genomes agree with both
reference tools, 1 is ambiguous, and 7 differ** — but only **two** of those
differing classes are gspec defects, and both are the same defect: **a formula
that names more than one serovar.**

### The defect in one sentence

A WHO-2007 antigenic formula is not injective: **57 of the panel's 1478
formula keys (3.9 %) name 2–3 serovars**, and gspec resolves every one of them
by **declaration order** (first row wins, the twins ride in `RUNNER_UP` +
`formula shared with …`). On this corpus **18/100 genomes** land on such a key.
Four of them disagree with the references only because declaration order
happened to pick the other twin — the other 14 are "right" by luck, not by
evidence:

| shared formula | rows, in declaration order | corpus genomes hitting it |
| --- | --- | --- |
| `C2-C3:d:1,2` | Muenchen, Virginia | 5 |
| `C2-C3:k:1,5` | Blockley, Haardt | 1 |
| `C2-C3:r:1,5` vs `C2-C3:r,[i]:1,5` | Hindmarsh, Bovismorbificans | 2 (**wrong**) |
| `C2-C3:r:l,w` vs `C2-C3:r,[i]:l,w` | Goldcoast, Brikama | 2 (**wrong**) |
| `B:b:1,2` | Paratyphi B, Paratyphi B var. Java | 1 |
| `C2-C3:z10:e,n,x` | Hadar, Istanbul | 1 |
| (57 keys in total) | … | 18 |

**The `r` vs `r,[i]` pairs are one antigen written two ways** (`[i]` is a
bracketed optional sub-factor), so Hindmarsh ≡ Bovismorbificans and
Goldcoast ≡ Brikama *by formula alone*. Which twin is reported is decided by
which spelling won the H1 ranking by a few hundredths of a percent:

```
GCF_002061825.1  formula C2-C3:r:1,5      -> Hindmarsh   (H1 near-tie: r,[i])
GCF_004681245.1  formula C2-C3:r,[i]:l,w  -> Brikama     (H1 near-tie: r)
```

### How the two references resolve it (neither uses antigen evidence)

| tool | mechanism | evidence it is not antigenic |
| --- | --- | --- |
| SeqSero2 | a hard-coded 30-entry `rename_dict` (`Hindmarsh→Bovismorbificans`, `Bardo→Newport`, …) | the map contains no sequence data at all; for `Goldcoast`/`Brikama` it has **no** entry and SeqSero2 prints `Goldcoast or Brikama` |
| SISTR | cgMLST nearest reference genome | `cgMLST_distance` 0.0152–0.0273 against a named reference genome |

So both agree the antigenic formula cannot decide these pairs. **An honest
`A or B` is the correct output; picking one is not.**

### What was tried and did NOT work (do not repeat these)

**1. O-factor markers (comparative genomics, 13 complete genomes).**
O:6 distinguishes both pairs (Hindmarsh `8,20` vs Bovismorbificans `6,8,20`;
Goldcoast `6,8` vs Brikama `8,20`). Four complete Goldcoast (`6,8`), five
`6,8,20` (Bovismorbificans ×3 + Newport ×2) and four `8,20` (Albany) genomes
were compared over the rfb neighbourhood (±20 kb around `wzx`, 24 CDS each):

- 21 of 28 gene clusters are present in all 13 genomes (conserved backbone);
- the "clean" split by presence/absence (`wzxC`/`wcaK`/`wcaL` vs
  `wzzB`/`ugd`/`gndA`) is a **position artefact**: whole-genome BLAST finds
  those genes in *every* genome at 98.4–100 %, they merely sit at different
  distances from the rfb cluster;
- at allele level the best candidate (a `glycosyltransferase family 4`) has
  **97.55 % within-group vs 95.84 % between-group** identity. Every marker that
  actually works in this panel has **100 % within-group** (14/14 Enteritidis
  genomes are identical over `wbaV`), so this is an order of magnitude too
  noisy.

**2. Differential k-mers (mash/SalmID style).** Two granularities, k=31:

- whole genome (4.5–5.0 M k-mers each): **0** k-mers present in all of the
  9 "factor-6" genomes and absent from all 7 `8,20` genomes;
- rfb-targeted (41 kb window around `wzx`, ~41 k k-mers each): still **0**, and
  the containment distributions do not separate — **within-group 0.8429
  (0.734–1.000) vs between-group 0.8212 (0.746–0.999)**, a 0.022 gap with full
  overlap. `BOV_016454185` vs `HIN_001478005` (different O) reach 0.9993, while
  `HIN_002031165` vs `BRI_021104055` (same O) sit at 0.7343.

Why the same method works for subspecies but not here: the six subspecies
differ by 2–5 % of the whole genome (hundreds of thousands of k-mers, measured
100/100 on this corpus), while O:6/O:20 sit in one small rfb segment and are
drowned by strain-level variation.

**3. Data caveat — public serotype labels are not verified.** Comparative work
must use **agglutination-confirmed** strains. Two examples from this pull:
`GCF_004194605.3` is labelled *Albany* (`8,20`) yet SeqSero2 calls its O as `8`
and its rfb cluster groups with the `6,8` genomes; and the only public *Brikama*
assembly (`GCA_021104055.1`) is a divergent contig-level genome that shares only
0.73–0.75 k-mer containment with everything else, so it cannot anchor a marker.
A serovar label in an INSDC record is a submitter assertion, not a phenotype.

### Recommended resolution (not implemented — deliberately deferred)

1. **`A or B` output for one-to-many formulas** — report both names (the
   `formula shared with …` note already carries the data; today the first name
   is promoted to `PHENOTYPE` and the twin demoted to `RUNNER_UP`). This is the
   only change that is correct without new evidence, and it alone turns the 4
   "wrong" calls into 4 honest candidate pairs.
2. **A nomenclature merge table** (SeqSero2's `rename_dict` idea) for the pairs
   WHO has formally merged (`Hindmarsh→Bovismorbificans`, `Bardo→Newport`).
   Cheap, but it papers over the ambiguity rather than resolving it.
3. **cgMLST** (SISTR's answer) breaks every pair — but it needs a reference
   dataset that this bundle does not carry, and it is a new subsystem.
4. **An O:6/O:20 marker** would be the real fix. It requires verified
   serotype data (see the caveat above) and either a hand-curated rfb
   annotation or ORF prediction across a larger, cleanly labelled panel.

### The rest of the remaining disagreement (not naming ties)

| genome | gspec | SISTR | SeqSero2 | class |
| --- | --- | --- | --- | --- |
| `GCF_002053225.1` | *ambiguous* (`Chichester`, `Ouagadougou`, `Stratford`) | `I E1:g,[s],t:1,2` (a formula) | `I 1,3,19:g,s,t:1,2` (a formula) | **all three fail to name it**: H1 ties `g,[s],t`/`g,s,t`/`g,t`/`i` and H2 ties `1,2`/`1,2,7`/`1,6`/`1,5`, and the widened H1 class reaches SISTR's over-merged `g`-family class. The O group *is* resolved (E4 over E1, `wbaP` gap 1.18 pp) |
| `GCF_001831865.1` | **`38:k:-`** (no serovar) | `I P:k:-` | `I 38:k:-` | **fixed 2026-11-17, twice**: the `h2 k` hit was 87.15 % identity, so the old code invented `P:k:k` — the duplicate-phase rule (no WHO row has H1 == H2) now drops the weaker phase, and the serogroup→numeric-O map renders `P` as `38`. All three tools now state the **same formula** |
| `GCF_019710395.2` | `B:-:-` (no serovar) | `Fulica` + 18 candidates | `I 4:-:-` | **genome has no H gene at all**; no tool can name a serovar from O alone. gspec's answer is the honest one |
| `GCF_002046265.1` | `Blockley` | `Galiema` | `Blockley` | **not a gspec defect** — 2 of 3 agree; `C2-C3:k:1,5` is Blockley/Haardt, and SISTR's cgMLST lands elsewhere |

### Systemic quality signals worth fixing independently

- **No call reaches `high` confidence**: on the 100-genome corpus the panel
  emits 99 × `low` + 1 × `ambiguous` and **zero `high`**. `perfect` requires a
  100 %-identity/100 %-coverage hit *and* a single name for that gene, and
  real assemblies almost never clear both against 1400+ references. The
  confidence ladder is therefore not discriminating anything today.
- **220 formula rows were dropped as unreachable** at build time (serovars
  whose alleles no reference FASTA covers), so the panel cannot name them at
  all — a coverage gap recorded in `bundled.json`.

## Provenance

Sequences are extracted from RefSeq/GenBank annotated CDS (`NC_003197.2`,
`NC_011294.1`, `NC_006511.1`, `NC_003198.1`) via NCBI E-utilities. INSDC
sequence data carries no copyright restrictions; the selection, naming and
threshold calibration are gspec's own, which is what the bundled-database
curation gate requires. `bundled.json` carries the accessions, the retrieval
date and the caveats so `gspec db info -d sen_sero` shows them to a user.
