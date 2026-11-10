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

The port was run against a local SISTR install
(`…/site-packages/sistr/data`, SISTR db 1.1.3) and ships as `sen_sero_sistr`:

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
  exact first and only widens a token to its class when nothing matches, which
  is SISTR's own order; a widened key that fits several rows is reported as an
  ambiguous call with the candidates, not pipe-joined. That took reachability
  from 95.4 % to 97.4 %.

The remaining 41 rows name alleles that SISTR's reference FASTAs themselves do
not contain (`z67`, `z91`, `g,m,t`, …). `gspec db info -d sen_sero_sistr` shows
the Apache-2.0 attribution; `tests/test_sen_sero_sistr.py` pins the numbers and
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
   form of SISTR's pipe-joined `Typhimurium|Lagos`.
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
python scripts/salmonella/regress_serovar.py --db sen_sero_sistr \
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
h1: g,m,q == g,q     g,p == g,p,s == g,p,u     g,t == i
h2: z35 == z52
```

These are folded into the document's `similarity_groups`
(`import_upstream.py --cross-labels`), so both names' rows stay reachable and the
answer is an honest ambiguity instead of whichever label won on bitscore. The audit
also flagged a record whose antigen name is empty (`h1__f_g_t`), which the merge
skips.

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

## Bundle size

`sen_sero_sistr` ships 2.9 MiB (sequences + `serotype.json` + `kmer.json` +
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
    --db-dir src/gspec/data/dbs/sen_sero_sistr --fastq reads.fq
#   o: B (1.00 … runner-up 0.00)   h1: i (1.00 … 0.11)   h2: 1,2 (1.00 … 0.07)
```

Caveats, measured: 128/133 antigens have a non-empty differential set — the five
that do not (`f,g,[s]`, `f,g,s`, `g,m,s`, `g,q`, `z4,z23,z32`) are wholly
contained in another antigen of their family, so k-mers cannot separate them
either. The subtraction uses one representative per other antigen, so the sets
are slightly permissive.

**The CLI uses it.** `bundled.py` now ships `kmer.json` inside the bundle and
materializes the sets on first use, and `run_screen_reads` routes a database with
a k-mer substrate away from minimap2: `gspec screen -d sen_sero_sistr --r1
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

## Provenance

Sequences are extracted from RefSeq/GenBank annotated CDS (`NC_003197.2`,
`NC_011294.1`, `NC_006511.1`, `NC_003198.1`) via NCBI E-utilities. INSDC
sequence data carries no copyright restrictions; the selection, naming and
threshold calibration are gspec's own, which is what the bundled-database
curation gate requires. `bundled.json` carries the accessions, the retrieval
date and the caveats so `gspec db info -d sen_sero` shows them to a user.
