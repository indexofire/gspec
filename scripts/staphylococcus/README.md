# sau_ent curation — S. aureus enterotoxin panel

`sau_ent` is the bundled enterotoxin gene panel: 25 references converted from
the MIT-licensed `indexofire/enterotoxin` abricate database, eleven curated
additions cut from public-domain INSDC records, and the `nuc` species control
(37 records total).
Build with `build_sau_ent.py` (offline: assets are committed; upstream is
pinned to commit `30c022eeb712`).

## Literature review (2026-11-17): polymorphism and missing types

The review that motivated the additions, kept here for the next curation
pass.

**Polymorphism is the family norm, not the exception.** Lefebvre et al.
(Toxins 14:249, 2022) quantified it while building their 24-toxin LC-HRMS
panel: eight genes carry six to ten variants each — sec (10), seu (10),
sem (9), sen (9), sek (8), seg (6), sey (6), tsst1 (6). Single-gene studies
agree: SEC has ≥6 variants incl. bovine/ovine forms whose ~3 amino-acid
differences shift host tropism (Marr 1993; Etter 2020; Johler 2016); SEB
variants alter Vβ specificity (Kohler 2012); seg/sei (Blaiotta 2004/2006)
and sek (Aguilar 2014) vary; and the egc locus is a recombination "nursery"
continuously minting new genes and pseudogenes (Jarraud 2001; Thomas 2006;
Dicks 2021). Practical consequence: one reference per type under-detects —
a divergent allele can fall below the default 80 %/80 % floors.

**Types the upstream panel missed**: the audited addition table lives in the
third-pass section below (the first-pass rows for sel27/selw predate the
2026-11-18 nomenclature audit and were superseded by it).

**Third pass (2026-11-17 evening → 2026-11-18 full audit): the
sel26/sel27/sel28/SElW nomenclature, audited against every reachable
primary source.** The audit was prompted by contradictory anchors — the
2025 review table says SEl28 = MN450304, 250 aa (the 753 bp gene), while the
2026 nomenclature-clarification review (Franjus et al., "Staphylococcal
enterotoxins and food poisoning: new players in the game", HAL hal-05698176)
says SEl28 = 726 nt. Every anchor, weighed:

| anchor | sel26 | sel27 | sel28 |
| --- | --- | --- | --- |
| Zhang IJMM 2018 (PMID 29574061, abstract verbatim: "designated SEl26 and 27") + the MF370874/78 CDS product labels | 744 bp (MF370874:2651-3394) | 753 bp (MF370874:3420-4172) | — |
| Franjus 2026 (dedicated nomenclature ruling, table adapted from Hu/Zhu) | SEl26-1 = 744 nt | SEl27-1 = 753 nt | **SEl26-2 = 726 nt** |
| MF370874/MF370878 description lines ("SEl27, and SEl28 genes") — no paper behind them | — | 744 bp | 753 bp |
| MN450304 CDS labels (unpublished Fanelli deposit) | — | 726 bp | 753 bp |
| 2025 review (PMC12733317, follows MN450304's labels) | — | 726 bp? | 753 bp (250 aa) |

Ruling adopted: the **primary designation paper + the dedicated 2026
clarification** — `sel26`/`sel27` keep Zhang's MF370874 sequences and
`sel28` = the 726 bp MN450304:2472-3197 gene (a "SEl26-family variant",
TTG start, ≤11 % aa to every other panel gene). The description-line and
MN450304 numberings follow no peer-reviewed paper and are not adopted. An
earlier pass of this panel had adopted the MN450304 numbering and was
reverted after this audit.

**SElW / SElU2, same audit.** EF030428 (Thomas 2006, 771 bp) is the gene
Collery & Smyth 2007 renamed "selw" — the name most WGS tools use. Franjus
2026 restores **SElU2** for it: the 2012 "true SElW" (UniProt A0A0H2WYW8,
705 bp, chromosomal, Wilson-era) is a DIFFERENT protein with no emetic
activity, which their review removes from the SE family altogether.
Measured here: seu (AY205305) vs EF030428 = 52 % anchored identity — a
distinct gene, not a seu variant (Franjus's "99 % SElU/SElU2" claim
concerns other, ψent-recombinant seu alleles). The panel label is now
`selu2` with the synonymy in the product string; the 2012 SElW is
deliberately not added (excluded from the family by the ruling review).

Also confirmed during the audit: the MF370874/75/76 "SEl28" description
mention has no CDS behind it (the three hypotheticals translate to 8-10 %
aa vs any SE — not family members), and MF370878 is internally inconsistent
the same way (description "SEl27 and SEl28", features SEl26/SEl27).
Franjus's other rulings recorded for the future: SEl02 = SEZ-2 (variant);
SEl01 (Hisatsune) and SEl29p–SEl33 (Dicks) "should not be designated as
SE" pending validation — consistent with this panel's exclusions.

| gap | source added | citation |
| --- | --- | --- |
| `seu` (egc2/3 core; outbreaks with no classical SE) | AY205305 (allele diversity 99.1–99.9 % ⇒ one reference detects all) | Letertre 2003; Umeda 2017 |
| `selu2` (egc4; widely called selw) | EF030428 | Thomas 2006; naming per Franjus 2026 |
| `sec` SEC3 allele (96.6 % to SEC1) | X51661 | Hovde 1990 |
| `sec` bovine allele (SaPIbov) | AF217235 complement(1229..2044) | Fitzgerald 2001 |
| `tsst1` bovine-SaPI allele (97.6 % to prototype) | AF217235 4242..4946 | Fitzgerald 2001 |
| egc3 variant alleles for `sei`/`seu`/`sen` (96.7–97.7 %) | AY158703 82..810, 964..1749, 1788..2543 (strain 382F; SEG variant on the same record is partial, skipped) | Letertre et al. |
| `sel28` (SEl26-family variant; TTG start) | MN450304 2472..3197 | Fanelli deposit; naming per Zhang 2018 + Franjus 2026 |
| `sel30` | PZ361155 32..784 (strain PR-14) | record annotation; Dicks 2021 |

**Second pass (2026-11-17, same day): what was checked and rejected or
deferred.**

- `sec2`: **no complete standalone INSDC record exists** — the ENA
  description search surfaces only partial cds entries plus LC074874
  ("complete cds" in the title but a 602 bp record holding a 360 bp CDS;
  rejected). SEC1↔SEC2 protein identity is 97.37 %, so the SEC1 nucleotide
  reference detects SEC2 genomes at ~96–97 % and reports them as `sec`
  (type-level). A dedicated allele waits on the original 1980s–90s
  submission or an annotated carrier genome.
- KF729631/AB084256/KT223833 (complete "enterotoxin C" single-gene
  records): extracted and compared — KF729631 is 100 % the SaPIbov allele
  we already ship, AB084256 is 100 % SEC3, KT223833 (Mu50) is 99.4 % SEC3.
  All redundant; skipped and recorded.
- MH737695 (seu "mRNA, complete cds"): 99.62 % to our seu — redundant.
- `sel28`–`sel33` (Dicks 2021): no UniProt entries, and the carrying NCTC
  strains are listed only in the paper's PDF (unparseable in the dev
  environment). Deferred: needs the strain table, then per-genome homology
  extraction. Modern-isolate prevalence is low (sel29p absent from any
  contemporary strain), so the marginal detection value does not justify
  rushed curation.
- `SE01`/`SE02` (Suzuki 2020 — SE02 has demonstrated emetic activity and a
  real outbreak): sequences live in the paywalled full text / the outbreak
  strains' DDBJ WGS records, unreachable this session (ScienceDirect and
  NCBI both blocked). Deferred with a precise recovery plan: find the DDBJ
  accession in the paper's data-availability statement, fetch the genome,
  locate `se2` by annotation or SEB-group homology, cut and validate the
  CDS.

**Deliberately not added (yet):**

- The remaining sec/seu/sem/sen/sek/seg/sey variant alleles (Lefebvre's
  counts) — each is a candidate for the same INSDC-extraction route when
  divergence data justify it.
- `sel34`/`sel35`: a 2026 bioRxiv preprint's machine-learning calls, not
  peer-reviewed — watch, do not bundle.
- egc pseudogenes (`psient1/2`, `seh-2p`, `ses-3p`): screening pseudogenes
  invites false calls; excluded on purpose.

**Licensing note:** NAuRA (ANSES) is the reference allele-aware screening
tool for this family, but it is GPL-3.0 — its sequences can never ride the
MIT wheel. Everything in this panel is re-sourced from public-domain INSDC
records; the `aatxaerandio/SE_database` UniProt catalogue (no license
declared) was used as a variant inventory only, no sequences taken.

## Phenotype backtest (2026-11-18): Fanelli et al. 2022, Toxins 14:731

The dataset the panel was missing: **12 Algerian food isolates with
whole-genome sequences (ENA BioProject PRJNA712985, WGS masters
JAHLTR01/JAFNJM01/JAFNJN01/JAHLTS01/JAHLTT01/JAHLTU01/JAFNJO01/JAFNJP01/
JAFRED01/JAFNJQ01/JAFNJR01/JAHLTV01) AND ELISA-quantified SE production**
(SA04 SEA 156-345 ng/mL + SEH; SA18 SEB 2871-14739; SA20 SED 322-399 +
SER 36720-63176; SA46 SEH n.d.; SA51 SEA 64-75). One ENA request per WGS
master (`ena/browser/api/fasta/<PREFIX>0100000001`, gzipped) fetches the
whole contig set.

**Result: 63/63 expected gene calls detected (recall 100 %), 4 extra calls
(precision 94.0 %), all four explainable and none a panel error:**

| call | evidence | verdict |
| --- | --- | --- |
| `selu2` in SA02/SA51 (egc1 per the paper) | 99.5-99.6 % id, **100 % cov** vs EF030428 | full-length selu2-homologous sequence is present; the paper tables the locus as psi-ent pseudogenes — a table-vs-genome discrepancy, and nucleotide screening cannot judge functional status either way |
| `sed` in SA20 | 99.6 % id, 99.9 % cov | the paper's "pseudogene" is a full-length nucleotide with an internal stop — invisible to blastn |
| `seh` in SA46 | 87.4 % id, 97.9 % cov | divergent pseudogene allele; above the 80/80 floors |

The backtest also **found a real gap**: the upstream panel never carried
`seo` at all (all 5 egc strains missed it). Fixed with the canonical
AF285760:370-1155 (Jarraud 2001's original egc cluster deposit, 786 bp,
ATT start; 100 % identity to SA02/SA82). LC074825 ("enterotoxin O,
complete cds") was checked and rejected — its CDS is 618 bp, a truncated
seo. After the fix: **recall 100 %**.

Truth table (paper Table 1, panel gene ids): SA01 selX; SA02
seg/sei/sem/sen/seo/selX; SA04 sea/seh/sek/seq/selX; SA07 selX; SA08
seg/sei/sem/sen/seo/selu2/selX; SA10 selX; SA18 seb/selX; SA20
sej/ser/selX (+sed pseudogene); SA24 = SA08; SA46 sey/selX (+seh
pseudogene); SA51 sea/seg/sei/sem/sen/seo/selX/sel30; SA82
seg/sei/sem/sen/seo/selu2/selX/tsst1. The paper's "selw" (Okumura 2012,
chromosomal, all 12 strains) is deliberately NOT a panel target — it is
the non-emetic 2012 protein the 2026 review excludes (see the audit above).
