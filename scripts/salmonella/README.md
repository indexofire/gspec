# sen_sero curation — Salmonella enterica subsp. enterica serotyping panel

`sen_sero` is a **gspec-curated gene panel** for *Salmonella enterica* subsp.
*enterica* molecular serotyping, built from annotated INSDC CDS records.

## Scope (deliberate)

Only **subsp. enterica** is targeted. The antigenic formula of the other five
subspecies (salamae, arizonae, diarizonae, houtenae, indica) can be read the
same way, but *which* subspecies a genome belongs to cannot be decided from
antigen genes — it needs biochemical tests. A panel that claimed to type them
would be wrong about the subspecies it printed, so nothing here does.

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
| `rfbJ` (group B, CDP-abequose synthase) vs everything else | 13–36 % | group B (O:4) *is* separable at the gene level |
| `fliC_i` (H:i) vs `fliC_r` (H:r) | 87.8–90.4 % | closest H1 pair in the panel; still far below the 98 % floor |
| `fljB` H:1,2 of Typhimurium vs Heidelberg | 99.9 % | one shared H2 marker covers both group-B serovars |
| `fliC_gm` (H:g,m) vs `fliC_gp` (H:g,p) | **99.80 %** (3 nt in 1518) | **inseparable by template identity at any sane floor** — the reason Dublin/Gallinarum/Pullorum are out of scope |
| `invA` vs everything | ≤25 % | control gene is unambiguous |

The panel therefore ships **five serovars**: Typhimurium and Heidelberg (group B,
`rfbJ`), Enteritidis and Typhi (group D1, `rfbE`), and the provisional
Paratyphi A (group A, `rfbS`, pending the rfb cluster stage).

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

## Known limitations of this first slice

1. **O:2 (group A) vs O:9 (group D1) is not resolved by gene presence.** The
   two groups differ in whether the paratose 2-epimerase is functional, not in
   which genes are present (Paratyphi A carries a pseudogene), so a presence
   panel cannot separate them. The next increment is an `rfb`-locus **cluster**
   stage (minimap2 asm20 against per-O-group locus references), the SeqSero2
   shape, fused with this panel through a composite database.
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
