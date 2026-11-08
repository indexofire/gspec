# Composite databases: several panels and methods in one database

A composite database is a database directory whose `composite.json` plans an
analysis instead of describing one panel: several databases, several methods,
fused into one typing result. It is discovered by `gspec db list` (`KIND=composite`),
addressed by `gspec screen -d NAME`, and built with `gspec db build NAME -C
composite.json --bundle`. Real serotyping workflows need several
evidence layers whose tools, thresholds and databases differ: a marker panel (blastn, strict
identity, tiny references), a full-locus alignment (blastn with union coverage, or
minimap2), an optional PCR/SNP layer. Until now those layers were stitched together by
scripts, which put the *fusion policy* — not just its parameters — outside configuration.

A **composite document** (`gspec.composite/1`) makes the fusion declarative: an ordered list of
stages, each naming a database (and optionally one of its schemes) plus the acceptance rule
for that stage's call. Stages evaluate in order and the **first accepted call wins**; every
stage's raw call is kept as evidence, so a two-layer call carries its audit trail.

```json
{
  "schema": "gspec.composite/1",
  "name": "cje_penner",
  "fallback": "untypeable (atypical/novel CPS locus or non-jejuni)",
  "stages": [
    {"name": "marker",  "db": "cje_hs_marker", "scheme": "hs_marker_vote",
     "tier": "B", "min_score": 0.30},
    {"name": "locus_a", "db": "cje_cps_locus",  "scheme": "tier_a",
     "tier": "A", "min_score": 0.85},
    {"name": "locus_b", "db": "cje_cps_locus",  "scheme": "tier_b",
     "tier": "B", "min_score": 0.50},
    {"name": "locus_c", "db": "cje_cps_locus",  "scheme": "tier_c",
     "tier": "C", "min_score": 0.30}
  ]
}
```

```console
$ gspec db build cje_penner -C composite.json --bundle
$ gspec screen -d cje_penner -t 16 -j 8 *.fna
FILE                  COMPOSITE    PHENOTYPE  TIER  STAGE     SCORE   EVIDENCE
GCF_000009085.1.fna   cje_penner   HS2        B     marker    1.0000  marker=HS2(1.0000)|locus_a=HS4c (HS4/13/43/50/62/65 complex)(0.9957)|…
GCF_000163975.1.fna   cje_penner   HS37       B     locus_b   0.7042  marker=-|locus_a=-|locus_b=HS37(0.7042)|locus_c=HS37(0.7042)
```

`--format json` emits `gspec.composite_result/1`: the fused call plus a per-stage
`evidence` map (`phenotype`, `score`, `confidence`, `scheme`), so downstream consumers can
audit or re-fuse without re-aligning.

## Stage fields

| Field | Meaning |
|---|---|
| `name` | Evidence label; must be unique across the composite |
| `db` | Installed database to screen and type |
| `scheme` | Scheme name inside that database; omit when the document declares exactly one |
| `tier` | Label that flows to the result's `TIER` column (`A`, `B`, `C`, …) |
| `min_score` | Reject a stage's call below this score (0..1) |
| `require_confidence` | `any` (default) or `high` (accept only high-confidence calls) |
| `min_cluster_cov` / `min_gene_cov` | Cluster databases only: the `-C` / `-g` gate for this stage, so one locus database can serve several tiers |

## Installing a composite (so `-p NAME` works from anywhere)

A composite **is** a database — one whose document plans several panels — so it is
created, listed and addressed like any other:

```console
$ gspec db build NAME -C cje_penner.json          # validate + copy (+ primer sidecars)
~/.local/share/gspec/db/<composite>/cje_penner/composite.json
$ gspec db list
cje_penner	marker_a(db:cje_hs_marker), …, pcr(pcr:primers)
$ gspec screen -d cje_penner -t 16 -j 8 *.fna
```

Layout (one directory per composite, so relative primer files travel with it):

```text
<datadir>/<composite>/<name>/composite.json
<datadir>/<composite>/<name>/<primer files>
```

`--install` validates the document, copies it as `composite.json` plus every relative
primer file it references, and checks that every `kind: "db"` stage's database is
installed (a PCR stage's `db` is a label, not a database). `--force` replaces an
install; `--name` renames. `-p` accepts a path or an installed name — a path that
exists always wins, so older invocations keep working.

### Self-contained bundle (`--bundle`)

A composite is often handed to someone else, and the panels it names are the
bulk of what has to travel with it. `--bundle` copies them in and pins `db_root`, so
the installed directory runs on a machine that has none of those databases:

```console
$ gspec db build NAME -C cje_penner.json --bundle
$ find ~/.local/share/gspec/db/<composite>/cje_penner -maxdepth 2 -type d
…/<composite>/cje_penner
…/<composite>/cje_penner/dbs
…/<composite>/cje_penner/dbs/cje_hs_marker
…/<composite>/cje_penner/dbs/cje_cps_locus
```

```json
{"schema": "gspec.composite/1", "name": "cje_penner", "db_root": "dbs", "stages": [ … ]}
```

A composite database directory therefore holds `composite.json`, its primer
sidecars, its manifest and (with `--bundle`) `dbs/`. Copying the directory into a
datadir — or pointing `-D` at a directory that contains it — is all a colleague
needs.

Resolution order for a stage's `db`: the bundle's `dbs/<name>` first, then the
datadir. Copy the directory anywhere and run it with an empty datadir — the result
is identical (verified byte for byte against the installed databases). `db_root`
must be relative to the document and may not traverse upwards.

Bundling is **not** the same as merging the databases into one: the layers
(tier schemes, PCR stage, conflict rules, per-stage evidence) are preserved, only
the distribution unit becomes a single directory.

## How a run executes

1. **Group by database.** Stages naming the same database share one screen run, so a locus
   database carrying `tier_a`/`tier_b`/`tier_c` schemes is aligned once and typed once
   (three stages, one alignment). Cluster stages re-screen only when their gate overrides
   differ.
2. **Screen.** Each stage's database is screened through the same path `gspec screen` uses,
   so its own `screen.json` defaults apply — including the `gspec.screen/2` evidence model
   (`coverage: "union"`, `dedup: "subject"`, `min_aligned_len`), which is what makes a
   divergent or repeat-laden locus measurable at all.
3. **Type.** Gene databases designate per file through their typing document; cluster
   databases are already designated by their screen (`gspec.cluster/1` or `/2`), and the
   composite reads the scheme's call from it.
4. **Fuse.** Stages are visited in document order; the first accepted call wins. Nothing
   about thresholds or phenotype logic lives in the composite — it only decides which stage's
   call wins and keeps every stage's evidence.

## Why not one database?

Three reasons, all of them structural rather than stylistic:

- **A database is one engine, one input kind.** A gene database is blastn over a FASTA
  panel; a cluster database is minimap2 over annotated loci. Fusing them inside one database
  would mean one engine or two incompatible evidence models.
- **Cluster rules have no locus scope.** `cluster_match` scores *the best locus* of the
  report for every rule, so a rule cannot say "only evaluate the HS4c reference". Per-tier
  specificity belongs to a gene database plus `best_hit`, not to a cluster rule set.
- **Fallback needs an order.** The primitives (`cutoff`, `ambiguity_margin`, `compose`) rank
  *inside* one scheme. "Marker first, then tiers, else untypeable" is an ordering across
  schemes and databases — that is what a composite declares.

## PCR stages: refine a call with an assay

### Stage kinds: the methods a scheme can mix

A composite document is an **analysis scheme**: ordered evidence stages, each one
method, fused into one call. Four method families are available today, and a
document may combine them freely:

| `kind` | Method | Engine | Evidence it contributes |
|---|---|---|---|
| `db` (gene) | marker / allele panel | blastn (+ `gspec.screen/2` coverage model) | per-gene coverage/identity → `weighted_genes` / `exact_set` / `best_hit` |
| `db` (cluster) | whole-locus typing | minimap2 (asm20) | locus coverage/identity + key genes/k-mers → `cluster_match` |
| `kmer` | alignment-free k-mer substrate | none (pure Python, process-parallel) | differential k-mer set hit fraction → phenotype |
| `pcr` | in-silico PCR | `gseek.ispcr` (pure stdlib) | product size per pair → override / subtype bracket |

A `kmer` stage reads the database's `gspec.kmerdb/1` sidecar (the same substrate
`spn_sero` uses for SNP-level pairs) and maps sets to phenotypes declaratively:

```json
{"name": "snps", "kind": "kmer", "db": "cje_cps_kmer", "tier": "K",
 "sets": [{"set": "HS4.diff", "phenotype": "HS4c (HS4/13/43/50/62/65 complex)",
           "min_fraction": 0.5}]}
```

The same scheme can therefore read blastn markers, blastn loci, minimap2 loci,
k-mer sets and PCR products, in the order that matches the assay it implements.
Mixing is measured, not assumed: on the *C. jejuni* Penner set the four method
families agree with the CjejuniTyper reference at 100% (marker tier A), 94.2%
(minimap2 locus db) and 94.8% (k-mer sets) independently, and the fused answer is
98.6% on the 600-genome check set.

A `kind: "pcr"` stage is the refinement half of a composite. It does not designate
from scratch: it runs an in-silico PCR (`gseek.ispcr`, pure stdlib — no exonerate
or ipcress) and either overrides the call the cascade made (the assay is
authoritative) or brackets a complex type into a subtype. Everything is declared:
the primer set, the mismatch budget, and the ordered rules.

```json
{"name": "pcr", "kind": "pcr", "db": "primers", "tier": "PCR",
 "primers": "cje_penner_primers.json", "max_mismatches": 2,
 "rules": [
   {"products": [{"name": "HS55_target", "min_size": 300, "max_size": 400}],
    "set": "HS55"},
   {"when": ["HS4c"],
    "products": [{"name": "HS4B_target", "min_size": 625, "max_size": 680}],
    "set": "HS4c (HS16/64)"}
 ]}
```

| Field | Meaning |
|---|---|
| `primers` | `gseek.primers/1` JSON or a TSV table (`name`, `forward`, `reverse`, optional `min_size`, `max_size`, `max_mismatches`); a relative path resolves next to the composite document |
| `max_mismatches` | Substitution budget per primer (a pair may override it in the table) |
| `rules` | Ordered rules; the first applicable one wins, mirroring a published if/elif cascade |
| `rules[].products` | Every listed product must amplify inside its own size window |
| `rules[].when` | Restrict the rule to calls whose base type matches (`HS4c (HS4 … complex)` matches `"HS4c"`) — the bracket case |
| `rules[].set` | The literal new phenotype, brackets included |
| `report_all_products` | Scan the whole primer set for evidence instead of only the pairs the rules name (default off: fewer genome scans) |

Semantics worth knowing:

- **Only the rule-named primer pairs are scanned** by default, so a 30-pair
  reference panel costs the runtime of the 7 pairs the rules actually use.
- **A PCR stage never invents a call**: without a matching product the cascade's
  call, tier and stage stay untouched, and `EVIDENCE` shows `pcr=-`.
- **The call's evidence is auditable**: `EVIDENCE` renders `pcr=HS55_target(334bp)`.
- **`when` gates on the base type**, so a bracket rule cannot fire on an
  unrelated call, and a bracketed call keeps the assay's stage/tier (`PCR`).
- **Substitutions only**: primer matching allows `max_mismatches` substitutions
  per primer; indel-bearing primer sites (which ipcress aligns) are not modelled.

## Running a scheme through `screen`

An installed scheme is also addressable as a database name, because that is how
most people look for one:

```console
$ gspec screen -d cje_penner_quad -t 16 -j 8 *.fna     # fused typing result
```

`screen` refuses database-only flags for a scheme (each stage carries its own
thresholds in the document) and emits the same table as `gspec screen`
(`gspec.composite_result/1` under `--format json`). `gspec screen` stays for
path-based runs, `--install/--list`, and iteration on a document that is not
installed yet.

## Relation to `gspec typing`

`gspec typing` designates one database's screen table and renders one row per (file, scheme).
`gspec screen` orchestrates several of those and renders one row per file. Use `typing`
when one panel answers your question, and a composite when a call needs several layers.

A composite names a *plan*, not a panel: `screen -d <composite>` runs it end to end.
For a single-panel alternative with the same content, build the panels into one database
instead (see the `cje_penner_single` example) — that trades the layering for a plain
`screen`/`typing` pair.
