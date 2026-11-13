# Campylobacter database builders

Developer scripts behind the bundled `cam_vf` database.

## `cam_vf` (gene)

`build_cam_vf.py` cuts the sixteen virulence markers from public INSDC
records (the license-hygiene rebuild of the earlier locally-built panel —
sequences byte-identical except `cgtB`, strand-corrected to the coding
orientation; annotations re-derived from INSDC, no third-party database
terms):

- NCTC 11168 chromosome 1 `NC_002163.1`: cdtA/B/C, ciaB, wlaN, and the
  cgtA/neuA1 fusion (annotated `neuA1`);
- 81-176 plasmid pVir `AF226280.2` (RefSeq NC_005012): virB4, virB11;
- 81-176 plasmid pVir `CP000550.1` (RefSeq NC_008770): virB8/9/10;
- NCTC 13255 `LR134499.1`: ggt;
- strain 108 T6SS locus `JX436460.1`: hcp (annotated TssD);
- `AF400048.1`: cstII; `AF130984.1` (OH4384): cgtA, cgtB.

```bash
for id in NC_002163.1gbwithparts AF226280.2 CP000550.1 LR134499.1 JX436460.1 AF400048.1 AF130984.1; do
  ret=gb; [ "$id" = NC_002163.1gbwithparts ] && ret=gbwithparts
  curl -o "${id%%.*}.gb" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${id%%.*}&rettype=$ret&retmode=text"
done
python scripts/campylobacter/build_cam_vf.py \
    --nctc11168 NC_002163.gb --pvir-af AF226280.gb --pvir-cp CP000550.gb \
    --nctc13255 LR134499.gb --t6ss JX436460.gb --cstii AF400048.gb \
    --los AF130984.gb --dest src/gspec/data/dbs/cam_vf
```

The builder emits only public-domain INSDC content; no GPL or
non-commercial database files are copied.
