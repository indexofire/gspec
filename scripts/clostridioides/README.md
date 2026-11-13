# Clostridioides database builders

Developer scripts behind the bundled `cdi_vf` database.

## `cdi_vf` (gene)

`build_cdi_vf.py` cuts the four toxin markers from public INSDC records
(the license-hygiene rebuild of the earlier locally-built panel —
sequences byte-identical; annotations re-derived from INSDC, no
third-party database terms):

- strain 630 chromosome `NC_009089.1`: tcdA (WP_011860904.1), tcdB
  (WP_009902069.1), the PaLoc toxins;
- the CCUG 20309 CDT locus `AF271719.1`: cdtA (AAF81760.1), cdtB
  (AAF81761.1).

```bash
curl -o NC_009089.gb "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_009089.1&rettype=gbwithparts&retmode=text"
curl -o AF271719.gb "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AF271719.1&rettype=gb&retmode=text"
python scripts/clostridioides/build_cdi_vf.py \
    --six30 NC_009089.gb --cdt AF271719.gb \
    --dest src/gspec/data/dbs/cdi_vf
```

The builder emits only public-domain INSDC content; no GPL or
non-commercial database files are copied.
