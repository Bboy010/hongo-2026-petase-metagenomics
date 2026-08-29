# Statistical analyses and figures — ACMI-D-26-00103

Reproducible code for the statistical analyses and the R-generated figures of:

> *Novel putative PETase candidates from metagenomic mining of Ébrié lagoon and
> Kassembié lake, Côte d'Ivoire.* — Access Microbiology, manuscript
> ACMI-D-26-00103.

Everything is in a **single script**, [`unique_analysis_file.R`](unique_analysis_file.R).

## Layout

```
hongo-2026-petase-metagenomics/   # repository root
├── unique_analysis_file.R   # all analyses and all figures
├── data/                    # inputs (see below)
├── Pictures/                # figures written by the script
└── results/                 # tables written by the script
```

## Running it

Requires **R ≥ 4.4** (tested with 4.5.2) and the packages `readr`, `dplyr`,
`tidyr`, `tibble`, `stringr`, `ggplot2`, `ggrepel`, `RColorBrewer`, `pheatmap`,
`readxl`.

```r
install.packages(c("readr", "dplyr", "tidyr", "tibble", "stringr",
                   "ggplot2", "ggrepel", "RColorBrewer", "pheatmap", "readxl"))
setwd("hongo-2026-petase-metagenomics")   # the repository root
source("unique_analysis_file.R")
```

or, from a terminal:

```sh
git clone https://github.com/Bboy010/hongo-2026-petase-metagenomics.git
cd hongo-2026-petase-metagenomics
Rscript unique_analysis_file.R
```

The script creates `Pictures/` and `results/` if they do not exist, prints the
quality-control checks as it goes, and stops with an explicit error if an input
or a package is missing. A full run takes about two minutes, most of it spent
reading the eggNOG-mapper annotations.

## Which figure comes from where

| Manuscript | Content | Produced by |
|---|---|---|
| Figure 1 | Pipeline schematic | not R (vector drawing) |
| **Figure 2** | Differential-abundance heatmap | script, section 2.4 |
| **Figure 3** | Relative abundance, 20 most abundant species | script, section 2.2 |
| Figure 4 | Circular phylogenomic tree | iTOL v7 — the cross-sample detection check behind it is in section 3 |
| **Figure 5** | KEGG xenobiotic biodegradation pathways | script, section 4.1 |
| **Figure 6** | CAZyme classes | script, section 5.1 |
| **Figure 7** | MAG × CAZyme-class heatmap | script, section 5.2 |
| Figure 8 | ColabFold model, coverage and pLDDT | ColabFold output |
| Figure 9 | Multiple sequence alignment | ESPript |
| Figure S1 | Krona charts | Krona HTML export |
| **Figure S2** | Relative abundance of the *uncharacterized* SGBs, 5 per site | script, section 2.3 |
| **Figure S3** | Volcano plot, log2 fold change | script, section 2.5 |
| **Figure S4** | eggNOG COG functional categories | script, section 4.2 |
| **Figure S5** | Enzyme Commission classes | script, section 4.3 |
| **Figure S6** | 100 most abundant CAZyme subclasses | script, section 5.3 |

Figures 3, 5, 6, S5 and S6 regenerate **byte-for-byte identical** to the images
published in the article; Figures 2, 7, S3 and S4 are re-rendered at higher
resolution than the copies embedded in the manuscript.

**Figure S2 differs from the published image on purpose.** The caption, the
Results section and the Supplementary Materials list all describe it as "the ten
most prevalent *uncharacterized* SGBs", but the published image ranks all 730
SGBs and therefore shows eight named species. The script restricts the ranking
to the 131 uncharacterized SGBs (76 detected only at Biétry, 55 only at
Kassembié — the figures quoted in the Results), and takes the five most abundant
at each site: those SGBs are strictly site-exclusive and the six most abundant
ones all occur in Kassembie_1, so a pooled ranking would collapse the figure
into a single bar and hide the Biétry taxa.

## Input data

| File | Content |
|---|---|
| `merged_abundance_table_species.csv` | MetaPhlAn 4 merged profile (`mpa_vJan25_CHOCOPhlAnSGB_202503`), `;`-separated, 5 samples, species and SGB levels |
| `figure4_data_circular_phylogenetic_data_overlap_check.csv` | per-bin cross-sample mapping coverage (support for Figure 4) |
| `kegg_xenobiotic_pathways.csv` | KofamScan + KEGG Mapper, gene count per xenobiotic-degradation map |
| `eggnog_COG_category_counts.tsv` | COG category counts derived from the eggNOG-mapper output |
| `cazyme_classes.csv` | dbCAN2, pooled non-redundant catalogue, 6 CAZyme classes |
| `subclass_CAZyme.xlsx` | dbCAN2, pooled non-redundant catalogue, all CAZyme subclasses |
| `heatmap_data.tsv` | dbCAN2, gene count per MAG and CAZyme class (49 MAGs) |
| `taxo.xlsx` | GTDB-Tk (release r226) classification of the 49 MAGs |

**Not tracked by git:** `data/eggnog_results.emapper.annotations`, the raw
eggNOG-mapper output (117 MB, above GitHub's 100 MB file limit). The script uses
it when it is present and falls back on `data/eggnog_COG_category_counts.tsv`
otherwise; both routes give the same Figure S4. Request the raw file from the
corresponding author if you need to recompute the counts from scratch.

## Related records

- Raw reads: NCBI BioProject **PRJNA1444035** (SRA SRS28551350–SRS28551354),
  Zenodo mirror **10.5281/zenodo.19146671**.
- 49 MAGs: NCBI BioSamples **SAMN61982482–SAMN61982530**; assemblies at
  DDBJ/ENA/GenBank under **JCBOTC000000000–JCBOUY000000000**.
- Reference database, supplementary tables and structures: supplementary data
  archive of the article.

## Notes on two analytical choices

- **No inferential test.** With three samples at the polluted site and two at
  the reference site, the between-site comparison is descriptive: a log2 fold
  change of the mean relative abundances with a 0.001 pseudo-count. This is what
  the manuscript reports, and no *p*-value is derived from it.
- **Row set of Figure 2.** The heatmap keeps every row carrying a species name,
  that is both the species-level and the SGB-level entries, which are then
  summed; each taxon is therefore counted twice. The row z-score removes that
  constant factor, so the colours are unaffected — the only consequence is that
  the prevalence filter acts at 0.025 % rather than 0.05 %. The behaviour is
  kept as published so that the figure is reproduced exactly, and it is
  commented as such in the script.

## Licence

Code released under the MIT licence. Data files are covered by the licence of
the corresponding deposit (CC-BY-4.0).
