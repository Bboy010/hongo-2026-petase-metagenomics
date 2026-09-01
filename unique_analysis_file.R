################################################################################
##  Metagenomic survey of plastic-polluted (Bietry) and reference (Kassembie)
##  sites - Abidjan lagoon, Cote d'Ivoire.
##
##  Manuscript: ACMI-D-26-00103 (Access Microbiology)
##
##  This single script reproduces every statistical analysis and every
##  R-generated figure of the manuscript and of its supplementary material.
##
##  HOW TO RUN
##  ----------
##      setwd("<repository root>")         # working directory = this folder
##      source("unique_analysis_file.R")
##  or, from a terminal:
##      Rscript unique_analysis_file.R
##
##  Inputs  : ./data/      (see the table below)
##  Figures : ./Pictures/  (created if absent)
##  Tables  : ./results/   (created if absent)
##
##  FIGURE MAP - which manuscript figure is produced where
##  ------------------------------------------------------
##   Figure 1  Pipeline schematic ............... not R (Mermaid/vector drawing)
##   Figure 2  Differential taxa heatmap ........ section 2.4   <- this script
##   Figure 3  Top-20 species, relative abund. .. section 2.2   <- this script
##   Figure 4  Circular phylogenomic tree ....... iTOL v7; the cross-sample
##             detection check behind it is reproduced in section 3
##   Figure 5  KEGG xenobiotic pathways ......... section 4.1   <- this script
##   Figure 6  CAZyme classes ................... section 5.1   <- this script
##   Figure 7  MAG x CAZyme heatmap ............. section 5.2   <- this script
##   Figure 8  ColabFold structure / pLDDT ...... not R (ColabFold output)
##   Figure 9  Multiple sequence alignment ...... not R (ESPript)
##   Figure S1 Krona charts ..................... not R (Krona HTML export)
##   Figure S2 Top-10 SGB, relative abundance ... section 2.3   <- this script
##   Figure S3 Volcano plot (log2 fold change) .. section 2.5   <- this script
##   Figure S4 eggNOG COG categories ............ section 4.2   <- this script
##   Figure S5 EC classes ....................... section 4.3   <- this script
##   Figure S6 Top-100 CAZyme subclasses ........ section 5.3   <- this script
##
##  INPUT FILES
##  -----------
##   data/merged_abundance_table_species.csv ....... MetaPhlAn 4 merged profile
##                                                   (";"-separated, 5 samples)
##   data/figure4_data_circular_phylogenetic_data_overlap_check.csv
##                                                   per-bin cross-sample
##                                                   mapping coverage
##   data/kegg_xenobiotic_pathways.csv ............. KofamScan + KEGG Mapper
##                                                   gene counts per map
##   data/eggnog_results.emapper.annotations ....... eggNOG-mapper raw output
##                                                   (117 MB - NOT in the repo,
##                                                   see README.md; optional)
##   data/eggnog_COG_category_counts.tsv ........... COG counts derived from the
##                                                   file above (used as input
##                                                   when it is absent)
##   data/cazyme_classes.csv ....................... dbCAN2, pooled catalogue,
##                                                   6 CAZyme classes
##   data/subclass_CAZyme.xlsx ..................... dbCAN2, pooled catalogue,
##                                                   all CAZyme subclasses
##   data/heatmap_data.tsv ......................... dbCAN2 counts per MAG and
##                                                   CAZyme class (49 MAGs)
##   data/taxo.xlsx ................................ GTDB-Tk r226 classification
##                                                   of the 49 MAGs
##
##  Tested with R 4.5.2.
################################################################################


## =============================================================================
## 0. SET-UP
## =============================================================================

## ---- 0.1 Packages -----------------------------------------------------------
## Install once, if needed:
##   install.packages(c("readr", "dplyr", "tidyr", "tibble", "stringr",
##                      "ggplot2", "ggrepel", "RColorBrewer", "pheatmap",
##                      "readxl"))

required_packages <- c("readr", "dplyr", "tidyr", "tibble", "stringr",
                       "ggplot2", "ggrepel", "RColorBrewer", "pheatmap",
                       "readxl")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
       "\nInstall them with install.packages(c(\"",
       paste(missing_packages, collapse = "\", \""), "\"))", call. = FALSE)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(pheatmap)
  library(readxl)
})

## ---- 0.2 Paths --------------------------------------------------------------
data_dir    <- "data"
figure_dir  <- "Pictures"
results_dir <- "results"

if (!dir.exists(data_dir)) {
  stop("Directory '", data_dir, "' not found. Set the working directory to the ",
       "folder that contains this script before running it.", call. = FALSE)
}
dir.create(figure_dir,  showWarnings = FALSE)
dir.create(results_dir, showWarnings = FALSE)

## ---- 0.3 Sample metadata ----------------------------------------------------
## Bietry  = plastic-polluted site;  Kassembie = reference ("control") site.
## Matrix wording follows the manuscript: "Sediment" (not "Soil") and "Water".
sample_ids <- c("Bietry_1", "Bietry_2", "Bietry_3", "Kassembie_1", "Kassembie_2")

metadata <- data.frame(
  sample = sample_ids,
  site   = c("Bietry", "Bietry", "Bietry", "Kassembie", "Kassembie"),
  status = c("Polluted", "Polluted", "Polluted", "Control", "Control"),
  matrix = c("Sediment", "Sediment", "Water", "Sediment", "Water"),
  stringsAsFactors = FALSE
)

## ---- 0.4 Shared graphical settings ------------------------------------------
theme_set(theme_minimal())

## Colours reused across figures.
col_polluted <- "#D55E00"
col_control  <- "#009E73"
col_sediment <- "#8B5E3C"
col_water    <- "#56B4E9"
col_highlight <- "#e74c3c"   # class highlighted in the EC / CAZyme bar charts
col_neutral   <- "#95a5a6"
col_blue      <- "#3498db"

## ggrepel places labels by simulated annealing: fix the seed so that the
## volcano plot is byte-reproducible from one run to the next.
set.seed(20260731)

## ---- 0.5 Small helpers ------------------------------------------------------

## The MetaPhlAn table was exported from a spreadsheet under a French locale:
## about thirty cells carry a comma decimal mark ("5,00E-05"). Reading the
## sample columns as text and converting them here avoids the silent NAs that
## read_delim() would otherwise produce.
as_abundance <- function(x) {
  x <- sub(",", ".", x, fixed = TRUE)
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x
}

## Drop a MetaPhlAn rank prefix ("s__Psychrobacter_celer" -> "Psychrobacter_celer").
strip_rank <- function(x) str_remove(x, "^[kpcofgst]__")

## Replace the single underscores of a taxon name by spaces, leaving the "x__"
## rank prefixes untouched ("s__Psychrobacter_celer" -> "s__Psychrobacter celer").
clean_taxon <- function(x) str_replace_all(x, "(?<=[^_])_(?=[^_])", " ")

## The MetaPhlAn 4 database still ships a few pre-2021 phylum names. The
## manuscript uses the GTDB nomenclature throughout, so they are renamed here
## (this is what Figure 2 shows: Bacillota, not Firmicutes).
gtdb_phylum <- c(Firmicutes     = "Bacillota",
                 Proteobacteria = "Pseudomonadota",
                 Actinobacteria = "Actinomycetota",
                 Bacteroidetes  = "Bacteroidota")

to_gtdb_phylum <- function(x) if_else(x %in% names(gtdb_phylum),
                                      unname(gtdb_phylum[x]), x)

## Site names as they are spelled on the figures. They are built from a Unicode
## code point so that this script stays pure ASCII and renders identically
## whatever the locale of the machine running it.
e_acute <- intToUtf8(0x00E9)                        # "e" with an acute accent
site_polluted_label <- paste0("Bi", e_acute, "try")      # Bietry, polluted
site_control_label  <- paste0("Kassembi", e_acute)       # Kassembie, reference

message("Working directory : ", normalizePath("."))
message("Figures will be written to: ", file.path(normalizePath("."), figure_dir))


## =============================================================================
## 1. READ-BASED TAXONOMIC PROFILE (MetaPhlAn 4) - LOADING AND QUALITY CONTROL
## =============================================================================

message("\n[1] Loading the MetaPhlAn 4 profile ...")

metaphlan <- read_delim(
  file.path(data_dir, "merged_abundance_table_species.csv"),
  delim          = ";",
  escape_double  = FALSE,
  trim_ws        = TRUE,
  col_types      = cols(.default = col_character()),
  show_col_types = FALSE
) %>%
  mutate(across(all_of(sample_ids), as_abundance))

## The table stacks every rank: rows with a species (s__) but no SGB (t__) are
## the species level, rows with an SGB are the strain level. Each of the two
## sets sums to 100 % per sample; mixing them is what made the submitted
## Figure 3 add up to 200 % (reviewer comment M5).
species_level <- metaphlan %>%
  filter(!is.na(species), species != "", is.na(SGB) | SGB == "")

sgb_level <- metaphlan %>%
  filter(!is.na(SGB), SGB != "")

qc_totals <- tibble(
  sample        = sample_ids,
  species_total = vapply(sample_ids,
                         function(s) sum(species_level[[s]]), numeric(1)),
  sgb_total     = vapply(sample_ids,
                         function(s) sum(sgb_level[[s]]), numeric(1))
)

message("    species-level rows : ", nrow(species_level))
message("    SGB-level rows     : ", nrow(sgb_level))
message("    per-sample totals (must be 100 %):")
print(as.data.frame(qc_totals), row.names = FALSE)

stopifnot(all(abs(qc_totals$species_total - 100) < 0.01),
          all(abs(qc_totals$sgb_total     - 100) < 0.01))

write_csv(qc_totals, file.path(results_dir, "Table_QC_metaphlan_totals.csv"))


## -----------------------------------------------------------------------------
## 1.1 Taxonomic richness per rank, and share of undescribed SGBs
## -----------------------------------------------------------------------------

taxonomic_summary <- sgb_level %>%
  summarise(
    Classes  = n_distinct(class),
    Orders   = n_distinct(order),
    Families = n_distinct(family),
    Genera   = n_distinct(genus),
    Species  = n_distinct(species),
    SGB      = n_distinct(SGB)
  )

message("\n    Number of distinct taxa per rank (SGB-level rows):")
print(as.data.frame(taxonomic_summary), row.names = FALSE)
write_csv(taxonomic_summary,
          file.path(results_dir, "Table_taxonomic_richness.csv"))

## An SGB whose species name still carries an SGB/GGB code, or is flagged
## "unclassified", has no described representative in the reference database.
sgb_novelty <- sgb_level %>%
  mutate(status = if_else(grepl("SGB|GGB|unclassified", species),
                          "Undescribed", "Named species"))

novelty_summary <- sgb_novelty %>%
  group_by(status) %>%
  summarise(
    n_SGB                     = n(),
    mean_abundance_Bietry     = mean((Bietry_1 + Bietry_2 + Bietry_3) / 3),
    total_abundance_Bietry    = sum((Bietry_1 + Bietry_2 + Bietry_3) / 3),
    mean_abundance_Kassembie  = mean((Kassembie_1 + Kassembie_2) / 2),
    total_abundance_Kassembie = sum((Kassembie_1 + Kassembie_2) / 2),
    .groups = "drop"
  )

message("\n    Described vs undescribed SGBs:")
print(as.data.frame(novelty_summary), row.names = FALSE)

write_csv(sgb_novelty %>% select(kingdom:SGB, all_of(sample_ids), status),
          file.path(results_dir, "Table_SGB_novelty_full.csv"))
write_csv(novelty_summary,
          file.path(results_dir, "Table_SGB_novelty_summary.csv"))


## =============================================================================
## 2. FIGURES BUILT ON THE MetaPhlAn PROFILE
## =============================================================================

## -----------------------------------------------------------------------------
## 2.1 Differential abundance between sites (used by sections 2.4 and 2.5)
## -----------------------------------------------------------------------------
## Three replicates at the polluted site, two at the reference site, so the
## comparison is descriptive: a log2 fold change of the mean relative
## abundances, with a 0.001 pseudo-count to keep zeros finite. No inferential
## test is performed and none is reported in the manuscript (n = 1-3 per site).

add_fold_change <- function(df) {
  df %>%
    mutate(
      mean_polluted = (Bietry_1 + Bietry_2 + Bietry_3) / 3,
      mean_control  = (Kassembie_1 + Kassembie_2) / 2,
      LFC = log2((mean_polluted + 0.001) / (mean_control + 0.001))
    )
}

## -----------------------------------------------------------------------------
## 2.2 FIGURE 3 - relative abundance of the 20 most abundant species
## -----------------------------------------------------------------------------
message("\n[2.2] Figure 3 - relative abundance of the 20 most abundant species")

species_abundance <- species_level %>%
  mutate(taxon = str_replace_all(strip_rank(species), "_", " ")) %>%
  select(taxon, all_of(sample_ids)) %>%
  group_by(taxon) %>%
  summarise(across(all_of(sample_ids), sum), .groups = "drop") %>%
  mutate(total = rowSums(across(all_of(sample_ids))))

top20_species <- species_abundance %>%
  arrange(desc(total)) %>%
  slice_head(n = 20) %>%
  pull(taxon)

fig3_data <- species_abundance %>%
  mutate(taxon = if_else(taxon %in% top20_species, taxon, "Other")) %>%
  pivot_longer(all_of(sample_ids), names_to = "sample", values_to = "abundance") %>%
  group_by(sample, taxon) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  mutate(taxon  = factor(taxon,  levels = c(top20_species, "Other")),
         sample = factor(sample, levels = sample_ids))

fig3_palette <- c(colorRampPalette(brewer.pal(12, "Paired"))(20), "grey80")

fig3 <- ggplot(fig3_data, aes(x = sample, y = abundance, fill = taxon)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = fig3_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                     limits = c(0, 100.5)) +
  labs(x = "", y = "Relative abundance (%)", fill = "Species (SGB)") +
  theme_minimal() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.text  = element_text(size = 7, face = "italic"),
        legend.key.size = unit(0.35, "cm")) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(figure_dir, "Figure3_top20_species_relative_abundance.png"),
       fig3, width = 11, height = 7, dpi = 300)

message("    per-sample totals: ",
        paste(round(tapply(fig3_data$abundance, fig3_data$sample, sum), 1),
              collapse = " "))
write_csv(species_abundance %>% arrange(desc(total)),
          file.path(results_dir, "Table_species_relative_abundance.csv"))

## -----------------------------------------------------------------------------
## 2.3 FIGURE S2 - the most abundant UNCHARACTERIZED SGBs, five per site
## -----------------------------------------------------------------------------
## The caption, the Results and the Supplementary Materials list all describe
## this figure as "the ten most prevalent uncharacterized SGBs", so the ranking
## is restricted to the 131 SGBs that have no described representative in the
## MetaPhlAn 4 reference (section 1.1): 76 detected only at Bietry, 55 only at
## Kassembie. Ranking the 730 SGBs without that filter returns eight named
## species and contradicts the caption.
##
## The ten are taken as the five most abundant at each site rather than the ten
## most abundant overall: the uncharacterized SGBs are strictly site-exclusive
## and the six most abundant ones all occur in Kassembie_1, so a pooled ranking
## would collapse the figure into a single bar and hide the Bietry taxa.
message("[2.3] Figure S2 - most abundant uncharacterized SGBs, five per site")

sgb_clean <- sgb_level %>%
  mutate(uncharacterized = grepl("SGB|GGB|unclassified", species)) %>%
  mutate(across(kingdom:SGB, strip_rank))

uncharacterized_sgb <- sgb_clean %>%
  filter(uncharacterized) %>%
  mutate(mean_bietry    = (Bietry_1 + Bietry_2 + Bietry_3) / 3,
         mean_kassembie = (Kassembie_1 + Kassembie_2) / 2,
         taxon          = paste(genus, species, sep = " "))

message("    uncharacterized SGBs available for the ranking: ",
        nrow(uncharacterized_sgb))

top5_bietry <- uncharacterized_sgb %>%
  arrange(desc(mean_bietry)) %>% slice_head(n = 5)
top5_kassembie <- uncharacterized_sgb %>%
  arrange(desc(mean_kassembie)) %>% slice_head(n = 5)

## The two sets are disjoint because no uncharacterized SGB occurs at both sites.
stopifnot(length(intersect(top5_bietry$SGB, top5_kassembie$SGB)) == 0)

figS2_taxa <- c(top5_bietry$taxon, top5_kassembie$taxon)

figS2_data <- uncharacterized_sgb %>%
  filter(SGB %in% c(top5_bietry$SGB, top5_kassembie$SGB)) %>%
  pivot_longer(all_of(sample_ids), names_to = "sample", values_to = "abundance") %>%
  mutate(sample = factor(sample, levels = sample_ids),
         taxon  = factor(taxon,  levels = figS2_taxa))

figS2 <- ggplot(figS2_data, aes(x = sample, y = abundance, fill = taxon)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "", y = "Relative Abundance (%)", fill = "Species (SGB)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(figure_dir, "FigureS2_top10_uncharacterized_SGB.png"),
       figS2, width = 12, height = 8, dpi = 300)

message("    Bietry top 5   : ", paste(top5_bietry$SGB,    collapse = ", "))
message("    Kassembie top 5: ", paste(top5_kassembie$SGB, collapse = ", "))

write_csv(uncharacterized_sgb %>%
            select(phylum, family, genus, species, SGB,
                   all_of(sample_ids), mean_bietry, mean_kassembie) %>%
            arrange(desc(pmax(mean_bietry, mean_kassembie))),
          file.path(results_dir, "Table_uncharacterized_SGBs.csv"))

## -----------------------------------------------------------------------------
## 2.4 FIGURE 2 - heatmap of the taxa that differ most between the two sites
## -----------------------------------------------------------------------------
message("[2.4] Figure 2 - differential-abundance heatmap")

## NOTE ON THE ROW SET (kept identical to the published figure)
## Every row carrying a species name is retained here, i.e. the species-level
## row AND the SGB-level rows of that species, which are then summed. Because
## the SGBs of a species add up to that species, each taxon is counted twice.
## The row z-score removes that constant factor, so the colours are unaffected;
## the only consequence is on the prevalence filter below, which therefore acts
## at 0.025 % rather than 0.05 %. This is the exact row selection published as
## Figure 2 and is left unchanged for reproducibility.
differential <- metaphlan %>%
  filter(!is.na(species), species != "") %>%
  group_by(phylum, family, genus, species) %>%
  summarise(across(all_of(sample_ids), sum), .groups = "drop") %>%
  add_fold_change() %>%
  filter(mean_polluted + mean_control > 0.05)   # drop ultra-rare taxa

## 25 taxa most enriched at the polluted site + 25 most enriched at the reference
top_polluted <- differential %>% filter(LFC > 0) %>% arrange(desc(LFC)) %>%
  slice_head(n = 25)
top_control  <- differential %>% filter(LFC < 0) %>% arrange(LFC) %>%
  slice_head(n = 25)

selected <- bind_rows(top_polluted, top_control) %>%
  mutate(row_label = paste0(clean_taxon(genus), ";", clean_taxon(species)))

heatmap_matrix <- selected %>%
  select(row_label, all_of(sample_ids)) %>%
  column_to_rownames("row_label") %>%
  as.matrix()

## Row z-score: the colour encodes relative enrichment across samples, not the
## raw abundance (a taxon at 20 % and one at 0.2 % are put on the same scale).
heatmap_z <- t(scale(t(heatmap_matrix)))

## Column order matters: pheatmap stacks the annotation tracks from the bottom
## up, so "Site" first puts the Matrix track on top, as in the published figure.
column_annotation <- data.frame(
  Site   = metadata$status,
  Matrix = metadata$matrix,
  row.names = metadata$sample
)

row_annotation <- selected %>%
  mutate(Phylum = to_gtdb_phylum(strip_rank(phylum))) %>%
  select(row_label, Phylum) %>%
  column_to_rownames("row_label")

phylum_levels <- unique(row_annotation$Phylum)
phylum_colours <- setNames(
  colorRampPalette(brewer.pal(min(length(phylum_levels), 12), "Set3"))(
    length(phylum_levels)),
  phylum_levels
)

annotation_colours <- list(
  Site   = c(Polluted = col_polluted, Control = col_control),
  Matrix = c(Sediment = col_sediment, Water   = col_water),
  Phylum = phylum_colours
)

fig2_args <- list(
  mat    = heatmap_z,
  color  = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-3, 3, length.out = 101),

  cluster_rows             = TRUE,
  cluster_cols             = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D2",

  annotation_col       = column_annotation,
  annotation_row       = row_annotation,
  annotation_colors    = annotation_colours,
  annotation_names_row = TRUE,
  annotation_names_col = TRUE,

  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row  = 7,
  fontsize_col  = 9,
  angle_col     = 45,

  border_color = NA,
  cellwidth    = 38,
  cellheight   = 10,

  legend_breaks = c(-3, -2, -1, 0, 1, 2, 3),
  legend_labels = c("-3", "-2", "-1", "0", "+1", "+2", "+3 SD"),

  main = NA   # the published figure carries its title in the caption only
)

## Written twice: PNG for the manuscript, PDF for a vector version.
invisible(do.call(pheatmap, c(fig2_args, list(
  filename = file.path(figure_dir, "Figure2_differential_taxa_heatmap.png"),
  width = 10, height = 10))))
invisible(do.call(pheatmap, c(fig2_args, list(
  filename = file.path(figure_dir, "Figure2_differential_taxa_heatmap.pdf"),
  width = 10, height = 10))))

write_csv(selected %>% select(phylum, family, genus, species,
                              all_of(sample_ids),
                              mean_polluted, mean_control, LFC),
          file.path(results_dir, "Table_Figure2_selected_taxa.csv"))

## -----------------------------------------------------------------------------
## 2.5 FIGURE S3 - volcano plot of the SGB log2 fold change
## -----------------------------------------------------------------------------
message("[2.5] Figure S3 - volcano plot (log2 fold change)")

volcano_data <- sgb_level %>%
  add_fold_change() %>%
  mutate(genus_label = str_replace_all(strip_rank(genus), "_", " "))

## Candidate plastizyme-associated SGBs: abundant at the polluted site and at
## least four times more abundant there than at the reference site.
candidates <- volcano_data %>%
  filter(mean_polluted > 0.1, LFC > 2) %>%
  arrange(desc(LFC))

figS3 <- ggplot(volcano_data, aes(x = LFC, y = mean_polluted)) +
  geom_point(aes(color = LFC > 2), alpha = 0.5, show.legend = FALSE) +
  scale_color_manual(values = c("darkgrey", "firebrick")) +
  geom_vline(xintercept = 2, linetype = 2, colour = "grey40") +
  geom_text_repel(data = head(candidates, 10), aes(label = genus_label),
                  size = 3, max.overlaps = 15) +
  labs(x = paste0("Log2 fold change (", site_polluted_label, " vs ",
                  site_control_label, ")"),
       y = paste0("Mean relative abundance, ", site_polluted_label, " (%)")) +
  theme_minimal() +
  theme(text = element_text(size = 12))

ggsave(file.path(figure_dir, "FigureS3_volcano_LFC.png"),
       figS3, width = 10, height = 6, dpi = 300)

message("    candidate SGBs (mean_polluted > 0.1 and LFC > 2): ", nrow(candidates),
        " | LFC range: ", paste(round(range(candidates$LFC), 2), collapse = " to "))

write_csv(candidates %>% select(phylum, family, genus, species, SGB,
                                mean_polluted, mean_control, LFC),
          file.path(results_dir, "Table_candidate_SGBs_LFC2.csv"))


## =============================================================================
## 3. MAG CROSS-SAMPLE DETECTION - the check behind FIGURE 4
## =============================================================================
## Figure 4 itself (circular phylogenomic tree of the MAGs assigned to phyla
## reported to degrade plastics) was drawn in iTOL v7 from the GTDB-Tk tree.
## What is reproduced here is the check that supports its reading: are the bins
## recovered at one site also detected at the other?
##
## Coverage was computed within two mapping batches, because the two matrices
## were co-assembled separately: sediment (Bietry_1, Bietry_2, Kassembie_1) and
## water (Bietry_3, Kassembie_2). An empty cell means "outside the batch", a
## zero means "in the batch but not detected".

message("\n[3] MAG cross-sample detection check (support for Figure 4)")

bin_coverage <- read_csv(
  file.path(data_dir,
            "figure4_data_circular_phylogenetic_data_overlap_check.csv"),
  show_col_types = FALSE
)

bin_detection <- bin_coverage %>%
  mutate(bin  = str_remove(bin, "\\.fa$"),
         site = if_else(str_starts(bin, "Bietry"), "Bietry", "Kassembie")) %>%
  pivot_longer(all_of(sample_ids), names_to = "sample", values_to = "coverage") %>%
  filter(!is.na(coverage)) %>%          # drop the samples outside the batch
  left_join(metadata %>% select(sample, sample_site = site), by = "sample") %>%
  group_by(bin, site) %>%
  summarise(
    samples_tested   = n(),
    detected_own     = sum(coverage > 0 & sample_site == first(site)),
    detected_other   = sum(coverage > 0 & sample_site != first(site)),
    tested_other     = sum(sample_site != first(site)),
    .groups = "drop"
  ) %>%
  mutate(site_exclusive = tested_other > 0 & detected_other == 0)

detection_summary <- bin_detection %>%
  group_by(site) %>%
  summarise(
    bins                 = n(),
    tested_at_other_site = sum(tested_other > 0),
    also_detected_there  = sum(detected_other > 0),
    site_exclusive       = sum(site_exclusive),
    .groups = "drop"
  )

message("    bins analysed: ", nrow(bin_detection))
print(as.data.frame(detection_summary), row.names = FALSE)

write_csv(bin_detection,
          file.path(results_dir, "Table_MAG_cross_sample_detection.csv"))
write_csv(detection_summary,
          file.path(results_dir, "Table_MAG_cross_sample_summary.csv"))


## =============================================================================
## 4. FUNCTIONAL ANNOTATION OF THE NON-REDUNDANT GENE CATALOGUE
## =============================================================================

## -----------------------------------------------------------------------------
## 4.1 FIGURE 5 - KEGG xenobiotic biodegradation and metabolism pathways
## -----------------------------------------------------------------------------
## KofamScan assigns the KEGG Orthology identifiers, KEGG Mapper reconstructs
## the pathways. Following Supplementary Table S4, the three cytochrome-P450
## maps (00980, 00982, 00983) describe one conceptual pathway - the initial
## oxidation of polyethylene alkane chains - and are pooled: 21 KEGG maps
## therefore become 19 pathways, for a total of 381 genes.

message("\n[4.1] Figure 5 - KEGG xenobiotic pathways")

kegg <- read_csv(file.path(data_dir, "kegg_xenobiotic_pathways.csv"),
                 show_col_types = FALSE) %>%
  mutate(id = sprintf("%05d", as.integer(id)))

p450_maps <- c("00980", "00982", "00983")

kegg_grouped <- bind_rows(
  kegg %>% filter(!id %in% p450_maps),
  tibble(id      = "00980",
         pathway = "Metab. of xenobiotics by cyt. P450 (980/982/983)",
         genes   = sum(kegg$genes[kegg$id %in% p450_maps]))
) %>%
  mutate(label = if_else(id == "00980", pathway,
                         paste0(pathway, " (", id, ")"))) %>%
  arrange(genes) %>%
  mutate(label = factor(label, levels = label))

fig5 <- ggplot(kegg_grouped, aes(x = genes, y = label)) +
  geom_col(fill = col_blue, width = 0.72) +
  geom_text(aes(label = genes), hjust = -0.3, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Number of genes", y = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.text.y = element_text(size = 8))

ggsave(file.path(figure_dir, "Figure5_KEGG_xenobiotic_pathways.png"),
       fig5, width = 9, height = 7, dpi = 300)

message("    pathways: ", nrow(kegg_grouped),
        " | total genes: ", sum(kegg_grouped$genes))
stopifnot(sum(kegg_grouped$genes) == 381)

write_csv(kegg_grouped %>% select(id, pathway, genes),
          file.path(results_dir, "Table_KEGG_xenobiotic_pathways.csv"))

## -----------------------------------------------------------------------------
## 4.2 FIGURE S4 - eggNOG-mapper COG functional categories
## -----------------------------------------------------------------------------
## A protein may carry several one-letter COG categories; the string is split so
## that each assignment is counted once. The raw eggNOG-mapper output is 117 MB
## and is therefore distributed through the data repository rather than through
## git (see data/README.md). When it is absent, the pre-computed category counts
## shipped with the repository are used instead: both routes give the same plot.

message("[4.2] Figure S4 - eggNOG COG categories")

eggnog_file <- file.path(data_dir, "eggnog_results.emapper.annotations")
cog_counts_file <- file.path(data_dir, "eggnog_COG_category_counts.tsv")

if (file.exists(eggnog_file)) {

  message("    reading the raw eggNOG-mapper annotations (this takes a minute) ...")
  eggnog <- read_delim(
    eggnog_file,
    delim     = "\t",
    comment   = "#",
    col_names = c("query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
                  "max_annot_lvl", "COG_category", "Description",
                  "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway",
                  "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
                  "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"),
    col_types      = cols(.default = col_character()),
    show_col_types = FALSE
  )
  message("    annotated proteins: ", nrow(eggnog))

  cog_counts <- eggnog %>%
    filter(!is.na(COG_category), COG_category != "-") %>%
    mutate(COG_category = strsplit(COG_category, "")) %>%
    unnest(COG_category) %>%
    count(COG_category, name = "proteins", sort = TRUE)

  ## Refresh the small table that stands in for the 117 MB file.
  write_tsv(cog_counts, cog_counts_file)

} else {

  message("    raw annotations not found, using ", cog_counts_file)
  cog_counts <- read_tsv(cog_counts_file, show_col_types = FALSE) %>%
    arrange(desc(proteins))

}

figS4_data <- cog_counts %>% slice_head(n = 20)

figS4 <- ggplot(figS4_data,
                aes(x = reorder(COG_category, -proteins), y = proteins,
                    fill = proteins)) +
  geom_col(show.legend = FALSE) +
  scale_fill_viridis_c(option = "plasma") +
  labs(x = "COG Category", y = "Number of proteins") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

ggsave(file.path(figure_dir, "FigureS4_COG_categories.png"),
       figS4, width = 10, height = 6, dpi = 300)

message("    COG categories shown: ", nrow(figS4_data),
        " | assignments counted: ", sum(cog_counts$proteins))

write_csv(cog_counts, file.path(results_dir, "Table_COG_categories.csv"))

## -----------------------------------------------------------------------------
## 4.3 FIGURE S5 - Enzyme Commission classes
## -----------------------------------------------------------------------------
## Counts of the six main EC classes in the pooled non-redundant catalogue, as
## reported in the manuscript. Hydrolases (EC 3) are highlighted because they
## contain the ester-bond-cleaving activities relevant to polyester hydrolysis;
## oxidoreductases (EC 1) are highlighted as the class initiating polyolefin
## oxidation. EC 7 (translocases) is not shown, following the manuscript.

message("[4.3] Figure S5 - EC classes")

ec_classes <- data.frame(
  enzyme_class = c("Oxidoreductases (EC 1)", "Transferases (EC 2)",
                   "Hydrolases (EC 3)", "Lyases (EC 4)",
                   "Isomerases (EC 5)", "Ligases (EC 6)"),
  genes = c(13085, 23735, 15860, 6909, 5151, 5614),
  stringsAsFactors = FALSE
)

figS5 <- ggplot(ec_classes,
                aes(x = reorder(enzyme_class, genes), y = genes,
                    fill = enzyme_class)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Hydrolases (EC 3)"      = col_highlight,
    "Oxidoreductases (EC 1)" = col_blue,
    "Transferases (EC 2)"    = col_neutral,
    "Lyases (EC 4)"          = col_neutral,
    "Isomerases (EC 5)"      = col_neutral,
    "Ligases (EC 6)"         = col_neutral)) +
  labs(x = "", y = "Number of genes") +
  theme_minimal() +
  theme(text = element_text(size = 12))

ggsave(file.path(figure_dir, "FigureS5_EC_classes.png"),
       figS5, width = 8, height = 5, dpi = 300)

message("    total genes in EC 1-6: ", sum(ec_classes$genes))
write_csv(ec_classes, file.path(results_dir, "Table_EC_classes.csv"))


## =============================================================================
## 5. CARBOHYDRATE-ACTIVE ENZYMES (dbCAN2)
## =============================================================================

cazyme_classes_order <- c("AA", "CBM", "CE", "GH", "GT", "PL")

## -----------------------------------------------------------------------------
## 5.1 FIGURE 6 - CAZyme classes in the pooled gene catalogue
## -----------------------------------------------------------------------------
message("\n[5.1] Figure 6 - CAZyme classes")

## The percentages printed on the bars are recomputed from the gene counts
## rather than read from the pct column of the file, so that the figure can
## never drift from the data. The values are identical to the stored ones
## (43.4, 37.5, 7.7, 6.6, 4.1, 0.8). Each is correctly rounded; because five of
## the six round up, they sum to 100.1 rather than 100.0.
cazyme_classes <- read_csv(file.path(data_dir, "cazyme_classes.csv"),
                           show_col_types = FALSE) %>%
  mutate(pct   = round(genes / sum(genes) * 100, 1),
         label = paste0(class, " (", name, ")")) %>%
  arrange(genes) %>%
  mutate(label     = factor(label, levels = label),
         highlight = if_else(class == "CE", "CE", "other"))

fig6 <- ggplot(cazyme_classes, aes(x = genes, y = label, fill = highlight)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = paste0(genes, " (", pct, "%)")), hjust = -0.08,
            size = 3) +
  scale_fill_manual(values = c(CE = col_highlight, other = col_neutral),
                    guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = "Number of genes", y = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.text.y = element_text(size = 9))

ggsave(file.path(figure_dir, "Figure6_CAZyme_classes.png"),
       fig6, width = 9, height = 5, dpi = 300)

message("    total CAZyme genes in the pooled catalogue: ",
        sum(cazyme_classes$genes))

## -----------------------------------------------------------------------------
## 5.2 FIGURE 7 - CAZyme classes across the 49 dereplicated MAGs
## -----------------------------------------------------------------------------
message("[5.2] Figure 7 - MAG x CAZyme-class heatmap")

cazyme_by_mag <- read_tsv(file.path(data_dir, "heatmap_data.tsv"),
                          col_names = c("Bin", "Class", "Count"),
                          show_col_types = FALSE)

## GTDB-Tk classification of the same MAGs. The genus is read from the "g__"
## field; MAGs whose genus field is empty (g__;) keep the label "Unknown"
## rather than falling back on a higher rank.
mag_taxonomy <- read_excel(file.path(data_dir, "taxo.xlsx")) %>%
  mutate(
    Genus   = if_else(grepl("g__[^;]+", classification),
                      sub(".*g__([^;]+).*", "\\1", classification),
                      NA_character_),
    Species = if_else(grepl("s__[^;]+", classification),
                      sub(".*s__([^;]+).*", "\\1", classification),
                      NA_character_),
    Label   = paste0(user_genome, " (", coalesce(Genus, "Unknown"), ")")
  )

fig7_data <- cazyme_by_mag %>%
  left_join(mag_taxonomy, by = c("Bin" = "user_genome")) %>%
  mutate(Class = factor(Class, levels = cazyme_classes_order))

fig7 <- ggplot(fig7_data, aes(x = Class, y = reorder(Label, Count),
                              fill = Count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  theme_minimal() +
  ## Parentheses in the row labels hold the GTDB-Tk genus (46 of 49 MAGs;
  ## "Unknown" for the 3 MAGs with no genus assignment), hence "Genera".
  labs(x = "CAZyme Class", y = "MAG(Genera)")

ggsave(file.path(figure_dir, "Figure7_MAG_CAZyme_heatmap.png"),
       fig7, width = 12.5, height = 9.6, dpi = 300)
ggsave(file.path(figure_dir, "Figure7_MAG_CAZyme_heatmap.pdf"),
       fig7, width = 12.5, height = 9.6)

message("    MAGs: ", n_distinct(fig7_data$Bin),
        " | CAZyme genes in MAGs: ", sum(fig7_data$Count))

write_csv(fig7_data %>% select(Bin, Label, Genus, Species, Class, Count),
          file.path(results_dir, "Table_CAZyme_by_MAG.csv"))

## -----------------------------------------------------------------------------
## 5.3 FIGURE S6 - the 100 most abundant CAZyme subclasses
## -----------------------------------------------------------------------------
## Counts are the union of the three dbCAN2 tools (HMMER, eCAMI, DIAMOND)
## deduplicated per gene, which is the method used in the manuscript
## (GT2 = 1650, GT4 = 1415, GH23 = 406, CBM50 = 401).
## Bars, not a connected line: subclasses are categories, not a series
## (reviewer comment R1). The y axis states what is actually plotted, the
## number of genes (reviewer comment R2).

message("[5.3] Figure S6 - top-100 CAZyme subclasses")

## The dbCAN2 export carries five rows whose label is an EC number rather than
## a CAZy family (2.4.1.12, 2.4.1.-, 3.2.1.23, 3.2.1.86, 3.2.1.133; 6 genes in
## total). They are annotations that leaked out of the "family|EC" strings of
## the dbCAN2 overview file, not CAZyme families, so they are dropped here.
## They all fall far outside the top 100 (ranks 339 and 420-423, 1-2 genes
## each), so Figure S6 is unchanged; only the family count is affected: the
## catalogue holds 418 CAZy families and subfamilies, not the 423 rows of the
## spreadsheet.
cazyme_subclasses_all <- read_excel(file.path(data_dir, "subclass_CAZyme.xlsx")) %>%
  rename(subclass = Sub_CAZymes, genes = number) %>%
  filter(!str_detect(subclass, "^[0-9]+[.][0-9]+[.][0-9]+[.]"))

message("    CAZy families and subfamilies: ", nrow(cazyme_subclasses_all),
        " | genes: ", sum(cazyme_subclasses_all$genes))

cazyme_subclasses <- cazyme_subclasses_all %>%
  arrange(desc(genes)) %>%
  slice_head(n = 100) %>%
  mutate(class = str_extract(subclass, "^[A-Z]+")) %>%
  mutate(class = factor(class, levels = cazyme_classes_order)) %>%
  arrange(class, subclass) %>%
  mutate(subclass = factor(subclass, levels = subclass))

stopifnot(!any(is.na(cazyme_subclasses$class)))

cazyme_palette <- c(AA = "#4E79A7", CBM = "#F28E2B", CE = "#59A14F",
                    GH = "#E15759", GT = "#B07AA1", PL = "#9C755F")

figS6 <- ggplot(cazyme_subclasses, aes(x = subclass, y = genes, fill = class)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = cazyme_palette, name = "CAZyme class") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(x = "CAZyme subclass", y = "Number of genes") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                      size = 6.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "top"
  )

ggsave(file.path(figure_dir, "FigureS6_CAZyme_subclasses_top100.png"),
       figS6, width = 13, height = 6, dpi = 300)

message("    subclasses: ", nrow(cazyme_subclasses),
        " | most abundant: ", max(cazyme_subclasses$genes), " genes")
message("    control GT2/GT4/GH23/CBM50: ",
        paste(cazyme_subclasses$genes[match(c("GT2", "GT4", "GH23", "CBM50"),
                                            cazyme_subclasses$subclass)],
              collapse = " "))

write_csv(cazyme_subclasses,
          file.path(results_dir, "Table_CAZyme_subclasses_top100.csv"))


## =============================================================================
## 6. SESSION INFORMATION
## =============================================================================

message("\nAll figures written to ", figure_dir, "/ and all tables to ",
        results_dir, "/")

writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "sessionInfo.txt"))
print(sessionInfo())

################################################################################
## End of script
################################################################################
