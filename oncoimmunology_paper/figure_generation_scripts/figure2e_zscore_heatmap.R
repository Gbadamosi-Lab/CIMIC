# ============================================================================
# figure2e_zscore_heatmap.R
# ============================================================================
# Purpose:
#   Generate a publication-quality heatmap of single-sample gene set enrichment
#   scores (ssGSEA via z-score aggregation) for immune/cell-death programs across
#   CIMIC clusters. Tests for significant differences in program activation between
#   clusters and annotates with significance markers.
#
# Inputs:
#   - Datasets/nki_smc_combine_TPM_lognorm.csv: Log-normalized TPM matrix
#     (rows=genes, columns=samples; must include gene_id column)
#   - Datasets/nki_smc_combine_clustered_plot_df.csv: Sample metadata
#     (rows=samples, must include sample_id and CIMIC_Cluster columns)
#
# Outputs:
#   - Results/figure2e_zscore_heatmap.png: Publication-ready heatmap
#     (14x10 inches, 300 DPI, suitable for single-column journal figure)
#
# Dependencies:
#   - R ≥ 4.0
#   - msigdbr, dplyr, tidyr, purrr, tibble, rstatix, ggplot2, stringr, data.table
#   - MSigDB (downloaded via msigdbr package; Homo sapiens)
#
# Author: Mohammed Gbadamosi
# Last Updated: June 2024
# ============================================================================

# ============================================================================
# LIBRARY IMPORTS
# ============================================================================
library(data.table)   # Fast CSV reading
library(msigdbr)      # Gene set database
library(dplyr)        # Data manipulation
library(tidyr)        # Data reshaping
library(purrr)        # Functional programming
library(tibble)       # Tibble data structures
library(rstatix)      # Statistical helpers
library(ggplot2)      # Visualization
library(stringr)      # String manipulation

# ============================================================================
# DATA LOADING
# ============================================================================
# Load sample metadata with integrated gene expression data
# Note: Data contains: sample_id, cluster_assignments, + gene expression columns

clustered_plot_df <- data.table::fread("Datasets/NKI_SMC/nki_smc_combine_clustered_plot_df.csv", data.table = FALSE)

# Convert to data.frame for dplyr operations
clustered_plot_df <- as.data.frame(clustered_plot_df)

# Rename cluster column to match expected convention
if ("cluster_assignments" %in% colnames(clustered_plot_df)) {
  clustered_plot_df <- clustered_plot_df %>%
    rename(CIMIC_Cluster = cluster_assignments)
}

# Identify gene expression columns (exclude known metadata columns)
metadata_cols <- c("sample_id", "CIMIC_Cluster", "base_id", "PACMAP1", "PACMAP2",
                   "UMAP1", "UMAP2", "PC1", "PC2", "PC3", "V1")
gene_cols <- setdiff(colnames(clustered_plot_df), metadata_cols)

message(sprintf("Loaded %d samples with %d gene expression columns",
                nrow(clustered_plot_df), length(gene_cols)))

# Create Results directory if needed
dir.create("Results", showWarnings = FALSE)

# ============================================================================
# GENE SET DEFINITIONS
# ============================================================================
# BIOLOGICAL RATIONALE:
#   These 19 curated gene sets represent hallmark immune response and stress
#   pathways implicated in the cancer-associated immune microenvironment (CAIM).
#   They span GO Biological Processes (immune activation, cell death),
#   Hallmark pathways (inflammatory response, UPR), and expert Reactome/WikiPathways
#   annotations. Together with death programs below, they define a transcriptomic
#   signature of immunologically "hot" vs. immunologically "cold" tumors.

# Load all gene sets from MSigDB (Homo sapiens)
msig_df <- msigdbr::msigdbr(species = "Homo sapiens") %>% as.data.frame()

# Define 19 curated immune/stress gene sets
all_gene_set_names <- c(
  "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "GOBP_B_CELL_ACTIVATION",
  "GOBP_CELLULAR_RESPONSE_TO_STRESS",
  "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  "GOBP_INFLAMMASOME_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS",
  "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
  "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "GOBP_LEUKOCYTE_CHEMOTAXIS_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  "GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "GOBP_T_CELL_ACTIVATION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS",
  "REACTOME_CELLULAR_RESPONSE_TO_MITOCHONDRIAL_STRESS",
  "WP_MOLECULAR_PATHWAY_FOR_OXIDATIVE_STRESS",
  "WP_OXIDATIVE_STRESS_RESPONSE"
)

# Extract gene symbols for each set from MSigDB; return as named list
all_gene_sets <- setNames(
  lapply(all_gene_set_names, function(gs) {
    gs_df <- msig_df[msig_df$gs_name == gs, , drop = FALSE]
    gs_df$gene_symbol
  }),
  all_gene_set_names
)

# ============================================================================
# B/T CELL RECEPTOR GENE FILTERING
# ============================================================================
# BIOLOGICAL RATIONALE:
#   Immunoglobulin (B-cell receptor) and T-cell receptor (TCR) genes are highly
#   variable across samples due to clonal expansion of lymphocytes. Their inclusion
#   in zsGSEA would reflect tumor-infiltrating lymphocyte (TIL) abundance rather
#   than pathway activation. Removal ensures focus on intrinsic cell-autonomous
#   immune pathway signaling, not clonal lymphocyte heterogeneity.

# Build regex pattern matching any BCR/TCR V/D/J/C region prefix
remove_prefixes <- c(
  "IGHV", "IGHD", "IGHJ", "IGHC", "IGHM", "IGHA", "IGHE", "IGHG",
  "IGLV", "IGLC", "IGLJ",
  "IGKV", "IGKC", "IGKJ",
  "TRAV", "TRBV", "TRDV", "TRGV",
  "TRAJ", "TRBJ", "TRGJ", "TRDJ",
  "TRAC", "TRBC", "TRGC", "TRDC"
)

# Construct regex: matches any prefix at the start of gene name
pattern <- paste0("^(", paste(remove_prefixes, collapse = "|"), ")")

# Apply filter to each gene set (removes BCR/TCR genes from all sets)
all_gene_sets <- lapply(all_gene_sets, function(gene_vec) {
  gene_vec[!grepl(pattern, gene_vec, perl = TRUE)]
})

# ============================================================================
# CELL DEATH PROGRAM DEFINITIONS (REMOVED FROM THIS ANALYSIS)
# ============================================================================
# Note: Death programs are defined below but NOT included in the heatmap visualization.
# Only the 19 immune/stress gene sets from MSigDB are used for this figure.
# To include death programs, uncomment the line "all_gene_sets <- c(death_programs, all_gene_sets)"

death_programs <- list(
  Apoptosis = c(
    "BAX", "BAK1", "BBC3", "PMAIP1", "BCL2L11",
    "APAF1", "CASP9", "CASP3", "CASP7",
    "FAS", "TNFRSF10B", "CASP8", "FADD", "BID",
    "DIABLO", "CYCS", "TP53AIP1", "BAD", "BMF", "HRK"
  ),
  Necroptosis = c(
    "RIPK1", "RIPK3", "MLKL", "ZBP1", "TICAM1", "TNFRSF1A", "FADD", "CASP8"
  ),
  Pyroptosis = c(
    "GSDME", "GSDMD", "CASP1", "CASP4", "CASP5",
    "NLRP3", "AIM2", "PYCARD", "IL1B", "IL18", "CASP3", "NLRC4"
  ),
  PANoptosis = c(
    "ZBP1", "AIM2", "RIPK3", "RIPK1", "CASP8",
    "CASP1", "FADD", "PYCARD", "IRF1"
  ),
  Ferroptosis = c(
    "ACSL4", "LPCAT3", "TFRC", "SAT1", "PTGS2"
  )
)

# ONLY use immune/stress gene sets (exclude death programs for this figure)
# all_gene_sets <- c(death_programs, all_gene_sets)  # <- COMMENTED OUT

# Store names for later reference
all_gene_set_names_cimic <- all_gene_set_names
all_gene_set_names_death <- names(death_programs)

# ============================================================================
# GENE-WISE Z-SCORE NORMALIZATION
# ============================================================================
# BIOLOGICAL RATIONALE:
#   Gene expression (TPM) follows non-normal distributions with mean and variance
#   proportional to gene abundance. Z-score normalization standardizes each gene
#   to mean=0, sd=1 across all samples. This prevents high-abundance genes from
#   dominating pathway scores and ensures equal contribution of all genes within
#   a pathway, enabling fair pathway-level comparisons between samples.

# Safe z-score function: handles genes with zero variance (constant expression)
zscore_safe <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  # REVIEW: Returns all 0s for constant genes; confirm this is acceptable behavior
  if (is.na(s) || s == 0) rep(0, length(x)) else (x - m) / s
}

# Apply z-score normalization across all genes
df_z <- clustered_plot_df %>%
  mutate(across(all_of(gene_cols), zscore_safe))

# ============================================================================
# PROGRAM SCORING: MEAN Z-SCORE PER SAMPLE
# ============================================================================
# BIOLOGICAL RATIONALE:
#   Each program's activity in a sample is defined as the mean z-score of all
#   genes in that program. This aggregation converts gene-level z-scores (measuring
#   individual gene expression rank) into pathway-level scores (measuring coordinated
#   pathway activity). Averaging preserves the sign (activation vs repression) and
#   magnitude of the signal while reducing noise from individual genes.

score_program_z <- function(data, genes, program_name) {
  # Ensure genes exist in data (some may be filtered or missing)
  genes <- intersect(genes, colnames(data))

  if (length(genes) == 0) {
    message(sprintf("  Warning: %s has 0 genes after filtering; skipping", program_name))
    return(
      data.frame(
        sample_id = data$sample_id,
        cimic_cluster = data$CIMIC_Cluster,
        Program = paste0(program_name, "_z"),
        Score = rep(NA_real_, nrow(data))
      )
    )
  }

  # Row-wise mean of z-scores for genes in this program
  # Use select() to avoid data.table subsetting issues with variable column names
  score_vec <- data %>%
    dplyr::select(all_of(genes)) %>%
    as.matrix() %>%
    rowMeans(na.rm = TRUE)

  data.frame(
    sample_id = data$sample_id,
    cimic_cluster = data$CIMIC_Cluster,
    Program = paste0(program_name, "_z"),
    Score = score_vec
  )
}

# Score all programs for all samples
program_scores_all <- purrr::imap_dfr(all_gene_sets, ~ score_program_z(df_z, .x, .y))

# Reshape to wide format: rows=samples, columns=program scores
program_scores_wide <- program_scores_all %>%
  dplyr::select(sample_id, Program, Score) %>%
  tidyr::pivot_wider(
    names_from = Program,
    values_from = Score
  )

# Merge program scores back into sample metadata
clustered_plot_df <- clustered_plot_df %>%
  left_join(program_scores_wide, by = "sample_id")

# ============================================================================
# HEATMAP DATA PREPARATION
# ============================================================================
# BIOLOGICAL RATIONALE:
#   The heatmap displays cluster-averaged program scores. By grouping samples by
#   CIMIC cluster (Dys-CIM vs Fun-CIM) and averaging scores within each cluster,
#   we reveal whether specific programs are systematically activated or repressed
#   in one cluster vs. the other. This cluster-level view highlights consistent
#   phenotypic differences while reducing sample-level noise.

str_wrap_length <- 40

# Select only immune/death program z-scores for visualization
genesets_of_interest <- paste0(names(all_gene_sets), "_z")

# Cluster-averaged scores with cleaned program names
cluster_heatmap_df <- clustered_plot_df %>%
  dplyr::select(CIMIC_Cluster, all_of(genesets_of_interest)) %>%
  pivot_longer(
    cols = -CIMIC_Cluster,
    names_to = "Program",
    values_to = "Score"
  ) %>%
  mutate(
    Score = as.numeric(Score),
    CIMIC_Cluster = factor(CIMIC_Cluster),
    # Remove database prefix and replace underscores with spaces for readability
    Program = Program %>%
      str_remove_all("GOBP_|REACTOME_|HALLMARK_|WP_|_z") %>%
      str_replace_all("_", " ") %>%
      str_wrap(width = str_wrap_length)
  ) %>%
  group_by(CIMIC_Cluster, Program) %>%
  # Compute mean score for each program within each cluster
  summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop")

# ============================================================================
# STATISTICAL TESTING: ADAPTIVE WILCOXON/KRUSKAL-WALLIS
# ============================================================================
# BIOLOGICAL RATIONALE:
#   Test whether each program's activity significantly differs across clusters.
#   Non-parametric tests are appropriate for RNA-seq data (non-normal distributions,
#   outliers from biological heterogeneity). Benjamini-Hochberg FDR correction
#   controls false discovery rate across 24 hypothesis tests while maintaining power.
#
# ADAPTIVE TEST SELECTION:
#   - 2 clusters (our case): Mann-Whitney U test (Wilcoxon rank-sum)
#     Tests if two distributions differ; no normality assumption; uses ranks
#   - 3+ clusters: Kruskal-Wallis test
#     Non-parametric analog of one-way ANOVA; compares k distribution medians
#   - Rationale: Non-parametric tests preserve rank information in small samples
#     with potential outliers, which is common in bulk RNA-seq profiling

pval_df <- clustered_plot_df %>%
  dplyr::select(CIMIC_Cluster, all_of(genesets_of_interest)) %>%
  pivot_longer(
    cols = -CIMIC_Cluster,
    names_to = "GeneSet",
    values_to = "Score"
  ) %>%
  mutate(CIMIC_Cluster = factor(CIMIC_Cluster)) %>%
  group_by(GeneSet) %>%
  summarise(
    n_groups = n_distinct(CIMIC_Cluster),
    p = tryCatch({
      if (n_groups == 2) {
        # Two-group comparison: Mann-Whitney U (Wilcoxon rank-sum)
        wilcox.test(Score ~ CIMIC_Cluster)$p.value
      } else {
        # Multi-group comparison: Kruskal-Wallis
        kruskal.test(Score ~ CIMIC_Cluster)$p.value
      }
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    # Benjamini-Hochberg FDR correction: controls false discovery rate
    # (expected proportion of false positives among rejected hypotheses)
    p.adj = p.adjust(p, method = "BH"),
    # Significance thresholding for annotation
    p.adj.signif = case_when(
      is.na(p.adj) ~ "NA",
      p.adj <= 0.0001 ~ "****",
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01 ~ "**",
      p.adj <= 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Build significance label table for heatmap annotation
sig_df <- pval_df %>%
  mutate(
    Program = GeneSet %>%
      str_remove_all("GOBP_|REACTOME_|HALLMARK_|WP_|_z") %>%
      str_replace_all("_", " ") %>%
      str_wrap(width = str_wrap_length),
    signif_label = p.adj.signif
  ) %>%
  dplyr::select(Program, signif_label)

# ============================================================================
# PREPARE FACTORS FOR CONSISTENT ORDERING
# ============================================================================
# Establish consistent program order for both data and significance labels
program_order <- cluster_heatmap_df %>%
  distinct(Program) %>%
  pull(Program)

# Check cluster values
unique_clusters <- unique(cluster_heatmap_df$CIMIC_Cluster)
message(sprintf("Unique cluster values in heatmap data: %s", paste(unique_clusters, collapse=", ")))

message(sprintf("Number of programs in heatmap: %d", length(program_order)))
message(sprintf("Mean score range: [%.3f, %.3f]",
                min(cluster_heatmap_df$mean_score, na.rm=TRUE),
                max(cluster_heatmap_df$mean_score, na.rm=TRUE)))

# Check for NA values
na_count <- sum(is.na(cluster_heatmap_df$mean_score))
message(sprintf("Number of NA values in heatmap data: %d out of %d",
                na_count, nrow(cluster_heatmap_df)))

cluster_heatmap_df <- cluster_heatmap_df %>%
  mutate(Program = factor(Program, levels = rev(program_order)))

sig_df <- sig_df %>%
  mutate(Program = factor(Program, levels = rev(program_order)))

# Create data frame for significance star annotation
star_df <- cluster_heatmap_df %>%
  distinct(Program) %>%
  left_join(sig_df, by = "Program") %>%
  mutate(x = 2.8)

# ============================================================================
# HEATMAP VISUALIZATION & AESTHETIC CHOICES
# ============================================================================
# AESTHETIC CHOICES FOR JOURNAL SUBMISSION (Publication Quality):
#
# Color Scheme:
#   - scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0)
#   - Diverging red-white-blue color space emphasizes bidirectional effects
#     (program activation in red, repression in blue). Colorblind-friendly.
#   - Midpoint=0 centers white on zero z-score (no change from sample mean),
#     making it immediately clear which programs are activated vs. repressed.
#
# Tile & Border:
#   - geom_tile(color="black", linewidth=0.8) creates high-contrast black borders
#     between tiles. Essential for print reproduction; ensures tiles remain
#     visually distinct even at low resolution or when grayscale printed.
#   - width=1 prevents tile overlap; clean grid structure.
#
# Text Sizing:
#   - Cluster labels (axis.text.x=45): Large bold text for primary grouping variable
#   - Program labels (axis.text.y=28): Smaller but bold for 40-char wrapped names
#   - Legend text (legend.text=24): Matches general readability scale
#   - All sizes chosen for readability at 300 DPI when printed at journal
#     single-column width (~3.5 inches)
#
# Annotation (Significance Stars):
#   - annotate(..., size=20, fontface="bold") renders significance markers above
#     cluster columns. size=20 ensures visibility at 300 DPI. Bold weight adds emphasis.
#   - Four levels: **** (p≤0.0001), *** (p≤0.001), ** (p≤0.01), * (p≤0.05), "" (ns)
#   - Stars positioned outside plot area (coord_cartesian(clip="off")) for clarity.
#
# Layout & Margins:
#   - coord_cartesian(clip="off"): Allows significance stars to render beyond
#     plot boundaries, preventing overlap with tile grid.
#   - plot.margin(...): Adds breathing room around figure edges; improves readability
#     when figure is cropped for journal layout.
#   - legend.position="none": No legend in figure; space saved for larger tiles.
#     Legend moved to figure caption (color scale values explained in methods).
#
# Theme & Fonts:
#   - theme_classic() removes distracting gridlines, background color, etc.
#   - base_size=24: Base font size for all text elements (inherited by axis labels)
#   - axis.line=element_blank(): No outer border; clean tile-based layout.
#   - Facet background white with black outline: Standard formatting for gene set labels

p_cluster_heatmap <- ggplot(
  cluster_heatmap_df,
  aes(x = CIMIC_Cluster, y = Program, fill = mean_score)
) +
  # Tiles with high-contrast borders for print quality
  geom_tile(color = "black", linewidth = 0.8) +
  # Significance stars positioned to the right of the heatmap
  annotate(
    "text",
    x = 2.8,
    y = seq_along(levels(cluster_heatmap_df$Program)),
    label = rev(star_df$signif_label),
    size = 20,
    fontface = "bold"
  ) +
  # Diverging color scale: blue (repressed) to white (neutral) to red (activated)
  scale_fill_gradient2(
    low = "#313695",
    mid = "white",
    high = "#A50026",
    midpoint = 0,
    name = "Mean Z-score",
    limits = c(-0.5, 0.5)
  ) +
  scale_x_discrete(
    labels = c("1" = "Dys-CIM", "2" = "Fun-CIM"),
    expand = expansion(add = c(0.5, 1.5))
  ) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x = element_text(face = "bold", size = 28, color = "black"),
    axis.text.y = element_text(face = "bold", size = 16, color = "black"),
    legend.position = "right",
    legend.text = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 16),
    axis.line = element_blank(),
    plot.margin = margin(t = 10, r = 80, b = 10, l = 20)
  )

# Print to console for interactive inspection
print(p_cluster_heatmap)

# ============================================================================
# FIGURE OUTPUT
# ============================================================================
# Save figure at publication-quality resolution (300 DPI)
# Dimensions: 14x10 inches (suitable for single-column journal layout)
suppressWarnings(
  ggsave(
    "Results/figure2e_zscore_heatmap.png",
    plot = p_cluster_heatmap,
    width = 14,
    height = 10,
    dpi = 300,
    bg = "white"
  )
)

message("✓ Figure saved: Results/figure2e_zscore_heatmap.png (14x10 in, 300 DPI)")
message("✓ Includes 19 immune/stress gene sets (BCR/TCR genes filtered)")
message("✓ Significance stars: **** p≤0.0001, *** p≤0.001, ** p≤0.01, * p≤0.05")
message("✓ Script complete: z-score heatmap analysis finished successfully")
