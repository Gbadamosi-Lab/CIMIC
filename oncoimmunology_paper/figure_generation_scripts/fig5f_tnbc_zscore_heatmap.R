# ============================================================================
# fig5f_tnbc_zscore_heatmap.R  
# ============================================================================
# Purpose:
#   Generate a publication-quality heatmap of single-sample gene set enrichment
#   scores (ssGSEA via z-score aggregation) for immune/cell-death programs across
#   CIMIC clusters in the TNBC_CL_Epirubicin dataset. Tests for significant
#   differences in program activation between clusters and annotates with
#   significance markers.
#
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
#
# Inputs:
#   - Datasets/TNBC_CL_Epirubicin/tnbc_cl_epi_clustered_plot_df.csv: Sample
#     metadata with gene expression columns and cluster assignments.
#
# Outputs:
#   - Results/fig5f_tnbc_zscore_heatmap.png: Publication-ready heatmap
#     (18x14 inches, 300 DPI)
#   - Results/fig5f_tnbc_zscore_heatmap_fry_camera_statistics.csv: fry +
#     block-adjusted cameraPR statistics plus the mean Z-score per cluster
#     shown on the heatmap tiles
#
# Dependencies:
#   - R >= 4.0
#   - msigdbr, dplyr, tidyr, purrr, tibble, rstatix, limma, ggplot2, stringr,
#     data.table
#
# Author: Mohammed Gbadamosi (adapted by GitHub Copilot)
# ============================================================================

# ============================================================================
# LIBRARY IMPORTS
# ============================================================================
library(data.table)
library(msigdbr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(rstatix)
library(limma)
library(ggplot2)
library(stringr)


# ============================================================================
# FRY + CAMERA RESULTS EXPORT
# ============================================================================
# Human cohorts:
#   - fry() for self-contained gene-set testing
#   - camera() for competitive gene-set testing
#
# Repeated cell-line cohorts:
#   - fry() with block + duplicateCorrelation
#   - blocked limma model followed by cameraPR() on moderated t-statistics
#
# Positive Direction / "Up" means enrichment in Dys-CIM because the contrast is:
#   Dys_CIM - Fun_CIM
# ============================================================================

make_fry_camera_table <- function(
    expr_matrix,
    gene_sets,
    design,
    contrast,
    cohort_type = c("human", "cell_line"),
    block = NULL,
    correlation = NULL,
    analysis_pool = "Main",
    standalone = FALSE
) {
  cohort_type <- match.arg(cohort_type)

  # --------------------------------------------------------------------------
  # Validate and filter gene sets
  # --------------------------------------------------------------------------
  if (is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
    stop("gene_sets must be a named list.")
  }

  gene_sets <- lapply(
    gene_sets,
    function(gs) intersect(unique(gs), rownames(expr_matrix))
  )

  gene_sets <- gene_sets[lengths(gene_sets) >= 5]

  if (length(gene_sets) == 0) {
    stop("No gene sets retained at the minimum threshold of 5 measured genes.")
  }

  # --------------------------------------------------------------------------
  # FRY: self-contained gene-set test
  # --------------------------------------------------------------------------
  if (cohort_type == "human") {

    fry_result <- limma::fry(
      y        = expr_matrix,
      index    = gene_sets,
      design   = design,
      contrast = contrast,
      sort     = "none"
    )

  } else {

    if (is.null(block) || is.null(correlation)) {
      stop(
        "For cohort_type = 'cell_line', both block and correlation are required."
      )
    }

    if (length(block) != ncol(expr_matrix)) {
      stop("Length of block must equal the number of expression-matrix samples.")
    }

    fry_result <- limma::fry(
      y           = expr_matrix,
      index       = gene_sets,
      design      = design,
      contrast    = contrast,
      block       = block,
      correlation = correlation,
      sort        = "none"
    )
  }

  fry_table <- data.frame(
    GeneSet = rownames(fry_result),
    fry_result,
    row.names = NULL,
    check.names = FALSE
  )

  names(fry_table)[-1] <- paste0("Fry_", names(fry_table)[-1])

  # --------------------------------------------------------------------------
  # CAMERA: competitive gene-set test
  # --------------------------------------------------------------------------
  if (cohort_type == "human") {

    # Direct CAMERA is appropriate for independent patient samples.
    camera_result <- limma::camera(
      y              = expr_matrix,
      index          = gene_sets,
      design         = design,
      contrast       = contrast,
      inter.gene.cor = 0.01,
      use.ranks      = FALSE,
      trend.var      = FALSE,
      sort           = FALSE
    )

    camera_implementation <- "camera"

  } else {

    # camera() does not accept block/correlation directly.
    # First fit the blocked limma model, then pass the resulting moderated
    # t-statistics to cameraPR().
    blocked_fit <- limma::lmFit(
      object      = expr_matrix,
      design      = design,
      block       = block,
      correlation = correlation
    )

    blocked_fit <- limma::contrasts.fit(
      fit       = blocked_fit,
      contrasts = contrast
    )

    blocked_fit <- limma::eBayes(
      blocked_fit,
      trend = FALSE
    )

    if (ncol(blocked_fit$t) != 1) {
      stop("The CAMERA export block expects exactly one contrast.")
    }

    moderated_t <- blocked_fit$t[, 1]
    names(moderated_t) <- rownames(expr_matrix)

    camera_result <- limma::cameraPR(
      statistic      = moderated_t,
      index          = gene_sets,
      use.ranks      = FALSE,
      inter.gene.cor = 0.01,
      sort           = FALSE
    )

    camera_implementation <- "cameraPR_block_adjusted"
  }

  camera_table <- data.frame(
    GeneSet = rownames(camera_result),
    camera_result,
    row.names = NULL,
    check.names = FALSE
  )

  names(camera_table)[-1] <- paste0(
    "Camera_",
    names(camera_table)[-1]
  )

  # --------------------------------------------------------------------------
  # Combine FRY and CAMERA results
  # --------------------------------------------------------------------------
  combined_table <- merge(
    fry_table,
    camera_table,
    by = "GeneSet",
    all = TRUE,
    sort = FALSE
  )

  # Restore the original gene-set ordering.
  combined_table$GeneSet <- factor(
    combined_table$GeneSet,
    levels = names(gene_sets)
  )

  combined_table <- combined_table[
    order(combined_table$GeneSet),
    ,
    drop = FALSE
  ]

  combined_table$GeneSet <- as.character(combined_table$GeneSet)

  combined_table <- combined_table %>%
    dplyr::mutate(
      AnalysisPool = analysis_pool,
      CohortType = cohort_type,
      Contrast = "Dys-CIM - Fun-CIM",
      FryTest = ifelse(
        cohort_type == "cell_line",
        "fry with cell-line block and duplicateCorrelation",
        "fry with independent patient samples"
      ),
      CameraTest = camera_implementation,
      PrimaryInference = "Nominal PValue"
    ) %>%
    dplyr::relocate(
      GeneSet,
      AnalysisPool,
      CohortType,
      Contrast,
      FryTest,
      CameraTest,
      PrimaryInference
    )

  # A separately tested prespecified pathway is not part of a pooled
  # multiple-testing family. Set its reported FDR columns to NA to avoid
  # implying that it was adjusted together with the main pathway collection.
  if (standalone) {
    fdr_columns <- grep(
      "FDR|adj",
      colnames(combined_table),
      ignore.case = TRUE,
      value = TRUE
    )

    for (column_name in fdr_columns) {
      combined_table[[column_name]] <- NA_real_
    }
  }

  combined_table
}



# ============================================================================
# DATA LOADING
# ============================================================================
clustered_plot_df <- data.table::fread(
  "oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epi_clustered_plot_df.csv",
  data.table = FALSE
)

clustered_plot_df <- as.data.frame(clustered_plot_df)

if ("cluster_assignments" %in% colnames(clustered_plot_df)) {
  clustered_plot_df <- clustered_plot_df %>%
    dplyr::rename(CIMIC_Cluster = cluster_assignments)
} else if (!"CIMIC_Cluster" %in% colnames(clustered_plot_df)) {
  stop("Error: Neither 'cluster_assignments' nor 'CIMIC_Cluster' column found in the metadata CSV.")
}

# FIX (#4): fail loudly rather than silently coercing to NA later.
clustered_plot_df$CIMIC_Cluster <- as.character(clustered_plot_df$CIMIC_Cluster)
if (!all(unique(clustered_plot_df$CIMIC_Cluster) %in% c("1", "2"))) {
  stop(sprintf(
    "CIMIC_Cluster contains unexpected values: %s. Expected only '1' (Fun-CIM) and '2' (Dys-CIM).",
    paste(setdiff(unique(clustered_plot_df$CIMIC_Cluster), c("1", "2")), collapse = ", ")
  ))
}

metadata_cols <- c(
  "sample_id", "CIMIC_Cluster", "base_id", "PACMAP1", "PACMAP2",
  "UMAP1", "UMAP2", "PC1", "PC2", "PC3", "V1"
)

gene_cols <- setdiff(colnames(clustered_plot_df), metadata_cols)

message(sprintf(
  "Loaded %d samples with %d gene expression columns",
  nrow(clustered_plot_df), length(gene_cols)
))

dir.create("oncoimmunology_paper/Results", showWarnings = FALSE)

# ============================================================================
# GENE SET DEFINITIONS
# ============================================================================

msig_df <- msigdbr::msigdbr(species = "Homo sapiens") %>% as.data.frame()

all_gene_set_names <- c(
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "GOBP_CELLULAR_RESPONSE_TO_STRESS",
  "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  "GOBP_INFLAMMASOME_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS",
  "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
  "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS",
  "REACTOME_CELLULAR_RESPONSE_TO_MITOCHONDRIAL_STRESS",
  "WP_MOLECULAR_PATHWAY_FOR_OXIDATIVE_STRESS",
  "WP_OXIDATIVE_STRESS_RESPONSE"
)

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
remove_prefixes <- c(
  "IGHV", "IGHD", "IGHJ", "IGHC", "IGHM", "IGHA", "IGHE", "IGHG",
  "IGLV", "IGLC", "IGLJ",
  "IGKV", "IGKC", "IGKJ",
  "TRAV", "TRBV", "TRDV", "TRGV",
  "TRAJ", "TRBJ", "TRGJ", "TRDJ",
  "TRAC", "TRBC", "TRGC", "TRDC"
)

pattern <- paste0("^(", paste(remove_prefixes, collapse = "|"), ")")

all_gene_sets <- lapply(all_gene_sets, function(gene_vec) {
  gene_vec[!grepl(pattern, gene_vec, perl = TRUE)]
})

all_gene_set_names_cimic <- all_gene_set_names


# ============================================================================
# GENE-WISE Z-SCORE NORMALIZATION (FOR HEATMAP DISPLAY ONLY)
# ============================================================================
# DISPLAY-ONLY. It must never be handed to fry()/limma downstream -- only to
# score_program_z() for the heatmap tile values.

zscore_safe <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(x)) else (x - m) / s
}

df_z <- clustered_plot_df %>%
  mutate(across(all_of(gene_cols), zscore_safe))

# ============================================================================
# PROGRAM SCORING: MEAN Z-SCORE PER SAMPLE (DISPLAY ONLY)
# ============================================================================
score_program_z <- function(data, genes, program_name) {
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

program_scores_all <- purrr::imap_dfr(all_gene_sets, ~ score_program_z(df_z, .x, .y))

program_scores_wide <- program_scores_all %>%
  dplyr::select(sample_id, Program, Score) %>%
  tidyr::pivot_wider(
    names_from = Program,
    values_from = Score
  )

# This left_join only adds the *_z summary columns; it does NOT overwrite the
# original gene_cols columns in clustered_plot_df, so clustered_plot_df keeps
# the raw expression values used below for the actual statistical test.
clustered_plot_df <- clustered_plot_df %>%
  left_join(program_scores_wide, by = "sample_id")

# ============================================================================
# HEATMAP DATA PREPARATION
# ============================================================================
str_wrap_length <- 40

genesets_of_interest <- paste0(names(all_gene_sets), "_z")

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
    Program = Program %>%
      str_remove_all("GOBP_|REACTOME_|HALLMARK_|WP_|_z") %>%
      str_replace_all("_", " ") %>%
      str_wrap(width = str_wrap_length)
  ) %>%
  group_by(CIMIC_Cluster, Program) %>%
  summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop")

# ============================================================================
# STATISTICAL TESTING: fry on gene-level expression matrix
# ============================================================================
#
# C1 = Fun-CIM, C2 = Dys-CIM
# Contrast: Dys-CIM - Fun-CIM
# Positive direction: enriched in Dys-CIM
# Negative direction: enriched in Fun-CIM
#
# ============================================================================

expr_matrix <- clustered_plot_df %>%
  dplyr::select(sample_id, all_of(gene_cols)) %>%
  tibble::column_to_rownames("sample_id") %>%
  t()

expr_matrix <- as.matrix(expr_matrix)
storage.mode(expr_matrix) <- "numeric"


fry_gene_sets <- lapply(
  all_gene_sets,
  function(gs) intersect(gs, rownames(expr_matrix))
)

fry_gene_sets <- fry_gene_sets[lengths(fry_gene_sets) >= 5]

message(sprintf("Running fry on %d gene sets", length(fry_gene_sets)))

# ---------------------------------------------------------------------------
# Design matrix
# ---------------------------------------------------------------------------
cluster_factor <- factor(
  clustered_plot_df$CIMIC_Cluster,
  levels = c("1", "2")
)

design <- model.matrix(~ 0 + cluster_factor)
colnames(design) <- c("Fun_CIM", "Dys_CIM")

# ---------------------------------------------------------------------------
# Cell-line blocking
# ---------------------------------------------------------------------------
if (!"cell_line" %in% colnames(clustered_plot_df)) {
  clustered_plot_df$cell_line <- sub(
    "^delta_([^_]+)_.*$", "\\1", clustered_plot_df$sample_id
  )
}

# FIX (#4): fail loudly if the regex didn't actually match (common source of
# a silently-broken block factor -> meaningless duplicateCorrelation).
if (any(clustered_plot_df$cell_line == clustered_plot_df$sample_id) ||
    any(is.na(clustered_plot_df$cell_line))) {
  stop("cell_line extraction failed for one or more sample_id values -- check the regex against actual sample_id formats.")
}

block_factor <- factor(clustered_plot_df$cell_line)

# ---------------------------------------------------------------------------
# Estimate within-cell-line correlation
# ---------------------------------------------------------------------------
dc <- limma::duplicateCorrelation(
  expr_matrix,
  design,
  block = block_factor
)

message(sprintf("Estimated consensus within-cell-line correlation: %.3f", dc$consensus.correlation))

# ---------------------------------------------------------------------------
# Dys-CIM versus Fun-CIM contrast
# ---------------------------------------------------------------------------
contrast <- limma::makeContrasts(
  Dys_CIM - Fun_CIM,
  levels = design
)

# ---------------------------------------------------------------------------
# Run fry
# FIX (#3): pass block + correlation through so the cell-line structure
# estimated above actually informs the test (previously computed but unused).
# ---------------------------------------------------------------------------
fry_res <- limma::fry(
  y = expr_matrix,
  index = fry_gene_sets,
  design = design,
  contrast = contrast,
  block = block_factor,
  correlation = dc$consensus.correlation
)

# ---------------------------------------------------------------------------
# Extract significance
# ---------------------------------------------------------------------------
sig_df <- fry_res %>%
  as.data.frame() %>%
  rownames_to_column("Program") %>%
  mutate(
    Program = Program %>%
      str_remove_all("GOBP_|REACTOME_|HALLMARK_|WP_") %>%
      str_replace_all("_", " ") %>%
      str_wrap(width = str_wrap_length),
    signif_label = case_when(
      PValue <= 0.0001 ~ "****",
      PValue <= 0.001 ~ "***",
      PValue <= 0.01 ~ "**",
      PValue <= 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  dplyr::select(Program, signif_label)

# ============================================================================
# PREPARE FACTORS FOR CONSISTENT ORDERING
# ============================================================================
program_order <- cluster_heatmap_df %>%
  distinct(Program) %>%
  pull(Program)

cluster_heatmap_df <- cluster_heatmap_df %>%
  mutate(Program = factor(Program, levels = rev(program_order)))

sig_df <- sig_df %>%
  mutate(Program = factor(Program, levels = rev(program_order)))

star_df <- cluster_heatmap_df %>%
  distinct(Program) %>%
  left_join(sig_df, by = "Program") %>%
  mutate(x = 2.8)

# ============================================================================
# HEATMAP VISUALIZATION
# ============================================================================
p_cluster_heatmap <- ggplot(
  cluster_heatmap_df,
  aes(x = CIMIC_Cluster, y = Program, fill = mean_score)
) +
  geom_tile(color = "black", linewidth = 0.8) +
  annotate(
    "text",
    x = 2.55,
    y = seq_along(levels(cluster_heatmap_df$Program)),
    label = rev(star_df$signif_label),
    size = 18,
    fontface = "bold",
    hjust = 0
  ) +
  scale_fill_gradient2(
    low = "#313695",
    mid = "white",
    high = "#A50026",
    midpoint = 0,
    name = "Mean\nZ-score",
    limits = c(-0.2, 0.2)
  ) +
  scale_x_discrete(
    labels = c("1" = "Fun\nCIM", "2" = "Dys\nCIM"),
    expand = expansion(add = c(0.5, 1.5))
  ) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x = element_text(face = "bold", size = 34, color = "black"),
    axis.text.y = element_text(face = "bold", size = 30, color = "black"),
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(face = "bold", size = 30),
    legend.title = element_text(face = "bold", size = 30),
    axis.line = element_blank()
  )

# ============================================================================
# FIGURE OUTPUT
# ============================================================================
ggsave(
  "oncoimmunology_paper/Results/fig5f_tnbc_zscore_heatmap.png",
  plot = p_cluster_heatmap,
  width = 18,
  height = 14,
  dpi = 300,
  bg = "white"
)

message("Figure saved: Results/fig5f_tnbc_zscore_heatmap.png (18x14 in, 300 DPI)")
message("Includes 13 immune/stress gene sets (BCR/TCR genes filtered)")
message("Significance stars: **** FDR<=0.0001, *** FDR<=0.001, ** FDR<=0.01, * FDR<=0.05")
message("Statistical test (fry) uses raw/log-scale expression with cell-line blocking applied")
message("Script complete: heatmap + fry analysis finished successfully")



# ============================================================================
# EXPORT FRY + BLOCK-ADJUSTED CAMERA STATISTICS
# FIGURE 4F: TNBC CELL-LINE Z-SCORE HEATMAP
# ============================================================================
# Contrast:
#   Dys_CIM - Fun_CIM
#
# Direction interpretation:
#   Up   = enriched in Dys-CIM
#   Down = enriched in Fun-CIM
#
# fry:
#   Self-contained gene-set test using cell-line blocking.
#
# cameraPR:
#   Competitive gene-set test based on moderated t-statistics from a limma
#   model fitted with the same cell-line block and duplicate correlation.
#
# The raw/log-scale expression matrix is used for both statistical tests.
# Gene-wise z-scores remain display-only.
# ============================================================================

# Ensure every tested gene set contains measured genes and retains its name
tested_gene_sets <- lapply(
  all_gene_sets,
  function(gs) intersect(unique(gs), rownames(expr_matrix))
)

tested_gene_sets <- tested_gene_sets[lengths(tested_gene_sets) >= 5]

if (length(tested_gene_sets) == 0) {
  stop("No gene sets contain at least 5 measured genes.")
}

# Convert gene symbols to row indices for consistent use by limma
tested_gene_indices <- limma::ids2indices(
  gene.sets = tested_gene_sets,
  identifiers = rownames(expr_matrix),
  remove.empty = TRUE
)

# ---------------------------------------------------------------------------
# FRY: self-contained test with cell-line blocking
# ---------------------------------------------------------------------------
# Re-run here with sort = "none" so output order matches the input pathway list.
fry_export_res <- limma::fry(
  y           = expr_matrix,
  index       = tested_gene_indices,
  design      = design,
  contrast    = contrast,
  block       = block_factor,
  correlation = dc$consensus.correlation,
  sort        = "none"
)

fry_export_df <- fry_export_res %>%
  as.data.frame() %>%
  tibble::rownames_to_column("GeneSet")

names(fry_export_df)[-1] <- paste0(
  "Fry_",
  names(fry_export_df)[-1]
)

# ---------------------------------------------------------------------------
# CAMERA: competitive test using the same blocked cell-line model
# ---------------------------------------------------------------------------
# camera() does not directly receive block/correlation arguments. Therefore:
#   1. Fit the blocked limma model.
#   2. Apply the Dys-CIM - Fun-CIM contrast.
#   3. Calculate moderated t-statistics.
#   4. Run cameraPR() on those block-adjusted statistics.

camera_fit <- limma::lmFit(
  object      = expr_matrix,
  design      = design,
  block       = block_factor,
  correlation = dc$consensus.correlation
)

camera_fit <- limma::contrasts.fit(
  fit       = camera_fit,
  contrasts = contrast
)

camera_fit <- limma::eBayes(
  camera_fit,
  trend = FALSE
)

if (ncol(camera_fit$t) != 1) {
  stop("Expected exactly one contrast for the CAMERA analysis.")
}

moderated_t <- camera_fit$t[, 1]
names(moderated_t) <- rownames(expr_matrix)

camera_res <- limma::cameraPR(
  statistic      = moderated_t,
  index          = tested_gene_indices,
  use.ranks      = FALSE,
  inter.gene.cor = 0.01,
  sort           = FALSE
)

camera_df <- camera_res %>%
  as.data.frame() %>%
  tibble::rownames_to_column("GeneSet")

names(camera_df)[-1] <- paste0(
  "Camera_",
  names(camera_df)[-1]
)

# ---------------------------------------------------------------------------
# Combine FRY and CAMERA statistics
# ---------------------------------------------------------------------------
fig5f_fry_camera_stats <- dplyr::full_join(
  fry_export_df,
  camera_df,
  by = "GeneSet"
)

# Restore the original pathway order
fig5f_fry_camera_stats$GeneSet <- factor(
  fig5f_fry_camera_stats$GeneSet,
  levels = names(tested_gene_indices)
)

fig5f_fry_camera_stats <- fig5f_fry_camera_stats %>%
  dplyr::arrange(GeneSet) %>%
  dplyr::mutate(
    GeneSet = as.character(GeneSet),
    Analysis = "Figure 4F TNBC cell-line immune/stress pathways",
    Contrast = "Dys-CIM - Fun-CIM",
    Fry_Test = paste0(
      "fry with cell-line block; consensus correlation = ",
      signif(dc$consensus.correlation, 4)
    ),
    Camera_Test = paste0(
      "cameraPR using moderated t-statistics from blocked limma model; ",
      "consensus correlation = ",
      signif(dc$consensus.correlation, 4)
    ),
    Direction_Interpretation = paste0(
      "Up = Dys-CIM; Down = Fun-CIM"
    ),
    Primary_Inference = "Nominal PValue; FDR retained for transparency"
  ) %>%
  dplyr::relocate(
    GeneSet,
    Analysis,
    Contrast,
    Direction_Interpretation,
    Fry_Test,
    Camera_Test,
    Primary_Inference
  )

# ---------------------------------------------------------------------------
# Append the mean Z-score per cluster (the values drawn as heatmap tiles).
#
# These are recomputed here directly from the per-sample `<gene set>_z` columns
# rather than joined off `cluster_heatmap_df`, because that data frame's
# `Program` labels are display strings (prefixes stripped, underscores replaced,
# and str_wrap()-ed with embedded newlines) and so cannot be matched reliably
# back to the raw MSigDB `GeneSet` identifiers used by the statistics table.
# The arithmetic is identical to the heatmap summarise() — same columns, same
# grouping, same mean(na.rm = TRUE) — only the join key differs.
#
# This must happen BEFORE fwrite() so the columns actually reach the CSV.
# ---------------------------------------------------------------------------
zscore_summary_f5 <- clustered_plot_df %>%
  dplyr::select(CIMIC_Cluster, all_of(genesets_of_interest)) %>%
  tidyr::pivot_longer(
    cols      = -CIMIC_Cluster,
    names_to  = "GeneSet",
    values_to = "Score"
  ) %>%
  dplyr::mutate(
    GeneSet = stringr::str_remove(GeneSet, "_z$"),
    Score   = as.numeric(Score)
  ) %>%
  dplyr::group_by(GeneSet, CIMIC_Cluster) %>%
  dplyr::summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop") %>%
  # TNBC_CL_Epirubicin coding: "1" = Fun-CIM, "2" = Dys-CIM
  # (note this is the reverse of the NKI_SMC and NEO cohorts)
  dplyr::mutate(
    ClusterColumn = dplyr::recode(
      as.character(CIMIC_Cluster),
      "1" = "MeanZscore_Cluster1_FunCIM",
      "2" = "MeanZscore_Cluster2_DysCIM"
    )
  ) %>%
  dplyr::select(GeneSet, ClusterColumn, mean_score) %>%
  tidyr::pivot_wider(names_from = ClusterColumn, values_from = mean_score) %>%
  dplyr::mutate(
    MeanZscore_Delta_DysCIM_minus_FunCIM =
      MeanZscore_Cluster2_DysCIM - MeanZscore_Cluster1_FunCIM
  )

fig5f_fry_camera_stats <- fig5f_fry_camera_stats %>%
  dplyr::left_join(zscore_summary_f5, by = "GeneSet")

if (anyNA(fig5f_fry_camera_stats$MeanZscore_Cluster2_DysCIM)) {
  warning(
    "Z-score columns are NA for: ",
    paste(
      fig5f_fry_camera_stats$GeneSet[
        is.na(fig5f_fry_camera_stats$MeanZscore_Cluster2_DysCIM)
      ],
      collapse = ", "
    )
  )
}

# ---------------------------------------------------------------------------
# Export CSV
# ---------------------------------------------------------------------------
fig5f_statistics_file <- paste0(
  "oncoimmunology_paper/Results/",
  "fig5f_tnbc_zscore_heatmap_fry_camera_statistics.csv"
)

data.table::fwrite(
  fig5f_fry_camera_stats,
  file = fig5f_statistics_file,
  na = "NA"
)

message(sprintf(
  "Saved FRY + CAMERA statistics: %s",
  fig5f_statistics_file
))

message(sprintf(
  "Exported %d pathways. Up = Dys-CIM; Down = Fun-CIM.",
  nrow(fig5f_fry_camera_stats)
))