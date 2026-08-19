# ============================================================================
# fig2e_zscore_heatmap_fry.R
# ============================================================================
# Purpose:
#   Z-score heatmap of 19 immune/stress gene sets across CIMIC clusters
#   (NKI_SMC cohort). Statistical testing via limma::fry with NOMINAL
#   p-values (no FDR) — justified by a priori pathway selection.
#   GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS is already in the 19-set list
#   and is tested inside the main fry call.
#   No cell-line blocking (one sample per patient).
#
# Inputs:
#   - oncoimmunology_paper/Datasets/NKI_SMC/nki_smc_combine_clustered_plot_df.csv
#
# Outputs:
#   - oncoimmunology_paper/Results/fig2e_zscore_heatmap_fry.png
#   - oncoimmunology_paper/Results/fig2e_zscore_heatmap_fry_camera_statistics.csv
#     (fry + camera statistics plus the mean Z-score per cluster shown on the
#      heatmap tiles)
#
# Dependencies: R ≥ 4.0
#   data.table, msigdbr, dplyr, tidyr, purrr, tibble, limma, ggplot2, stringr
#
# Author: Mohammed Gbadamosi
# ============================================================================

library(data.table)
library(msigdbr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
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
  "oncoimmunology_paper/Datasets/NKI_SMC/nki_smc_combine_clustered_plot_df.csv",
  data.table = FALSE
) %>% as.data.frame()

if ("cluster_assignments" %in% colnames(clustered_plot_df)) {
  clustered_plot_df <- clustered_plot_df %>%
    dplyr::rename(CIMIC_Cluster = cluster_assignments)
} else if (!"CIMIC_Cluster" %in% colnames(clustered_plot_df)) {
  stop("Neither 'cluster_assignments' nor 'CIMIC_Cluster' found.")
}

clustered_plot_df$CIMIC_Cluster <- as.character(clustered_plot_df$CIMIC_Cluster)
if (!all(unique(clustered_plot_df$CIMIC_Cluster) %in% c("1", "2"))) {
  stop(sprintf(
    "Unexpected CIMIC_Cluster values: %s",
    paste(setdiff(unique(clustered_plot_df$CIMIC_Cluster), c("1", "2")), collapse = ", ")
  ))
}

# NKI_SMC: "1" = Dys-CIM, "2" = Fun-CIM
metadata_cols <- c(
  "sample_id", "CIMIC_Cluster", "base_id", "PACMAP1", "PACMAP2",
  "UMAP1", "UMAP2", "PC1", "PC2", "PC3", "V1"
)
gene_cols <- setdiff(colnames(clustered_plot_df), metadata_cols)

message(sprintf("Loaded %d samples | %d gene columns", nrow(clustered_plot_df), length(gene_cols)))
dir.create("oncoimmunology_paper/Results", showWarnings = FALSE)

# ============================================================================
# GENE SET DEFINITIONS (19 sets — includes GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS)
# ============================================================================
msig_df <- msigdbr::msigdbr(species = "Homo sapiens") %>% as.data.frame()

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

all_gene_sets <- setNames(
  lapply(all_gene_set_names, function(gs) {
    msig_df[msig_df$gs_name == gs, "gene_symbol"]
  }),
  all_gene_set_names
)

# B/T cell receptor gene filtering
remove_prefixes <- c(
  "IGHV", "IGHD", "IGHJ", "IGHC", "IGHM", "IGHA", "IGHE", "IGHG",
  "IGLV", "IGLC", "IGLJ", "IGKV", "IGKC", "IGKJ",
  "TRAV", "TRBV", "TRDV", "TRGV",
  "TRAJ", "TRBJ", "TRGJ", "TRDJ",
  "TRAC", "TRBC", "TRGC", "TRDC"
)
bcr_pattern <- paste0("^(", paste(remove_prefixes, collapse = "|"), ")")
all_gene_sets <- lapply(all_gene_sets, function(g) g[!grepl(bcr_pattern, g, perl = TRUE)])

# ============================================================================
# Z-SCORE NORMALIZATION (DISPLAY ONLY — never fed to fry)
# ============================================================================
zscore_safe <- function(x) {
  s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
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
    return(data.frame(
      sample_id    = data$sample_id,
      cimic_cluster = data$CIMIC_Cluster,
      Program      = paste0(program_name, "_z"),
      Score        = rep(NA_real_, nrow(data))
    ))
  }
  score_vec <- data %>%
    dplyr::select(all_of(genes)) %>%
    as.matrix() %>%
    rowMeans(na.rm = TRUE)
  data.frame(
    sample_id    = data$sample_id,
    cimic_cluster = data$CIMIC_Cluster,
    Program      = paste0(program_name, "_z"),
    Score        = score_vec
  )
}

program_scores_all <- purrr::imap_dfr(all_gene_sets, ~ score_program_z(df_z, .x, .y))

program_scores_wide <- program_scores_all %>%
  dplyr::select(sample_id, Program, Score) %>%
  tidyr::pivot_wider(names_from = Program, values_from = Score)

clustered_plot_df <- clustered_plot_df %>%
  left_join(program_scores_wide, by = "sample_id")

# ============================================================================
# HEATMAP DATA PREPARATION
# ============================================================================
str_wrap_length <- 40
genesets_of_interest <- paste0(names(all_gene_sets), "_z")

cluster_heatmap_df <- clustered_plot_df %>%
  dplyr::select(CIMIC_Cluster, all_of(genesets_of_interest)) %>%
  pivot_longer(cols = -CIMIC_Cluster, names_to = "Program", values_to = "Score") %>%
  mutate(
    Score        = as.numeric(Score),
    CIMIC_Cluster = factor(CIMIC_Cluster),
    Program      = Program %>%
      str_remove_all("GOBP_|REACTOME_|HALLMARK_|WP_|_z") %>%
      str_replace_all("_", " ") %>%
      str_wrap(width = str_wrap_length)
  ) %>%
  group_by(CIMIC_Cluster, Program) %>%
  summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop")

# ============================================================================
# STATISTICAL TESTING: limma::fry — NOMINAL p-values, NO blocking
# (one patient per sample — no repeated measures structure)
# ============================================================================
expr_matrix <- clustered_plot_df %>%
  dplyr::select(sample_id, all_of(gene_cols)) %>%
  tibble::column_to_rownames("sample_id") %>%
  t() %>%
  as.matrix()
storage.mode(expr_matrix) <- "numeric"

fry_gene_sets <- lapply(all_gene_sets, function(gs) intersect(gs, rownames(expr_matrix)))
fry_gene_sets <- fry_gene_sets[lengths(fry_gene_sets) >= 5]
message(sprintf("Running fry on %d gene sets", length(fry_gene_sets)))

# NKI_SMC: "1" = Dys-CIM, "2" = Fun-CIM
cluster_factor <- factor(clustered_plot_df$CIMIC_Cluster, levels = c("2", "1"))
design <- model.matrix(~ 0 + cluster_factor)
colnames(design) <- c("Fun_CIM", "Dys_CIM")

contrast <- limma::makeContrasts(Dys_CIM - Fun_CIM, levels = design)

# No block — singular patients, no duplicateCorrelation needed
fry_res <- limma::fry(
  y        = expr_matrix,
  index    = fry_gene_sets,
  design   = design,
  contrast = contrast
)

# Extract NOMINAL p-values (no FDR — a priori selected pathways)
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
      PValue <= 0.001  ~ "***",
      PValue <= 0.01   ~ "**",
      PValue <= 0.05   ~ "*",
      TRUE             ~ ""
    )
  ) %>%
  dplyr::select(Program, signif_label)

# ============================================================================
# FACTOR ORDERING
# ============================================================================
program_order <- cluster_heatmap_df %>% distinct(Program) %>% pull(Program)

cluster_heatmap_df <- cluster_heatmap_df %>%
  mutate(Program = factor(Program, levels = rev(program_order)))

sig_df <- sig_df %>%
  mutate(Program = factor(Program, levels = rev(program_order)))

star_df <- cluster_heatmap_df %>%
  distinct(Program) %>%
  left_join(sig_df, by = "Program") %>%
  mutate(x = 2.8)

# ============================================================================
# HEATMAP PLOT
# ============================================================================
p_heatmap <- ggplot(
  cluster_heatmap_df,
  aes(x = CIMIC_Cluster, y = Program, fill = mean_score)
) +
  geom_tile(color = "black", linewidth = 0.8) +
  annotate(
    "text",
    x        = 2.55,
    y        = seq_along(levels(cluster_heatmap_df$Program)),
    label    = rev(star_df$signif_label),
    size     = 15,
    fontface = "bold",
    hjust    = 0
  ) +
  scale_fill_gradient2(
    low      = "#313695",
    mid      = "white",
    high     = "#A50026",
    midpoint = 0,
    name     = "Mean\nZ-score",
    limits   = c(-0.5, 0.5)
  ) +
  scale_x_discrete(
    labels = c("1" = "Dys\nCIM", "2" = "Fun\nCIM"),
    expand = expansion(add = c(0.5, 1.5))
  ) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x   = element_text(face = "bold", size = 28, color = "black"),
    axis.text.y   = element_text(face = "bold", size = 20, color = "black"),
    legend.position      = "right",
    legend.justification = "center",
    legend.text   = element_text(face = "bold", size = 20),
    legend.title  = element_text(face = "bold", size = 20),
    axis.line     = element_blank()
  )

# ============================================================================
# SAVE
# ============================================================================
suppressWarnings(
  ggsave(
    "oncoimmunology_paper/Results/fig2e_zscore_heatmap_fry.png",
    plot   = p_heatmap,
    width  = 14,
    height = 14,
    dpi    = 300,
    bg     = "white"
  )
)

message("✓ Saved: Results/fig2e_zscore_heatmap_fry.png")
message("  19 gene sets | NKI_SMC cohort | fry | nominal p | no block")


# ============================================================================
# EXPORT FRY + CAMERA STATISTICS: NKI/SMC FIGURE 2E
# ============================================================================

# ---------------------------------------------------------------------------
# Export FRY + CAMERA statistics with mean Z-score per cluster
# ---------------------------------------------------------------------------
fig2e_stats <- make_fry_camera_table(
  expr_matrix  = expr_matrix,
  gene_sets    = fry_gene_sets,
  design       = design,
  contrast     = contrast,
  cohort_type  = "human",
  analysis_pool = "19 a priori immune/stress pathways",
  standalone   = FALSE
)

# ---------------------------------------------------------------------------
# Mean Z-score per gene set per cluster (the values drawn as heatmap tiles).
#
# These are recomputed here directly from the per-sample `<gene set>_z` columns
# rather than joined off `cluster_heatmap_df`, because that data frame's
# `Program` labels are display strings (prefixes stripped, underscores replaced,
# and str_wrap()-ed with embedded newlines) and so cannot be matched reliably
# back to the raw MSigDB `GeneSet` identifiers used by the statistics table.
# The arithmetic is identical to the heatmap summarise() — same columns, same
# grouping, same mean(na.rm = TRUE) — only the join key differs.
# ---------------------------------------------------------------------------
zscore_summary <- clustered_plot_df %>%
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
  # NKI_SMC coding: "1" = Dys-CIM, "2" = Fun-CIM
  dplyr::mutate(
    ClusterColumn = dplyr::recode(
      as.character(CIMIC_Cluster),
      "1" = "MeanZscore_Cluster1_DysCIM",
      "2" = "MeanZscore_Cluster2_FunCIM"
    )
  ) %>%
  dplyr::select(GeneSet, ClusterColumn, mean_score) %>%
  tidyr::pivot_wider(names_from = ClusterColumn, values_from = mean_score) %>%
  dplyr::mutate(
    MeanZscore_Delta_DysCIM_minus_FunCIM =
      MeanZscore_Cluster1_DysCIM - MeanZscore_Cluster2_FunCIM
  )

fig2e_stats <- fig2e_stats %>%
  dplyr::left_join(zscore_summary, by = "GeneSet")

if (anyNA(fig2e_stats$MeanZscore_Cluster1_DysCIM)) {
  warning(
    "Z-score columns are NA for: ",
    paste(
      fig2e_stats$GeneSet[is.na(fig2e_stats$MeanZscore_Cluster1_DysCIM)],
      collapse = ", "
    )
  )
}

data.table::fwrite(
  fig2e_stats,
  file = paste0(
    "oncoimmunology_paper/Results/",
    "fig2e_zscore_heatmap_fry_camera_statistics.csv"
  ),
  na = "NA"
)

message(
  "Saved: Results/fig2e_zscore_heatmap_fry_camera_statistics.csv"
)