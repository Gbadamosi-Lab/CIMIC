# ============================================================================
# fig2f_2k_icd_programs_fry.R
# ============================================================================
# Purpose:
#   Box-plots of the five ICD programs (Apoptosis, Necroptosis, Pyroptosis,
#   PANoptosis, Ferroptosis) + standalone GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS
#   across CIMIC clusters for the NKI_SMC cohort.
#
#   Statistical approach: limma::fry, NOMINAL p-values (no FDR).
#   No blocking — one patient per sample.
#   NKI_SMC cluster coding: "1" = Dys-CIM, "2" = Fun-CIM.
#
#   The five ICD programs are tested together in one fry call.
#   GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS is tested in a SEPARATE
#   standalone fry call and is NOT included in the main pool.
#
# Inputs:
#   - oncoimmunology_paper/Datasets/NKI_SMC/nki_smc_combine_clustered_plot_df.csv
#
# Outputs:
#   - Results/fig2f_2k_icd_programs_fry.png           (faceted, 5 programs)
#   - Results/fig2f_2k_<Program>_fry.png              (individual, 5 programs)
#   - Results/fig2f_2k_Inflammatory_Cell_Apoptosis_fry.png  (standalone)
#   - Results/fig2f_2k_icd_programs_fry_camera_statistics.csv
#     (fry + camera statistics plus the mean Z-score per cluster plotted in
#      the box-plots)
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

metadata_cols <- c(
  "sample_id", "CIMIC_Cluster", "base_id", "PACMAP1", "PACMAP2",
  "UMAP1", "UMAP2", "PC1", "PC2", "PC3", "V1"
)
gene_cols <- setdiff(colnames(clustered_plot_df), metadata_cols)

message(sprintf("Loaded %d samples | %d gene columns", nrow(clustered_plot_df), length(gene_cols)))
dir.create("oncoimmunology_paper/Results", showWarnings = FALSE)

# ============================================================================
# GENE SET DEFINITIONS
# ============================================================================

# --- Five ICD programs (main pool) ---
death_programs <- list(
  Apoptosis = c(
    "BAX", "BAK1", "BBC3", "PMAIP1", "BCL2L11",
    "APAF1", "CASP9", "CASP3", "CASP7",
    "FAS", "TNFRSF10B", "CASP8", "FADD", "BID",
    "DIABLO", "CYCS", "TP53AIP1", "BAD", "BMF", "HRK"
  ),
  Necroptosis = c(
    "RIPK1", "RIPK3", "MLKL", "ZBP1", "TICAM1",
    "TNFRSF1A", "FADD", "CASP8"
  ),
  Pyroptosis = c(
    "GSDME", "GSDMD", "CASP1", "CASP4", "CASP5",
    "NLRP3", "AIM2", "PYCARD", "IL1B", "IL18",
    "CASP3", "NLRC4"
  ),
  PANoptosis = c(
    "ZBP1", "AIM2", "RIPK3", "RIPK1", "CASP8",
    "CASP1", "FADD", "PYCARD", "IRF1"
  ),
  Ferroptosis = c(
    "ACSL4", "LPCAT3", "TFRC", "SAT1", "PTGS2"
  )
)

# --- Standalone gene set (NOT in main fry pool) ---
msig_df <- msigdbr::msigdbr(species = "Homo sapiens") %>% as.data.frame()

inflam_apop_genes <- msig_df[
  msig_df$gs_name == "GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS",
  "gene_symbol"
]
inflam_apop_label <- "Inflammatory Cell\nApoptotic Process"

# ============================================================================
# Z-SCORE NORMALIZATION (DISPLAY ONLY)
# ============================================================================
zscore_safe <- function(x) {
  s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(x)) else (x - m) / s
}

df_z <- clustered_plot_df %>%
  dplyr::mutate(across(all_of(gene_cols), zscore_safe))

# ============================================================================
# PROGRAM SCORING
# ============================================================================
score_program_z <- function(data, genes, prog_name) {
  genes <- intersect(genes, colnames(data))
  if (length(genes) == 0) {
    return(data.frame(
      sample_id    = data$sample_id,
      CIMIC_Cluster = data$CIMIC_Cluster,
      Program      = prog_name,
      Score        = rep(NA_real_, nrow(data))
    ))
  }
  score_vec <- data %>%
    dplyr::select(all_of(genes)) %>%
    as.matrix() %>%
    rowMeans(na.rm = TRUE)
  data.frame(
    sample_id    = data$sample_id,
    CIMIC_Cluster = data$CIMIC_Cluster,
    Program      = prog_name,
    Score        = score_vec
  )
}

program_scores_all <- purrr::imap_dfr(death_programs, ~ score_program_z(df_z, .x, .y))

# Score for the standalone set
inflam_apop_scores <- score_program_z(df_z, inflam_apop_genes, inflam_apop_label)

# ============================================================================
# BUILD EXPRESSION MATRIX
# ============================================================================
expr_matrix <- clustered_plot_df %>%
  dplyr::select(sample_id, all_of(gene_cols)) %>%
  tibble::column_to_rownames("sample_id") %>%
  t() %>%
  as.matrix()
storage.mode(expr_matrix) <- "numeric"

# NKI_SMC: "1" = Dys-CIM, "2" = Fun-CIM
cluster_factor <- factor(clustered_plot_df$CIMIC_Cluster, levels = c("2", "1"))
design <- model.matrix(~ 0 + cluster_factor)
colnames(design) <- c("Fun_CIM", "Dys_CIM")

contrast <- limma::makeContrasts(Dys_CIM - Fun_CIM, levels = design)

# ============================================================================
# MAIN FRY: five ICD programs — NOMINAL p-values, no block
# ============================================================================
fry_gene_sets <- lapply(death_programs, function(gs) intersect(gs, rownames(expr_matrix)))
fry_gene_sets <- fry_gene_sets[lengths(fry_gene_sets) >= 5]
message(sprintf("Main fry: %d gene sets", length(fry_gene_sets)))

fry_res <- limma::fry(
  y        = expr_matrix,
  index    = fry_gene_sets,
  design   = design,
  contrast = contrast
)

str_wrap_length <- 40

# Use FDR-adjusted p-values for significance labeling in the main ICD program analysis.
sig_df <- fry_res %>%
  as.data.frame() %>%
  rownames_to_column("Program") %>%
  mutate(
    signif_label = case_when(
      FDR <= 0.0001 ~ "****",
      FDR <= 0.001  ~ "***",
      FDR <= 0.01   ~ "**",
      FDR <= 0.05   ~ "*",
      FDR <  0.15   ~ sprintf("p = %.3g", FDR),
      TRUE          ~ "ns"
    )
  ) %>%
  dplyr::select(Program, signif_label)

# ============================================================================
# STANDALONE FRY: GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS
# ============================================================================
inflam_idx <- list(
  inflam_apop_label = intersect(inflam_apop_genes, rownames(expr_matrix))
)

if (length(inflam_idx[[1]]) >= 5) {
  fry_inflam <- limma::fry(
    y        = expr_matrix,
    index    = inflam_idx,
    design   = design,
    contrast = contrast
  )
  inflam_p     <- fry_inflam["inflam_apop_label", "PValue"]
  inflam_signif <- case_when(
    inflam_p <= 0.0001 ~ "****",
    inflam_p <= 0.001  ~ "***",
    inflam_p <= 0.01   ~ "**",
    inflam_p <= 0.05   ~ "*",
    inflam_p <  0.15   ~ sprintf("p = %.3g", inflam_p),
    TRUE               ~ "ns"
  )
  message(sprintf("Standalone fry [%s]: PValue = %.4f (%s)",
                  inflam_apop_label, inflam_p, inflam_signif))
} else {
  inflam_signif <- "n/a"
  message("Standalone gene set has < 5 genes in expr_matrix — skipped.")
}

inflam_sig_df <- data.frame(
  Program      = inflam_apop_label,
  signif_label = inflam_signif,
  stringsAsFactors = FALSE
)

# ============================================================================
# SHARED PLOT HELPERS
# ============================================================================
cluster_labels <- c("1" = "Dys-CIM", "2" = "Fun-CIM")
fill_colors    <- c("1" = "#F8766D", "2" = "#00BFC4")

# ============================================================================
# FACETED PLOT (5 ICD programs) — style matches individual plots
# ============================================================================
plot_df <- program_scores_all %>%
  dplyr::left_join(sig_df, by = "Program") %>%
  dplyr::mutate(CIMIC_Cluster = factor(CIMIC_Cluster, levels = c("1", "2")))

star_df_facet <- plot_df %>%
  group_by(Program) %>%
  summarise(y = max(Score, na.rm = TRUE) * 0.90, .groups = "drop") %>%
  left_join(sig_df, by = "Program")

p_facet <- ggplot(
  plot_df,
  aes(x = CIMIC_Cluster, y = Score, fill = CIMIC_Cluster)
) +
  geom_boxplot(
    outlier.shape = NA,
    colour = "black",
    width  = 0.6,
    alpha  = 0.75
  ) +
  geom_jitter(
    width  = 0.15,
    size   = 6,
    shape  = 16,
    colour = "black",
    alpha  = 0.75
  ) +
  geom_text(
    data        = star_df_facet,
    aes(x = 1.5, y = y, label = signif_label),
    inherit.aes = FALSE,
    size        = 10,
    fontface    = "bold"
  ) +
  facet_wrap(~Program, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = fill_colors) +
  scale_x_discrete(labels = cluster_labels, expand = expansion(add = c(0.5, 0.5))) +
  labs(
    title = "ICD Programs — NKI_SMC (fry, nominal p)",
    y     = "Z-score of\nProgram Activation",
    x     = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 16),
    axis.text.x      = element_text(face = "bold", size = 24, colour = "black"),
    axis.text.y      = element_text(face = "bold", size = 14, colour = "black"),
    axis.title.y     = element_text(face = "bold", size = 16, colour = "black"),
    legend.position  = "none",
    panel.grid       = element_blank(),
    panel.border     = element_rect(linewidth = 1, fill = NA)
  )

ggsave(
  "oncoimmunology_paper/Results/fig2f_2k_icd_programs_fry.png",
  plot   = p_facet,
  width  = 12,
  height = 14,
  dpi    = 300,
  bg     = "white"
)
message("✓ Saved: Results/fig2f_2k_icd_programs_fry.png")

# ============================================================================
# INDIVIDUAL PLOTS (5 ICD programs)
# ============================================================================
make_individual_boxplots <- function(sig_df) {
  for (prog in unique(program_scores_all$Program)) {

    prog_df <- program_scores_all %>%
      dplyr::filter(Program == prog) %>%
      dplyr::left_join(sig_df %>% dplyr::filter(Program == prog), by = "Program") %>%
      dplyr::mutate(CIMIC_Cluster = factor(CIMIC_Cluster, levels = c("1", "2")))

    star_y <- max(prog_df$Score, na.rm = TRUE) * 0.90

    p_ind <- ggplot(prog_df, aes(x = CIMIC_Cluster, y = Score, fill = CIMIC_Cluster)) +
      geom_boxplot(outlier.shape = NA, colour = "black", width = 0.6, alpha = 0.75) +
      geom_jitter(width = 0.15, size = 6, shape = 16, colour = "black", alpha = 0.75) +
      geom_text(
        aes(x = 1.5, y = star_y, label = signif_label),
        inherit.aes = FALSE, size = 16, fontface = "bold"
      ) +
      labs(
        title = prog,
        y     = "Z-score of\nProgram Activation",
        x     = NULL
      ) +
      scale_fill_manual(values = fill_colors) +
      scale_x_discrete(labels = cluster_labels, expand = expansion(add = c(0.5, 0.5))) +
      theme_classic(base_size = 18) +
      theme(
        axis.text.x  = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y  = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 30, colour = "black"),
        plot.title   = element_text(face = "bold", size = 36, hjust = 0.5),
        legend.position = "none",
        panel.grid   = element_blank(),
        panel.border = element_rect(linewidth = 1, fill = NA)
      )

    out_file <- sprintf(
      "oncoimmunology_paper/Results/fig2f_2k_%s_fry.png",
      gsub(" ", "_", prog)
    )
    ggsave(out_file, plot = p_ind, width = 7, height = 7, dpi = 300, bg = "white")
    message(sprintf("  ✓ Saved: %s", out_file))
  }
}

make_individual_boxplots(sig_df)

# ============================================================================
# STANDALONE PLOT: GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS
# ============================================================================
inflam_plot_df <- inflam_apop_scores %>%
  dplyr::mutate(
    CIMIC_Cluster = factor(CIMIC_Cluster, levels = c("1", "2")),
    signif_label  = inflam_signif
  )

star_y_inflam <- max(inflam_plot_df$Score, na.rm = TRUE) * 0.90

p_inflam <- ggplot(inflam_plot_df, aes(x = CIMIC_Cluster, y = Score, fill = CIMIC_Cluster)) +
  geom_boxplot(outlier.shape = NA, colour = "black", width = 0.6, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 6, shape = 16, colour = "black", alpha = 0.75) +
  geom_text(
    aes(x = 1.5, y = star_y_inflam, label = signif_label),
    inherit.aes = FALSE, size = 16, fontface = "bold"
  ) +
  labs(
    title = "Inflammatory Cell\nApoptotic Process",
    y     = "Z-score of\nProgram Activation",
    x     = NULL
  ) +
  scale_fill_manual(values = fill_colors) +
  scale_x_discrete(labels = cluster_labels, expand = expansion(add = c(0.5, 0.5))) +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x   = element_text(face = "bold", size = 30, colour = "black"),
    axis.text.y   = element_text(face = "bold", size = 30, colour = "black"),
    axis.title.y  = element_text(face = "bold", size = 30, colour = "black"),
    plot.title    = element_text(face = "bold", size = 36, hjust = 0.5),
    plot.caption  = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "none",
    panel.grid    = element_blank(),
    panel.border  = element_rect(linewidth = 1, fill = NA)
  )

ggsave(
  "oncoimmunology_paper/Results/fig2f_2k_Inflammatory_Cell_Apoptosis_fry.png",
  plot   = p_inflam,
  width  = 7,
  height = 7,
  dpi    = 300,
  bg     = "white"
)
message("✓ Saved: Results/fig2f_2k_Inflammatory_Cell_Apoptosis_fry.png")
message("  [Standalone fry — excluded from main pool]")


# ============================================================================
# EXPORT FRY + CAMERA STATISTICS: NKI/SMC ICD PROGRAMS
# ============================================================================
# The five ICD programs form one multiple-testing family; the prespecified
# GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS pathway is tested standalone and its
# FDR columns are therefore blanked by make_fry_camera_table(standalone = TRUE).
# ============================================================================

fig2f_main_stats <- make_fry_camera_table(
  expr_matrix   = expr_matrix,
  gene_sets     = fry_gene_sets,
  design        = design,
  contrast      = contrast,
  cohort_type   = "human",
  analysis_pool = "Five a priori ICD programs",
  standalone    = FALSE
)

fig2f_inflammatory_set <- list(
  GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS =
    intersect(inflam_apop_genes, rownames(expr_matrix))
)

fig2f_inflammatory_stats <- make_fry_camera_table(
  expr_matrix   = expr_matrix,
  gene_sets     = fig2f_inflammatory_set,
  design        = design,
  contrast      = contrast,
  cohort_type   = "human",
  analysis_pool = "Standalone prespecified pathway",
  standalone    = TRUE
)

fig2f_all_stats <- dplyr::bind_rows(
  fig2f_main_stats,
  fig2f_inflammatory_stats
)

# ---------------------------------------------------------------------------
# Mean Z-score per program per cluster (the values shown in the box-plots).
#
# Same source as the plots: per-sample program scores produced by
# score_program_z() from the gene-wise z-scored matrix `df_z`. The summary is
# keyed on the raw gene-set identifier used by the statistics table rather than
# on the display label, so the join cannot silently drop rows.
# ---------------------------------------------------------------------------
icd_zscore_long <- dplyr::bind_rows(
  program_scores_all %>%
    dplyr::transmute(
      GeneSet       = Program,
      CIMIC_Cluster = as.character(CIMIC_Cluster),
      Score
    ),
  inflam_apop_scores %>%
    dplyr::transmute(
      GeneSet       = "GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS",
      CIMIC_Cluster = as.character(CIMIC_Cluster),
      Score
    )
)

zscore_summary <- icd_zscore_long %>%
  dplyr::group_by(GeneSet, CIMIC_Cluster) %>%
  dplyr::summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop") %>%
  # NKI_SMC coding: "1" = Dys-CIM, "2" = Fun-CIM
  dplyr::mutate(
    ClusterColumn = dplyr::recode(
      CIMIC_Cluster,
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

fig2f_all_stats <- fig2f_all_stats %>%
  dplyr::left_join(zscore_summary, by = "GeneSet")

if (anyNA(fig2f_all_stats$MeanZscore_Cluster1_DysCIM)) {
  warning(
    "Z-score columns are NA for: ",
    paste(
      fig2f_all_stats$GeneSet[is.na(fig2f_all_stats$MeanZscore_Cluster1_DysCIM)],
      collapse = ", "
    )
  )
}

data.table::fwrite(
  fig2f_all_stats,
  file = paste0(
    "oncoimmunology_paper/Results/",
    "fig2f_2k_icd_programs_fry_camera_statistics.csv"
  ),
  na = "NA"
)

message(
  "Saved: Results/fig2f_2k_icd_programs_fry_camera_statistics.csv"
)