# ============================================================================
# fig5i_5l_tnbc_cl_icd_programs_fry.R   [UPDATED]
# ============================================================================
# Purpose:
#   Box-plots of the five ICD programs across CIMIC clusters for the
#   TNBC_CL_Epirubicin dataset. Statistical approach: limma::fry with
#   duplicateCorrelation (cell-line blocking). NOMINAL p-values (no FDR).
#
#   Changes from previous version:
#     1. Faceted plot now matches individual plot style (alpha=0.75, size=6,
#        shape=16 dots, panel.border, axis text sizes consistent).
#     2. Points are plain dots (shape=16) — no shape-by-cell-line mapping.
#     3. GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS added as a STANDALONE
#        separate fry call, NOT included in the main pool.
#
#   TNBC_CL: "1" = Fun-CIM, "2" = Dys-CIM
#
# Inputs:
#   - oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epi_clustered_plot_df.csv
#
# Outputs:
#   - Results/fig5i_5l_tnbc_cl_icd_programs_fry.png           (faceted)
#   - Results/fig5i_5l_tnbc_cl_icd_<Program>_fry.png          (individual x5)
#   - Results/fig5i_5l_tnbc_cl_Inflammatory_Cell_Apoptosis_fry.png (standalone)
#   - Results/fig5i_5l_tnbc_cl_icd_programs_fry_camera_statistics.csv
#     (fry + block-adjusted cameraPR statistics plus the mean Z-score per
#      cluster plotted in the box-plots)
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
  "oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epi_clustered_plot_df.csv",
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

# Five ICD programs (main pool)
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

# Standalone gene set
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
      cimic_cluster = data$CIMIC_Cluster,
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
    cimic_cluster = data$CIMIC_Cluster,
    Program      = prog_name,
    Score        = score_vec
  )
}

# Cell-line blocking
if (!"cell_line" %in% colnames(clustered_plot_df)) {
  clustered_plot_df$cell_line <- sub(
    "^delta_([^_]+)_.*$", "\\1", clustered_plot_df$sample_id
  )
}

if (any(clustered_plot_df$cell_line == clustered_plot_df$sample_id) ||
    any(is.na(clustered_plot_df$cell_line))) {
  stop("cell_line extraction failed — check regex against sample_id format.")
}


# Compute program scores and retain cell line information for shape mapping
program_scores_all <- purrr::imap_dfr(death_programs, ~ score_program_z(df_z, .x, .y)) %>%
  dplyr::left_join(
    clustered_plot_df %>% dplyr::select(sample_id, cell_line),
    by = "sample_id"
  )
inflam_apop_scores <- score_program_z(df_z, inflam_apop_genes, inflam_apop_label) %>%
  dplyr::left_join(
    clustered_plot_df %>% dplyr::select(sample_id, cell_line),
    by = "sample_id"
  )

# ============================================================================
# EXPRESSION MATRIX FOR FRY
# ============================================================================
expr_matrix <- clustered_plot_df %>%
  dplyr::select(sample_id, all_of(gene_cols)) %>%
  tibble::column_to_rownames("sample_id") %>%
  t() %>%
  as.matrix()
storage.mode(expr_matrix) <- "numeric"

# TNBC_CL: "1" = Fun-CIM, "2" = Dys-CIM
cluster_factor <- factor(clustered_plot_df$CIMIC_Cluster, levels = c("1", "2"))
design <- model.matrix(~ 0 + cluster_factor)
colnames(design) <- c("Fun_CIM", "Dys_CIM")


block_factor <- factor(clustered_plot_df$cell_line)

dc <- limma::duplicateCorrelation(expr_matrix, design, block = block_factor)
message(sprintf("Consensus within-cell-line correlation: %.3f", dc$consensus.correlation))

contrast <- limma::makeContrasts(Dys_CIM - Fun_CIM, levels = design)

# ============================================================================
# MAIN FRY: five ICD programs — NOMINAL p-values
# ============================================================================
fry_gene_sets <- lapply(death_programs, function(gs) intersect(gs, rownames(expr_matrix)))
fry_gene_sets <- fry_gene_sets[lengths(fry_gene_sets) >= 5]
message(sprintf("Main fry: %d gene sets", length(fry_gene_sets)))

fry_res <- limma::fry(
  y           = expr_matrix,
  index       = fry_gene_sets,
  design      = design,
  contrast    = contrast,
  block       = block_factor,
  correlation = dc$consensus.correlation
)

str_wrap_length <- 40

sig_df <- fry_res %>%
  as.data.frame() %>%
  rownames_to_column("Program") %>%
  mutate(
    signif_label = case_when(
      PValue <= 0.0001 ~ "****",
      PValue <= 0.001  ~ "***",
      PValue <= 0.01   ~ "**",
      PValue <= 0.05   ~ "*",
      PValue <  0.15   ~ sprintf("p = %.3g", PValue),
      TRUE             ~ "ns"
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
    y           = expr_matrix,
    index       = inflam_idx,
    design      = design,
    contrast    = contrast,
    block       = block_factor,
    correlation = dc$consensus.correlation
  )
  inflam_p      <- fry_inflam["inflam_apop_label", "PValue"]
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
  message("Standalone gene set < 5 genes — skipped.")
}

# ============================================================================
# SHARED AESTHETICS
# ============================================================================
cluster_labels <- c("2" = "Dys-CIM", "1" = "Fun-CIM")
fill_colors    <- c("2" = "#F8766D", "1" = "#00BFC4")

shared_theme <- theme_classic(base_size = 18) +
  theme(
    axis.text.x   = element_text(face = "bold", size = 24, colour = "black"),
    axis.text.y   = element_text(face = "bold", size = 24, colour = "black"),
    axis.title.y  = element_text(face = "bold", size = 24, colour = "black"),
    plot.title    = element_text(face = "bold", size = 24, hjust = 0.5),
    legend.position = "none",
    panel.grid    = element_blank(),
    panel.border  = element_rect(linewidth = 1, fill = NA)
  )


# ============================================================================
# CELL LINE SHAPE DEFINITIONS
# ============================================================================
# Define a vector of shape codes
# From original script: fig5i_5l_tnbc_cl_icd_programs_fry.R
# Using a specific set of shapes to ensure visual distinctiveness
shape_values <- c(
  16, 17, 15, 18, 8, 3, 4, 7, 10
)

# Create a named vector mapping cell lines to shapes.
# This ensures that "BT549" always gets the same shape, "DU4475" the next, etc.
cell_line_shapes <- setNames(
  rep(
    shape_values,
    length.out = length(unique(clustered_plot_df$cell_line))
  ),
  sort(unique(clustered_plot_df$cell_line))
)

# ============================================================================
# FACETED PLOT — now matches individual plot style
# ============================================================================
plot_df_facet <- program_scores_all %>%
  dplyr::left_join(sig_df, by = "Program") %>%
  dplyr::mutate(CIMIC_Cluster = factor(cimic_cluster, levels = c("2", "1")))

star_df_facet <- plot_df_facet %>%
  group_by(Program) %>%
  summarise(y = max(Score, na.rm = TRUE) * 0.90, .groups = "drop") %>%
  left_join(sig_df, by = "Program")

p_facet <- ggplot(
  plot_df_facet,
  aes(x = CIMIC_Cluster, y = Score, fill = CIMIC_Cluster)
) +
  geom_boxplot(
    outlier.shape = NA,
    colour = "black",
    width  = 0.6,
    alpha  = 0.75
  ) +
  # Use the cell_line column for the shape aesthetic
  geom_jitter(
    aes(shape = cell_line),
    width  = 0.15,
    size   = 6,
    colour = "black",
    alpha  = 0.75
  ) +
  # Add the manual shape scale
  scale_shape_manual(name = "Cell line", values = cell_line_shapes) +
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
    title = "ICD Programs — TNBC_CL (fry, nominal p)",
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
    # Show legend on the right
    legend.position  = "right",
    panel.grid       = element_blank(),
    panel.border     = element_rect(linewidth = 1, fill = NA)
  )

ggsave(
  "oncoimmunology_paper/Results/fig5i_5l_tnbc_cl_icd_programs_fry.png",
  plot   = p_facet,
  width  = 12,
  height = 14,
  dpi    = 300,
  bg     = "white"
)
message("✓ Saved: Results/fig5i_5l_tnbc_cl_icd_programs_fry.png")

# ============================================================================
# INDIVIDUAL PLOTS (5 ICD programs)
# ============================================================================
make_individual_boxplots <- function(sig_df) {
  for (prog in unique(program_scores_all$Program)) {

    prog_df <- program_scores_all %>%
      dplyr::filter(Program == prog) %>%
      dplyr::left_join(sig_df %>% dplyr::filter(Program == prog), by = "Program") %>%
      dplyr::mutate(CIMIC_Cluster = factor(cimic_cluster, levels = c("2", "1")))

    star_y <- max(prog_df$Score, na.rm = TRUE) * 0.90

    p_ind <- ggplot(prog_df, aes(x = CIMIC_Cluster, y = Score, fill = CIMIC_Cluster)) +
      geom_boxplot(outlier.shape = NA, colour = "black", width = 0.6, alpha = 0.75) +
      # Map shape to cell_line
      geom_jitter(
        aes(shape = cell_line),
        width = 0.15,
        size = 6,
        colour = "black",
        alpha = 0.75
      ) +
      geom_text(
        aes(x = 1.5, y = star_y, label = signif_label),
        inherit.aes = FALSE, size = 10, fontface = "bold"
      ) +
      labs(
        title = prog,
        y     = "Z-score of\nProgram Activation",
        x     = NULL,
        shape = "Cell line" # This will be the legend title
      ) +
      scale_fill_manual(values = fill_colors, guide = "none") + # Hide fill legend
      # Apply the manual shape mapping
      scale_shape_manual(name = "Cell line", values = cell_line_shapes) +
      scale_x_discrete(labels = cluster_labels, expand = expansion(add = c(0.5, 0.5))) +
      theme_classic(base_size = 18) +
      theme(
        axis.text.x  = element_text(face = "bold", size = 24, colour = "black"),
        axis.text.y  = element_text(face = "bold", size = 24, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 24, colour = "black"),
        plot.title   = element_text(face = "bold", size = 24, hjust = 0.5),
        # Show legend on the right for individual plots
        legend.position = "right",
        legend.title    = element_text(face = "bold", size = 14),
        legend.text     = element_text(size = 12),
        panel.grid   = element_blank(),
        panel.border = element_rect(linewidth = 1, fill = NA)
      )

    out_file <- sprintf(
      "oncoimmunology_paper/Results/fig5i_5l_tnbc_cl_icd_%s_fry.png",
      gsub(" ", "_", prog)
    )
    ggsave(out_file, plot = p_ind, width = 7, height = 5, dpi = 300, bg = "white")
    message(sprintf("  ✓ Saved: %s", out_file))
  }
}


make_individual_boxplots(sig_df)

# ============================================================================
# STANDALONE PLOT: GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS
# ============================================================================
inflam_plot_df <- inflam_apop_scores %>%
  dplyr::mutate(
    CIMIC_Cluster = factor(cimic_cluster, levels = c("2", "1")),
    signif_label  = inflam_signif
  )
star_y_inflam <- max(inflam_plot_df$Score, na.rm = TRUE) * 0.90

p_inflam <- ggplot(inflam_plot_df, aes(x = CIMIC_Cluster, y = Score, fill = CIMIC_Cluster)) +
  geom_boxplot(outlier.shape = NA, colour = "black", width = 0.6, alpha = 0.75) +
  geom_jitter(aes(shape = cell_line), width = 0.15, size = 6, colour = "black", alpha = 0.75) +
  geom_text(
    aes(x = 1.5, y = star_y_inflam, label = signif_label),
    inherit.aes = FALSE, size = 10, fontface = "bold"
  ) +
  labs(
    title   = "Inflammatory Cell\nApoptotic Process",
    y       = "Z-score of\nProgram Activation",
    x       = NULL
  ) +
   scale_fill_manual(values = fill_colors, guide = "none") + # Hide fill legend
      # Apply the manual shape mapping
      scale_shape_manual(name = "Cell line", values = cell_line_shapes) +
  scale_x_discrete(labels = cluster_labels, expand = expansion(add = c(0.5, 0.5))) +
  shared_theme +
  theme(
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic"),
    plot.title   = element_text(face = "bold", size = 20, hjust = 0.5),
     # Show legend on the right for individual plots
        legend.position = "right",
        legend.text     = element_text(size = 12),
        panel.grid   = element_blank(),
        panel.border = element_rect(linewidth = 1, fill = NA)
  )

ggsave(
  "oncoimmunology_paper/Results/fig5i_5l_tnbc_cl_Inflammatory_Cell_Apoptosis_fry.png",
  plot   = p_inflam,
  width  = 7,
  height = 5,
  dpi    = 300,
  bg     = "white"
)
message("✓ Saved: Results/fig5i_5l_tnbc_cl_Inflammatory_Cell_Apoptosis_fry.png")
message("  [Standalone fry — excluded from main pool]")


# ============================================================================
# EXPORT FRY + BLOCK-ADJUSTED CAMERA STATISTICS: TNBC CELL LINES
# ============================================================================

fig5i_main_stats <- make_fry_camera_table(
  expr_matrix   = expr_matrix,
  gene_sets     = fry_gene_sets,
  design        = design,
  contrast      = contrast,
  cohort_type   = "cell_line",
  block         = block_factor,
  correlation   = dc$consensus.correlation,
  analysis_pool = "Five a priori ICD programs",
  standalone    = FALSE
)

fig5i_inflammatory_set <- list(
  GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS =
    intersect(inflam_apop_genes, rownames(expr_matrix))
)

fig5i_inflammatory_stats <- make_fry_camera_table(
  expr_matrix   = expr_matrix,
  gene_sets     = fig5i_inflammatory_set,
  design        = design,
  contrast      = contrast,
  cohort_type   = "cell_line",
  block         = block_factor,
  correlation   = dc$consensus.correlation,
  analysis_pool = "Standalone prespecified pathway",
  standalone    = TRUE
)

fig5i_all_stats <- dplyr::bind_rows(
  fig5i_main_stats,
  fig5i_inflammatory_stats
)

# ---------------------------------------------------------------------------
# Mean Z-score per program per cluster (the values shown in the box-plots).
#
# Same source as the plots: per-sample program scores produced by
# score_program_z() from the gene-wise z-scored matrix `df_z`. The summary is
# keyed on the raw gene-set identifier used by the statistics table rather than
# on the display label, so the join cannot silently drop rows.
#
# Note: score_program_z() names its cluster column `cimic_cluster` in this
# script (lower case), unlike the NKI_SMC and NEO ICD scripts.
# ---------------------------------------------------------------------------
icd_zscore_long <- dplyr::bind_rows(
  program_scores_all %>%
    dplyr::transmute(
      GeneSet       = Program,
      CIMIC_Cluster = as.character(cimic_cluster),
      Score
    ),
  inflam_apop_scores %>%
    dplyr::transmute(
      GeneSet       = "GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS",
      CIMIC_Cluster = as.character(cimic_cluster),
      Score
    )
)

zscore_summary <- icd_zscore_long %>%
  dplyr::group_by(GeneSet, CIMIC_Cluster) %>%
  dplyr::summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop") %>%
  # TNBC_CL coding: "1" = Fun-CIM, "2" = Dys-CIM
  # (note this is the reverse of the NKI_SMC and NEO cohorts)
  dplyr::mutate(
    ClusterColumn = dplyr::recode(
      CIMIC_Cluster,
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

fig5i_all_stats <- fig5i_all_stats %>%
  dplyr::left_join(zscore_summary, by = "GeneSet")

if (anyNA(fig5i_all_stats$MeanZscore_Cluster2_DysCIM)) {
  warning(
    "Z-score columns are NA for: ",
    paste(
      fig5i_all_stats$GeneSet[
        is.na(fig5i_all_stats$MeanZscore_Cluster2_DysCIM)
      ],
      collapse = ", "
    )
  )
}

data.table::fwrite(
  fig5i_all_stats,
  file = paste0(
    "oncoimmunology_paper/Results/",
    "fig5i_5l_tnbc_cl_icd_programs_fry_camera_statistics.csv"
  ),
  na = "NA"
)

message(
  paste0(
    "Saved: Results/",
    "fig5i_5l_tnbc_cl_icd_programs_fry_camera_statistics.csv"
  )
)