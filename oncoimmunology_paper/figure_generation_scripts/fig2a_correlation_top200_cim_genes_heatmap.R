# ============================================================================
# fig2a_correlation_top100_cim_genes_heatmap.R
# ============================================================================
# Purpose: Generate a heatmap of the top 100 genes most strongly correlated
#          with each CIM trajectory (Dys‑CIM and Fun‑CIM) using the results of
#          the rank‑biserial correlation analysis (correlation_analysis.R).
#          The heatmap visualises the expression patterns of these genes across
#          samples, split by trajectory.
#
# Input:
#   oncoimmunology_paper/Datasets/NKI_SMC/nki_smc_combine_clustered_plot_df.csv
#   oncoimmunology_paper/Results/nki_smc_correlation_analysis.csv
# Output:
#   oncoimmunology_paper/Results/fig2a_correlation_top100_cim_genes_heatmap.png
# Author: Mohammed Gbadamosi (adapted from fig2a_induction_top100_cim_genes_heatmap.R)
# ============================================================================

library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(fastcluster)

# ---------------------------------------------------------------------------
# Define explicit repository‑relative paths
# ---------------------------------------------------------------------------
clustered_path <- "oncoimmunology_paper/Datasets/NKI_SMC/nki_smc_combine_clustered_plot_df.csv"
correlation_path <- "oncoimmunology_paper/Results/nki_smc_correlation_analysis.csv"
# Gene‑set definitions (used for annotation labeling)
gene_set_path  <- "oncoimmunology_paper/helper_code/human_full_cimic_gene_sets.R"

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
clustered_plot_df <- fread(clustered_path) |> as.data.frame(check.names = FALSE)
cor_df <- fread(correlation_path) |> as.data.frame(check.names = FALSE)
# Load gene set script and create bundle for downstream use
source(gene_set_path)
bundle <- list(msig_gene_sets = msig_gene_sets)
# Genes present in the CIM gene‑set collection
cimic_genes <- unique(unlist(bundle$msig_gene_sets))

cutoff <- 0.6



# ============================================================================
# Helper functions
# ============================================================================

clustered_gene_cols <- function(df) {
  exclude_cols <- c(
    "sample_id",
    "cluster_assignments",
    "trajectory",
    "cell_line",
    "cell_line_id"
  )
  intersect(
    setdiff(colnames(df), exclude_cols),
    colnames(df)
  )
}


# ---------------------------------------------------------------------------
# Select top correlated genes per trajectory
# ---------------------------------------------------------------------------
# Keep only genes that are significantly associated (padj < 0.05) and have a
# correlation magnitude >= 0.6 (as defined in correlation_analysis.R).
sig_corr <- cor_df %>%
  dplyr::filter(padj < 0.05, abs(r_rb) >= cutoff)

n_top <- 200
# Dys‑CIM corresponds to trajectory "Dys-CIM" (negative correlation)
# Fun‑CIM corresponds to trajectory "Fun-CIM" (positive correlation)

# Total number of significant genes per trajectory (before selecting top N)
total_dys_genes <- sum(sig_corr$trajectory == "Dys-CIM")
total_fun_genes <- sum(sig_corr$trajectory == "Fun-CIM")

dys_genes <- sig_corr %>%
  dplyr::filter(trajectory == "Dys-CIM") %>%
  dplyr::arrange(r_rb) %>%               # most negative first
  dplyr::slice_head(n = n_top) %>%
  dplyr::pull(gene_id)

fun_genes <- sig_corr %>%
  dplyr::filter(trajectory == "Fun-CIM") %>%
  dplyr::arrange(desc(r_rb)) %>%        # most positive first
  dplyr::slice_head(n = n_top) %>%
  dplyr::pull(gene_id)

# ---------------------------------------------------------------------------
# Define specific genes to label for each trajectory (provided by user)
# ---------------------------------------------------------------------------
fun_hits <- c(
  "CIITA", "CCR2", "TLR4", "JAK3", "CD69", "CD28", "CD40", "BTN3A3",
  "IL1R1", "LCP2", "PTGER4", "CARD8", "CASP10", "NLRP1", "ATM",
  "STAT5B", "LYZ", "IRAK3", "SCARF1", "FCRL5", "IGLL5", "SLAMF1",
  "SLAMF7", "PTPN22", "TNFRSF17", "IRF4"
)

dys_hits <- c(
  # Mitochondrial proteostasis
  "CLPP",
  "PHB1",
  "PHB2",
  "TIMM23",
  "TIMM50",
  "PAM16",
  "C1QBP",

  # Mitochondrial translation / mitoribosome
  "MALSU1",
  "MRPL14",
  "MRPS34",

  # Oxidative phosphorylation
  "ATP5IF1",
  "ATP5F1C",

  # Proteasome / protein quality control
  "PSMC3",
  "PSMA7",

  # Oxidative stress adaptation
  "PRDX5",

  # DNA damage tolerance
  "RAD51",
  "FEN1"
)


heatmap_genes <- unique(c(dys_genes, fun_genes))
# Identify which of these genes are present in the clustered expression matrix
gene_cols <- clustered_gene_cols(clustered_plot_df)
genes_present <- intersect(heatmap_genes, gene_cols)

# ---------------------------------------------------------------------------
# Prepare expression matrix (genes x samples)
# ---------------------------------------------------------------------------
# Order samples so that Dys‑CIM (cluster 2) appears first, matching the visual
# convention used in other heatmaps.
ordered_samples <- clustered_plot_df %>%
  dplyr::arrange(desc(cluster_assignments)) %>%
  dplyr::pull(sample_id)

mat <- clustered_plot_df %>%
  dplyr::select(sample_id, dplyr::all_of(genes_present)) %>%
  tibble::column_to_rownames("sample_id")

mat <- mat[ordered_samples, , drop = FALSE]
# Transpose and z‑score standardize per gene across samples
heatmap_mat <- t(as.matrix(mat))

# Mapping clusters to trajectory display labels
cluster_map <- c("1" = "Dys-CIM", "2" = "Fun-CIM")
# Determine cluster assignment for each ordered sample
sample_clusters <- clustered_plot_df$cluster_assignments[match(ordered_samples, clustered_plot_df$sample_id)]

# Compute number of genes selected for each trajectory (for annotation)
count_dys_genes <- length(dys_genes)
count_fun_genes <- length(fun_genes)

# Create column split factor with gene counts in titles
col_split <- factor(
  ifelse(
    as.character(sample_clusters) == "1",
    paste0("Dys-CIM\n(N = ", total_dys_genes, ")"),
    paste0("Fun-CIM\n(N = ", total_fun_genes, ")")
  ),
  levels = c(
    paste0("Dys-CIM\n(N = ", total_dys_genes, ")"),
    paste0("Fun-CIM\n(N = ", total_fun_genes, ")")
  )
)

# Row split categories for annotation
# Row split categories now include gene counts for each trajectory
# Reverse the order of row splits so that Dys‑CIM specific genes appear on the right side of the heatmap.
gene_group <- factor(
  dplyr::case_when(
    rownames(heatmap_mat) %in% dys_genes ~ paste0("Dys-CIM specific\n(N = ", total_dys_genes, ")"),
    rownames(heatmap_mat) %in% fun_genes ~ paste0("Fun-CIM specific\n(N = ", total_fun_genes, ")"),
    TRUE                                 ~ "Other"
  ),
  levels = c(
    paste0("Fun-CIM specific\n(N = ", total_fun_genes, ")"),
    paste0("Dys-CIM specific\n(N = ", total_dys_genes, ")")
  )
)

# ---------------------------------------------------------------------------
# Gene Annotations & Labeling
# Dys-CIM labels on LEFT
# Fun-CIM labels on RIGHT
# ---------------------------------------------------------------------------

dys_highlight <- intersect(dys_hits, rownames(heatmap_mat))
fun_highlight <- intersect(fun_hits, rownames(heatmap_mat))

dys_idx <- which(rownames(heatmap_mat) %in% dys_highlight)
fun_idx <- which(rownames(heatmap_mat) %in% fun_highlight)

left_ha <- NULL
right_ha <- NULL

if (length(dys_idx) > 0) {
  left_ha <- rowAnnotation(
    DysCIM = anno_mark(
      at = dys_idx,
      labels = rownames(heatmap_mat)[dys_idx],
      side = "left",
      labels_gp = gpar(
        fontsize = 14,
        fontface = "bold",
        col = "black"
      ),
      labels_rot = 0,
      link_width = unit(6, "mm")
    )
  )
}

if (length(fun_idx) > 0) {
  right_ha <- rowAnnotation(
    FunCIM = anno_mark(
      at = fun_idx,
      labels = rownames(heatmap_mat)[fun_idx],
      side = "right",
      labels_gp = gpar(
        fontsize = 14,
        fontface = "bold",
        col = "black"
      ),
      labels_rot = 0,
      link_width = unit(6, "mm")
    )
  )
}

# ---------------------------------------------------------------------------
# Render and Save Heatmap
# ---------------------------------------------------------------------------
png(
  filename = "oncoimmunology_paper/Results/fig2a_correlation_top200_cim_genes_heatmap.png",
  width = 9,
  height = 14,
  units = "in",
  res = 300
)

color_fn <- circlize::colorRamp2(c(-2, 0, 2), c("#0000FF", "#FFFFFF", "#FF0000"))

ht <- ComplexHeatmap::Heatmap(
  heatmap_mat,
  name = "Delta GE",
  col = color_fn,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  row_split = gene_group,
  column_split = col_split,
  row_title = NULL,
  row_gap = grid::unit(3, "mm"),
  column_gap = grid::unit(3, "mm"),
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_row_dend = FALSE,
  left_annotation = left_ha,
right_annotation = right_ha,
  use_raster = TRUE,
  raster_quality = 3,
  rect_gp = grid::gpar(col = NA),
  column_title_gp = grid::gpar(fontsize = 25, fontface = "bold"),
  heatmap_legend_param = list(
    title = expression(bold(paste(Delta, "GE"))),
    title_gp = grid::gpar(fontsize = 18, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 15),
    at = c(-2, -1, 0, 1, 2),
    legend_height = grid::unit(4, "cm")
  )
)

ComplexHeatmap::draw(ht, merge_legend = TRUE)

dev.off()

# End of script