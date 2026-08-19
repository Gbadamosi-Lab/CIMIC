# ============================================================================
# fig5e_correlation_top200_cim_genes_heatmap.R
# ============================================================================
# Purpose:
# Generate a heatmap of the top 200 genes most strongly correlated with each
# CIM trajectory (Dys-CIM and Fun-CIM) for the TNBC_CL_Epirubicin dataset.
#
# Uses replicate-aware correlation results from:
#   tnbc_cl_correlation_analysis.csv
#
# Output:
#   oncoimmunology_paper/Results/fig5e_correlation_top200_cim_genes_heatmap.png
# ============================================================================


library(data.table)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(fastcluster)


# ============================================================================
# Paths
# ============================================================================

clustered_path <-
  "oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epi_clustered_plot_df.csv"

correlation_path <-
  "oncoimmunology_paper/Results/tnbc_cl_correlation_analysis.csv"

gene_set_path <-
  "oncoimmunology_paper/helper_code/cell_line_full_cimic_gene_sets.R"

output_path <-
  "oncoimmunology_paper/Results/fig5e_correlation_top200_cim_genes_heatmap.png"


# ============================================================================
# Load data
# ============================================================================

if (!file.exists(clustered_path)) {
  stop("Missing clustered dataset: ", clustered_path)
}

if (!file.exists(correlation_path)) {
  stop("Missing correlation file: ", correlation_path)
}

if (!file.exists(gene_set_path)) {
  stop("Missing gene set file: ", gene_set_path)
}


clustered_plot_df <-
  fread(clustered_path) |>
  as.data.frame(check.names = FALSE)


cor_df <-
  fread(correlation_path) |>
  as.data.frame(check.names = FALSE)


source(gene_set_path)


if (!exists("cell_line_cimic_gene_sets")) {
  stop("cell_line_cimic_gene_sets not found after sourcing gene set file")
}


cimic_genes <- unique(
  unlist(cell_line_cimic_gene_sets)
)


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
    c(colnames(df))
  )
}


# ============================================================================
# Extract cluster information
# ============================================================================


if (!"cluster_assignments" %in% colnames(clustered_plot_df)) {
  stop("cluster_assignments column missing from clustered_plot_df")
}


cluster_assignments <-
  as.character(clustered_plot_df$cluster_assignments)


# Confirm expected clusters

if (!all(unique(cluster_assignments) %in% c("1", "2"))) {
  stop(
    "Unexpected cluster labels detected: ",
    paste(unique(cluster_assignments), collapse = ", ")
  )
}


# ============================================================================
# Select top correlated genes
# ============================================================================


cutoff <- 0.5
n_top <- 200


sig_corr <- cor_df %>%
  filter(
    padj < 0.05,
    abs(r_rb) >= cutoff
  )


if (nrow(sig_corr) == 0) {
  stop("No significant correlations after filtering")
}


total_dys_genes <-
  sum(sig_corr$trajectory == "Dys-CIM")


total_fun_genes <-
  sum(sig_corr$trajectory == "Fun-CIM")


# ============================================================================
# Correlation direction:
#
# C1 = Fun-CIM  -> negative r_rb
# C2 = Dys-CIM  -> positive r_rb
#
# Therefore:
# Dys-CIM genes = strongest positive correlations
# Fun-CIM genes = strongest negative correlations
# ============================================================================


dys_genes <- sig_corr %>%
  filter(trajectory == "Dys-CIM") %>%
  arrange(desc(r_rb)) %>%
  slice_head(n = n_top) %>%
  pull(gene_id)


fun_genes <- sig_corr %>%
  filter(trajectory == "Fun-CIM") %>%
  arrange(r_rb) %>%
  slice_head(n = n_top) %>%
  pull(gene_id)


if (length(dys_genes) == 0) {
  warning("No Dys-CIM genes identified")
}

if (length(fun_genes) == 0) {
  warning("No Fun-CIM genes identified")
}


# ============================================================================
# User-defined annotation genes
# ============================================================================


fun_hits <- c(
  # Immediate early response
  "FOS",
  "FOSB",
  "EGR1",
  "BTG2",
  "KLF5",

  # IFN / viral mimicry
  "IDO1",
  "IFI35",
  "STING1",
  "STAT2",
  "RTP4",
  "LINC00513",
  "APOBEC3G",

  # Antigen processing
  "PSMB8",
  "SYK",
  "PLCG2",
  "MICB-DT"
)


dys_hits <- c(
  # Proteostasis
  "PSMD14",
  "UBQLN1",
  "UBE3C",
  "USP7",
  "OLA1",

  # Mitochondrial stress
  "TRAP1",
  "VDAC1",
  "CHCHD3",
  "TAMM41",
  "NARS2",
  "ATP5F1A",

  # Metabolism
  "GLUD1",
  "ME2",
  "ALDH18A1",
  "BCAT1",
  "NAPRT",
  "IMPDH1",

  # Autophagy
  "ATG7",
  "VMP1",
  "FAF1",

  # Plasticity
  "WWTR1",
  "ZEB1",
  "ZEB2",
  "NT5E",

  # DNA damage
  "DDB1",
  "SMURF2",
  "PARP8"
)


# ============================================================================
# Build heatmap matrix
# ============================================================================


heatmap_genes <- unique(
  c(dys_genes, fun_genes)
)


gene_cols <- clustered_gene_cols(clustered_plot_df)


genes_present <- intersect(
  heatmap_genes,
  gene_cols
)


if (length(genes_present) == 0) {
  stop("No selected genes found in expression matrix")
}


mat <- clustered_plot_df %>%
  dplyr::select(
    sample_id,
    all_of(genes_present)
  ) %>%
  column_to_rownames("sample_id")


# Remove zero variance genes

mat <- mat[, apply(mat, 2, sd, na.rm = TRUE) > 0, drop = FALSE]


if (ncol(mat) == 0) {
  stop("All genes removed due to zero variance")
}


# ============================================================================
# Sample ordering:
#
# C2 = Dys-CIM (left)
# C1 = Fun-CIM (right)
# ============================================================================


ordered_samples <- clustered_plot_df %>%
  mutate(
    cluster_assignments =
      factor(
        cluster_assignments,
        levels = c("2", "1")
      )
  ) %>%
  arrange(cluster_assignments) %>%
  pull(sample_id)


mat <- mat[
  ordered_samples, ,
  drop = FALSE
]


# Gene x sample matrix

heatmap_mat <- t(as.matrix(mat))


# Remove any remaining NA rows

heatmap_mat <- heatmap_mat[
  apply(
    heatmap_mat,
    1,
    function(x) all(is.finite(x))
  ), ,
  drop = FALSE
]


# END PART 1

# ============================================================================
# Row grouping
# ============================================================================

# Determine trajectory membership
# Priority:
#   Dys-CIM > Fun-CIM
# to avoid ambiguous overlap

gene_group <- case_when(
  rownames(heatmap_mat) %in% dys_genes ~
    paste0(
      "Dys-CIM specific\n(N = ",
      total_dys_genes,
      ")"
    ),
  rownames(heatmap_mat) %in% fun_genes ~
    paste0(
      "Fun-CIM specific\n(N = ",
      total_fun_genes,
      ")"
    ),
  TRUE ~ "Other"
)


gene_group <- factor(
  gene_group,
  levels = c(
    paste0(
      "Fun-CIM specific\n(N = ",
      total_fun_genes,
      ")"
    ),
    paste0(
      "Dys-CIM specific\n(N = ",
      total_dys_genes,
      ")"
    ),
    "Other"
  )
)


# ============================================================================
# Column splitting
# ============================================================================


sample_cluster_df <- clustered_plot_df %>%
  dplyr::select(
    sample_id,
    cluster_assignments
  ) %>%
  dplyr::mutate(
    cluster_assignments =
      as.character(cluster_assignments)
  )


sample_clusters <- sample_cluster_df$cluster_assignments[
  match(
    ordered_samples,
    sample_cluster_df$sample_id
  )
]


col_split <- factor(
  ifelse(
    sample_clusters == "2",
    paste0(
      "Dys-CIM\n(N = ",
      total_dys_genes,
      ")"
    ),
    paste0(
      "Fun-CIM\n(N = ",
      total_fun_genes,
      ")"
    )
  ),
  levels = c(
    paste0(
      "Dys-CIM\n(N = ",
      total_dys_genes,
      ")"
    ),
    paste0(
      "Fun-CIM\n(N = ",
      total_fun_genes,
      ")"
    )
  )
)


# ============================================================================
# Gene highlighting annotations
# ============================================================================


dys_highlight <- intersect(
  dys_hits,
  rownames(heatmap_mat)
)


fun_highlight <- intersect(
  fun_hits,
  rownames(heatmap_mat)
)


dys_idx <- which(
  rownames(heatmap_mat) %in% dys_highlight
)


fun_idx <- which(
  rownames(heatmap_mat) %in% fun_highlight
)


left_ha <- NULL
right_ha <- NULL


# ---------------------------------------------------------------------------
# Dys-CIM labels LEFT (previously on right)
# ---------------------------------------------------------------------------

if (length(dys_idx) > 0) {
  left_ha <- rowAnnotation(
    DysCIM = anno_mark(
      at = dys_idx,
      labels = rownames(heatmap_mat)[dys_idx],
      side = "left",
      labels_gp = gpar(
        fontsize = 24,
        fontface = "bold",
        col = "black"
      ),
      labels_rot = 0,
      link_width = unit(6, "mm")
    )
  )
}


# ---------------------------------------------------------------------------
# Fun-CIM labels RIGHT (previously on left)
# ---------------------------------------------------------------------------

if (length(fun_idx) > 0) {
  right_ha <- rowAnnotation(
    FunCIM = anno_mark(
      at = fun_idx,
      labels = rownames(heatmap_mat)[fun_idx],
      side = "right",
      labels_gp = gpar(
        fontsize = 24,
        fontface = "bold",
        col = "black"
      ),
      labels_rot = 0,
      link_width = unit(6, "mm")
    )
  )
}


# ============================================================================
# Render heatmap
# ============================================================================


png(
  filename = output_path,
  width = 9,
  height = 14,
  units = "in",
  res = 300
)


color_fn <- circlize::colorRamp2(
  c(-2, 0, 2),
  c(
    "#0000FF",
    "#FFFFFF",
    "#FF0000"
  )
)


ht <- Heatmap(
  heatmap_mat,
  name = "Delta GE",
  col = color_fn,


  # clustering
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",


  # splitting

  row_split = gene_group,
  column_split = col_split,
  row_title = NULL,
  row_gap =
    unit(
      3,
      "mm"
    ),
  column_gap =
    unit(
      3,
      "mm"
    ),


  # labels

  show_column_names = FALSE,
  show_row_names = FALSE,
  show_row_dend = FALSE,


  # annotations

  left_annotation = left_ha,
  right_annotation = right_ha,


  # rasterization

  use_raster = TRUE,
  raster_quality = 3,
  rect_gp =
    gpar(
      col = NA
    ),
    column_title_gp =
      gpar(
        fontsize = 30,
        fontface = "bold"
      ),
  heatmap_legend_param = list(
    title =
      expression(
        bold(
          paste(
            Delta,
            "GE"
          )
        )
      ),
      title_gp =
        gpar(
          fontsize = 24,
          fontface = "bold"
        ),
      labels_gp =
        gpar(
          fontsize = 24
        ),
    at = c(-2, -1, 0, 1, 2),
    legend_height =
      unit(
        4,
        "cm"
      )
  )
)


draw(
  ht,
  merge_legend = TRUE
)


dev.off()


message(
  "Heatmap successfully written to: ",
  output_path
)


# ============================================================================
# End script
# ============================================================================