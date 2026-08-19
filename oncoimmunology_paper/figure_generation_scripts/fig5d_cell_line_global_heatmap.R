# ============================================================================
# fig5d_cell_line_global_heatmap.R
# ============================================================================
# Purpose: Create a heatmap visualising global induction differences for the
#          TNBC_CL_Epirubicin cell line dataset.
# Input:  oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epirubicin_master_induction_df.csv
# Output: oncoimmunology_paper/Results/fig5d_cell_line_global_heatmap.png
# ============================================================================

# Load necessary libraries ----------------------------------------------------
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(fastcluster) # Ultra-fast hclust replacement

# Create output directory
output_dir <- "oncoimmunology_paper/Results"
if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE)

# Load the master induction data frame for the cell line dataset
master_induction_df <- fread(
  "oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epirubicin_master_induction_df.csv",
  data.table = FALSE
)

# -------------------------------------------------------------------------
# Filter to protein‑coding genes only (using the pre‑saved list)
# -------------------------------------------------------------------------
protein_genes <- readRDS("oncoimmunology_paper/helper_code/protein_coding_genes.rds")
# Keep only rows where gene_id is in the protein‑coding list
master_induction_df <- master_induction_df %>%
  dplyr::filter(gene_id %in% protein_genes)

# 1. Prepare data matrix (use the same columns as the original heatmap)
mat <- master_induction_df %>%
  dplyr::select(gene_id, cluster1_mean, cluster2_mean) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

colnames(mat) <- c("Dys-CIM\nN = 4", "Fun-CIM\nN = 5")

# 2. Define color palette (same as original)
color_fn <- colorRamp2(
  c(-0.6, 0, 0.6),
  c("#3182bd", "#ffffff", "#e6550d")
)

# 3. Compute fast hierarchical clustering on rows
row_dend <- as.dendrogram(fastcluster::hclust(dist(mat)))

# 4. Define Bottom Annotation Track (Pushes labels cleanly below border)
column_ha <- HeatmapAnnotation(
  labels = anno_text(
    c("Dys-CIM\nN = 4", "Fun-CIM\nN = 5"),
    location = 0.5,
    just = "center",
    rot = 0,
    gp = gpar(fontsize = 20, fontface = "bold")
  ),
  annotation_height = unit(1.5, "cm"),
  show_annotation_name = FALSE
)

# 5. Generate Heatmap and save to PNG
png(
  filename = file.path(output_dir, "fig5d_cell_line_global_heatmap.png"),
  width = 6,
  height = 12,
  units = "in",
  res = 300
)

ht <- Heatmap(
  mat,
  name = "Delta GE",
  col = color_fn,
  width = unit(9.5, "cm"),
  height = unit(18, "cm"),
  cluster_rows = row_dend,
  cluster_columns = FALSE,
  show_row_dend = TRUE,
  row_dend_side = "left",
  row_dend_width = unit(2, "cm"),
  show_column_names = FALSE,
  bottom_annotation = column_ha,
  use_raster = TRUE,
  raster_quality = 3,
  rect_gp = gpar(col = NA),
  border = TRUE,
  border_gp = gpar(col = "black", lwd = 1.2),
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 14, fontface = "plain"),
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_max_height = unit(2.5, "cm"),
  heatmap_legend_param = list(
    title = expression(bold(paste(Delta, " GE"))),
    title_gp = gpar(fontsize = 18, fontface = "bold"),
    labels_gp = gpar(fontsize = 15),
    at = c(-0.5, -0.25, 0, 0.25, 0.5),
    legend_height = unit(4, "cm")
  )
)

draw(ht)

dev.off()
