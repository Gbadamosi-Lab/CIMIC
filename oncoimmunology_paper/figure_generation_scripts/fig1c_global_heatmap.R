# ============================================================================
# fig1c_global_heatmap.R
# ============================================================================
# Purpose: Create a heatmap that visualises the global difference in induction
# patterns across the two CIM trajectories (Dys‑CIM vs Fun‑CIM).
# The underlying data are the per‑gene induction statistics from the master
# induction matrix (Δ sign) – i.e. the sign of the differential expression between
# the two clusters.  The script loads the master induction CSV, extracts the
# differential column, converts it to a sign matrix (‑1, 0, 1) and plots a heatmap.
#
# Input:
#   oncoimmunology_paper/Datasets/NKI_SMC/nki_smc_master_induction_df.csv
# Output:
#   oncoimmunology_paper/Results/fig1c_global_heatmap.png
# Author: Mohammed Gbadamosi (adapted from fig1b_pacmap_projection.R)
# ============================================================================
# Load necessary libraries ----------------------------------------------------
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(fastcluster) # Ultra-fast hclust replacement

# Create output directory
dir.create("oncoimmunology_paper/Results", showWarnings = FALSE)

# Load the master induction data frame
master_induction_df <- fread(
  "oncoimmunology_paper/Datasets/NKI_SMC/nki_smc_master_induction_df.csv",
  data.table = FALSE
)

# 1. Prepare data matrix
mat <- master_induction_df %>%
  dplyr::select(gene_id, cluster1_mean, cluster2_mean) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

colnames(mat) <- c("Dys-CIM\nN = 11", "Fun-CIM\nN = 25")

# 2. Define color palette
color_fn <- colorRamp2(
  c(-0.6, 0, 0.6), 
  c("#3182bd", "#ffffff", "#e6550d")
)

# 3. Compute fast hierarchical clustering on rows
row_dend <- as.dendrogram(fastcluster::hclust(dist(mat)))


# 4. Define Bottom Annotation Track (Pushes labels cleanly below border)
column_ha <- HeatmapAnnotation(
  labels = anno_text(
    c("Dys-CIM\nN = 11", "Fun-CIM\nN = 25"),
    location = 0.5,
    just = "center",
    rot = 0,
    gp = gpar(fontsize = 20, fontface = "bold")
  ),
  annotation_height = unit(1.5, "cm"),
  show_annotation_name = FALSE
)

# 5. Generate Heatmap
png(
  filename = "oncoimmunology_paper/Results/fig1c_global_heatmap.png",
  width = 6,
  height = 12,
  units = "in",
  res = 300
)

ht <- Heatmap(
  mat,
  name = "Delta GE",
  col = color_fn,
  
# Set explicit physical dimensions for a tall vertical rectangle
  width = unit(9.5, "cm"),
  height = unit(18, "cm"),


  # Dendrogram & Fast Clustering
  cluster_rows = row_dend,
  cluster_columns = FALSE,
  show_row_dend = TRUE,
  row_dend_side = "left",
  row_dend_width = unit(2, "cm"),

  # Hide built-in column names (replaced by bottom_annotation)
  show_column_names = FALSE,
  bottom_annotation = column_ha,
  
  # Performance optimizations
  use_raster = TRUE,            # Rasterizes dense rows to avoid SVG/vector lag
  raster_quality = 3,           # Keeps image crisp at 300 DPI
  rect_gp = gpar(col = NA),     # Removes individual grid borders (huge speedup)

  # Heatmap Outer Border
  border = TRUE,
  border_gp = gpar(col = "black", lwd = 1.2),
  
  # Row & Column Labels
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 14, fontface = "plain"),
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_max_height = unit(2.5, "cm"), # Reserves dedicated bottom clearance
  
  # Legend customization
  heatmap_legend_param = list(
    title = expression(bold(paste(Delta, " GE"))),
    title_gp = gpar(fontsize = 18, fontface = "bold"),
    labels_gp = gpar(fontsize = 15),
    at = c(-0.5, -0.25, 0, 0.25, 0.5),
    legend_height = unit(4, "cm")
  )
)

draw(
  ht, 
  # padding = unit(c(5, 5, 8, 5), "mm") # Top, Right, Bottom, Left outer margins
)

dev.off()