# ============================================================================
# fig5c_cell_line_pacmap_projection.R
# ============================================================================
# Purpose: Perform PaCMAP projection plot for the TNBC_CL_Epirubicin cell line dataset.
# Inputs: Datasets/TNBC_CL_Epirubicin/tnbc_cl_epi_clustered_plot_df.csv
# Outputs: PaCMAP projection plot saved to Results folder.
# ============================================================================

# Load necessary libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ggforce)
library(stringr)
library(scales)
library(ggrepel)

# Load the dataset using fread for faster reading
# Use explicit relative path from the repository root
cell_line_data <- fread("oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epi_clustered_plot_df.csv")

# Create output directory if it doesn't exist
# Use explicit relative path from the repository root
dir.create("oncoimmunology_paper/Results", showWarnings = FALSE)

# Check if the required columns exist
required_columns <- c("sample_id", "PACMAP1", "PACMAP2", "cluster_assignments")
if (!all(required_columns %in% colnames(cell_line_data))) {
  stop("The required columns do not exist in the dataset.")
}
# Extract only the cell line name from the sample ID (e.g. "delta_CELLLINE_EPIRUBICIN_1" -> "CELLLINE")
cell_line_data$base_id <- sub("^delta_([^_]+)_.*$", "\\1", cell_line_data$sample_id)

# Define explicit colors: C1 = blue, C2 = salmon (standard manuscript colors)
# Ensure the vector is named by cluster label for correct mapping
cols <- c(
  "1" = "#00bfc4",   # blue (default Matplotlib blue)
  "2" = "#f8766d"    # orange/salmon tone
)

# Plot using ggplot, color by cluster
final_cluster_plot_pacmap <- ggplot(
  cell_line_data,
  aes(
    x = PACMAP1,
    y = PACMAP2,
    label = base_id,  # Patient/sample label column
    color = as.factor(cluster_assignments)
  )
) +
  scale_color_manual(values = cols, name = "Cluster") +
  scale_fill_manual(values = cols, name = "Cluster") +
  # Use an ellipse instead of a hull for a more circular appearance around each cluster
  ggforce::geom_mark_hull(
    aes(group = as.factor(cluster_assignments),
        fill = as.factor(cluster_assignments), label = NULL),
        concavity = 100,               # Forces convex shape
    expand = unit(5, "mm"),        # Adds outer padding
    radius = unit(5, "mm"),        # Smoothes the sharp corners
    , alpha = 0.2
  ) +
  ggrepel::geom_label_repel(
    size = 8, fontface = "bold", label.size = 0.3, fill = "white",
    alpha = 1, box.padding = 0.25, max.overlaps = 100
  ) +
  labs(
    x = "PacMAP 1", y = "PacMAP 2",
    fill = "Cluster",
    color = "Cluster"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.text = element_text(colour = "black", size = 30),
    axis.text.x = element_text(colour = "black", size = 30),
    axis.text.y = element_text(colour = "black", size = 30),
    axis.title = element_text(colour = "black", size = 36),
    axis.ticks = element_line(linewidth = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

# Save the plot
ggsave(
  "oncoimmunology_paper/Results/fig5c_cell_line_pacmap_projection.png",
  plot = final_cluster_plot_pacmap,
  width = 9.6, height = 6, dpi = 300
)
