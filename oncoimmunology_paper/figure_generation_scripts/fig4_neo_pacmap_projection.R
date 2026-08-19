# ============================================================================
# fig4_neo_pacmap_projection.R
# ============================================================================
# Purpose: Perform PaCMAP projection and dimensionality reduction on the NEO dataset.
# Inputs: Datasets/NEO/neo_clustered_plot_df.csv
# Outputs: PaCMAP projection plot saved to Results folder.
# Author: Mohammed Gbadamosi

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
neo_data <- fread("oncoimmunology_paper/Datasets/NEO/neo_clustered_plot_df.csv")


# Create output directory if it doesn't exist
# Use explicit relative path from the repository root
dir.create("oncoimmunology_paper/Results", showWarnings = FALSE)

# Check if the required columns exist
required_columns <- c("sample_id", "PACMAP1", "PACMAP2", "cluster_assignments")
if (!all(required_columns %in% colnames(neo_data))) {
  stop("The required columns do not exist in the dataset.")
}
# Extract the cell line identifier (e.g., "NEO42_patient_42" -> "NEO42")
neo_data$base_id <- sub("_.*", "", neo_data$sample_id)

# Generate a color palette using hue_pal from scales
num_clusters <- length(unique(neo_data$cluster_assignments))
cols <- scales::hue_pal()(num_clusters)

# Plot using ggplot, color by cluster
final_cluster_plot_pacmap <- ggplot(
  neo_data,
  aes(
    x = PACMAP1,
    y = PACMAP2,
    label = base_id,  # Patient/sample label column
    color = as.factor(cluster_assignments)
  )
) +
  scale_color_manual(values = cols, name = "Cluster") +
  scale_fill_manual(values = cols, name = "Cluster") +
   ggforce::geom_mark_hull(
    aes(group = as.factor(cluster_assignments),
        fill = as.factor(cluster_assignments), label = NULL),
        concavity = 100,               # Forces convex shape
    expand = unit(5, "mm"),        # Adds outer padding
    radius = unit(5, "mm"),        # Smoothes the sharp corners
    , alpha = 0.2
  )  +
  ggrepel::geom_label_repel(size = 10, fontface = "bold", label.size = 0.3, fill = "white",
                   alpha = 1, box.padding = 0.25, max.overlaps = 100) +
  labs(x = "PacMAP 1", y = "PacMAP 2",
       fill = "Cluster",
       color = "Cluster"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.text = element_text(colour = "black", size = 40),
    axis.title = element_text(colour = "black", size = 45),
    axis.ticks = element_line(linewidth = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

# Save the plot
ggsave("oncoimmunology_paper/Results/fig4_neo_pacmap_projection.png", plot = final_cluster_plot_pacmap, width = 10, height = 9.2, dpi = 300)
