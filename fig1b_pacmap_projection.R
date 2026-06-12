# ============================================================================
# fig1b_pacmap_projection.R
# ============================================================================
# Purpose: Perform PaCMAP projection and dimensionality reduction on the NKI SMC dataset.
# Inputs: Datasets/NKI_SMC/nki_smc_combine_clustered_plot_df.csv
# Outputs: PaCMAP projection plot
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
nki_smc_data <- fread("Datasets/NKI_SMC/nki_smc_combine_clustered_plot_df.csv")

# Check if the required columns exist
required_columns <- c("sample_id", "PACMAP1", "PACMAP2", "cluster_assignments")
if (!all(required_columns %in% colnames(nki_smc_data))) {
  stop("The required columns do not exist in the dataset.")
}
# Remove last "_delta" from sample ID
nki_smc_data$base_id <- gsub("_delta$", "", nki_smc_data$sample_id)

# Generate a color palette using hue_pal from scales
num_clusters <- length(unique(nki_smc_data$cluster_assignments))
cols <- scales::hue_pal()(num_clusters)

# Plot using GGPlot, color by cluster
final_cluster_plot_pacmap <- ggplot(
  nki_smc_data,
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
    concavity = 2, expand = unit(2, "mm"), alpha = 0.2
  ) +
  ggrepel::geom_label_repel(size = 4, fontface = "bold", label.size = 0.3, fill = "white",
                   alpha = 1, box.padding = 0.25, max.overlaps = 100) +
  labs(x = "PacMAP 1", y = "PacMAP 2",
       fill = "Cluster",
       color = "Cluster"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    axis.text = element_text(colour = "black", size = 24),
    axis.text.x = element_text(colour = "black", size = 24),
    axis.text.y = element_text(colour = "black", size = 24),
    axis.title = element_text(colour = "black", size = 30),
    axis.ticks = element_line(linewidth = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

# Save the plot
ggsave("fig1b_pacmap_projection.png", plot = final_cluster_plot_pacmap, width = 8, height = 8, dpi = 300)