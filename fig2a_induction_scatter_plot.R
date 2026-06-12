# ============================================================================
# fig2a_induction_scatter_plot.R
# ============================================================================
# Purpose: Create a differential induction scatter plot
# Inputs: Datasets/NKI_SMC/nki_smc_master_induction_df.csv
# Outputs: Scatter plot
# Author: Mohammed Gbadamosi

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggrepel)
library(data.table)

# Load data
master_induction_df <- fread("Datasets/NKI_SMC/nki_smc_master_induction_df.csv")

# Define the two clusters of interest and other parameters
clusters_oi <- c("cluster1", "cluster2")
slope_cutoff <- log2(1.5)
pval_cutoff  <- 0.05
ind_cutoff <- log2(1.50)
penetrance_target_cutoff <- 0.20  # ≥20% induced in target cluster
penetrance_other_cutoff  <- 0.20  # ≤20% induced in other cluster

# Build column names dynamically
diff_col     <- paste0(clusters_oi[1], "_", clusters_oi[2], "_differential")
diff_col_rev <- paste0(clusters_oi[2], "_", clusters_oi[1], "_differential")
pval_col     <- "stat_test_adj_p"
mean_col     <- paste0(clusters_oi[1], "_mean")
mean_col_rev <- paste0(clusters_oi[2], "_mean")
pen_target <- paste0(clusters_oi[1], "_penetrance")
pen_other  <- paste0(clusters_oi[2], "_penetrance")

C1_class <- "Fun-CIM"
C2_class <- "Dys-CIM"

ind_plot_category_1 <- paste0("Induced in ", C1_class, " only")
ind_plot_category_2 <- paste0("Induced in ", C2_class, " only")
ind_plot_category_3 <- paste0("Highly induced in ", C1_class)
ind_plot_category_4 <- paste0("Highly induced in ", C2_class)
ind_plot_category_5 <- paste0("Highly repressed in ", C1_class)
ind_plot_category_6 <- paste0("Highly repressed in ", C2_class)
ind_plot_category_NS <- "NS"

# Prepare plot data
plot_data <- master_induction_df %>%
  dplyr::select(
    gene_id,
    C1 = all_of(mean_col),
    C2 = all_of(mean_col_rev),
    pval = all_of(pval_col),
    diff = all_of(diff_col)
  ) %>%
  dplyr::mutate(
    significant = pval < pval_cutoff & abs(diff) > slope_cutoff,
    category = dplyr::case_when(
      !significant ~ ind_plot_category_NS,
      C1 > 0 & C2 < 0 ~ ind_plot_category_1,
      C1 < 0 & C2 > 0 ~ ind_plot_category_2,
      C1 > 0 & C2 > 0 & (C1 - C2) > slope_cutoff ~ ind_plot_category_3,
      C1 > 0 & C2 > 0 & (C2 - C1) > slope_cutoff ~ ind_plot_category_4,
      C1 < 0 & C2 < 0 & (C2 - C1) > slope_cutoff ~ ind_plot_category_5,
      C1 < 0 & C2 < 0 & (C1 - C2) > slope_cutoff ~ ind_plot_category_6,
      TRUE ~ ind_plot_category_NS
    )
  )

# Count categories
category_counts <- plot_data %>%
  dplyr::count(category) %>%
  dplyr::mutate(label_with_n = paste0(category, " (N = ", n, ")"))

# Build named vector to recode legend labels
n_labels <- setNames(category_counts$label_with_n, category_counts$category)

# Select top genes
n_genes <- 15
available_genes <- unique(master_induction_df$gene_id)
top_up <- plot_data %>%
  dplyr::filter(gene_id %in% available_genes, significant,
                category %in% c(ind_plot_category_1, ind_plot_category_3)) %>%
  dplyr::mutate(dist_from_diag = C1 - C2) %>%
  dplyr::arrange(dplyr::desc(dist_from_diag)) %>%
  dplyr::slice_head(n = n_genes) %>%
  dplyr::pull(gene_id)
top_down <- plot_data %>%
  dplyr::filter(gene_id %in% available_genes, significant,
                category %in% c(ind_plot_category_2, ind_plot_category_4)) %>%
  dplyr::mutate(dist_from_diag = C2 - C1) %>%
  dplyr::arrange(dplyr::desc(dist_from_diag)) %>%
  dplyr::slice_head(n = n_genes) %>%
  dplyr::pull(gene_id)

# Define top genes for labeling
c1 <- c("FOS", "FOSB", "EGR1", "ISG15", "ISG20", "MX1", "IFI35", "DHX58", "XAF1", "OAS1", "OAS2", "OAS3", "P2Y11", "ATF4")
c2 <- c("VIM", "ZEB1", "CDK6", "HMGA1", "LIFE", "XRCC4", "EIF2AK4", "EIF3E", "EIF3H", "HSPA5", "UVRAG", "ATF6")
top_genes <- c(c1, c2)

# Add label column to plot data
plot_data <- plot_data %>%
  dplyr::mutate(label = dplyr::if_else(gene_id %in% top_genes, gene_id, NA_character_))

# Create scatter plot
size_n <- 18
ind_graph <- ggplot(plot_data, aes(x = C2, y = C1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1.5) +
  geom_abline(slope = 1, intercept = slope_cutoff, linetype = "dashed", color = "steelblue", linewidth = 1.5) +
  geom_abline(slope = 1, intercept = -slope_cutoff, linetype = "dashed", color = "firebrick", linewidth = 1.5) +
  geom_point(aes(color = category, shape = significant), size = 1.5, alpha = 0.9) +
  ggrepel::geom_label_repel(
    data = subset(plot_data, !is.na(label)),
    aes(label = label),
    size = 4.5,
    fontface = "bold",
    color = "black",
    fill = "white",
    label.size = 0.5,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.4,
    max.overlaps = Inf
  ) +
  scale_shape_manual(
    name = "Significance",
    values = c("TRUE" = 16, "FALSE" = 17),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "NS")
  ) +
  scale_color_manual(
    name = "Category",
    values = setNames(
      c("black", "#E41A1C", "#377EB8", "#FF7F00", "#984EA3", "#4DAF4A", "#A65628"),
      c(ind_plot_category_NS, ind_plot_category_1, ind_plot_category_2,
        ind_plot_category_3, ind_plot_category_4, ind_plot_category_5,
        ind_plot_category_6)
    ),
    labels = n_labels
  ) +
  labs(
    x = paste0("Dys-CIM mean induction\n(Δlog2[TPM+1])"),
    y = paste0("Fun-CIM mean induction\n(Δlog2[TPM+1])")
  ) +
  ggplot2::coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3)) +
  ggplot2::scale_x_continuous(breaks = seq(-3, 3, 1)) +
  ggplot2::scale_y_continuous(breaks = seq(-3, 3, 1)) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(face = "bold", size = size_n, color = "black"),
    axis.title.y = element_text(face = "bold", size = size_n, color = "black"),
    axis.text.x = element_text(face = "bold", size = size_n, color = "black"),
    axis.text.y = element_text(face = "bold", size = size_n, color = "black"),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0),
    axis.line = element_blank()
  )

# Print the plot
print(ind_graph)

# Save the plot
ggsave("fig2a_induction_scatter_plot.png", plot = ind_graph, width = 8, height = 8, dpi = 300)