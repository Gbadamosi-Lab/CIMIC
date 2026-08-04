# figure2c_new.R
# ============================================================================
# Purpose: Generate overrepresentation analysis barplots for cluster-specific genes
# Inputs: Datasets/nki_smc_master_induction_df.csv
# Outputs: ORA simplified plots saved to Results/ directory
# Author: Mohammed Gbadamosi

library(data.table)
library(msigdbr)
library(clusterProfiler)
library(openxlsx)
library(dplyr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(paletteer)
library(enrichplot)
library(patchwork)
library(tidyverse)
library(stringr)

# Load data
master_induction_df <- fread("Datasets/NKI_SMC/nki_smc_master_induction_df.csv")

# Parameters
optimal_k <- 2
padj_cutoff <- 0.05
differential_threshold <- log2(1.50)
induction_threshold <- log2(1.15)
penetrance_threshold <- 0.20
pen_target_cutoff <- 0.50
pen_other_cutoff <- 0.50
universe_genes <- master_induction_df$gene_id
n_count <- 15

# Create output directory if it doesn't exist
dir.create("Results", showWarnings = FALSE)

# Generate all non-empty, non-full subsets
cluster_subsets <- unlist(
  lapply(1:(optimal_k - 1), function(m) combn(1:optimal_k, m, simplify = FALSE)),
  recursive = FALSE
)

ORA_induction_pic_list <- list()

for (subset in cluster_subsets) {
  cluster_name <- paste0("cluster_", paste(subset, collapse = "_"))
  cluster_gene_set_name <- paste0(cluster_name, "_specific_induced_genes")

  # Work on a temp copy
  df_tmp <- master_induction_df

  # Combined mean for subset
  target_means <- paste0("cluster", subset, "_mean")
  penetrance_name <- paste0("cluster", subset, "_penetrance")
  df_tmp <- df_tmp %>%
    rowwise() %>%
    mutate(combined_mean = mean(c_across(all_of(target_means)), na.rm = TRUE)) %>%
    ungroup()

  # Remaining clusters
  other_clusters <- setdiff(1:optimal_k, subset)
  other_means <- paste0("cluster", other_clusters, "_mean")
  other_pen_cols <- paste0("cluster", other_clusters, "_penetrance")

  # Differential columns: only subset vs other
  cluster_diff_cols <- c()
  for (pos in subset) {
    cluster_diff_cols <- c(
      cluster_diff_cols,
      paste0("cluster", pos, "_cluster", other_clusters, "_differential")
    )
  }
  cluster_diff_cols <- unique(cluster_diff_cols)

  # Filtering
  filtered_genes <- df_tmp %>%
    filter(
      stat_test_adj_p < padj_cutoff,
      combined_mean > induction_threshold,
      if_all(all_of(other_means), ~ . < 0),
      if_all(all_of(cluster_diff_cols), ~ . > differential_threshold),
      if_all(all_of(paste0("cluster", subset, "_penetrance")), ~ . >= pen_target_cutoff),
      if_all(all_of(other_pen_cols), ~ . < pen_other_cutoff)
    )

  assign(cluster_gene_set_name, filtered_genes, envir = .GlobalEnv)
  message(cluster_gene_set_name, " in progress...")

  # ORA analysis
  ORA_induction <- enrichGO(
    gene          = filtered_genes$gene_id,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = padj_cutoff,
    qvalueCutoff  = padj_cutoff,
    universe      = universe_genes,
    keyType       = "SYMBOL",
    minGSSize     = 10,
    maxGSSize     = 500
  )

  ORA_name <- paste0(cluster_gene_set_name, "_ORA_analysis")
  ORA_name_simplified <- paste0("simplified_", cluster_gene_set_name, "_ORA_analysis")

  # Simplify ORA results
  ORA_simplified <- clusterProfiler::simplify(
    ORA_induction,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang"
  )

  if (!is.null(ORA_induction) && nrow(ORA_induction@result) > 0) {
    Check_ORA_induction <- rownames_to_column(ORA_induction@result)
    assign(ORA_name, Check_ORA_induction, envir = .GlobalEnv)

    simp_Check_ORA_induction <- rownames_to_column(ORA_simplified@result)
    assign(ORA_name_simplified, simp_Check_ORA_induction, envir = .GlobalEnv)

    # Top terms for plotting
    top_per_ontology <- ORA_simplified@result %>%
      filter(Count >= 2) %>%
      slice_min(order_by = p.adjust, n = n_count, with_ties = FALSE) %>%
      mutate(Description = str_wrap(Description, width = 20))

    # Create plot
    p <- ggplot(top_per_ontology, aes(x = Description, y = Count)) +
      geom_col(aes(fill = p.adjust), color = "black") +
      coord_flip() +
      labs(
        title = paste0(cluster_name, "\nTop ", n_count, " Enriched Terms"),
        x = NULL,
        y = "Gene Count",
        fill = "Adj_p"
      ) +
      scale_fill_gradient(
        low = "red",
        high = "blue",
        trans = "reverse",
        labels = scales::label_scientific(digits = 2)
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text.y = element_text(size = 18, angle = 0, face = "bold"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18)
      )

    ORA_induction_pic_list[[ORA_name]] <- p

    # Save plot
    plot_filename <- paste0("Results/figure2c_simplified_ORA_", cluster_name, ".png")
    ggsave(plot_filename, plot = p, width = 14, height = 10, dpi = 300)
    message("Saved: ", plot_filename)

    # Save ORA results as CSV
    csv_filename <- paste0("Results/figure2c_simplified_ORA_", cluster_name, ".csv")
    fwrite(simp_Check_ORA_induction, csv_filename)
    message("Saved: ", csv_filename)

  } else if (is.null(ORA_induction)) {
    assign(ORA_name, NULL, envir = .GlobalEnv)
    message(paste0(ORA_name, " is empty - skipping plot generation"))
  }
}

message("All ORA plots generated and saved to Results/ directory")
