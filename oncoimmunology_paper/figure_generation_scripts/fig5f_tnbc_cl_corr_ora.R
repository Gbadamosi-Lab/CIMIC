# ============================================================================
# fig5f_tnbc_cl_corr_ora.R
# ============================================================================
# Purpose: Generate over‑representation analysis (ORA) barplots for genes
#          significantly correlated with each CIM trajectory (Dys‑CIM and
#          Fun‑CIM) based on the replicate‑aware correlation results for the
#          TNBC_CL_Epirubicin cell‑line dataset.
#
# Input:  oncoimmunology_paper/Results/tnbc_cl_correlation_analysis.csv
# Output: ORA plots and CSV files saved to oncoimmunology_paper/Results/
# ============================================================================

library(data.table)
library(msigdbr)
library(clusterProfiler)
library(openxlsx)
library(dplyr)
library(org.Hs.eg.db)
library(paletteer)
library(enrichplot)
library(patchwork)
library(tidyverse)
library(stringr)

# ---------------------------------------------------------------------------
# Load correlation results
# ---------------------------------------------------------------------------
correlation_path <- "oncoimmunology_paper/Results/tnbc_cl_correlation_analysis.csv"
cor_df <- fread(correlation_path) |> as.data.frame(check.names = FALSE)

# Parameters
padj_cutoff <- 0.15
r_rb_cutoff <- 0.5

# Universe of genes for enrichment (all genes tested)
universe_genes <- cor_df$gene_id

# Create output directory if it does not exist
output_dir <- "oncoimmunology_paper/Results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# Define gene sets for each trajectory (C1 = Fun‑CIM, C2 = Dys‑CIM)
# ---------------------------------------------------------------------------
# Fun‑CIM corresponds to negative correlation (C1)
fun_genes_df <- cor_df %>%
  filter(padj < padj_cutoff, r_rb <= -r_rb_cutoff)

# Dys‑CIM corresponds to positive correlation (C2)
dys_genes_df <- cor_df %>%
  filter(padj < padj_cutoff, r_rb >= r_rb_cutoff)

# Extract gene symbols
fun_genes <- fun_genes_df$gene_id
dys_genes <- dys_genes_df$gene_id

# ---------------------------------------------------------------------------
# Helper function to run ORA and save results
# ---------------------------------------------------------------------------
run_ora <- function(gene_list, set_name) {
  if (length(gene_list) == 0) {
    message(set_name, " has no genes after filtering – skipping.")
    return(NULL)
  }

  ora_res <- enrichGO(
    gene          = gene_list,
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

  if (is.null(ora_res) || nrow(ora_res@result) == 0) {
    message(set_name, " ORA returned no results.")
    return(NULL)
  }

  # Simplify redundant GO terms
  ora_simplified <- clusterProfiler::simplify(
    ora_res,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang"
  )

  # Prepare data for plotting – top 15 terms with at least 2 genes
  top_terms <- ora_simplified@result %>%
    filter(Count >= 2) %>%
    slice_min(order_by = p.adjust, n = 15, with_ties = FALSE) %>%
    mutate(Description = str_wrap(Description, width = 40))

  p <- ggplot(top_terms, aes(x = Description, y = Count)) +
    geom_col(aes(fill = p.adjust), colour = "black") +
    coord_flip() +
    labs(
      title = paste0(set_name, " – Top Enriched GO BP Terms"),
      x = NULL,
      y = "Gene Count",
      fill = "Adj p"
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
      axis.text.y = element_text(size = 20, colour = "black", face = "bold"),
      axis.text.x = element_text(size = 20, colour = "black", face = "bold"),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18)
    )

  # Save plot
  plot_file <- file.path(output_dir, paste0("fig5f_", set_name, "_ORA.png"))
  ggsave(plot_file, plot = p, width = 14, height = 10, dpi = 300)
  message("Saved plot: ", plot_file)

  # Save ORA table
  csv_file <- file.path(output_dir, paste0("fig5f_", set_name, "_ORA.csv"))
  fwrite(as.data.frame(ora_simplified@result), csv_file)
  message("Saved ORA table: ", csv_file)

  return(p)
}

# ---------------------------------------------------------------------------
# Run ORA for each trajectory
# ---------------------------------------------------------------------------
plot_fun <- run_ora(fun_genes, "FunCIM")
plot_dys <- run_ora(dys_genes, "DysCIM")

message("ORA analysis for TNBC_CL completed.")
