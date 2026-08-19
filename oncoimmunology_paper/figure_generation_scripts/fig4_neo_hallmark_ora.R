# fig4_neo_hallmark_ora.R
# ============================================================================
# Purpose: Run Over‑Representation Analysis (ORA) on the Hallmark gene sets
#          from MSigDB using the Dys‑CIM and Fun‑CIM gene lists derived from
#          the NEO dataset correlation analysis (neo_correlation_analysis.R).
#          The script produces bar‑plots of the top enriched Hallmark pathways
#          and saves both the plots and the enrichment tables.
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
# Load correlation results (NEO dataset)
# ---------------------------------------------------------------------------
correlation_path <- "oncoimmunology_paper/Results/neo_correlation_analysis.csv"
cor_df <- fread(correlation_path) |> as.data.frame(check.names = FALSE)

# Parameters for selecting significant genes
padj_cutoff <- 0.15
r_rb_cutoff <- 0.5

# Universe of genes (all tested genes)
universe_genes <- cor_df$gene_id

# ---------------------------------------------------------------------------
# Define gene sets for each trajectory (identical to the GO analysis)
# ---------------------------------------------------------------------------
# Dys‑CIM (C1) – negative correlation
dys_genes_df <- cor_df %>%
  dplyr::filter(trajectory == "Dys-CIM", padj < padj_cutoff, r_rb <= -r_rb_cutoff)

# Fun‑CIM (C2) – positive correlation
fun_genes_df <- cor_df %>%
  dplyr::filter(trajectory == "Fun-CIM", padj < padj_cutoff, r_rb >= r_rb_cutoff)

# Extract gene symbols
dys_genes <- dys_genes_df$gene_id
fun_genes <- fun_genes_df$gene_id

# ---------------------------------------------------------------------------
# Retrieve Hallmark gene sets from MSigDB (category "H")
# ---------------------------------------------------------------------------
hallmark_msigdb <- msigdbr::msigdbr(
  species = "Homo sapiens",
  category = "H"
) %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::distinct()

# Create TERM2GENE data frame required by clusterProfiler::enricher
term2gene <- hallmark_msigdb %>%
  dplyr::rename(term = gs_name, gene = gene_symbol)

# Output directory (same as other analyses)
output_dir <- "oncoimmunology_paper/Results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# Helper function to run ORA on Hallmark sets and save results
# ---------------------------------------------------------------------------
run_hallmark_ora <- function(gene_list, set_name) {
  if (length(gene_list) == 0) {
    message(set_name, " has no genes after filtering – skipping.")
    return(NULL)
  }

  ora_res <- clusterProfiler::enricher(
    gene         = gene_list,
    TERM2GENE    = term2gene,
    pvalueCutoff = padj_cutoff,
    qvalueCutoff = padj_cutoff,
    universe     = universe_genes,
    minGSSize    = 5,
    maxGSSize    = 500,
    pAdjustMethod = "BH"
  )

  if (is.null(ora_res) || nrow(ora_res@result) == 0) {
    message(set_name, " ORA returned no results.")
    return(NULL)
  }

  # Filter to statistically significant pathways (adjusted p <= 0.05)
  sig_terms <- ora_res@result %>%
    dplyr::filter(p.adjust <= padj_cutoff)

  if (nrow(sig_terms) == 0) {
    message(set_name, " ORA returned no significant pathways.")
    return(NULL)
  }

  # Keep top 15 significant pathways (or fewer if not enough) and wrap long names
  top_terms <- sig_terms %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 15) %>%
    dplyr::mutate(Description = stringr::str_wrap(Description, width = 40))

  p <- ggplot(top_terms, aes(x = Description, y = Count)) +
    geom_col(aes(fill = p.adjust), colour = "black") +
    coord_flip() +
    labs(
      title = paste0(set_name, " – Top Enriched Hallmark Pathways"),
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
  plot_file <- file.path(output_dir, paste0("fig4_neo_", set_name, "_Hallmark_ORA.png"))
  ggsave(plot_file, plot = p, width = 14, height = 10, dpi = 300)
  message("Saved plot: ", plot_file)

  # Save enrichment table
  csv_file <- file.path(output_dir, paste0("fig4_neo_", set_name, "_Hallmark_ORA.csv"))
  data.table::fwrite(as.data.frame(ora_res@result), csv_file)
  message("Saved ORA table: ", csv_file)

  return(p)
}

# ---------------------------------------------------------------------------
# Run Hallmark ORA for each trajectory
# ---------------------------------------------------------------------------
plot_dys_hallmark <- run_hallmark_ora(dys_genes, "DysCIM")
plot_fun_hallmark <- run_hallmark_ora(fun_genes, "FunCIM")

message("Hallmark ORA analysis for NEO completed.")
