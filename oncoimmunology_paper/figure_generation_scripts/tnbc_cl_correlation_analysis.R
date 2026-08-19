# ============================================================================
# tnbc_cl_correlation_analysis.R
# ============================================================================
# Purpose: Compute rank‑biserial correlation (r_rb) between gene expression
#          and the two CIM trajectories (Dys‑CIM vs Fun‑CIM) for the
#          TNBC_CL_Epirubicin cell‑line dataset, accounting for biological
#          replicates using limma's duplicateCorrelation.
#
# Input:
#   oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epi_clustered_plot_df.csv
# Output:
#   oncoimmunology_paper/Results/tnbc_cl_correlation_analysis.csv
# ============================================================================

# Load necessary libraries ----------------------------------------------------
library(data.table)   # fast CSV reading
library(dplyr)        # data manipulation
library(limma)        # differential expression with duplicateCorrelation

# ---------------------------------------------------------------------------
# Define explicit repository‑relative paths
# ---------------------------------------------------------------------------
clustered_path <- "oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epi_clustered_plot_df.csv"

# ---------------------------------------------------------------------------
# Load the clustered expression matrix (samples x genes)
# ---------------------------------------------------------------------------
clustered_plot_df <- fread(clustered_path) |> as.data.frame(check.names = FALSE)

# ---------------------------------------------------------------------------
# Prepare vectors for analysis
# ---------------------------------------------------------------------------
cutoff <- 0.5

# Ensure the column containing cluster assignments is present
if (!"cluster_assignments" %in% colnames(clustered_plot_df)) {
  stop("cluster_assignments column not found in clustered data.")
}

# Identify gene columns (numeric expression columns). Exclude known metadata.
metadata_cols <- c(
  "sample_id", "UMAP1", "UMAP2", "PACMAP1", "PACMAP2",
  "PC1", "PC2", "PC3", "base_id", "cluster_assignments"
)

gene_cols <- setdiff(colnames(clustered_plot_df), metadata_cols)

# ---------------------------------------------------------------------------
# Derive block (cell‑line) information from sample IDs
# Sample IDs look like "delta_CELLLINE_EPIRUBICIN_1"
# We keep only the cell‑line part ("CELLLINE") as the block factor.
# ---------------------------------------------------------------------------
clustered_plot_df$cell_line <- sub("^delta_([^_]+)_.*$", "\\1", clustered_plot_df$sample_id)
block_factor <- factor(clustered_plot_df$cell_line)

# ---------------------------------------------------------------------------
# Prepare expression matrix for limma (genes x samples)
# ---------------------------------------------------------------------------
M <- as.matrix(clustered_plot_df[, gene_cols, drop = FALSE])
storage.mode(M) <- "double"
Y <- t(M)  # genes x samples

# ---------------------------------------------------------------------------
# Design matrix for the two CIM trajectories (cluster_assignments)
# ---------------------------------------------------------------------------
clusterf <- factor(clustered_plot_df$cluster_assignments)
if (nlevels(clusterf) != 2) {
  stop("point‑biserial correlation requires exactly 2 clusters.")
}

design <- model.matrix(~ clusterf)

# ---------------------------------------------------------------------------
# Fit replicate‑aware limma model using duplicateCorrelation
# ---------------------------------------------------------------------------
# Estimate the intra‑cell‑line correlation
dc <- duplicateCorrelation(Y, design, block = block_factor)

# Fit the linear model with the estimated correlation
fit <- lmFit(Y, design, block = block_factor, correlation = dc$consensus.correlation)
fit <- eBayes(fit, robust = TRUE, trend = TRUE)

# ---------------------------------------------------------------------------
# Extract statistics (point‑biserial correlation, p‑value)
# ---------------------------------------------------------------------------
coef <- 2L  # Cluster 2 vs Cluster 1

t_stat <- fit$t[, coef]
# Effective degrees of freedom (residual + prior)
dfree <- fit$df.residual + fit$df.prior
r_pb <- t_stat / sqrt(t_stat^2 + dfree)   # point‑biserial correlation
pval <- fit$p.value[, coef]

cor_df <- data.frame(
  gene_id = gene_cols,
  r_rb = as.numeric(r_pb),
  p_value = as.numeric(pval),
  stringsAsFactors = FALSE
)

# Adjust p‑values using Benjamini‑Hochberg (BH) method
cor_df$padj <- p.adjust(cor_df$p_value, method = "BH")

# ---------------------------------------------------------------------------
# Assign trajectory labels based on significance and correlation direction
# NOTE: In this dataset, C1 corresponds to **Fun‑CIM** and C2 to **Dys‑CIM**.
# The coefficient (cluster 2 vs cluster 1) therefore reflects Dys‑CIM
# relative to Fun‑CIM. A positive r_rb indicates higher expression in
# Dys‑CIM, and a negative r_rb indicates higher expression in Fun‑CIM.
# ---------------------------------------------------------------------------
cor_df <- cor_df %>%
  dplyr::mutate(
    trajectory = dplyr::case_when(
      padj < 0.05 & abs(r_rb) >= cutoff & r_rb > 0  ~ "Dys-CIM",
      padj < 0.05 & abs(r_rb) >= cutoff & r_rb < 0  ~ "Fun-CIM",
      TRUE                                          ~ "NS"
    )
  )

# ---------------------------------------------------------------------------
# Save results to CSV in the Results folder
# ---------------------------------------------------------------------------
output_path <- "oncoimmunology_paper/Results/tnbc_cl_correlation_analysis.csv"
write.csv(cor_df, file = output_path, row.names = FALSE)

# End of script
