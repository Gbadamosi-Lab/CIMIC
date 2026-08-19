# ============================================================================
# neo_correlation_analysis.R
# ============================================================================
# Purpose: Compute rank‑biserial correlation (r_rb) between gene expression
#          and the two CIM trajectories (Dys‑CIM vs Fun‑CIM) using the
#          clustered expression matrix for the NEO dataset.
#
# Input:
#   oncoimmunology_paper/Datasets/NEO/neo_clustered_plot_df.csv
# Output:
#   oncoimmunology_paper/Results/neo_correlation_analysis.csv
# Author: Mohammed Gbadamosi (adapted from nki_smc_correlation_analysis.R)
# ============================================================================

# Load necessary libraries ----------------------------------------------------
library(data.table)   # fast CSV reading
library(dplyr)        # data manipulation

# ---------------------------------------------------------------------------
# Define explicit repository‑relative paths
# ---------------------------------------------------------------------------
clustered_path <- "oncoimmunology_paper/Datasets/NEO/neo_clustered_plot_df.csv"

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
metadata_cols <- c("sample_id", "UMAP1", "UMAP2", "PACMAP1", "PACMAP2",
                   "PC1", "PC2", "PC3", "base_id", "cluster_assignments")

gene_cols <- setdiff(colnames(clustered_plot_df), metadata_cols)

# Split samples by cluster (1 = Dys‑CIM, 2 = Fun‑CIM)
cluster_vec <- clustered_plot_df$cluster_assignments
group1_idx <- which(cluster_vec == 1)
group2_idx <- which(cluster_vec == 2)

n1 <- length(group1_idx)
n2 <- length(group2_idx)

# ---------------------------------------------------------------------------
# Compute point‑biserial correlation using limma (moderated t statistic)
# ---------------------------------------------------------------------------
library(limma)

# Prepare expression matrix (genes x samples) for limma
M <- as.matrix(clustered_plot_df[, gene_cols, drop = FALSE])
storage.mode(M) <- "double"
Y <- t(M)  # genes x samples

clusterf <- factor(clustered_plot_df$cluster_assignments)
stopifnot("point-biserial needs exactly 2 clusters" = nlevels(clusterf) == 2)
design <- stats::model.matrix(~ clusterf)

fit <- limma::lmFit(Y, design)
fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)
coef <- 2L  # Cluster 2 vs Cluster 1
t_stat <- fit$t[, coef]
dfree <- fit$df.residual + fit$df.prior  # effective degrees of freedom
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
# ---------------------------------------------------------------------------
cor_df <- cor_df %>%
  dplyr::mutate(
    trajectory = dplyr::case_when(
      padj < 0.05 & abs(r_rb) >= cutoff & r_rb < 0 ~ "Dys-CIM",
      padj < 0.05 & abs(r_rb) >= cutoff & r_rb > 0 ~ "Fun-CIM",
      TRUE                                      ~ "NS"
    )
  )

# ---------------------------------------------------------------------------
# Save results to CSV in the Results folder
# ---------------------------------------------------------------------------
output_path <- "oncoimmunology_paper/Results/neo_correlation_analysis.csv"
write.csv(cor_df, file = output_path, row.names = FALSE)

# End of script