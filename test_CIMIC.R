
---

## 2. test_CIMIC.R

```r
# =============================================================================
# CIMIC Test Script
# =============================================================================
# Purpose:  Demonstrate CIMIC pipeline usage on a small simulated dataset
# Input:    Simulated ΔGE matrix (50 samples x 200 genes)
# Expected: Two stable clusters with discriminatory gene sets identified
# Runtime:  ~5-10 minutes (depends on CCP_iter and system)
#
# Usage: Rscript test_CIMIC.R
# =============================================================================

# ---- Load pipeline ----------------------------------------------------------
source("CIMIC_Release_1.0.0.R")

set.seed(2024L)

# ---- Simulate a ΔGE matrix --------------------------------------------------
# 50 paired samples (post - pre log2[TPM+1])
# 200 genes, two underlying biological groups with different response patterns

n_samples <- 50
n_genes   <- 200

# Create two groups with distinct transcriptional response signatures
group1_idx <- 1:25
group2_idx <- 26:50

delta_ge <- matrix(rnorm(n_samples * n_genes, mean = 0, sd = 0.5),
                   nrow = n_samples, ncol = n_genes)

# Inject signal: genes 1-50 upregulated in group 1, downregulated in group 2
delta_ge[group1_idx, 1:50]  <- delta_ge[group1_idx, 1:50]  + 1.5
delta_ge[group2_idx, 1:50]  <- delta_ge[group2_idx, 1:50]  - 1.5

# Genes 51-100: opposite pattern
delta_ge[group1_idx, 51:100] <- delta_ge[group1_idx, 51:100] - 1.0
delta_ge[group2_idx, 51:100] <- delta_ge[group2_idx, 51:100] + 1.0

# Assign row and column names
rownames(delta_ge) <- paste0("delta_PT", sprintf("%03d", 1:n_samples), "_sample")
colnames(delta_ge) <- paste0("GENE", sprintf("%03d", 1:n_genes))

cat("Input matrix dimensions:", nrow(delta_ge), "samples x", ncol(delta_ge), "genes\n")

# ---- Define a small custom gene set -----------------------------------------
# In real usage, leave all_gene_sets = NULL to use the default 19 CIM pathways
# Here we define minimal custom sets to keep runtime short

custom_gene_sets <- list(
  SET_A = paste0("GENE", sprintf("%03d", 1:50)),    # signal genes - group 1 up
  SET_B = paste0("GENE", sprintf("%03d", 51:100)),  # signal genes - group 2 up
  SET_C = paste0("GENE", sprintf("%03d", 101:150))  # noise genes
)

# ---- Create output directory ------------------------------------------------
output_dir <- file.path(getwd(), "CIMIC_test_output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- Run CIMIC --------------------------------------------------------------
cat("\n Running CIMIC pipeline on simulated data...\n")
cat("This may take several minutes.\n\n")

results <- CIM_feature_selection_by_gene_set_pacmap(
  clustering_matrix  = delta_ge,
  all_gene_sets      = custom_gene_sets,   # NULL to use default CIM pathways used in original paper
  clustering_alg     = "hc",
  max_k              = 3,                  # test k = 2 and 3
  CCP_iter           = 100,                # reduced for speed; use 5000 in practice
  adj_pval_thresh    = 0.05,
  max_pipeline_iter  = 5,                  # reduced for speed; use 30 in practice
  seed               = 2024L,
  clustering_metrics = c("pac", "silhouette_combined_avg", "item_cluster_consensus"),
  filter_approach    = c("app_one"),        # run approach 1 only for test
  working_dir        = output_dir,
  verbose            = TRUE
)

# ---- Inspect results --------------------------------------------------------
cat("\n========== CIMIC TEST RESULTS ==========\n")

cat("\nFinal discriminatory gene set (Approach 1):\n")
cat(length(results$iterated_by_gene_sets), "genes retained\n")
cat(head(results$iterated_by_gene_sets, 10), "\n")

cat("\nIteration log:\n")
print(results$iter_log)

cat("\nFinal cluster assignments:\n")
if (!is.null(results$final_df_app_one)) {
  cluster_table <- table(results$final_df_app_one$cluster_assignments)
  print(cluster_table)
  cat("\nExpected: 2 clusters with ~25 samples each\n")
} else {
  cat("No cluster assignments produced — consider relaxing adj_pval_thresh\n")
}

cat("\nOutput files written to:", output_dir, "\n")
cat("\n CIMIC test complete.\n")