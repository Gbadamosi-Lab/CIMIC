# =============================================================================
# MAKE MASTER INDUCTION DATAFRAME DIRECTLY (NEO dataset)
# =============================================================================
# Purpose:
#   Create the NEO master‑induction dataframe directly from the final
#   clustered‑sample dataframe without relying on master_induction_limma_lib.R.
#
# Input:
#   clustered_df (NEO)
#
# Expected structure:
#   - One row per sample
#   - gene_id columns = genes
#   - cluster_assignments = final CIMIC cluster
#   - Optional metadata columns are excluded from the gene matrix
#
# Statistical approach:
#   - LIMMA robust + trend
#   - No duplicateCorrelation because NEO is a patient cohort
#   - H_value = original Wilcoxon W for 2 clusters
#   - p_value = LIMMA moderated t-test
#   - stat_test_adj_p = BH‑adjusted LIMMA p‑value
#
# Output:
#   master_induction_df
# =============================================================================

# -----------------------------------------------------------------------------
# 1. REQUIRED PACKAGES
# -----------------------------------------------------------------------------

library(dplyr)
library(limma)
library(data.table)

# -----------------------------------------------------------------------------
# 2. DEFINE YOUR INPUT DATAFRAME (NEO)
# -----------------------------------------------------------------------------
# Change this ONLY if your clustered dataframe has a different name.

# file location (NEO clustered matrix)
clustered_csv_location <- "D:/OneDrive - University of Florida/Gbadamosi Lab/Mohammed Gbadamosi/Gbadamosi_Lab_GitHub/CIMIC/oncoimmunology_paper/Results/neo_cimic_results/CIM_states_results_hc_pacmap_app_three/clustered_samples_app_three.csv"


# Load the clustered dataframe
clustered_df <- fread(clustered_csv_location) |> as.data.frame(check.names = FALSE)

# The dataframe must contain:
#   cluster_assignments
#   gene columns

# -----------------------------------------------------------------------------
# 4b. DEFINE OUTPUT DIRECTORY
# -----------------------------------------------------------------------------
# Directory where the master induction CSV files will be saved.
# Mirrors the NKI script but points to the NEO dataset folder.
output_dir <- "D:/OneDrive - University of Florida/Gbadamosi Lab/Mohammed Gbadamosi/Gbadamosi_Lab_GitHub/CIMIC/oncoimmunology_paper/Datasets/NEO"

# -----------------------------------------------------------------------------
# 3. DEFINE METADATA COLUMNS
# -----------------------------------------------------------------------------

meta_cols <- c(
    "sample_id",
    "UMAP1",
    "UMAP2",
    "PACMAP1",
    "PACMAP2",
    "PC1",
    "PC2",
    "PC3",
    "base_id",
    "cluster_assignments"
)

# -----------------------------------------------------------------------------
# 4. IDENTIFY GENE COLUMNS
# -----------------------------------------------------------------------------

gene_cols <- setdiff(
    colnames(clustered_df),
    meta_cols
)

cat("\n============================================================\n")
cat("MASTER INDUCTION DATAFRAME (NEO)\n")
cat("============================================================\n")

cat("Samples:", nrow(clustered_df), "\n")
cat("Genes:", length(gene_cols), "\n")
cat("Clusters:", paste(sort(unique(clustered_df$cluster_assignments)), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 5. CREATE CLUSTER FACTOR
# -----------------------------------------------------------------------------

clusters <- factor(
    clustered_df$cluster_assignments
)

cluster_levels <- levels(clusters)
 n_clusters <- length(cluster_levels)

stopifnot(
    n_clusters >= 2
)

# -----------------------------------------------------------------------------
# 6. CREATE GENE MATRIX
# -----------------------------------------------------------------------------
# LIMMA requires rows = genes, columns = samples

M <- as.matrix(
    clustered_df[, gene_cols, drop = FALSE]
)

storage.mode(M) <- "double"

Y <- t(M)
rownames(Y) <- gene_cols

# -----------------------------------------------------------------------------
# 7. RUN LIMMA
# -----------------------------------------------------------------------------

design <- model.matrix(
    ~ clusters
)

fit <- lmFit(
    Y,
    design
)

fit <- eBayes(
    fit,
    robust = TRUE,
    trend = TRUE
)

# -----------------------------------------------------------------------------
# 8. EXTRACT LIMMA P‑VALUES
# -----------------------------------------------------------------------------

coef_cols <- grep(
    "^clusters",
    colnames(design)
)

tt <- topTable(
    fit,
    coef = coef_cols,
    number = Inf,
    sort.by = "none"
)

# Preserve original gene order
tt <- tt[
    match(gene_cols, rownames(tt)),
    ,
    drop = FALSE
]

limma_p <- tt$P.Value

# -----------------------------------------------------------------------------
# 9. CALCULATE ORIGINAL H_VALUE
# -----------------------------------------------------------------------------

H_value <- numeric(length(gene_cols))

for (i in seq_along(gene_cols)) {
    gene_values <- M[, i]
    test <- wilcox.test(
        gene_values ~ clusters,
        exact = FALSE
    )
    H_value[i] <- as.numeric(test$statistic)
}

# -----------------------------------------------------------------------------
# 10. BUILD BASIC MASTER TABLE
# -----------------------------------------------------------------------------

master_induction_df <- data.frame(
    gene_id = gene_cols,
    H_value = H_value,
    p_value = limma_p,
    stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 11. ADD BH‑ADJUSTED LIMMA P‑VALUE
# -----------------------------------------------------------------------------

master_induction_df$stat_test_adj_p <- p.adjust(
    master_induction_df$p_value,
    method = "BH"
)

# -----------------------------------------------------------------------------
# 12. ADD PER‑CLUSTER SUMMARY STATISTICS
# -----------------------------------------------------------------------------

for (k in seq_along(cluster_levels)) {
    idx <- clusters == cluster_levels[k]
    cluster_values <- M[idx, , drop = FALSE]

    # Mean
    master_induction_df[[paste0("cluster", k, "_mean")]] <-
        colMeans(cluster_values, na.rm = TRUE)

    # SD
    master_induction_df[[paste0("cluster", k, "_sd")]] <-
        apply(cluster_values, 2, sd, na.rm = TRUE)

    # Number of non‑missing observations
    cluster_n <- apply(cluster_values, 2, function(x) sum(!is.na(x)))

    # SEM
    master_induction_df[[paste0("cluster", k, "_sem")]] <-
        master_induction_df[[paste0("cluster", k, "_sd")]] / sqrt(cluster_n)

    # Penetrance
    master_induction_df[[paste0("cluster", k, "_penetrance")]] <-
        apply(cluster_values, 2, function(x) mean(x > 0, na.rm = TRUE))
}

# -----------------------------------------------------------------------------
# 13. ADD PAIRWISE DIFFERENTIALS
# -----------------------------------------------------------------------------

if (n_clusters >= 2) {
    for (a in 1:(n_clusters - 1)) {
        for (b in (a + 1):n_clusters) {
            mean_a <- master_induction_df[[paste0("cluster", a, "_mean")]]
            mean_b <- master_induction_df[[paste0("cluster", b, "_mean")]]
            sd_a   <- master_induction_df[[paste0("cluster", a, "_sd")]]
            sd_b   <- master_induction_df[[paste0("cluster", b, "_sd")]]

            # A - B
            master_induction_df[[paste0("cluster", a, "_cluster", b, "_differential")]] <- mean_a - mean_b
            # B - A
            master_induction_df[[paste0("cluster", b, "_cluster", a, "_differential")]] <- mean_b - mean_a

            # Propagated SD
            sd_ab <- sqrt(sd_a^2 + sd_b^2)
            master_induction_df[[paste0("cluster", a, "_cluster", b, "_sd")]] <- sd_ab
            master_induction_df[[paste0("cluster", b, "_cluster", a, "_sd")]] <- sd_ab
        }
    }
}

# -----------------------------------------------------------------------------
# 14. ADD LIMMA TEST PROVENANCE
# -----------------------------------------------------------------------------

if (n_clusters == 2) {
    master_induction_df$limma_test <- "limma_2grp(mod_t)"
} else {
    master_induction_df$limma_test <- "limma_multi(mod_F)"
}

# -----------------------------------------------------------------------------
# 15. REORDER COLUMNS
# -----------------------------------------------------------------------------

master_induction_df <- master_induction_df %>%
    select(
        gene_id,
        H_value,
        p_value,
        stat_test_adj_p,
        everything()
    )

# -----------------------------------------------------------------------------
# 16. BASIC SANITY CHECKS
# -----------------------------------------------------------------------------

stopifnot(
    nrow(master_induction_df) == length(gene_cols)
)
stopifnot(
    !anyDuplicated(master_induction_df$gene_id)
)

# -----------------------------------------------------------------------------
# 17. REPORT RESULTS
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("MASTER INDUCTION COMPLETE (NEO)\n")
cat("============================================================\n")

cat("Genes:", nrow(master_induction_df), "\n")
cat("Columns:", ncol(master_induction_df), "\n")

cat(
    "LIMMA significant genes (BH < 0.05):",
    sum(master_induction_df$stat_test_adj_p < 0.05, na.rm = TRUE),
    "\n"
)

cat(
    "LIMMA significant genes (BH < 0.15):",
    sum(master_induction_df$stat_test_adj_p < 0.15, na.rm = TRUE),
    "\n"
)

# -----------------------------------------------------------------------------
# 19. OPTIONAL: WRITE TO CSV
# -----------------------------------------------------------------------------
# Uncomment when you are satisfied with the result.


    # Write the master induction dataframe and the clustered dataframe to CSV files.
    # Use filenames that reflect the NEO dataset.
    csv_file_master_induction <- file.path(output_dir, "neo_master_induction_df.csv")
    fwrite(master_induction_df, csv_file_master_induction)

    csv_file_cluster_induction <- file.path(output_dir, "neo_clustered_plot_df.csv")
    fwrite(clustered_df, csv_file_cluster_induction)

