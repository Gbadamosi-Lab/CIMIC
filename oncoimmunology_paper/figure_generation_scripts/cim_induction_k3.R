
# =============================================================================
# MAKE MASTER INDUCTION DATAFRAME DIRECTLY  —  k = 3
# =============================================================================
# Purpose:
#   Create the NKI master-induction dataframe by merging the cluster
#   assignment file with the expression delta matrix, then running
#   LIMMA + Kruskal-Wallis for k=3 clusters.
#
# Statistical approach:
#   - LIMMA robust + trend  |  omnibus moderated F-test (k=3)
#   - H_value = Kruskal-Wallis H statistic
#   - stat_test_adj_p = BH-adjusted LIMMA p-value
#
# CIM state annotation:
#   Cluster 1 → Fun-CIM
#   Cluster 2 → Dys-CIM
#   Cluster 3 → Fun-CIM
#
# Inputs:
#   1. nki_smc_combine_k3.csv                  — sample_id + cluster_assignment
#   2. nki_smc_combine_initial_clustering_mat.csv — sample_id + gene deltas
#
# All outputs → Results/nki_smc_k3/
# =============================================================================

# -----------------------------------------------------------------------------
# 1. REQUIRED PACKAGES
# -----------------------------------------------------------------------------
library(dplyr)
library(limma)
library(data.table)
library(ggplot2)
library(ggvenn)        # install.packages("ggvenn")
library(ggalluvial)    # install.packages("ggalluvial")

# Explicitly resolve namespace conflicts from data.table
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
# -----------------------------------------------------------------------------
# 2. DEFINE YOUR INPUT FILES AND MERGE
# -----------------------------------------------------------------------------
# File 1: cluster assignments (sample_id + cluster_assignment only)
clustered_csv_location <- paste0(
  "C:/Users/lopesdelima.i/OneDrive - University of Florida/",
  "Gbadamosi Lab/Mohammed Gbadamosi/Gbadamosi_Lab_GitHub/CIMIC/",
  "oncoimmunology_paper/Datasets/NKI_SMC/nki_smc_combine_k3.csv"
)

# File 2: expression delta matrix (sample_id + gene columns)
expression_csv_location <- paste0(
  "C:/Users/lopesdelima.i/OneDrive - University of Florida/",
  "Gbadamosi Lab/Mohammed Gbadamosi/Gbadamosi_Lab_GitHub/CIMIC/",
  "oncoimmunology_paper/Datasets/NKI_SMC/",
  "nki_smc_combine_initial_clustering_mat.csv"
)

# Load both files
cluster_df    <- fread(clustered_csv_location)  |> as.data.frame(check.names = FALSE)
expression_df <- fread(expression_csv_location) |> as.data.frame(check.names = FALSE)

# Ensure the first column of the expression matrix is named sample_id
colnames(expression_df)[1] <- "sample_id"

# Merge on sample_id — inner join keeps only samples present in both files [1]
clustered_df <- inner_join(
  cluster_df,
  expression_df,
  by = "sample_id"
)

cat("Samples after merge:", nrow(clustered_df), "\n")
cat("Columns after merge:", ncol(clustered_df), "\n")

# Verify no samples were lost
if (nrow(clustered_df) != nrow(cluster_df)) {
  warning(
    "Sample count changed after merge! ",
    "cluster file: ", nrow(cluster_df), " | ",
    "merged: ",       nrow(clustered_df),
    " — check for sample_id mismatches."
  )
}

# -----------------------------------------------------------------------------
# 3. OUTPUT DIRECTORY  —  create nki_smc_k3 subfolder if it does not exist
# -----------------------------------------------------------------------------
output_dir <- paste0(
  "C:/Users/lopesdelima.i/OneDrive - University of Florida/",
  "Gbadamosi Lab/Mohammed Gbadamosi/Gbadamosi_Lab_GitHub/CIMIC/",
  "oncoimmunology_paper/Results/nki_smc_k3"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
} else {
  cat("Output directory already exists:", output_dir, "\n")
}

# -----------------------------------------------------------------------------
# 4. DEFINE METADATA COLUMNS
# -----------------------------------------------------------------------------
# Only sample_id and cluster_assignment are metadata.
# All remaining columns after the merge are gene expression deltas.
meta_cols <- c(
  "sample_id",
  "cluster_assignment"
)

# -----------------------------------------------------------------------------
# 5. IDENTIFY GENE COLUMNS
# -----------------------------------------------------------------------------
gene_cols <- setdiff(
  colnames(clustered_df),
  meta_cols
)

cat("\n============================================================\n")
cat("MASTER INDUCTION DATAFRAME\n")
cat("============================================================\n")
cat("Samples:", nrow(clustered_df), "\n")
cat("Genes:",   length(gene_cols),  "\n")
cat("Clusters:",
    paste(sort(unique(clustered_df$cluster_assignment)), collapse = ", "),
    "\n")

# -----------------------------------------------------------------------------
# 6. CREATE CLUSTER FACTOR
# -----------------------------------------------------------------------------
clusters <- factor(
  clustered_df$cluster_assignment
)
cluster_levels <- levels(clusters)
n_clusters     <- length(cluster_levels)
stopifnot(n_clusters >= 2)

# -----------------------------------------------------------------------------
# 7. CREATE GENE MATRIX
# -----------------------------------------------------------------------------
# LIMMA requires rows = genes, columns = samples
M <- as.matrix(
  clustered_df[, gene_cols, drop = FALSE]
)
storage.mode(M) <- "double"
Y <- t(M)
rownames(Y) <- gene_cols

# -----------------------------------------------------------------------------
# 8. RUN LIMMA
# -----------------------------------------------------------------------------
design <- model.matrix(~ clusters)
fit    <- lmFit(Y, design)
fit    <- eBayes(fit, robust = TRUE, trend = TRUE)

# -----------------------------------------------------------------------------
# 9. EXTRACT LIMMA P-VALUES  (omnibus moderated F-test for k=3)
# -----------------------------------------------------------------------------
# coef_cols captures both contrast columns for k=3,
# so topTable returns the omnibus moderated F-test p-value.
coef_cols <- grep("^clusters", colnames(design))

tt <- topTable(
  fit,
  coef    = coef_cols,
  number  = Inf,
  sort.by = "none"
)

tt      <- tt[match(gene_cols, rownames(tt)), , drop = FALSE]
limma_p <- tt$P.Value

# -----------------------------------------------------------------------------
# 10. CALCULATE H_VALUE  (Kruskal-Wallis H, appropriate for k=3)
# -----------------------------------------------------------------------------
# Replaces Wilcoxon W (2-group only) with Kruskal-Wallis H.
# Column name H_value is preserved for downstream compatibility.
H_value <- numeric(length(gene_cols))
for (i in seq_along(gene_cols)) {
  gene_values <- M[, i]
  test        <- kruskal.test(gene_values ~ clusters)
  H_value[i]  <- as.numeric(test$statistic)
}

# -----------------------------------------------------------------------------
# 11. BUILD BASIC MASTER TABLE
# -----------------------------------------------------------------------------
master_induction_df <- data.frame(
  gene_id = gene_cols,
  H_value = H_value,
  p_value = limma_p,
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 12. ADD BH-ADJUSTED LIMMA P-VALUE
# -----------------------------------------------------------------------------
master_induction_df$stat_test_adj_p <- p.adjust(
  master_induction_df$p_value,
  method = "BH"
)

# -----------------------------------------------------------------------------
# 13. ADD PER-CLUSTER SUMMARY STATISTICS
# -----------------------------------------------------------------------------
# For each cluster: mean, SD, SEM, penetrance (proportion > 0)
for (k in seq_along(cluster_levels)) {
  cluster_name   <- cluster_levels[k]
  idx            <- clusters == cluster_name
  cluster_values <- M[idx, , drop = FALSE]
  
  master_induction_df[[paste0("cluster", k, "_mean")]] <-
    colMeans(cluster_values, na.rm = TRUE)
  
  master_induction_df[[paste0("cluster", k, "_sd")]] <-
    apply(cluster_values, 2, sd, na.rm = TRUE)
  
  cluster_n <- apply(cluster_values, 2, function(x) sum(!is.na(x)))
  
  master_induction_df[[paste0("cluster", k, "_sem")]] <-
    master_induction_df[[paste0("cluster", k, "_sd")]] / sqrt(cluster_n)
  
  master_induction_df[[paste0("cluster", k, "_penetrance")]] <-
    apply(cluster_values, 2, function(x) mean(x > 0, na.rm = TRUE))
}

# -----------------------------------------------------------------------------
# 14. ADD PAIRWISE DIFFERENTIALS
# -----------------------------------------------------------------------------
# k=3 produces all three pairs + reverses:
#   C1vC2, C1vC3, C2vC3  (and C2vC1, C3vC1, C3vC2)
if (n_clusters >= 2) {
  for (a in 1:(n_clusters - 1)) {
    for (b in (a + 1):n_clusters) {
      mean_a <- master_induction_df[[paste0("cluster", a, "_mean")]]
      mean_b <- master_induction_df[[paste0("cluster", b, "_mean")]]
      sd_a   <- master_induction_df[[paste0("cluster", a, "_sd")]]
      sd_b   <- master_induction_df[[paste0("cluster", b, "_sd")]]
      
      master_induction_df[[paste0("cluster", a, "_cluster", b, "_differential")]] <-
        mean_a - mean_b
      master_induction_df[[paste0("cluster", b, "_cluster", a, "_differential")]] <-
        mean_b - mean_a
      
      sd_ab <- sqrt(sd_a^2 + sd_b^2)
      master_induction_df[[paste0("cluster", a, "_cluster", b, "_sd")]] <- sd_ab
      master_induction_df[[paste0("cluster", b, "_cluster", a, "_sd")]] <- sd_ab
    }
  }
}

# -----------------------------------------------------------------------------
# 15. ADD LIMMA TEST PROVENANCE
# -----------------------------------------------------------------------------
if (n_clusters == 2) {
  master_induction_df$limma_test <- "limma_2grp(mod_t)"
} else {
  master_induction_df$limma_test <- "limma_multi(mod_F)"
}

# -----------------------------------------------------------------------------
# 16. REORDER COLUMNS
# -----------------------------------------------------------------------------
master_induction_df <- master_induction_df %>%
  select(gene_id, H_value, p_value, stat_test_adj_p, everything())

# -----------------------------------------------------------------------------
# 17. BASIC SANITY CHECKS
# -----------------------------------------------------------------------------
stopifnot(nrow(master_induction_df) == length(gene_cols))
stopifnot(!anyDuplicated(master_induction_df$gene_id))

# -----------------------------------------------------------------------------
# 18. REPORT RESULTS
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("MASTER INDUCTION COMPLETE\n")
cat("============================================================\n")
cat("Genes:",   nrow(master_induction_df), "\n")
cat("Columns:", ncol(master_induction_df), "\n")
cat(
  "LIMMA significant genes (BH < 0.05):",
  sum(master_induction_df$stat_test_adj_p < 0.05, na.rm = TRUE), "\n"
)
cat(
  "LIMMA significant genes (BH < 0.15):",
  sum(master_induction_df$stat_test_adj_p < 0.15, na.rm = TRUE), "\n"
)

# -----------------------------------------------------------------------------
# 19. SAVE MASTER INDUCTION AND CLUSTERED DATA
# -----------------------------------------------------------------------------
csv_file_master_induction <- file.path(
  output_dir, "nki_smc_master_induction_df.csv"
)
fwrite(master_induction_df, csv_file_master_induction)
cat("Master induction saved:", csv_file_master_induction, "\n")

csv_file_cluster_induction <- file.path(
  output_dir, "nki_smc_combine_clustered_plot_df.csv"
)
fwrite(clustered_df, csv_file_cluster_induction)
cat("Clustered plot df saved:", csv_file_cluster_induction, "\n")

# =============================================================================
# 20. VENN DIAGRAMS
# =============================================================================
# Four Venn diagrams are produced:
#
#   Venn 1 — Global DE (direction-independent):
#             All genes significantly different between any pair,
#             regardless of which cluster is higher.
#             abs(differential) >= DIFF_THRESHOLD
#
#   Venn 2 — Fun-CIM positive induction:
#             Genes positively induced in C1 (Fun-CIM) vs C2 (Dys-CIM)
#             and in C3 (Fun-CIM) vs C2 (Dys-CIM).
#             Overlap = core Fun-CIM induction signature.
#
#   Venn 3 — Dys-CIM specificity:
#             Genes positively induced in C2 (Dys-CIM) vs C1 and vs C3.
#             Overlap = robustly Dys-CIM-specific genes.
#
#   Venn 4 — Overall positive differential induction by cluster:
#             For each cluster, genes significantly more induced than at least
#             one of the other two clusters (union of either comparison).
#
# Thresholds:
#   ADJ_P_THRESHOLD = 0.05  — BH-adjusted LIMMA omnibus F-test
#   DIFF_THRESHOLD  = 0.15  — adjust to your assay's meaningful delta
#
# CIM annotation [1]:
#   Cluster 1 → Fun-CIM
#   Cluster 2 → Dys-CIM
#   Cluster 3 → Fun-CIM
# =============================================================================

# -----------------------------------------------------------------------------
# 20a. THRESHOLDS  —  adjust here if needed
# -----------------------------------------------------------------------------
ADJ_P_THRESHOLD <- 0.05
DIFF_THRESHOLD  <- 0.15

# -----------------------------------------------------------------------------
# 20b. HELPER FUNCTIONS
# -----------------------------------------------------------------------------
# Build a membership dataframe: all statistics preserved + TRUE/FALSE per set
build_membership_df <- function(master_df, venn_sets) {
  all_genes  <- unique(unlist(venn_sets))
  membership <- master_df[
    master_df$gene_id %in% all_genes, , drop = FALSE
  ]
  for (set_label in names(venn_sets)) {
    col_name <- gsub("[ \n]", "_", set_label)
    membership[[col_name]] <-
      membership$gene_id %in% venn_sets[[set_label]]
  }
  return(membership)
}

# Draw a styled ggvenn and return the plot object
draw_venn <- function(venn_sets, title_label, subtitle_label,
                      caption_label, fill_colors) {
  venn_plot <- ggvenn(
    venn_sets,
    fill_color      = fill_colors,
    stroke_size     = 0.5,
    set_name_size   = 4.5,
    text_size       = 5,
    show_percentage = FALSE
  )
  
  # Make both set names and region counts bold.
  for (i in seq_along(venn_plot$layers)) {
    if (inherits(venn_plot$layers[[i]]$geom, "GeomText")) {
      venn_plot$layers[[i]]$aes_params$fontface <- "bold"
    }
  }
  
  venn_plot +
    labs(title = NULL, subtitle = NULL, caption = NULL) +
    theme_void() +
    theme(plot.margin = margin(4, 4, 4, 4))
}

# -----------------------------------------------------------------------------
# 20c. VENN 1 — GLOBAL DE  (direction-independent, all three pairs)
# -----------------------------------------------------------------------------
venn_sets_global <- list(
  "C1 vs C2" = master_induction_df$gene_id[
    master_induction_df$stat_test_adj_p < ADJ_P_THRESHOLD &
      abs(master_induction_df$cluster1_cluster2_differential) >= DIFF_THRESHOLD
  ],
  "C1 vs C3" = master_induction_df$gene_id[
    master_induction_df$stat_test_adj_p < ADJ_P_THRESHOLD &
      abs(master_induction_df$cluster1_cluster3_differential) >= DIFF_THRESHOLD
  ],
  "C2 vs C3" = master_induction_df$gene_id[
    master_induction_df$stat_test_adj_p < ADJ_P_THRESHOLD &
      abs(master_induction_df$cluster2_cluster3_differential) >= DIFF_THRESHOLD
  ]
)

cat("\n--- Venn 1: Global DE (direction-independent) ---\n")
cat("C1 vs C2 :", length(venn_sets_global[["C1 vs C2"]]), "genes\n")
cat("C1 vs C3 :", length(venn_sets_global[["C1 vs C3"]]), "genes\n")
cat("C2 vs C3 :", length(venn_sets_global[["C2 vs C3"]]), "genes\n")
cat("Unique   :", length(unique(unlist(venn_sets_global))), "genes\n")

venn_plot_global <- draw_venn(
  venn_sets  = venn_sets_global,
  title_label   = "All Differentially Expressed Genes — All Pairwise Comparisons",
  subtitle_label = paste0(
    "BH-adj p < ", ADJ_P_THRESHOLD,
    "  |  |differential| \u2265 ", DIFF_THRESHOLD,
    "  |  direction-independent"
  ),
  caption_label  = paste0(
    "k = ", n_clusters, " clusters  |  ",
    length(unique(unlist(venn_sets_global))), " unique DE genes"
  ),
  fill_colors = c("#7EB6D4", "#f9948c", "#6DBF81")   # C1 blue, C2 gold, C3 green
)

# -----------------------------------------------------------------------------
# 20d. VENN 2 — FUN-CIM POSITIVE INDUCTION vs DYS-CIM
# -----------------------------------------------------------------------------
# Genes higher in C1 (Fun-CIM) than C2 (Dys-CIM)
set_C1_pos_vs_C2 <- master_induction_df$gene_id[
  master_induction_df$stat_test_adj_p < ADJ_P_THRESHOLD &
    master_induction_df$cluster1_cluster2_differential >= DIFF_THRESHOLD
]

# Genes higher in C3 (Fun-CIM) than C2 (Dys-CIM)
set_C3_pos_vs_C2 <- master_induction_df$gene_id[
  master_induction_df$stat_test_adj_p < ADJ_P_THRESHOLD &
    master_induction_df$cluster3_cluster2_differential >= DIFF_THRESHOLD
]

venn_sets_funcim <- list(
  "C1 Fun-CIM\nvs Dys-CIM" = set_C1_pos_vs_C2,
  "C3 Fun-CIM\nvs Dys-CIM" = set_C3_pos_vs_C2
)

cat("\n--- Venn 2: Fun-CIM positive induction vs Dys-CIM ---\n")
cat("C1 Fun-CIM vs Dys-CIM :", length(set_C1_pos_vs_C2), "genes\n")
cat("C3 Fun-CIM vs Dys-CIM :", length(set_C3_pos_vs_C2), "genes\n")
cat("Overlap (core Fun-CIM) :", length(intersect(set_C1_pos_vs_C2, set_C3_pos_vs_C2)), "genes\n")

venn_plot_funcim <- draw_venn(
  venn_sets      = venn_sets_funcim,
  title_label    = "Fun-CIM Positively Induced Genes vs Dys-CIM",
  subtitle_label = paste0(
    "BH-adj p < ", ADJ_P_THRESHOLD,
    "  |  differential \u2265 ", DIFF_THRESHOLD,
    "  |  C1 & C3 vs C2"
  ),
  caption_label  = paste0(
    "Overlap = core Fun-CIM induction signature (shared C1 \u2229 C3)\n",
    "C1 only / C3 only = subcluster-specific Fun-CIM induction"
  ),
  fill_colors = c("#7EB6D4", "#6DBF81")               # blue = C1, green = C3
)

# -----------------------------------------------------------------------------
# 20e. VENN 3 — DYS-CIM SPECIFICITY vs EACH FUN-CIM ARM
# -----------------------------------------------------------------------------
# Genes higher in C2 (Dys-CIM) than C1 (Fun-CIM)
set_C2_pos_vs_C1 <- master_induction_df$gene_id[
  master_induction_df$stat_test_adj_p < ADJ_P_THRESHOLD &
    master_induction_df$cluster2_cluster1_differential >= DIFF_THRESHOLD
]

# Genes higher in C2 (Dys-CIM) than C3 (Fun-CIM)
set_C2_pos_vs_C3 <- master_induction_df$gene_id[
  master_induction_df$stat_test_adj_p < ADJ_P_THRESHOLD &
    master_induction_df$cluster2_cluster3_differential >= DIFF_THRESHOLD
]

venn_sets_dyscim <- list(
  "Dys-CIM (C2)\nvs C1 Fun-CIM" = set_C2_pos_vs_C1,
  "Dys-CIM (C2)\nvs C3 Fun-CIM" = set_C2_pos_vs_C3
)

cat("\n--- Venn 3: Dys-CIM specificity vs Fun-CIM arms ---\n")
cat("C2 Dys-CIM vs C1 Fun-CIM :", length(set_C2_pos_vs_C1), "genes\n")
cat("C2 Dys-CIM vs C3 Fun-CIM :", length(set_C2_pos_vs_C3), "genes\n")
cat("Overlap (core Dys-CIM)   :", length(intersect(set_C2_pos_vs_C1, set_C2_pos_vs_C3)), "genes\n")

venn_plot_dyscim <- draw_venn(
  venn_sets      = venn_sets_dyscim,
  title_label    = "Dys-CIM Positively Induced Genes vs Fun-CIM Clusters",
  subtitle_label = paste0(
    "BH-adj p < ", ADJ_P_THRESHOLD,
    "  |  differential \u2265 ", DIFF_THRESHOLD,
    "  |  C2 vs C1 & C3"
  ),
  caption_label  = paste0(
    "Overlap = robustly Dys-CIM-specific induction (C2 > C1 \u2229 C2 > C3)\n",
    "C2 vs C1 only / C2 vs C3 only = Dys-CIM induction relative to one arm"
  ),
  fill_colors = c("#f9948c", "#fe8e06")               # both gold = C2
)

# -----------------------------------------------------------------------------
# 20f. VENN 4 — OVERALL POSITIVE DIFFERENTIAL INDUCTION BY CLUSTER
# -----------------------------------------------------------------------------
# A gene is assigned to a cluster when it is significantly more induced in that
# cluster than in at least one of the other clusters. "Any comparison" therefore
# uses the union (OR), not the intersection, of that cluster's pairwise contrasts.
set_C1_pos_any <- master_induction_df$gene_id[
  master_induction_df$stat_test_adj_p < ADJ_P_THRESHOLD &
    (
      master_induction_df$cluster1_cluster2_differential >= DIFF_THRESHOLD |
        master_induction_df$cluster1_cluster3_differential >= DIFF_THRESHOLD
    )
]

set_C2_pos_any <- master_induction_df$gene_id[
  master_induction_df$stat_test_adj_p < ADJ_P_THRESHOLD &
    (
      master_induction_df$cluster2_cluster1_differential >= DIFF_THRESHOLD |
        master_induction_df$cluster2_cluster3_differential >= DIFF_THRESHOLD
    )
]

set_C3_pos_any <- master_induction_df$gene_id[
  master_induction_df$stat_test_adj_p < ADJ_P_THRESHOLD &
    (
      master_induction_df$cluster3_cluster1_differential >= DIFF_THRESHOLD |
        master_induction_df$cluster3_cluster2_differential >= DIFF_THRESHOLD
    )
]

venn_sets_cluster_positive <- list(
  "C1" = set_C1_pos_any,
  "C2" = set_C2_pos_any,
  "C3" = set_C3_pos_any
)

cat("\n--- Venn 4: Overall positive differential induction by cluster ---\n")
cat("C1 positive in any comparison:", length(set_C1_pos_any), "genes\n")
cat("C2 positive in any comparison:", length(set_C2_pos_any), "genes\n")
cat("C3 positive in any comparison:", length(set_C3_pos_any), "genes\n")
cat(
  "Unique genes across clusters     :",
  length(unique(unlist(venn_sets_cluster_positive))), "genes\n"
)

venn_plot_cluster_positive <- draw_venn(
  venn_sets      = venn_sets_cluster_positive,
  title_label    = "Overall Positive Differential Induction by Cluster",
  subtitle_label = paste0(
    "BH-adj p < ", ADJ_P_THRESHOLD,
    "  |  differential >= ", DIFF_THRESHOLD,
    "  |  positive in at least one pairwise comparison"
  ),
  caption_label  = paste0(
    "C1: C1 > C2 or C3  |  C2: C2 > C1 or C3  |  ",
    "C3: C3 > C1 or C2"
  ),
  fill_colors = c("#7EB6D4", "#f9948c", "#6DBF81")   # C1 blue, C2 gold, C3 green
)

# -----------------------------------------------------------------------------
# 20g. SAVE ALL FOUR VENNS AS SVGs AND HIGH-RESOLUTION PNGs
# -----------------------------------------------------------------------------
svg_venn_global <- file.path(output_dir, "nki_smc_venn_global_DE.svg")
ggsave(svg_venn_global, venn_plot_global,
       device = "svg", width = 7, height = 6, units = "in")
cat("Venn 1 (global DE) saved        :", svg_venn_global, "\n")

png_venn_global <- file.path(output_dir, "nki_smc_venn_global_DE.png")
ggsave(png_venn_global, venn_plot_global,
       device = "png", width = 7, height = 6, units = "in",
       dpi = 600, bg = "white")
cat("Venn 1 (global DE) PNG saved    :", png_venn_global, "\n")

svg_venn_funcim <- file.path(output_dir, "nki_smc_venn_funcim_vs_dyscim.svg")
ggsave(svg_venn_funcim, venn_plot_funcim,
       device = "svg", width = 7, height = 6, units = "in")
cat("Venn 2 (Fun-CIM induction) saved:", svg_venn_funcim, "\n")

png_venn_funcim <- file.path(output_dir, "nki_smc_venn_funcim_vs_dyscim.png")
ggsave(png_venn_funcim, venn_plot_funcim,
       device = "png", width = 7, height = 6, units = "in",
       dpi = 600, bg = "white")
cat("Venn 2 (Fun-CIM induction) PNG saved:", png_venn_funcim, "\n")

svg_venn_dyscim <- file.path(output_dir, "nki_smc_venn_dyscim_specificity.svg")
ggsave(svg_venn_dyscim, venn_plot_dyscim,
       device = "svg", width = 7, height = 6, units = "in")
cat("Venn 3 (Dys-CIM specificity) saved:", svg_venn_dyscim, "\n")

png_venn_dyscim <- file.path(output_dir, "nki_smc_venn_dyscim_specificity.png")
ggsave(png_venn_dyscim, venn_plot_dyscim,
       device = "png", width = 7, height = 6, units = "in",
       dpi = 600, bg = "white")
cat("Venn 3 (Dys-CIM specificity) PNG saved:", png_venn_dyscim, "\n")

svg_venn_cluster_positive <- file.path(
  output_dir, "nki_smc_venn_cluster_positive_any_comparison.svg"
)
ggsave(svg_venn_cluster_positive, venn_plot_cluster_positive,
       device = "svg", width = 7, height = 6, units = "in")
cat("Venn 4 (positive induction by cluster) saved:",
    svg_venn_cluster_positive, "\n")

png_venn_cluster_positive <- file.path(
  output_dir, "nki_smc_venn_cluster_positive_any_comparison.png"
)
ggsave(png_venn_cluster_positive, venn_plot_cluster_positive,
       device = "png", width = 7, height = 6, units = "in",
       dpi = 600, bg = "white")
cat("Venn 4 (positive induction by cluster) PNG saved:",
    png_venn_cluster_positive, "\n")

# -----------------------------------------------------------------------------
# 20h. SAVE ALL MEMBERSHIP CSVs  (full statistics preserved)
# -----------------------------------------------------------------------------
fwrite(
  build_membership_df(master_induction_df, venn_sets_global),
  file.path(output_dir, "venn_global_DE_membership.csv")
)
cat("Venn 1 membership CSV saved.\n")

fwrite(
  build_membership_df(master_induction_df, venn_sets_funcim),
  file.path(output_dir, "venn_funcim_vs_dyscim_membership.csv")
)
cat("Venn 2 membership CSV saved.\n")

fwrite(
  build_membership_df(master_induction_df, venn_sets_dyscim),
  file.path(output_dir, "venn_dyscim_specificity_membership.csv")
)
cat("Venn 3 membership CSV saved.\n")

fwrite(
  build_membership_df(master_induction_df, venn_sets_cluster_positive),
  file.path(output_dir, "venn_cluster_positive_any_comparison_membership.csv")
)
cat("Venn 4 membership CSV saved.\n")


# -----------------------------------------------------------------
# 21a. CIM STATE MAP  (define the annotation for each cluster)
# -----------------------------------------------------------------
# Cluster 1 → Fun‑CIM
# Cluster 2 → Dys‑CIM
# Cluster 3 → Fun‑CIM
cim_map <- c(
  "1" = "Fun-CIM",
  "2" = "Dys-CIM",
  "3" = "Fun-CIM"
)
# -----------------------------------------------------------------------------
# 21b. BUILD SANKEY DATA
# -----------------------------------------------------------------------------
sankey_df <- clustered_df %>%
  mutate(
    cluster_label = paste0("C", cluster_assignment),
    cim_state     = cim_map[as.character(cluster_assignment)]
  ) %>%
  group_by(cim_state, cluster_label) %>%
  summarise(n = n(), .groups = "drop")

cim_totals <- sankey_df %>%
  group_by(cim_state) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  mutate(label = paste0(cim_state, "\n(n = ", total, ")"))

cluster_totals <- sankey_df %>%
  group_by(cluster_label) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  mutate(label = paste0(cluster_label, "\n(n = ", total, ")"))

sankey_df <- sankey_df %>%
  left_join(
    cim_totals     %>% select(cim_state,     cim_label      = label),
    by = "cim_state"
  ) %>%
  left_join(
    cluster_totals %>% select(cluster_label, cluster_label2 = label),
    by = "cluster_label"
  )

# Explicit top-to-bottom order keeps both Fun-CIM rows together on each side:
# Dys-CIM (n = 11) -> C2 (n = 11)
# Fun-CIM (n = 25) -> C1 (n = 14), C3 (n = 11)
cim_order_top_to_bottom     <- c("Dys-CIM", "Fun-CIM")
cluster_order_top_to_bottom <- c("C2", "C1", "C3")

cim_label_levels <- cim_totals %>%
  mutate(cim_state = factor(cim_state, levels = cim_order_top_to_bottom)) %>%
  arrange(cim_state) %>%
  pull(label)

cluster_label_levels <- cluster_totals %>%
  mutate(
    cluster_label = factor(
      cluster_label,
      levels = cluster_order_top_to_bottom
    )
  ) %>%
  arrange(cluster_label) %>%
  pull(label)

# ggalluvial draws factor levels from bottom to top, so reverse the desired order.
sankey_df$cim_label <- factor(
  sankey_df$cim_label,
  levels = rev(cim_label_levels)
)
sankey_df$cluster_label2 <- factor(
  sankey_df$cluster_label2,
  levels = rev(cluster_label_levels)
)

# -----------------------------------------------------------------------------
# 21c. COLOUR PALETTE
# -----------------------------------------------------------------------------
cim_colours <- c(
  "Dys-CIM" = "#7EB6D4",    # salmon/pink
  "Fun-CIM" =   "#f9948c"   # muted blue
)

cluster_colours <- c(
  "C1" = "#7EB6D4",    # blue
  "C2" = "#f9948c",    # gold/amber
  "C3" = "#6DBF81"     # green
)

# Match the labelled strata (which include counts) to their corresponding colors.
stratum_colours <- c(
  setNames(cim_colours[cim_totals$cim_state], cim_totals$label),
  setNames(cluster_colours[cluster_totals$cluster_label], cluster_totals$label)
)

# -----------------------------------------------------------------------------
# 21d. DRAW SANKEY
# -----------------------------------------------------------------------------
sankey_plot <- ggplot(
  sankey_df,
  aes(axis1 = cim_label, axis2 = cluster_label2, y = n)
) +
  geom_alluvium(
    aes(fill = cim_state),
    width    = 0.12,
    alpha    = 0.85,
    knot.pos = 0.4
  ) +
  geom_stratum(
    aes(fill = after_stat(stratum)),
    width     = 0.12,
    color     = "white",
    linewidth = 0.3
  ) +
  geom_text(
    stat     = "stratum",
    aes(label = after_stat(stratum)),
    size     = 2.4,
    fontface = "bold",
    color    = "black"
  ) +
  scale_x_discrete(
    limits = c("CIM State", "Cluster"),
    labels = NULL,
    expand = c(0.03, 0.03)
  ) +
  scale_fill_manual(
    values = c(cim_colours, cluster_colours, stratum_colours),
    guide  = "none"
  ) +
  labs(title = NULL, subtitle = NULL, caption = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.margin   = margin(2, 2, 2, 2)
  )

# -----------------------------------------------------------------------------
# 21e. SAVE SANKEY OUTPUTS
# -----------------------------------------------------------------------------
svg_sankey <- file.path(output_dir, "nki_smc_sankey.svg")
ggsave(
  svg_sankey,
  sankey_plot,
  device = "svg",
  width = 4.34,
  height = 2.51,
  units = "in",
  bg = "white"
)
cat("Sankey SVG saved:", svg_sankey, "\n")

png_sankey <- file.path(output_dir, "nki_smc_sankey.png")
ggsave(
  png_sankey,
  sankey_plot,
  device = "png",
  width = 4.34,
  height = 2.51,
  units = "in",
  dpi = 600,
  bg = "white"
)
cat("Sankey high-resolution PNG saved:", png_sankey, "\n")

csv_sankey <- file.path(output_dir, "nki_smc_sankey_counts.csv")
fwrite(sankey_df %>% select(cim_state, cluster_label, n), csv_sankey)
cat("Sankey counts CSV saved:", csv_sankey, "\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("\n============================================================\n")
cat("ALL OUTPUTS SAVED TO:\n", output_dir, "\n")
cat("============================================================\n")
cat("  nki_smc_master_induction_df.csv\n")
cat("  nki_smc_combine_clustered_plot_df.csv\n")
cat("  nki_smc_venn_diagram.svg\n")
cat("  venn_membership_full.csv\n")
cat("  venn_gene_sets.csv\n")
cat("  nki_smc_venn_global_DE.png\n")
cat("  nki_smc_venn_funcim_vs_dyscim.png\n")
cat("  nki_smc_venn_dyscim_specificity.png\n")
cat("  nki_smc_venn_cluster_positive_any_comparison.svg\n")
cat("  nki_smc_venn_cluster_positive_any_comparison.png\n")
cat("  venn_cluster_positive_any_comparison_membership.csv\n")
cat("  nki_smc_sankey.svg\n")
cat("  nki_smc_sankey.png\n")
cat("  nki_smc_sankey_counts.csv\n")
cat("============================================================\n")