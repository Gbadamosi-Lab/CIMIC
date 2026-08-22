# ============================================================================
# pre_post_treatment_correlation_analysis.R
# ============================================================================
# Purpose: For three treatment-response datasets (NKI_SMC, NEO, and
#          TNBC_CL_Epirubicin), compute per patient / cell-line Spearman
#          correlations, restricted to CIM genes, between:
#            (1) pre-treatment and post-treatment expression profiles
#                ("Pre vs Post")
#            (2) the treatment-induced delta (post - pre) and the
#                pre-treatment profile ("Delta vs Pre")
#          and visualize |rho| for each comparison as paired boxplots.
#
#          CIM gene universe: NKI_SMC and NEO (patient datasets) are
#          restricted to the genes built by human_full_cimic_gene_sets.R;
#          TNBC_CL_Epirubicin (cell lines) is restricted to the genes built
#          by cell_line_full_cimic_gene_sets.R.
#
#          NKI_SMC: column suffix "_B" = pre-treatment, "_S" = post-
#                   treatment, paired by shared sample-id prefix
#                   (e.g. "2046_B" / "2046_S").
#          NEO:     "_T1_B_Rep1" = pre-treatment, "_T4_S_Rep1" = post-
#                   treatment, paired by shared patient prefix
#                   (e.g. "NEO11 _Patient_11").
#          TNBC_CL_Epirubicin: "_DMSO_R#" = pre-treatment (vehicle),
#                   "_EPI_R#" = post-treatment (epirubicin). The DMSO
#                   replicates are averaged per cell line into a single
#                   shared baseline; deltas and Pre-vs-Post correlations
#                   are then computed for each EPI replicate against that
#                   shared baseline.
#
# Inputs:
#   oncoimmunology_paper/Datasets/NKI_SMC/nki_smc_combine_TPM_lognorm.csv
#   oncoimmunology_paper/Datasets/NEO/neo_tpm_lognorm.csv
#   oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epi_TPM_lognorm.csv
#   oncoimmunology_paper/helper_code/human_full_cimic_gene_sets.R
#   oncoimmunology_paper/helper_code/cell_line_full_cimic_gene_sets.R
#
# Outputs:
#   oncoimmunology_paper/Results/pre_post_treatment_correlation_nki_smc.csv
#   oncoimmunology_paper/Results/pre_post_treatment_correlation_neo.csv
#   oncoimmunology_paper/Results/pre_post_treatment_correlation_tnbc_cl.csv
#   oncoimmunology_paper/Results/pre_post_treatment_correlation_nki_smc.png
#   oncoimmunology_paper/Results/pre_post_treatment_correlation_neo.png
#   oncoimmunology_paper/Results/pre_post_treatment_correlation_tnbc_cl.png
# ============================================================================

library(data.table)   # fast CSV reading
library(dplyr)        # data manipulation
library(ggplot2)      # plotting
library(rstatix)       # paired Wilcoxon test

results_dir <- "oncoimmunology_paper/Results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# CIM gene universes
# Source each helper script into its own environment (avoids polluting this
# script's environment) and pull out its `all_unique_gene_ids` CIM gene list.
# cell_line_full_cimic_gene_sets.R writes an .rds to CIMIC_BASE_DIR; point it
# at a temp directory so sourcing has no side effects on the repo.
# ---------------------------------------------------------------------------
source_gene_set_env <- function(path, base_dir_override = NULL) {
  env <- new.env()
  if (!is.null(base_dir_override)) {
    old_val <- Sys.getenv("CIMIC_BASE_DIR", unset = NA)
    Sys.setenv(CIMIC_BASE_DIR = base_dir_override)
    on.exit({
      if (is.na(old_val)) Sys.unsetenv("CIMIC_BASE_DIR") else Sys.setenv(CIMIC_BASE_DIR = old_val)
    })
  }
  sys.source(path, envir = env)
  env
}

human_gene_set_env <- source_gene_set_env("oncoimmunology_paper/helper_code/human_full_cimic_gene_sets.R")
human_cim_genes <- human_gene_set_env$all_unique_gene_ids

cell_line_gene_set_env <- source_gene_set_env(
  "oncoimmunology_paper/helper_code/cell_line_full_cimic_gene_sets.R",
  base_dir_override = tempdir()
)
cell_line_cim_genes <- cell_line_gene_set_env$all_unique_gene_ids

stopifnot(length(human_cim_genes) > 0, length(cell_line_cim_genes) > 0)

# ---------------------------------------------------------------------------
# Spearman correlation of two numeric vectors -> rho and p-value
# ---------------------------------------------------------------------------
spearman_stat <- function(x, y) {
  test <- suppressWarnings(stats::cor.test(x, y, method = "spearman"))
  c(rho = unname(test$estimate), p_value = test$p.value)
}

# ---------------------------------------------------------------------------
# Given matched pre/post gene-expression matrices (genes x samples, columns
# paired 1:1 and identified by `sample_ids`), compute for each sample:
#   - "Pre vs Post"  : Spearman(pre, post)
#   - "Delta vs Pre" : Spearman(post - pre, pre)
# Returns a long data.frame (two rows per sample) with BH-adjusted p-values
# computed separately within each comparison type.
# ---------------------------------------------------------------------------
compute_pre_post_correlations <- function(pre_mat, post_mat, sample_ids) {
  stopifnot(ncol(pre_mat) == ncol(post_mat), ncol(pre_mat) == length(sample_ids))

  rows <- lapply(seq_along(sample_ids), function(i) {
    pre_vec   <- pre_mat[, i]
    post_vec  <- post_mat[, i]
    delta_vec <- post_vec - pre_vec

    pp <- spearman_stat(pre_vec, post_vec)
    dp <- spearman_stat(delta_vec, pre_vec)

    data.frame(
      sample_id  = sample_ids[i],
      comparison = c("Pre vs Post", "Delta vs Pre"),
      rho        = c(pp[["rho"]], dp[["rho"]]),
      p_value    = c(pp[["p_value"]], dp[["p_value"]]),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out$abs_rho <- abs(out$rho)

  out %>%
    dplyr::group_by(comparison) %>%
    dplyr::mutate(padj = stats::p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup() %>%
    as.data.frame()
}

# Compute average Spearman rho for each comparison type
average_correlation <- function(corr_df) {
  corr_df %>%
    dplyr::group_by(comparison) %>%
    dplyr::summarise(mean_rho = mean(rho, na.rm = TRUE)) %>%
    dplyr::ungroup()
}

# ---------------------------------------------------------------------------
# Shared boxplot builder: |rho| for "Delta vs Pre" vs "Pre vs Post"
# If a `cluster` column exists, the plot is faceted by cluster (C1 / C2).
# Optional `shape_var` can be supplied (e.g., cell line) to map point shapes.
# ---------------------------------------------------------------------------
build_boxplot <- function(corr_df, title, shape_var = NULL) {
  corr_df$comparison <- factor(corr_df$comparison, levels = c("Delta vs Pre", "Pre vs Post"))

  stat_res <- suppressWarnings(
    rstatix::wilcox_test(corr_df, abs_rho ~ comparison, paired = TRUE)
  )
  p_val <- stat_res$p
  signif_label <- dplyr::case_when(
    p_val <= 0.0001 ~ "****",
    p_val <= 0.001  ~ "***",
    p_val <= 0.01   ~ "**",
    p_val <= 0.05   ~ "*",
    TRUE            ~ sprintf("p = %.3g", p_val)
  )

  fill_colors <- c("Delta vs Pre" = "#F8766D", "Pre vs Post" = "#00BFC4")
  star_y <- max(corr_df$abs_rho, na.rm = TRUE) * 1.05

  p <- ggplot(corr_df, aes(x = comparison, y = rho, fill = comparison)) +
    geom_boxplot(outlier.shape = NA, colour = "black", width = 0.6, alpha = 0.75)

  # Jitter points – optionally map shape
  if (!is.null(shape_var) && shape_var %in% colnames(corr_df)) {
    # Limit to 6 distinct shapes for readability. Any extra categories are
    # assigned the first shape in the palette.
    shape_palette <- c(16, 17, 15, 18, 8, 3)
    # Replace missing values with a placeholder so they are not dropped.
    corr_df[[shape_var]] <- ifelse(is.na(corr_df[[shape_var]]), "Other", corr_df[[shape_var]])
    vals <- sort(unique(corr_df[[shape_var]]))
    # Map each value to a shape, recycling the palette if more than 6 levels.
    shape_map <- setNames(rep(shape_palette, length.out = length(vals)), vals)
    p <- p +
      geom_jitter(aes(shape = .data[[shape_var]]), width = 0.15, size = 5, colour = "black", alpha = 0.75) +
      scale_shape_manual(name = "Cell line", values = shape_map)
  } else {
    p <- p + geom_jitter(width = 0.15, size = 2.5, colour = "black", alpha = 0.6)
  }

  # Average rho annotations
  avg_df <- corr_df %>%
    dplyr::group_by(comparison) %>%
    dplyr::summarise(mean_rho = mean(rho, na.rm = TRUE)) %>%
    dplyr::ungroup()

  p <- p +
    annotate("text", x = 1.5, y = star_y, label = signif_label, size = 8, fontface = "bold") +
    annotate("text", x = 1, y = -0.95,
             label = sprintf("avg \u03C1 = %.3f", avg_df$mean_rho[avg_df$comparison == "Delta vs Pre"]),
             size = 5) +
    annotate("text", x = 2, y = -0.95,
             label = sprintf("avg \u03C1 = %.3f", avg_df$mean_rho[avg_df$comparison == "Pre vs Post"]),
             size = 5) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "black") +
    scale_fill_manual(values = fill_colors, guide = "none") +
    labs(title = title, x = NULL, y = expression("Spearman " * rho)) +
    coord_cartesian(ylim = c(-1, 1))

  # Facet by cluster if present
  if ("cluster" %in% colnames(corr_df)) {
    p <- p + facet_wrap(~cluster)
  }

  p +
    theme_classic(base_size = 18) +
    theme(
      axis.text.x   = element_text(face = "bold", size = 20, colour = "black"),
      axis.text.y   = element_text(face = "bold", size = 20, colour = "black"),
      axis.title.y  = element_text(face = "bold", size = 22, colour = "black"),
      plot.title    = element_text(face = "bold", size = 26, hjust = 0.5),
      legend.position = if (!is.null(shape_var)) "right" else "none",
      panel.grid    = element_blank(),
      panel.border  = element_rect(linewidth = 1, fill = NA)
    )
}

# ============================================================================
# NKI_SMC: "_B" = pre-treatment, "_S" = post-treatment
# ============================================================================
nki_path <- "oncoimmunology_paper/Datasets/NKI_SMC/nki_smc_combine_TPM_lognorm.csv"
nki_dt <- fread(nki_path)

nki_sample_cols <- grep("_(B|S)$", colnames(nki_dt), value = TRUE)
nki_prefixes    <- unique(sub("_(B|S)$", "", nki_sample_cols))

nki_pre_cols  <- paste0(nki_prefixes, "_B")
nki_post_cols <- paste0(nki_prefixes, "_S")
stopifnot(all(nki_pre_cols %in% colnames(nki_dt)), all(nki_post_cols %in% colnames(nki_dt)))

nki_pre_mat  <- as.matrix(nki_dt[, ..nki_pre_cols]);  storage.mode(nki_pre_mat)  <- "double"
nki_post_mat <- as.matrix(nki_dt[, ..nki_post_cols]); storage.mode(nki_post_mat) <- "double"

# Restrict to CIM genes
nki_keep <- nki_dt$gene_id %in% human_cim_genes
stopifnot(sum(nki_keep) > 0)
nki_pre_mat  <- nki_pre_mat[nki_keep, , drop = FALSE]
nki_post_mat <- nki_post_mat[nki_keep, , drop = FALSE]

nki_corr <- compute_pre_post_correlations(nki_pre_mat, nki_post_mat, nki_prefixes)
# Add CIMIC cluster information (C1 = Dys‑CIM, C2 = Fun‑CIM)
# Load clustered data and align sample IDs (remove trailing "_delta" if present)
clustered_nki <- data.table::fread("oncoimmunology_paper/Datasets/NKI_SMC/nki_smc_combine_clustered_plot_df.csv")
if ("cluster_assignments" %in% colnames(clustered_nki)) {
  clustered_nki <- clustered_nki %>% dplyr::rename(CIMIC_Cluster = cluster_assignments)
}
# Ensure sample_id matches the correlation IDs (strip "_delta" suffix)
clustered_nki <- clustered_nki %>%
  dplyr::mutate(sample_id = sub("_delta$", "", sample_id))
cluster_map_nki <- c("1" = "C1", "2" = "C2")
nki_corr <- nki_corr %>%
  dplyr::left_join(clustered_nki %>% dplyr::select(sample_id, CIMIC_Cluster), by = "sample_id") %>%
  dplyr::mutate(cluster = cluster_map_nki[as.character(CIMIC_Cluster)])

write.csv(nki_corr, file.path(results_dir, "pre_post_treatment_correlation_nki_smc.csv"), row.names = FALSE)

# Save average correlation per comparison for NKI_SMC
nki_avg <- average_correlation(nki_corr)
write.csv(nki_avg, file.path(results_dir, "pre_post_treatment_correlation_nki_smc_avg.csv"), row.names = FALSE)

nki_plot <- build_boxplot(nki_corr, "NKI_SMC")
ggsave(file.path(results_dir, "pre_post_treatment_correlation_nki_smc.png"),
       plot = nki_plot, width = 6, height = 6, dpi = 300)

# ============================================================================
# NEO: "_T1_B_Rep1" = pre-treatment, "_T4_S_Rep1" = post-treatment
# ============================================================================
neo_path <- "oncoimmunology_paper/Datasets/NEO/neo_tpm_lognorm.csv"
neo_dt <- fread(neo_path)

neo_pre_cols  <- grep("_T[0-9]+_B_Rep[0-9]+$", colnames(neo_dt), value = TRUE)
neo_post_cols <- grep("_T[0-9]+_S_Rep[0-9]+$", colnames(neo_dt), value = TRUE)

neo_pre_prefix  <- sub("_T[0-9]+_B_Rep[0-9]+$", "", neo_pre_cols)
neo_post_prefix <- sub("_T[0-9]+_S_Rep[0-9]+$", "", neo_post_cols)
stopifnot(setequal(neo_pre_prefix, neo_post_prefix))

# Reorder post columns so they line up 1:1 with pre columns by patient prefix
neo_post_cols <- neo_post_cols[match(neo_pre_prefix, neo_post_prefix)]

neo_pre_mat  <- as.matrix(neo_dt[, ..neo_pre_cols]);  storage.mode(neo_pre_mat)  <- "double"
neo_post_mat <- as.matrix(neo_dt[, ..neo_post_cols]); storage.mode(neo_post_mat) <- "double"

# Restrict to CIM genes
neo_keep <- neo_dt$gene_id %in% human_cim_genes
stopifnot(sum(neo_keep) > 0)
neo_pre_mat  <- neo_pre_mat[neo_keep, , drop = FALSE]
neo_post_mat <- neo_post_mat[neo_keep, , drop = FALSE]

neo_corr <- compute_pre_post_correlations(neo_pre_mat, neo_post_mat, neo_pre_prefix)
# Add CIMIC cluster information (C1 = Dys‑CIM, C2 = Fun‑CIM)
clustered_neo <- data.table::fread("oncoimmunology_paper/Datasets/NEO/neo_clustered_plot_df.csv")
if ("cluster_assignments" %in% colnames(clustered_neo)) {
  clustered_neo <- clustered_neo %>% dplyr::rename(CIMIC_Cluster = cluster_assignments)
}
# Align sample IDs if they contain a suffix (e.g., "_delta")
if ("sample_id" %in% colnames(clustered_neo)) {
  clustered_neo <- clustered_neo %>%
    dplyr::mutate(sample_id = sub("_delta$", "", sample_id))
}
cluster_map_neo <- c("1" = "C1", "2" = "C2")
neo_corr <- neo_corr %>%
  dplyr::left_join(clustered_neo %>% dplyr::select(sample_id, CIMIC_Cluster), by = "sample_id") %>%
  dplyr::mutate(cluster = cluster_map_neo[as.character(CIMIC_Cluster)])

write.csv(neo_corr, file.path(results_dir, "pre_post_treatment_correlation_neo.csv"), row.names = FALSE)

# Save average correlation per comparison for NEO
neo_avg <- average_correlation(neo_corr)
write.csv(neo_avg, file.path(results_dir, "pre_post_treatment_correlation_neo_avg.csv"), row.names = FALSE)

neo_plot <- build_boxplot(neo_corr, "NEO")
ggsave(file.path(results_dir, "pre_post_treatment_correlation_neo.png"),
       plot = neo_plot, width = 6, height = 6, dpi = 300)

# ============================================================================
# TNBC_CL_Epirubicin: "_DMSO_R#" = pre-treatment, "_EPI_R#" = post-treatment
# DMSO replicates are averaged per cell line into one shared baseline; each
# EPI replicate is compared against that same baseline.
# ============================================================================
tnbc_path <- "oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epi_TPM_lognorm.csv"
tnbc_dt <- fread(tnbc_path)

tnbc_dmso_cols <- grep("_DMSO_R[0-9]+$", colnames(tnbc_dt), value = TRUE)
tnbc_epi_cols  <- grep("_EPI_R[0-9]+$",  colnames(tnbc_dt), value = TRUE)

tnbc_dmso_mat <- as.matrix(tnbc_dt[, ..tnbc_dmso_cols]); storage.mode(tnbc_dmso_mat) <- "double"
tnbc_epi_mat  <- as.matrix(tnbc_dt[, ..tnbc_epi_cols]);  storage.mode(tnbc_epi_mat)  <- "double"

# Restrict to CIM genes (use the first of the two duplicate "gene_id" columns)
tnbc_gene_id <- tnbc_dt[[which(colnames(tnbc_dt) == "gene_id")[1]]]
tnbc_keep <- tnbc_gene_id %in% cell_line_cim_genes
stopifnot(sum(tnbc_keep) > 0)
tnbc_dmso_mat <- tnbc_dmso_mat[tnbc_keep, , drop = FALSE]
tnbc_epi_mat  <- tnbc_epi_mat[tnbc_keep, , drop = FALSE]

# Shared DMSO baseline per cell line (mean across its 3 replicates)
tnbc_dmso_cell_line <- sub("_DMSO_R[0-9]+$", "", tnbc_dmso_cols)
tnbc_baseline_mat <- sapply(split(seq_along(tnbc_dmso_cols), tnbc_dmso_cell_line), function(idx) {
  rowMeans(tnbc_dmso_mat[, idx, drop = FALSE])
})

# Expand the per-cell-line baseline so it aligns 1:1 with each EPI replicate
tnbc_epi_cell_line <- sub("_EPI_R[0-9]+$", "", tnbc_epi_cols)
stopifnot(all(tnbc_epi_cell_line %in% colnames(tnbc_baseline_mat)))
tnbc_pre_mat_expanded <- tnbc_baseline_mat[, tnbc_epi_cell_line, drop = FALSE]

tnbc_corr <- compute_pre_post_correlations(tnbc_pre_mat_expanded, tnbc_epi_mat, tnbc_epi_cols)
tnbc_corr$cell_line <- sub("_EPI_R[0-9]+$", "", tnbc_corr$sample_id)
# Add CIMIC cluster information (C1 = Fun‑CIM, C2 = Dys‑CIM for TNBC)
clustered_tnbc <- data.table::fread("oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epi_clustered_plot_df.csv")
if ("cluster_assignments" %in% colnames(clustered_tnbc)) {
  clustered_tnbc <- clustered_tnbc %>% dplyr::rename(CIMIC_Cluster = cluster_assignments)
}
# Keep a single row per cell line (base_id) for cluster assignment
clustered_tnbc <- clustered_tnbc %>%
  dplyr::select(base_id, CIMIC_Cluster) %>%
  dplyr::distinct(base_id, .keep_all = TRUE)
cluster_map_tnbc <- c("1" = "C2", "2" = "C1")
tnbc_corr <- tnbc_corr %>%
  dplyr::left_join(clustered_tnbc, by = c("cell_line" = "base_id")) %>%
  dplyr::mutate(cluster = cluster_map_tnbc[as.character(CIMIC_Cluster)])
write.csv(tnbc_corr, file.path(results_dir, "pre_post_treatment_correlation_tnbc_cl.csv"), row.names = FALSE)

# Save average correlation per comparison for TNBC_CL_Epirubicin
tnbc_avg <- average_correlation(tnbc_corr)
write.csv(tnbc_avg, file.path(results_dir, "pre_post_treatment_correlation_tnbc_cl_avg.csv"), row.names = FALSE)

tnbc_plot <- build_boxplot(tnbc_corr, "TNBC_CL_Epirubicin", shape_var = "cell_line")
ggsave(file.path(results_dir, "pre_post_treatment_correlation_tnbc_cl.png"),
       plot = tnbc_plot, width = 9, height = 6, dpi = 300)

# End of script
