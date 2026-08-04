# ============================================================================
# figure3b_pam50.R
# ============================================================================
# Purpose:
#   Classify NKI/SMC tumor samples into intrinsic breast cancer subtypes using
#   the PAM50 signature (genefu), then compare the distribution of PAM50
#   subtypes across CIMIC clusters (Dys-CIM vs Fun-CIM), separately at baseline
#   (pre-treatment biopsy) and post-treatment (surgery). Significance of the
#   association is tested with Fisher's exact test.
#
# Inputs:
#   - Datasets/NKI_SMC/nki_smc_combine_TPM_lognorm.csv: Log-normalized TPM matrix
#       (rows = genes, first column 'gene_id'; columns = samples named
#        <patient>_B [biopsy/baseline] and <patient>_S [surgery/post-treatment])
#   - Datasets/NKI_SMC/nki_smc_combine_clustered_plot_df.csv: Sample metadata
#       (sample_id like '<patient>_delta', plus cluster_assignments in {1, 2})
#
# Outputs:
#   - Results/figure3b_pam50_baseline.png        (pre-treatment, 300 DPI)
#   - Results/figure3b_pam50_posttreatment.png   (post-treatment, 300 DPI)
#
# Cluster convention (per request):
#   cluster_assignments == 1  ->  Dys-CIM
#   cluster_assignments == 2  ->  Fun-CIM
#
# Dependencies:
#   - R >= 4.0; genefu, org.Hs.eg.db, AnnotationDbi, data.table, tidyverse
#
# Author: Mohammed Gbadamosi
# Last Updated: June 2026
# ============================================================================

# ============================================================================
# LIBRARY IMPORTS
# ============================================================================
library(data.table)     # Fast CSV reading
library(genefu)         # PAM50 molecular subtyping
library(org.Hs.eg.db)   # Human gene annotation
library(AnnotationDbi)  # Annotation helpers
library(tidyverse)      # dplyr, tidyr, stringr, tibble, ggplot2, purrr

# Create output directory if needed
dir.create("Results", showWarnings = FALSE)

# ============================================================================
# DATA LOADING
# ============================================================================
# Log-normalized TPM matrix (rows = genes via 'gene_id', columns = samples)
combined_TPM_lognorm <- data.table::fread(
  "Datasets/NKI_SMC/nki_smc_combine_TPM_lognorm.csv",
  data.table = FALSE
)

# Sample metadata with CIMIC cluster assignments (per-patient, on the delta)
clustered_plot_df <- data.table::fread(
  "Datasets/NKI_SMC/nki_smc_combine_clustered_plot_df.csv",
  data.table = FALSE
)

message(sprintf("Loaded TPM matrix: %d genes x %d samples",
                nrow(combined_TPM_lognorm), ncol(combined_TPM_lognorm) - 1))
message(sprintf("Loaded %d clustered (delta) samples", nrow(clustered_plot_df)))

# ============================================================================
# PAM50 MOLECULAR SUBTYPING
# ============================================================================
# Load PAM50 centroid models bundled with genefu
data(pam50)
data(pam50.robust)

# Build expression matrix: genes in rows, then transpose to samples x genes
expr_mat <- combined_TPM_lognorm %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

texp_mat <- t(expr_mat)                    # samples (rows) x genes (cols)
texp_mat_numeric <- texp_mat
class(texp_mat_numeric) <- "numeric"

if (any(is.na(texp_mat_numeric))) {
  message("Warning: NAs detected in expression matrix; replacing with 0.")
  texp_mat_numeric[is.na(texp_mat_numeric)] <- 0
}

# Gene annotation: match by gene symbol (do.mapping = FALSE)
gene_info <- data.frame(Genes = colnames(texp_mat))
rownames(gene_info) <- gene_info$Genes

# PAM50 classification per sample
pam50_predictions <- molecular.subtyping(
  sbt.model = "pam50",
  data      = texp_mat_numeric,
  annot     = gene_info,
  do.mapping = FALSE
)

# Per-sample subtype calls
out <- data.frame(
  Sample        = rownames(texp_mat),
  PAM50_Subtype = pam50_predictions$subtype,
  stringsAsFactors = FALSE
)

message("PAM50 subtype calls (all samples):")
print(table(out$PAM50_Subtype))

# ============================================================================
# ANNOTATE TREATMENT (Baseline vs Post-Treatment)
# ============================================================================
# TPM sample names: <patient>_B = biopsy (baseline), <patient>_S = surgery
# (post-treatment). 'base' is the patient identifier shared with the metadata.
pam50_df <- out %>%
  dplyr::filter(Sample != "gene_id") %>%
  dplyr::select(Sample, PAM50_Subtype) %>%
  dplyr::mutate(
    base = stringr::str_remove(Sample, "_[BS]$"),
    Treatment = dplyr::case_when(
      stringr::str_detect(Sample, "_B$") ~ "Baseline",
      stringr::str_detect(Sample, "_S$") ~ "Post-Treatment"
    )
  )

# ============================================================================
# ANNOTATE CIMIC CLUSTER
# ============================================================================
# Metadata sample_id like '<patient>_delta'; 'base' links it to PAM50 samples.
cluster_df <- clustered_plot_df %>%
  dplyr::select(sample_id, cluster_assignments) %>%
  dplyr::mutate(
    base = stringr::str_remove(sample_id, "_delta$"),
    # Cluster convention: 1 = Dys-CIM, 2 = Fun-CIM
    CIMIC_Cluster = dplyr::case_when(
      cluster_assignments == 1 ~ "Dys-CIM",
      cluster_assignments == 2 ~ "Fun-CIM"
    )
  ) %>%
  dplyr::select(base, CIMIC_Cluster)

# Join PAM50 calls to CIMIC clusters by patient base id
plot_df <- pam50_df %>%
  dplyr::left_join(cluster_df, by = "base") %>%
  dplyr::filter(!is.na(CIMIC_Cluster), !is.na(Treatment))

baseline_df_pam50 <- plot_df %>% dplyr::filter(Treatment == "Baseline")
post_df_pam50     <- plot_df %>% dplyr::filter(Treatment == "Post-Treatment")

# ============================================================================
# PLOTTING FUNCTION: STACKED % BARPLOT + FISHER TEST
# ============================================================================
create_pam50_plot <- function(plot_df, condition_name) {

  # --- Format data: per-cluster subtype percentages --------------------------
  plot_data <- plot_df %>%
    dplyr::group_by(CIMIC_Cluster, PAM50_Subtype) %>%
    dplyr::summarise(total_n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(CIMIC_Cluster) %>%
    dplyr::mutate(
      total_samples = sum(total_n),
      percentage = (total_n / total_samples) * 100
    ) %>%
    dplyr::ungroup()

  # --- Sample-size labels (N per cluster) ------------------------------------
  sample_sizes <- plot_data %>%
    dplyr::group_by(CIMIC_Cluster) %>%
    dplyr::summarise(N = sum(total_n), .groups = "drop")

  labels_df <- data.frame(
    x = sample_sizes$CIMIC_Cluster,
    y = -8,
    label = paste0("N = ", sample_sizes$N)
  )

  # --- Fisher's exact test on the cluster x subtype contingency table --------
  contingency_table <- plot_data %>%
    dplyr::select(CIMIC_Cluster, PAM50_Subtype, total_n) %>%
    tidyr::pivot_wider(
      names_from  = PAM50_Subtype,
      values_from = total_n,
      values_fill = 0
    )

  mat <- as.matrix(contingency_table[, -1])
  rownames(mat) <- contingency_table$CIMIC_Cluster

  fisher_p <- fisher.test(mat, simulate.p.value = TRUE)$p.value
  fish_stats_P <- paste0("P = ", signif(fisher_p, 3))

  # --- Subtype colors --------------------------------------------------------
  subtype_colors <- c(
    "Basal"  = "#C44E52",
    "Her2"   = "#DD8DA6",
    "LumA"   = "#0072B2",
    "LumB"   = "#5C8C45",
    "Normal" = "#9467BD"
  )

  # --- Plot ------------------------------------------------------------------
  p <- ggplot(
    plot_data,
    aes(
      x = factor(CIMIC_Cluster, levels = c("Dys-CIM", "Fun-CIM")),
      y = percentage,
      fill = PAM50_Subtype
    )
  ) +
    geom_bar(stat = "identity", width = 0.85, color = "black", linewidth = 1) +
    geom_text(
      aes(label = ifelse(percentage > 0, paste0(round(percentage, 1), "%"), "")),
      position = position_stack(vjust = 0.5),
      size = 6, fontface = "bold", color = "black"
    ) +
    geom_text(
      data = labels_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE, size = 6, fontface = "bold"
    ) +
    annotate("text", x = 1.5, y = 110, label = fish_stats_P,
             size = 6, fontface = "bold") +
    scale_fill_manual(values = subtype_colors) +
    scale_y_continuous(limits = c(-15, 117), breaks = seq(0, 100, 20)) +
    scale_x_discrete(expand = expansion(mult = c(0.55, 0.55))) +
    labs(
      title = paste0("PAM50 Subtypes — ", condition_name),
      x = "", y = "Percentage", fill = ""
    ) +
    theme_classic(base_size = 18) +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.text    = element_text(colour = "black", size = 18, face = "bold"),
      axis.title   = element_text(colour = "black", size = 20, face = "bold"),
      axis.ticks   = element_line(linewidth = 1.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
      legend.position = "right",
      legend.text  = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 16, face = "bold")
    ) +
    guides(
      fill = guide_legend(
        title = NULL,
        keywidth  = unit(1, "cm"),
        keyheight = unit(0.8, "cm")
      )
    )

  return(p)
}

# ============================================================================
# BUILD & SAVE FIGURES
# ============================================================================
p_baseline <- create_pam50_plot(baseline_df_pam50, "Pre-Treatment")
p_post     <- create_pam50_plot(post_df_pam50,     "Post-Treatment")

print(p_baseline)
print(p_post)

ggsave("Results/figure3b_pam50_baseline.png",
       plot = p_baseline, width = 8, height = 7, dpi = 300, bg = "white")
ggsave("Results/figure3b_pam50_posttreatment.png",
       plot = p_post, width = 8, height = 7, dpi = 300, bg = "white")

message("✓ Saved: Results/figure3b_pam50_baseline.png")
message("✓ Saved: Results/figure3b_pam50_posttreatment.png")
message("✓ Cluster convention: 1 = Dys-CIM, 2 = Fun-CIM")
message("✓ Script complete: PAM50 subtype analysis finished successfully")
