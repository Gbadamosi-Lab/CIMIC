# =============================================================================
# figure2f_2k.R
# =============================================================================
# Purpose:
#   Create faceted box-plots of the five cell-death programs (Apoptosis,
#   Ferroptosis, Necroptosis, Pyroptosis, PANoptosis) across CIMIC clusters.
#   The script follows the same data sources, z-score preprocessing and
#   aesthetic choices as figure2e_zscore_heatmap.R for consistency.
#
#   Two figures are produced, differing only in the significance test used for
#   the stars (significance is based on RAW, uncorrected p-values):
#     - Results/figure2f_2k_boxplots.png        (Wilcoxon rank-sum)
#     - Results/figure2f_2k_boxplots_welch.png  (Welch's t-test)
# =============================================================================

# -------------------------------------------------------------------------
# 1. Load required libraries (identical to figure2e)
# -------------------------------------------------------------------------
library(data.table)   # fast CSV reading
library(dplyr)        # data manipulation
library(tidyr)        # reshaping
library(purrr)        # functional programming
library(ggplot2)      # plotting
library(stringr)      # string handling
library(rstatix)      # statistical helpers

# -------------------------------------------------------------------------
# 2. Load sample metadata (same source as figure2e)
# -------------------------------------------------------------------------
clustered_plot_df <- data.table::fread(
  "Datasets/NKI_SMC/nki_smc_combine_clustered_plot_df.csv",
  data.table = FALSE
) %>% as.data.frame()

# Align column name with figure2e conventions
if ("cluster_assignments" %in% colnames(clustered_plot_df)) {
  clustered_plot_df <- clustered_plot_df %>%
    rename(CIMIC_Cluster = cluster_assignments)
}

# -------------------------------------------------------------------------
# 3. Identify gene-expression columns (same logic as figure2e)
# -------------------------------------------------------------------------
metadata_cols <- c(
  "sample_id", "CIMIC_Cluster", "base_id", "PACMAP1", "PACMAP2",
  "UMAP1", "UMAP2", "PC1", "PC2", "PC3", "V1"
)
gene_cols <- setdiff(colnames(clustered_plot_df), metadata_cols)

message(sprintf("Loaded %d samples with %d gene expression columns",
                nrow(clustered_plot_df), length(gene_cols)))

# Create Results directory early (mirrors figure2e)
dir.create("Results", showWarnings = FALSE)

# -------------------------------------------------------------------------
# 4. Z-score each gene across all samples (identical to figure2e)
# -------------------------------------------------------------------------
zscore_safe <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(x)) else (x - m) / s
}
df_z <- clustered_plot_df %>%
  mutate(across(all_of(gene_cols), zscore_safe))

# -------------------------------------------------------------------------
# 5. Define the five cell-death programs (exactly as in figure2e)
# -------------------------------------------------------------------------
death_programs <- list(
  Apoptosis = c(
    "BAX", "BAK1", "BBC3", "PMAIP1", "BCL2L11",
    "APAF1", "CASP9", "CASP3", "CASP7",
    "FAS", "TNFRSF10B", "CASP8", "FADD", "BID",
    "DIABLO", "CYCS", "TP53AIP1", "BAD", "BMF", "HRK"
  ),
  Necroptosis = c(
    "RIPK1", "RIPK3", "MLKL", "ZBP1", "TICAM1",
    "TNFRSF1A", "FADD", "CASP8"
  ),
  Pyroptosis = c(
    "GSDME", "GSDMD", "CASP1", "CASP4", "CASP5",
    "NLRP3", "AIM2", "PYCARD", "IL1B", "IL18",
    "CASP3", "NLRC4"
  ),
  PANoptosis = c(
    "ZBP1", "AIM2", "RIPK3", "RIPK1", "CASP8",
    "CASP1", "FADD", "PYCARD", "IRF1"
  ),
  Ferroptosis = c(
    "ACSL4", "LPCAT3", "TFRC", "SAT1", "PTGS2"
  )
)

# -------------------------------------------------------------------------
# 6. Compute program scores - mean Z-score per sample for each program
# -------------------------------------------------------------------------
score_program_z <- function(data, genes, prog_name) {
  genes <- intersect(genes, colnames(data))
  if (length(genes) == 0) {
    return(data.frame(
      sample_id = data$sample_id,
      CIMIC_Cluster = data$CIMIC_Cluster,
      Program = prog_name,
      Score = rep(NA_real_, nrow(data))
    ))
  }
  score_vec <- data %>%
    dplyr::select(all_of(genes)) %>%
    as.matrix() %>%
    rowMeans(na.rm = TRUE)
  data.frame(
    sample_id = data$sample_id,
    CIMIC_Cluster = data$CIMIC_Cluster,
    Program = prog_name,
    Score = score_vec
  )
}
program_scores_all <- purrr::imap_dfr(death_programs,
                                      ~ score_program_z(df_z, .x, .y))

# -------------------------------------------------------------------------
# 7. Statistical testing (RAW, uncorrected p-values)
# -------------------------------------------------------------------------
# Significance stars are based on RAW p-values per program (no FDR correction).
# Two test variants are produced:
#   - "wilcox": Wilcoxon rank-sum (2 clusters); Kruskal-Wallis if >2 clusters
#   - "welch" : Welch's t-test, unequal variance (2 clusters);
#               one-way Welch ANOVA if >2 clusters

# Map a raw p-value to a significance label.
signif_label <- function(p) {
  case_when(
    is.na(p)      ~ "",
    p <= 0.0001   ~ "****",
    p <= 0.001    ~ "***",
    p <= 0.01     ~ "**",
    p <= 0.05     ~ "*",
    TRUE          ~ ""
  )
}

# Compute per-program raw p-values for a given test type ("wilcox" or "welch").
compute_stats <- function(scores, test = c("wilcox", "welch")) {
  test <- match.arg(test)
  scores %>%
    group_by(Program) %>%
    summarise(
      n_clusters = n_distinct(CIMIC_Cluster),
      p_val = tryCatch({
        if (test == "wilcox") {
          if (n_clusters == 2) wilcox.test(Score ~ CIMIC_Cluster)$p.value
          else                 kruskal.test(Score ~ CIMIC_Cluster)$p.value
        } else {
          if (n_clusters == 2) t.test(Score ~ CIMIC_Cluster, var.equal = FALSE)$p.value
          else                 oneway.test(Score ~ CIMIC_Cluster, var.equal = FALSE)$p.value
        }
      }, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(signif = signif_label(p_val))   # raw p-values, no FDR correction
}

# Ensure consistent ordering/labeling of clusters
clustered_plot_df$CIMIC_Cluster <- factor(clustered_plot_df$CIMIC_Cluster)

# -------------------------------------------------------------------------
# 8. Faceted box-plot builder (publication-ready, matches Figure 2E style)
# -------------------------------------------------------------------------
make_boxplot <- function(stat_df, test_label) {
  plot_df <- program_scores_all %>% left_join(stat_df, by = "Program")
  plot_df$CIMIC_Cluster <- factor(plot_df$CIMIC_Cluster)

  # Per-program star position: facets use free_y, so place each star just above
  # that program's own maximum score (a single global max would push stars
  # off-screen for low-scoring panels).
  star_df <- program_scores_all %>%
    group_by(Program) %>%
    summarise(y = max(Score, na.rm = TRUE) * 1.05, .groups = "drop") %>%
    left_join(stat_df, by = "Program")

  ggplot(plot_df,
         aes(x = CIMIC_Cluster, y = Score, fill = CIMIC_Cluster)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 1) +
    geom_text(data = star_df,
              aes(x = 1.5, y = y, label = signif),
              inherit.aes = FALSE, size = 6, fontface = "bold") +
    facet_wrap(~ Program, scales = "free_y") +
    labs(
      title = sprintf("Cell-death program Z-scores by CIMIC cluster (%s, raw p)", test_label),
      y = "Z-score of Program Activation",
      x = NULL
    ) +
    # Salmon = Dys-CIM (cluster 1), Teal = Fun-CIM (cluster 2), matching figure panels
    scale_fill_manual(values = c("1" = "#F8766D", "2" = "#00BFC4")) +
    scale_x_discrete(labels = c("1" = "Dys-CIM", "2" = "Fun-CIM"),
                     expand = expansion(add = c(0.5, 0.5))) +
    theme_classic(base_size = 18) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 16),
      axis.text.x = element_text(face = "bold", size = 28, colour = "black"),
      axis.text.y = element_text(colour = "black", size = 14),
      axis.title = element_text(colour = "black", size = 16),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 1)
    )
}

# -------------------------------------------------------------------------
# 9. Build both figures
# -------------------------------------------------------------------------
stat_wilcox <- compute_stats(program_scores_all, "wilcox")
stat_welch  <- compute_stats(program_scores_all, "welch")

p_wilcox <- make_boxplot(stat_wilcox, "Wilcoxon")
p_welch  <- make_boxplot(stat_welch,  "Welch's t-test")

print(p_wilcox)
print(p_welch)

# -------------------------------------------------------------------------
# 10. Save figures (same folder as figure2e outputs)
# -------------------------------------------------------------------------
if (!dir.exists("Results")) dir.create("Results", showWarnings = FALSE)

ggsave(
  filename = "Results/figure2f_2k_boxplots.png",
  plot = p_wilcox,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "Results/figure2f_2k_boxplots_welch.png",
  plot = p_welch,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

message("Wilcoxon raw p-values:")
print(stat_wilcox)
message("Welch's t-test raw p-values:")
print(stat_welch)
