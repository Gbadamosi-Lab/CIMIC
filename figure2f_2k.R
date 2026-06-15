# figure2f_2k.R
# ------------------------------------------------------------
# This script reproduces the Z‑score calculation for the cell‑death
# programs (as defined in `all_gene_set_names_death`) and generates
# faceted box‑plots similar to Figure 2F‑2K of the manuscript.
# ------------------------------------------------------------

# Load required packages ------------------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
# (Add any additional packages needed for score calculation, e.g. readr, GSVA)

# ------------------------------------------------------------
# 1. Define cell‑death gene‑set names
# ------------------------------------------------------------
# Replace the placeholder vector with the actual gene‑set names used in the
# manuscript (e.g. "Apoptosis", "Ferroptosis", "Necroptosis", …).
all_gene_set_names_death <- c(
  "Apoptosis",
  "Ferroptosis",
  "Necroptosis",
  "Pyroptosis",
  "PANoptosis"
  # ... add any other cell‑death programs here
)

# ------------------------------------------------------------
# 2. Load data (expression matrix and sample metadata)
# ------------------------------------------------------------
# NOTE: Update the file paths to point to your own data files.
expr_matrix <- readRDS("path/to/expression_matrix.rds")   # genes x samples
metadata    <- readRDS("path/to/sample_metadata.rds")   # must contain `cimic_cluster`

# ------------------------------------------------------------
# 3. Compute program scores (placeholder implementation)
# ------------------------------------------------------------
# The original `figure_2e` script contains the logic for scoring each
# program. Here we use a simple mean‑expression approach as a stand‑in.
# Replace with the exact method (e.g., ssGSEA, AUCell) if required.

compute_program_score <- function(gene_set) {
  # `gene_set` – character vector of gene symbols belonging to the program
  # Returns a numeric vector of scores (one per sample)
  intersect_genes <- intersect(rownames(expr_matrix), gene_set)
  if (length(intersect_genes) == 0) return(rep(NA, ncol(expr_matrix)))
  colMeans(expr_matrix[intersect_genes, , drop = FALSE])
}

# Create a tidy data frame of raw scores for each program
program_scores_raw <- lapply(all_gene_set_names_death, function(prog) {
  # In a real workflow, you would retrieve the gene list for `prog`
  # For this template we assume a function `get_genes_for_program()` exists.
  gene_list <- get_genes_for_program(prog)  # <-- user must define
  tibble(
    Program = prog,
    Sample  = colnames(expr_matrix),
    Score   = compute_program_score(gene_list)
  )
}) %>% bind_rows()

# Merge with metadata to obtain the CIMiC cluster assignment
program_scores_all <- program_scores_raw %>%
  left_join(metadata %>% select(Sample, cimic_cluster), by = "Sample") %>%
  mutate(cimic_cluster = factor(cimic_cluster))

# ------------------------------------------------------------
# 4. Convert raw scores to Z‑scores per program
# ------------------------------------------------------------
program_scores_all <- program_scores_all %>%
  group_by(Program) %>%
  mutate(Score = scale(Score, center = TRUE, scale = TRUE)[, 1]) %>%
  ungroup()

# ------------------------------------------------------------
# 5. Summary statistics and p‑values (same as in figure_2e)
# ------------------------------------------------------------
# 1. Summary stats per cluster
program_summary <- program_scores_all %>%
  group_by(cimic_cluster, Program) %>%
  summarise(
    mean_score   = mean(Score, na.rm = TRUE),
    median_score = median(Score, na.rm = TRUE),
    sd_score     = sd(Score, na.rm = TRUE),
    n            = n(),
    .groups = "drop"
  )

# 2. P‑values per program (handles k = 2 or k >= 3)
program_pvals <- program_scores_all %>%
  group_by(Program) %>%
  summarise(
    n_clusters = n_distinct(cimic_cluster),
    p_value = tryCatch({
      if (n_clusters == 2) {
        wilcox.test(Score ~ cimic_cluster)$p.value
      } else if (n_clusters > 2) {
        kruskal.test(Score ~ cimic_cluster)$p.value
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(p_value, method = "fdr"))

# 3. Merge summary and p‑values
program_summary_final <- program_summary %>%
  left_join(program_pvals, by = "Program")

# ------------------------------------------------------------
# 6. Prepare data for faceted box‑plots (only cell‑death programs)
# ------------------------------------------------------------
program_scores_plot_df <- program_scores_all %>%
  filter(Program %in% paste0(all_gene_set_names_death, "_z")) %>%
  mutate(
    Program = Program %>%
      str_remove_all("GOBP_|REACTOME_|HALLMARK_|_z") %>%
      str_replace_all("_", " "),
    Program = factor(
      Program,
      levels = unique(
        paste0(all_gene_set_names_death, "_z") %>%
          str_remove_all("GOBP_|REACTOME_|HALLMARK_|_z") %>%
          str_replace_all("_", " ")
      )
    )
  )

# ------------------------------------------------------------
# 7. Compute significance annotations
# ------------------------------------------------------------
pval_df_all <- program_scores_plot_df %>%
  group_by(Program) %>%
  summarise(
    n_groups = n_distinct(cimic_cluster),
    y.position = max(Score, na.rm = TRUE) * 1.1,
    pval = tryCatch({
      if (n_groups == 2) {
        wilcox.test(Score ~ cimic_cluster)$p.value
      } else if (n_groups > 2) {
        kruskal.test(Score ~ cimic_cluster)$p.value
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(pval, method = "fdr"),
    FDR.signif = case_when(
      is.na(FDR)                ~ "FDR = NA",
      FDR < 0.001              ~ "FDR < 0.001",
      TRUE                     ~ paste0("FDR = ", signif(FDR, 3))
    )
  )

# ------------------------------------------------------------
# 8. Faceted box‑plot for ALL cell‑death programs
# ------------------------------------------------------------
p_all_programs_p <- ggplot(program_scores_plot_df, aes(x = cimic_cluster, y = Score, fill = cimic_cluster)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 1) +
  geom_text(
    data = pval_df_all,
    aes(x = 1.5, y = y.position, label = FDR.signif),
    inherit.aes = FALSE,
    size = 12,
    fontface = "bold"
  ) +
  facet_wrap(~ Program, scales = "free_y") +
  labs(
    title = "Program scores by CIMiC cluster",
    y = "Z-score of\nProgram Activation",
    x = NULL
  ) +
  theme_bw(base_line_size = 2) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 18),
    axis.text = element_text(colour = "black", size = 24),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 24),
    axis.title.y = element_text(colour = "black", size = 24),
    legend.text = element_text(colour = "black", size = 24),
    legend.title = element_text(colour = "black", size = 24),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 2)
  )
print(p_all_programs_p)

# ------------------------------------------------------------
# 9. Faceted box‑plot for SIGNIFICANT programs only
# ------------------------------------------------------------
sig_programs <- pval_df_all %>%
  filter(!is.na(FDR), FDR < 0.05) %>%
  pull(Program)

program_scores_plot_df_sig <- program_scores_plot_df %>%
  filter(Program %in% sig_programs)

pval_df_sig <- pval_df_all %>%
  filter(Program %in% sig_programs)

p_all_programs_sig <- ggplot(program_scores_plot_df_sig, aes(x = cimic_cluster, y = Score, fill = cimic_cluster)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 1) +
  geom_text(
    data = pval_df_sig,
    aes(x = 1.5, y = y.position, label = FDR.signif),
    inherit.aes = FALSE,
    size = 12,
    fontface = "bold"
  ) +
  facet_wrap(~ Program, scales = "free_y") +
  labs(
    title = "Significant program scores by CIMiC cluster",
    y = "Z-score of\nProgram Activation",
    x = NULL
  ) +
  theme_bw(base_line_size = 2) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 18),
    axis.text = element_text(colour = "black", size = 24),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 24),
    axis.title.y = element_text(colour = "black", size = 24),
    legend.text = element_text(colour = "black", size = 24),
    legend.title = element_text(colour = "black", size = 24),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 2)
  )
print(p_all_programs_sig)

# ------------------------------------------------------------
# 10. Example: single‑program plot (e.g., Ferroptosis)
# ------------------------------------------------------------
single_gene_set <- "Ferroptosis"

program_scores_plot_df_sig <- program_scores_plot_df %>%
  filter(Program %in% single_gene_set) %>%
  mutate(Program = str_wrap(Program, width = 20))

pval_df_sig <- pval_df_all %>%
  filter(Program %in% single_gene_set) %>%
  mutate(Program = str_wrap(Program, width = 20))

y_pos_one <- max(program_scores_plot_df_sig$Score, na.rm = TRUE) + 0.1

p_program_individual <- ggplot(program_scores_plot_df_sig, aes(x = cimic_cluster, y = Score, fill = cimic_cluster)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 1) +
  geom_text(
    data = pval_df_sig,
    aes(x = 1.5, y = y_pos_one, label = FDR.signif),
    inherit.aes = FALSE,
    size = 15,
    fontface = "bold"
  ) +
  facet_wrap(~ Program, scales = "free_y") +
  labs(
    y = "Z-score of\nProgram Activation",
    x = NULL
  ) +
  theme_bw(base_line_size = 2) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 32),
    axis.text = element_text(colour = "black", size = 24),
    axis.title.x = element_blank(),
    axis.text.x = element_text(colour = "black", size = 24, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 24, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 24, face = "bold"),
    legend.text = element_text(colour = "black", size = 24),
    legend.title = element_text(colour = "black", size = 24),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 2)
  )
print(p_program_individual)

# End of script ------------------------------------------------
