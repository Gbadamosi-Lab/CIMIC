# ============================================================================
# figure_individual_gene_panels.R
# ============================================================================
# Purpose:
#   Per-gene two-panel summary across CIMIC trajectories (Dys-CIM vs Fun-CIM):
#     (left)  Induction boxplot  - Δ(log2[TPM+1]) distribution, Wilcoxon P
#     (right) Penetrance barplot - % of samples with positive induction
#                                   (Δ > 0), Fisher/χ² P
#   Also provides a point-biserial scatter and 2-/3-panel combiners.
#
# Input:
#   Datasets/NKI_SMC/nki_smc_combine_clustered_plot_df.csv
#     - sample_id like '<patient>_delta'
#     - cluster_assignments in {1, 2}   (1 = Dys-CIM, 2 = Fun-CIM)
#     - per-gene Δ(log2[TPM+1]) columns (e.g. EIF2AK3)
#
# Output (one PNG per gene of interest):
#   Results/figure_gene_panel_<GENE>.png
#
# Dependencies: dplyr, tidyr, ggplot2, ggpubr, rstatix, stringr, purrr,
#               scales, data.table
# ============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)      # ggarrange / annotate_figure
library(rstatix)     # wilcox_test, kruskal_test, etc.
library(stringr)     # str_wrap
library(purrr)
library(scales)      # percent axis labels

dir.create("Results", showWarnings = FALSE)

# ----------------------------------------------------------------------------
# Load data and build the CIMIC_Cluster label column expected by the plotters
# ----------------------------------------------------------------------------
clustered_plot_df_for_merge <- data.table::fread(
  "Datasets/NKI_SMC/nki_smc_combine_clustered_plot_df.csv",
  data.table = FALSE
) %>%
  dplyr::mutate(
    # Cluster convention: 1 = Dys-CIM, 2 = Fun-CIM
    CIMIC_Cluster = factor(
      dplyr::case_when(
        cluster_assignments == 1 ~ "Dys-CIM",
        cluster_assignments == 2 ~ "Fun-CIM"
      ),
      levels = c("Dys-CIM", "Fun-CIM")
    )
  )

message(sprintf("Samples per cluster: %s",
                paste(names(table(clustered_plot_df_for_merge$CIMIC_Cluster)),
                      table(clustered_plot_df_for_merge$CIMIC_Cluster),
                      sep = "=", collapse = ", ")))

# ============================================================================
# 1. Induction boxplot  ------------------------------------------------------
# ============================================================================
compare_induction_boxplot <- function(df, gene) {

  gene_col <- df %>%
    dplyr::select(dplyr::matches(paste0("^", gene, "(\\.x|\\.y)?$"))) %>%
    names()
  if (length(gene_col) == 0) {
    stop("Gene column '", gene, "' not found in the data frame")
  }
  gene_col <- gene_col[[1]]   # use the first match

  df2 <- df %>%
    dplyr::mutate(
      Expression = as.numeric(as.character(.data[[gene_col]]))
    ) %>%
    dplyr::select(CIMIC_Cluster, Expression)   # keep only what we need

  pval <- df2 %>%
    rstatix::wilcox_test(Expression ~ CIMIC_Cluster) %>%
    dplyr::pull(p)

  ptxt <- paste0("P = ", format.pval(pval, digits = 4))

  ggplot(df2,
         aes(x = CIMIC_Cluster,
             y = Expression,
             fill = CIMIC_Cluster)) +
    geom_boxplot(width = 0.6,
                 outlier.shape = NA,
                 alpha = 0.8) +
    geom_jitter(width = 0.15,
                size = 1.5,
                alpha = 0.4) +
    annotate("text",
             x = 1.5,
             y = max(df2$Expression, na.rm = TRUE) * 1.1,
             label = ptxt,
             size = 8,
             fontface = "bold") +
    scale_fill_manual(values = c("Dys-CIM" = "#2980B9",
                                 "Fun-CIM" = "#C0392B")) +
    labs(title = paste0(gene),
         y = "Δ(log2[TPM+1])",
         x = NULL) +
    theme_classic(base_size = 14) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "none",
      axis.text.x     = element_text(face = "bold", size = 20),
      axis.text.y     = element_text(face = "bold", size = 20),
      axis.title.y    = element_text(face = "bold", size = 20),
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 1.5),
      axis.line       = element_line(colour = "black", linewidth = 1.5)
    )
}

# ============================================================================
# 2. Penetrance stacked bar  -------------------------------------------------
# ============================================================================
compare_penetrance_plot <- function(df, gene) {
  # 2.1  Find the exact gene column
  gene_col <- df %>%
    dplyr::select(dplyr::matches(paste0("^", gene, "(\\.x|\\.y)?$"))) %>%
    names()
  if (length(gene_col) == 0) {
    stop("Gene column '", gene, "' not found in the data frame")
  }
  gene_col <- gene_col[[1]]

  # 2.2  Binary induction flag (Δ > 0 = positive induction)
  df2 <- df %>%
    dplyr::mutate(
      Expression = as.numeric(as.character(.data[[gene_col]])),
      Induction  = ifelse(Expression > 0,
                          "Positive Induction",
                          "Negative/No Induction")
    ) %>%
    dplyr::select(CIMIC_Cluster, Induction)

  # 2.3  Fisher or χ² test on the 2×2 table
  tab <- df2 %>%
    dplyr::count(CIMIC_Cluster, Induction) %>%
    tidyr::pivot_wider(names_from = Induction, values_from = n, values_fill = 0)

  # Guard against an absent column if a cluster has no negatives/positives
  for (col in c("Negative/No Induction", "Positive Induction")) {
    if (!col %in% names(tab)) tab[[col]] <- 0
  }
  mat <- as.matrix(tab[ , c("Negative/No Induction", "Positive Induction")])
  pval <- if (any(mat < 5) || sum(mat) < 50) {
    stats::fisher.test(mat)$p.value
  } else {
    stats::chisq.test(mat)$p.value
  }
  ptxt <- paste0("P = ", format.pval(pval, digits = 4))

  # 2.4  Data for the stacked-bar plot
  bar_df <- df2 %>%
    dplyr::group_by(CIMIC_Cluster, Induction) %>%
    dplyr::summarise(n = n(), .groups = "drop") %>%
    dplyr::group_by(CIMIC_Cluster) %>%
    dplyr::mutate(Total = sum(n), Prop = n / Total) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(
      CIMIC_Cluster,
      factor(Induction, levels = c("Negative/No Induction", "Positive Induction"))
    )

  # 2.5  Plot
  ggplot(bar_df,
         aes(x = CIMIC_Cluster, y = Prop, fill = Induction)) +
    geom_bar(stat = "identity",
             position = "fill",
             width = 0.7,
             aes(fill = factor(Induction,
                               levels = c("Negative/No Induction",
                                          "Positive Induction")))) +
    geom_text(aes(y = Prop,
                  label = ifelse(Prop > 0, paste0(round(Prop * 100), "%"), "")),
              position = position_stack(vjust = 0.5),
              colour = "black", size = 8, fontface = "bold") +
    geom_text(aes(x = 1.5, y = 1.06, label = ptxt),
              inherit.aes = FALSE, size = 8, fontface = "bold") +
    geom_text(aes(x = CIMIC_Cluster, y = -0.05, label = paste0("N = ", Total)),
              inherit.aes = FALSE, size = 6, fontface = "bold") +
    labs(title = paste0(gene), y = "Percentage", fill = "Induction") +
    scale_fill_manual(values = c("Positive Induction"    = "#3fb949ff",
                                 "Negative/No Induction" = "#f76661ff")) +
    scale_y_continuous(labels = scales::percent,
                       limits = c(-0.1, 1.12),
                       breaks = seq(0, 1, 0.2)) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    theme_classic(base_size = 18) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 18),
      legend.position = "top",
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 1.5),
      axis.line.x     = element_line(colour = "black", linewidth = 1.5),
      axis.line.y     = element_line(colour = "black", linewidth = 1.5),
      axis.text.x     = element_text(face = "bold", size = 20),
      axis.text.y     = element_text(face = "bold", size = 20),
      axis.title.y    = element_text(face = "bold", size = 20),
      legend.text     = element_text(size = 18),
      legend.title    = element_text(size = 18, face = "bold")
    )
}

# ============================================================================
# 3. Point-biserial correlation scatter  ------------------------------------
# ============================================================================
compare_pointbiserial_plot <- function(df, gene) {

  gene_col <- df %>%
    dplyr::select(dplyr::matches(paste0("^", gene, "(\\.x|\\.y)?$"))) %>%
    names()
  if (length(gene_col) == 0) {
    stop("Gene column '", gene, "' not found in the data frame")
  }
  gene_col <- gene_col[[1]]

  df2 <- df %>%
    dplyr::mutate(
      Expression  = as.numeric(as.character(.data[[gene_col]])),
      Cluster_bin = ifelse(CIMIC_Cluster == "Fun-CIM", 1, 0)
    ) %>%
    dplyr::select(CIMIC_Cluster, Cluster_bin, Expression)

  cor_res <- cor.test(df2$Cluster_bin, df2$Expression, method = "pearson")
  rtxt <- sprintf("r = %.3f\nP = %s",
                  cor_res$estimate,
                  format.pval(cor_res$p.value, digits = 4))

  ggplot(df2,
         aes(x = Cluster_bin, y = Expression, colour = CIMIC_Cluster)) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 1.2) +
    scale_x_continuous(breaks = c(0, 1), labels = c("Dys-CIM", "Fun-CIM")) +
    scale_colour_manual(values = c("Dys-CIM" = "#C0392B",
                                   "Fun-CIM" = "#2980B9")) +
    annotate("text",
             x = 0.5,
             y = max(df2$Expression, na.rm = TRUE) * 1.1,
             label = rtxt, size = 5, fontface = "bold") +
    labs(title = paste0(gene), x = NULL, y = "Δ(log2[TPM+1])") +
    theme_classic(base_size = 14) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 18),
      legend.position = "none",
      axis.text.x     = element_text(face = "bold", size = 18),
      axis.text.y     = element_text(face = "bold", size = 18),
      axis.title.y    = element_text(face = "bold", size = 18),
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 1.5),
      axis.line       = element_line(colour = "black", linewidth = 1.5)
    )
}

# ============================================================================
# 4. Panel combiners  --------------------------------------------------------
# ============================================================================
make_three_panel <- function(gene) {
  ind <- compare_induction_boxplot(clustered_plot_df_for_merge, gene)
  pen <- compare_penetrance_plot(clustered_plot_df_for_merge, gene)
  pb  <- compare_pointbiserial_plot(clustered_plot_df_for_merge, gene)

  ggpubr::ggarrange(ind, pen, pb, ncol = 3, nrow = 1, common.legend = FALSE) %>%
    ggpubr::annotate_figure(
      top = ggpubr::text_grob(paste0(gene, " – Summary"), face = "bold", size = 18)
    )
}

make_two_panel <- function(gene) {
  ind <- compare_induction_boxplot(clustered_plot_df_for_merge, gene)
  pen <- compare_penetrance_plot(clustered_plot_df_for_merge, gene)

  # Drop the per-panel titles (a single shared title is added below) and the
  # penetrance legend, matching the reference two-panel layout.
  ind <- ind + theme(plot.title = element_blank())
  pen <- pen + theme(plot.title = element_blank(),
                     axis.title.x = element_blank(),
                     legend.position = "none")

  ggpubr::ggarrange(ind, pen, ncol = 2, nrow = 1, common.legend = FALSE) %>%
    ggpubr::annotate_figure(
      top = ggpubr::text_grob(paste0(gene), face = "bold", size = 18)
    )
}

# ============================================================================
# 5. Build & save figures
# ============================================================================
genes_of_interest <- c("EIF2AK3")   # add more genes here as needed

for (g in genes_of_interest) {
  fig <- make_two_panel(g)
  out <- sprintf("Results/figure_gene_panel_%s.png", g)
  ggsave(out, plot = fig, width = 9, height = 5, dpi = 300, bg = "white")
  message(sprintf("Saved: %s", out))
}

message("Done.")
