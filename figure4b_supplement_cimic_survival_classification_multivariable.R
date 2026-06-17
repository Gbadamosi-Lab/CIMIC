# ============================================================================
# figure_cimic_survival_classification_multivariable.R
# ============================================================================
# Purpose:
#   Classify each CIMIC-associated gene as Beneficial / Adverse / NS based on
#   cross-cohort concordant MULTIVARIABLE survival associations
#   (METABRIC + SCANB), then visualize the distribution across Fun-CIM vs
#   Dys-CIM trajectories as a stacked percentage bar plot.
#
#   Two versions are produced:
#     - raw p-value   (< 0.05)  -> Results/..._multivariable_pvalue.png
#     - BH-adjusted p (< 0.05)  -> Results/..._multivariable_padj_BH.png
#
# Input:
#   Datasets/Multivariable_CIMIC_Clinical_Analysis.csv
#     - Category : trajectory label ("Fun-CIM" / "Dys-CIM")
#     - gene     : gene symbol
#     - {COHORT}_{ENDPOINT}_HR / _CI_lower / _CI_upper / _p_value / _p_adj_BH
#       COHORT   in {METABRIC, SCANB}
#       METABRIC endpoints: OS, DSS, RFI
#       SCANB    endpoints: OS, RFI, DRFI
#
# Classification (per gene):
#   Beneficial : HR < 1 AND p < 0.05 in >= 1 endpoint in BOTH cohorts
#   Adverse    : HR > 1 AND p < 0.05 in >= 1 endpoint in BOTH cohorts
#   NS         : everything else (including discordant or single-cohort signals)
#
# Dependencies: ggplot2, dplyr
# ============================================================================

library(dplyr)
library(ggplot2)

dir.create("Results", showWarnings = FALSE)

# ----------------------------------------------------------------------------
# 1. Load data
# ----------------------------------------------------------------------------
df <- read.csv(
  "Datasets/Multivariable_CIMIC_Clinical_Analysis.csv",
  stringsAsFactors = FALSE,
  check.names      = FALSE
)
names(df)[1] <- "Category"

# Endpoints per cohort
mb_endpoints    <- c("OS", "DSS", "RFI")
scanb_endpoints <- c("OS", "RFI", "DRFI")

# Shared aesthetics
outcome_colors <- c(
  "Adverse"    = "#FF6B6B",
  "NS"         = "#BFBFBF",
  "Beneficial" = "#4CAF50"
)
class_levels <- c("Adverse", "NS", "Beneficial")
cat_levels   <- c("Dys-CIM", "Fun-CIM")

# ----------------------------------------------------------------------------
# 2. Helper: does a gene show directional significance in >= 1 endpoint
#            within a given cohort?
#    data      : full data frame (all rows passed at once via vectorised loop)
#    cohort    : "METABRIC" or "SCANB"
#    endpoints : character vector of endpoint codes
#    direction : "benef" (HR < 1) or "adv" (HR > 1)
#    p_suffix  : "_p_value" or "_p_adj_BH"
# ----------------------------------------------------------------------------
cohort_signal <- function(data, cohort, endpoints, direction, p_suffix) {
  signal <- rep(FALSE, nrow(data))
  for (ep in endpoints) {
    hr_col <- paste0(cohort, "_", ep, "_HR")
    p_col  <- paste0(cohort, "_", ep, p_suffix)
    hr <- data[[hr_col]]
    p  <- data[[p_col]]
    dir_ok <- if (direction == "benef") hr < 1 else hr > 1
    cond   <- dir_ok & (p < 0.05)
    cond[is.na(cond)] <- FALSE
    signal <- signal | cond
  }
  signal
}

# ----------------------------------------------------------------------------
# 3. Scientific notation label for Fisher p-value
# ----------------------------------------------------------------------------
fmt_sci_p <- function(p) {
  if (is.na(p) || p <= 0) return("P < 1 %*% 10^-300")
  e <- floor(log10(p))
  m <- p / 10^e
  sprintf("P == %.2f %%*%% 10^%d", m, e)
}

# ----------------------------------------------------------------------------
# 4. Build figure for a given significance column suffix
# ----------------------------------------------------------------------------
build_figure <- function(data, p_suffix, p_type_label, out_file) {

  # Store reference BEFORE mutate so it can be passed to cohort_signal
  df_input <- data
  # --- Classify genes -------------------------------------------------------
  data <- data %>%
    mutate(
      mb_benef    = cohort_signal(df_input, "METABRIC", mb_endpoints,    "benef", p_suffix),
      mb_adv      = cohort_signal(df_input, "METABRIC", mb_endpoints,    "adv",   p_suffix),
      scanb_benef = cohort_signal(df_input, "SCANB",    scanb_endpoints, "benef", p_suffix),
      scanb_adv   = cohort_signal(df_input, "SCANB",    scanb_endpoints, "adv",   p_suffix),
      Classification = dplyr::case_when(
        mb_benef & scanb_benef ~ "Beneficial",
        mb_adv   & scanb_adv   ~ "Adverse",
        TRUE                   ~ "NS"
      )
    )

  # --- Per-category percentages ---------------------------------------------
  counts <- data %>%
    dplyr::count(Category, Classification, name = "n")

  plot_data <- expand.grid(
    Category       = cat_levels,
    Classification = class_levels,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(counts, by = c("Category", "Classification")) %>%
    dplyr::mutate(n = dplyr::coalesce(n, 0L)) %>%
    dplyr::group_by(Category) %>%
    dplyr::mutate(total = sum(n), percentage = n / total * 100) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Classification = factor(Classification, levels = c("Beneficial", "NS", "Adverse")),
      Category       = factor(Category, levels = cat_levels)
    )

  # N labels
  n_labels <- plot_data %>%
    dplyr::distinct(Category, total) %>%
    dplyr::mutate(label = paste0("N = ", total))

  # --- Fisher exact test ----------------------------------------------------
  cont_tab <- table(
    factor(data$Category,        levels = cat_levels),
    factor(data$Classification,  levels = class_levels)
  )
  fisher_p <- tryCatch(
    fisher.test(cont_tab)$p.value,
    error = function(e) fisher.test(cont_tab, simulate.p.value = TRUE, B = 1e6)$p.value
  )
  p_label <- fmt_sci_p(fisher_p)

  # --- Plot -----------------------------------------------------------------
  p <- ggplot(plot_data, aes(x = Category, y = percentage, fill = Classification)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.6) +
    geom_text(
      aes(label = ifelse(percentage > 0, paste0(round(percentage, 1), "%"), "")),
      position = position_stack(vjust = 0.5),
      color    = "white", fontface = "bold", size = 5
    ) +
    geom_text(
      data        = n_labels,
      aes(x = Category, y = -5, label = label),
      inherit.aes = FALSE,
      fontface    = "bold", color = "black", size = 5
    ) +
    annotate(
      "text", x = 1.5, y = 108,
      label    = p_label,
      parse    = TRUE,
      fontface = "bold", size = 5.5
    ) +
    scale_fill_manual(
      name   = "Outcome",
      values = outcome_colors,
      breaks = c("Beneficial", "NS", "Adverse")
    ) +
    scale_y_continuous(
      breaks = c(0, 25, 50, 75, 100),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0.02, 0.15))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x        = NULL,
      y        = "Percentage",
      subtitle = paste0("Significance: ", p_type_label)
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position  = "top",
      legend.title     = element_text(face = "bold", size = 14),
      legend.text      = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(hjust = 0.5, size = 12),
      axis.title.y     = element_text(face = "bold", size = 16),
      axis.text.y      = element_text(face = "bold", size = 12, color = "black"),
      axis.text.x      = element_text(face = "bold", size = 14, color = "black"),
      axis.line.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      panel.border     = element_blank(),
      plot.margin      = margin(t = 45, r = 15, b = 25, l = 15)
    )

  ggsave(out_file, plot = p, width = 6, height = 7, dpi = 300, bg = "white")

  message(sprintf("\n[%s] classification counts:", p_type_label))
  print(addmargins(cont_tab))
  message(sprintf("[%s] Fisher exact P = %.3g  ->  %s", p_type_label, fisher_p, out_file))

  invisible(p)
}

# ----------------------------------------------------------------------------
# 5. Generate both versions
# ----------------------------------------------------------------------------
build_figure(df, "_p_value",  "raw p < 0.05",
             "Results/figure_cimic_survival_classification_multivariable_pvalue.png")

build_figure(df, "_p_adj_BH", "BH-adjusted p < 0.05",
             "Results/figure_cimic_survival_classification_multivariable_padj_BH.png")

message("\nDone.")
