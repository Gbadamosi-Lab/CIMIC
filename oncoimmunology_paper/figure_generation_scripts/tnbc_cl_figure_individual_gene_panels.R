# ============================================================================
# tnbc_cl_figure_individual_gene_panels.R
# ============================================================================
# Purpose:
#   Per-gene two-panel summary across CIMIC trajectories for the TNBC_CL_Epirubicin dataset:
#
#     Left:
#       Induction boxplot of Delta(log2[TPM+1])
#       Statistical test: limma empirical-Bayes moderated t-test
#
#     Right:
#       Penetrance plot showing the proportion of samples with Delta > 0
#       Statistical test:
#         - Fisher's exact test if any expected cell count is < 5
#         - Pearson chi-squared test otherwise
#
# Cluster mapping for TNBC_CL_Epirubicin:
#   Cluster 1 = Dys-CIM
#   Cluster 2 = Fun-CIM
#
# limma contrast:
#   Fun-CIM - Dys-CIM
#
# Positive logFC:
#   Higher induction in Fun-CIM
#
# Negative logFC:
#   Higher induction in Dys-CIM
#
# Input:
#   oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/
#     tnbc_cl_epi_clustered_plot_df.csv
#
# Output:
#   oncoimmunology_paper/Results/tnbc_cl_epi_individual_gene_panels/figure_gene_panel_<GENE>.png
#
# Dependencies:
#   data.table, dplyr, tidyr, ggplot2, ggpubr, scales, limma
# ============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(limma)

# ---------------------------------------------------------------------------
# 2a. Helper: calculate average MHC‑I expression (Avg_MHC1)
# ---------------------------------------------------------------------------
calculate_avg_mhc1 <- function(df) {
  required_cols <- c("HLA-A", "HLA-B", "HLA-C")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing required MHC‑I columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  mhc_matrix <- as.matrix(df[required_cols])
  avg <- rowMeans(mhc_matrix, na.rm = TRUE)
  avg[is.nan(avg)] <- NA_real_
  avg
}

# ============================================================================
# 1. Paths and output directory
# ============================================================================

input_file <- paste0(
  "oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/",
  "tnbc_cl_epi_clustered_plot_df.csv"
)

output_dir <- "oncoimmunology_paper/Results/tnbc_cl_individual_gene_panels"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================================
# 2. Helper functions (plot parameters, formatting, etc.)
# ============================================================================
annotation_text_size <- 6
annotation_text_size_pointbiserial <- 6
axis_text_size_x <- 20
axis_text_size_y <- 20
plot_title_size_induction <- 24
plot_title_size_penetrance <- 24
plot_title_size_pointbiserial <- 24
axis_title_size_y <- 20
legend_text_size <- 18
legend_title_size <- 18
sample_size_text_size <- 7
show_sample_size <- TRUE
base_size_induction <- 18
base_size_penetrance <- 18
base_size_pointbiserial <- 18
main_title_size <- 24

format_p4 <- function(p) {
  if (length(p) != 1L || is.na(p) || !is.finite(p)) {
    return("NA")
  }
  if (p < 0.001) {
    return(formatC(p, format = "e", digits = 2))
  }
  if (p < 0.01) {
    return(formatC(p, format = "f", digits = 5))
  }
  if (p < 0.1) {
    return(formatC(p, format = "f", digits = 5))
  }
  formatC(p, format = "f", digits = 3)
}

format_p_annotation <- function(p, prefix = "P") {
  fp <- format_p4(p)
  if (fp == "NA") {
    return(paste0(prefix, " = NA"))
  }
  if (startsWith(fp, "<")) {
    return(paste0(prefix, " ", fp))
  }
  paste0(prefix, " = ", fp)
}

annotation_y <- function(x, fraction = 0.10) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(0)
  }
  limits <- range(x)
  spread <- diff(limits)
  if (!is.finite(spread) || spread == 0) {
    spread <- max(abs(limits), 1)
  }
  limits[[2]] + fraction * spread
}

find_gene_column <- function(df, gene) {
  escaped_gene <- gsub(
    "([][{}()+*^$|\\?.])",
    "\\\\1",
    gene
  )
  pattern <- paste0(
    "^",
    escaped_gene,
    "(\\.x|\\.y)?$"
  )
  matches <- grep(
    pattern,
    names(df),
    value = TRUE
  )
  if (length(matches) == 0L) {
    stop(
      "Gene column '",
      gene,
      "' was not found in the data frame."
    )
  }
  if (length(matches) > 1L) {
    warning(
      "Multiple columns matched gene '",
      gene,
      "': ",
      paste(matches, collapse = ", "),
      ". Using ",
      matches[[1]],
      "."
    )
  }
  matches[[1]]
}

# ============================================================================
# 3. Load data and define CIMIC labels
# ============================================================================

clustered_plot_df_for_merge <- data.table::fread(
  input_file,
  data.table = FALSE
) %>%
  as.data.frame()

# Compute Avg_MHC1 and add to dataframe
clustered_plot_df_for_merge$Avg_MHC1 <- calculate_avg_mhc1(clustered_plot_df_for_merge)

required_columns <- c(
  "sample_id",
  "cluster_assignments"
)
missing_required <- setdiff(
  required_columns,
  colnames(clustered_plot_df_for_merge)
)
if (length(missing_required) > 0L) {
  stop(
    "Required columns are missing: ",
    paste(missing_required, collapse = ", ")
  )
}

clustered_plot_df_for_merge <- clustered_plot_df_for_merge %>%
  dplyr::mutate(
    cluster_assignments = as.character(cluster_assignments),
    CIMIC_Cluster = factor(
      dplyr::case_when(
        cluster_assignments == "1" ~ "Fun-CIM",
        cluster_assignments == "2" ~ "Dys-CIM",
        TRUE ~ NA_character_
      ),
      levels = c(
        "Dys-CIM",
        "Fun-CIM"
      )
    )
  )

if (anyNA(clustered_plot_df_for_merge$CIMIC_Cluster)) {
  unexpected_values <- unique(
    clustered_plot_df_for_merge$cluster_assignments[
      is.na(clustered_plot_df_for_merge$CIMIC_Cluster)
    ]
  )
  stop(
    "Unexpected or missing cluster assignments: ",
    paste(unexpected_values, collapse = ", "),
    ". Expected only 1 or 2."
  )
}

cluster_counts <- table(
  clustered_plot_df_for_merge$CIMIC_Cluster
)
if (any(cluster_counts == 0L)) {
  stop(
    "Both Dys-CIM and Fun-CIM must contain at least one sample."
  )
}
message(
  sprintf(
    "Samples per cluster: %s",
    paste(
      names(cluster_counts),
      as.integer(cluster_counts),
      sep = "=",
      collapse = ", "
    )
  )
)

# ============================================================================
# 4. Identify gene‑expression columns
# ============================================================================
metadata_columns <- c(
  "sample_id",
  "cluster_assignments",
  "CIMIC_Cluster",
  "base_id",
  "PACMAP1",
  "PACMAP2",
  "UMAP1",
  "UMAP2",
  "PC1",
  "PC2",
  "PC3",
  "V1"
)

candidate_gene_columns <- setdiff(
  colnames(clustered_plot_df_for_merge),
  metadata_columns
)
if (length(candidate_gene_columns) == 0L) {
  stop("No candidate gene‑expression columns were found.")
}

clustered_plot_df_for_merge[
  candidate_gene_columns
] <- lapply(
    clustered_plot_df_for_merge[candidate_gene_columns],
    function(x) {
      suppressWarnings(
        as.numeric(as.character(x))
      )
    }
  )

valid_gene_column <- vapply(
  clustered_plot_df_for_merge[candidate_gene_columns],
  function(x) {
    finite_x <- x[is.finite(x)]
    length(finite_x) >= 2L && length(unique(finite_x)) >= 2L
  },
  logical(1)
)

gene_columns <- candidate_gene_columns[valid_gene_column]
if (length(gene_columns) == 0L) {
  stop(
    "No variable numeric gene‑expression columns remained after filtering."
  )
}
message(
  sprintf(
    "Using %d gene‑expression columns for the limma model.",
    length(gene_columns)
  )
)
# ============================================================================
# 5. Extract cell-line IDs from sample IDs
# ============================================================================
# Sample IDs are structured as:
#
#   delta_BT549_EPIRUBICIN_1
#   delta_BT549_EPIRUBICIN_2
#   delta_BT549_EPIRUBICIN_3
#
# The cell line is everything between "delta_" and "_EPIRUBICIN_<replicate>".

clustered_plot_df_for_merge <- clustered_plot_df_for_merge %>%
  dplyr::mutate(
    cell_line_id = sub(
      "^delta_(.*)_EPIRUBICIN_[0-9]+$",
      "\\1",
      as.character(sample_id)
    )
  )

# Check extraction
message("Cell-line assignments:")

print(
  clustered_plot_df_for_merge %>%
    dplyr::count(cell_line_id, CIMIC_Cluster)
)

# Verify that all sample IDs produced a plausible cell-line ID
if (
  any(
    is.na(clustered_plot_df_for_merge$cell_line_id) |
    clustered_plot_df_for_merge$cell_line_id == clustered_plot_df_for_merge$sample_id
  )
) {
  warning(
    "One or more sample IDs may not have been parsed correctly. ",
    "Check the cell_line_id column."
  )

  print(
    clustered_plot_df_for_merge %>%
      dplyr::select(
        sample_id,
        cell_line_id
      )
  )
}


# ============================================================================
# 6. Fit replicate-aware limma model using duplicateCorrelation
# ============================================================================
#
# Model:
#
#   Expression ~ CIMIC_Cluster
#
# with:
#
#   block = cell_line_id
#
# This accounts for correlation among replicate measurements from the
# same cell line.
#
# The workflow is:
#
#   1. Initial lmFit()
#   2. Estimate within-cell-line correlation using duplicateCorrelation()
#   3. Refit lmFit() using the consensus correlation
#   4. Apply empirical Bayes moderation with robust = TRUE
#
# robust = TRUE protects the empirical-Bayes variance moderation from
# genes with unusually large or influential variances.
# ============================================================================

# limma expects:
#   rows    = genes
#   columns = samples

expression_matrix <- t(
  as.matrix(
    clustered_plot_df_for_merge[
      ,
      gene_columns,
      drop = FALSE
    ]
  )
)

storage.mode(expression_matrix) <- "double"

rownames(expression_matrix) <- gene_columns

colnames(expression_matrix) <- make.unique(
  as.character(
    clustered_plot_df_for_merge$sample_id
  )
)

# lmFit cannot use Inf/-Inf values.
expression_matrix[!is.finite(expression_matrix)] <- NA_real_


# ============================================================================
# 6A. Define CIMIC cluster design
# ============================================================================

# Ensure that the first column (C1) corresponds to Fun-CIM and the second (C2) to Dys-CIM
cluster_factor <- factor(
  clustered_plot_df_for_merge$CIMIC_Cluster,
  levels = c(
    "Fun-CIM",
    "Dys-CIM"
  )
)

design <- stats::model.matrix(
  ~ 0 + cluster_factor
)

colnames(design) <- c(
  "Fun_CIM",
  "Dys_CIM"
)


# ============================================================================
# 6B. Define Fun-CIM vs Dys-CIM contrast
# ============================================================================

contrast_matrix <- limma::makeContrasts(
  Fun_CIM_vs_Dys_CIM = Fun_CIM - Dys_CIM,
  levels = design
)


# ============================================================================
# 6C. Initial model for duplicateCorrelation
# ============================================================================

fit_initial <- limma::lmFit(
  expression_matrix,
  design
)


# ============================================================================
# 6D. Estimate within-cell-line correlation
# ============================================================================

cell_line_id <- factor(
  clustered_plot_df_for_merge$cell_line_id
)

corfit <- limma::duplicateCorrelation(
  expression_matrix,
  design = design,
  block = cell_line_id
)

message(
  sprintf(
    "Estimated duplicateCorrelation consensus correlation: %.4f",
    corfit$consensus
  )
)


# ============================================================================
# 6E. Refit model accounting for cell-line correlation
# ============================================================================

limma_fit <- limma::lmFit(
  expression_matrix,
  design,
  block = cell_line_id,
  correlation = corfit$consensus
)


# ============================================================================
# 6F. Apply Fun-CIM vs Dys-CIM contrast
# ============================================================================

limma_fit <- limma::contrasts.fit(
  limma_fit,
  contrast_matrix
)


# ============================================================================
# 6G. Empirical-Bayes moderation
# ============================================================================
#
# robust = TRUE:
#   Makes the empirical-Bayes variance moderation less sensitive to
#   genes with unusually large variances.
#
# trend = TRUE:
#   Allows the prior variance to depend on average expression.
#
# Both are appropriate for this expression-level analysis.
# ============================================================================

limma_fit <- limma::eBayes(
  limma_fit,
  trend = TRUE,
  robust = TRUE
)


# ============================================================================
# 6H. Extract complete limma results
# ============================================================================

limma_results <- limma::topTable(
  limma_fit,
  coef = "Fun_CIM_vs_Dys_CIM",
  number = Inf,
  sort.by = "none"
)

limma_results$gene_column <- rownames(
  limma_results
)


# Store results in original gene-column order
limma_results <- limma_results[
  match(
    gene_columns,
    limma_results$gene_column
  ),
  ,
  drop = FALSE
]

rownames(limma_results) <- limma_results$gene_column


# ============================================================================
# 6I. Report model information
# ============================================================================

message(
  paste0(
    "Fitted replicate-aware limma model with contrast: ",
    "Fun-CIM minus Dys-CIM."
  )
)

message(
  paste0(
    "Blocking variable: cell_line_id"
  )
)

message(
  sprintf(
    "Number of samples: %d",
    ncol(expression_matrix)
  )
)

message(
  sprintf(
    "Number of cell lines: %d",
    length(unique(cell_line_id))
  )
)

message(
  sprintf(
    "Consensus within-cell-line correlation: %.4f",
    corfit$consensus
  )
)

message(
  "Empirical-Bayes moderation: trend = TRUE, robust = TRUE"
)


# ============================================================================
# 7. Extract limma statistics for one gene
# ============================================================================

get_limma_statistics <- function(gene_column) {

  if (!gene_column %in% rownames(limma_results)) {
    stop(
      "Gene column '",
      gene_column,
      "' was not included in the limma model."
    )
  }

  result <- limma_results[
    gene_column,
    ,
    drop = FALSE
  ]

  list(
    logFC = as.numeric(
      result$logFC
    ),

    moderated_t = as.numeric(
      result$t
    ),

    p_value = as.numeric(
      result$P.Value
    ),

    adjusted_p_value = as.numeric(
      result$adj.P.Val
    )
  )
}

# ============================================================================
# 7. Induction boxplot using limma p‑value
# ============================================================================
compare_induction_boxplot <- function(df, gene) {
  gene_column <- find_gene_column(df, gene)
  limma_stat <- get_limma_statistics(gene_column)
  df2 <- df %>%
    dplyr::transmute(
      CIMIC_Cluster = CIMIC_Cluster,
      Expression = as.numeric(.data[[gene_column]])
    ) %>%
    dplyr::filter(
      !is.na(CIMIC_Cluster),
      is.finite(Expression)
    )
  if (nrow(df2) == 0L) {
    stop(
      "No finite observations were available for gene ",
      gene,
      "."
    )
  }
  if (dplyr::n_distinct(df2$CIMIC_Cluster) != 2L) {
    stop(
      "Both CIMIC clusters require finite observations for gene ",
      gene,
      "."
    )
  }
  p_text <- paste0(
    "P = ",
    format_p4(limma_stat$p_value)
  )
  y_position <- annotation_y(df2$Expression, fraction = 0.12)
  ggplot(
    df2,
    aes(
      x = CIMIC_Cluster,
      y = Expression,
      fill = CIMIC_Cluster
    )
  ) +
    geom_boxplot(
      width = 0.6,
      outlier.shape = NA,
      alpha = 0.8
    ) +
    geom_jitter(
      width = 0.15,
      height = 0,
      size = 1.5,
      alpha = 1
    ) +
    annotate(
      geom = "text",
      x = 1.5,
      y = y_position,
      label = p_text,
      size = annotation_text_size,
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = c(
        "Dys-CIM" = "#2980B9",
        "Fun-CIM" = "#C0392B"
      ),
      drop = FALSE
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(0.05, 0.20)
      )
    ) +
    labs(
      title = gene,
      y = expression(Delta * "(log"[2] * "[TPM+1])"),
      x = NULL
    ) +
    theme_classic(
      base_size = base_size_induction
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = plot_title_size_induction
      ),
      legend.position = "none",
      axis.text.x = element_text(
        face = "bold",
        size = axis_text_size_x
      ),
      axis.text.y = element_text(
        face = "bold",
        size = axis_text_size_y
      ),
      axis.title.y = element_text(
        face = "bold",
        size = axis_title_size_y
      ),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1.5
      ),
      axis.line = element_line(
        colour = "black",
        linewidth = 1.5
      )
    )
}

# ============================================================================
# 8. Penetrance test and stacked‑bar plot
# ============================================================================
compare_penetrance_plot <- function(df, gene) {
  gene_column <- find_gene_column(df, gene)
  df2 <- df %>%
    dplyr::transmute(
      CIMIC_Cluster = CIMIC_Cluster,
      Expression = as.numeric(.data[[gene_column]])
    ) %>%
    dplyr::filter(
      !is.na(CIMIC_Cluster),
      is.finite(Expression)
    ) %>%
    dplyr::mutate(
      Induction = ifelse(
        Expression > 0,
        "Positive Induction",
        "Negative/No Induction"
      ),
      Induction = factor(
        Induction,
        levels = c(
          "Negative/No Induction",
          "Positive Induction"
        )
      )
    )
  if (nrow(df2) == 0L) {
    stop(
      "No finite observations were available for gene ",
      gene,
      "."
    )
  }
  count_df <- df2 %>%
    dplyr::count(
      CIMIC_Cluster,
      Induction,
      name = "n"
    ) %>%
    tidyr::complete(
      CIMIC_Cluster = factor(
        c("Dys-CIM", "Fun-CIM"),
        levels = c("Dys-CIM", "Fun-CIM")
      ),
      Induction = factor(
        c(
          "Negative/No Induction",
          "Positive Induction"
        ),
        levels = c(
          "Negative/No Induction",
          "Positive Induction"
        )
      ),
      fill = list(n = 0L)
    )
  penetrance_matrix <- count_df %>%
    tidyr::pivot_wider(
      names_from = Induction,
      values_from = n,
      values_fill = 0
    ) %>%
    dplyr::arrange(CIMIC_Cluster)
  row_names <- as.character(penetrance_matrix$CIMIC_Cluster)
  penetrance_matrix <- as.matrix(
    penetrance_matrix[
      ,
      c(
        "Negative/No Induction",
        "Positive Induction"
      ),
      drop = FALSE
    ]
  )
  rownames(penetrance_matrix) <- row_names
  valid_table <- (
    nrow(penetrance_matrix) == 2L &&
      ncol(penetrance_matrix) == 2L &&
      all(rowSums(penetrance_matrix) > 0L) &&
      all(colSums(penetrance_matrix) > 0L)
  )
  if (!valid_table) {
    p_value <- NA_real_
    test_name <- "Not testable"
  } else {
    chi_check <- suppressWarnings(
      stats::chisq.test(
        penetrance_matrix,
        correct = FALSE
      )
    )
    expected_counts <- chi_check$expected
    if (any(expected_counts < 5)) {
      test_result <- stats::fisher.test(penetrance_matrix)
      p_value <- test_result$p.value
      test_name <- "Fisher"
    } else {
      test_result <- stats::chisq.test(penetrance_matrix, correct = FALSE)
      p_value <- test_result$p.value
      test_name <- expression(chi^2)
    }
  }
  if (identical(test_name, "Fisher")) {
    p_text <- paste0("Fisher P = ", format_p4(p_value))
  } else if (identical(test_name, "Not testable")) {
    p_text <- "P = NA"
  } else {
    p_text <- paste0("\u03C7\u00B2 P = ", format_p4(p_value))
  }
  bar_df <- count_df %>%
    dplyr::group_by(CIMIC_Cluster) %>%
    dplyr::mutate(
      Total = sum(n),
      Prop = ifelse(Total > 0, n / Total, NA_real_)
    ) %>%
    dplyr::ungroup()
  sample_size_df <- bar_df %>%
    dplyr::distinct(CIMIC_Cluster, Total)
  ggplot(
    bar_df,
    aes(
      x = CIMIC_Cluster,
      y = Prop,
      fill = Induction
    )
  ) +
    geom_col(
      width = 0.7,
      colour = "black",
      linewidth = 0.3
    ) +
    geom_text(
      aes(
        label = ifelse(
          Prop > 0,
          paste0(round(Prop * 100), "%"),
          ""
        )
      ),
      position = position_stack(vjust = 0.5),
      colour = "black",
      size = 7,
      fontface = "bold"
    ) +
    annotate(
      geom = "text",
      x = 1.5,
      y = 1.07,
      label = p_text,
      size = annotation_text_size,
      fontface = "bold"
    ) +
    { if (show_sample_size) {
        geom_text(
          data = sample_size_df,
          aes(
            x = CIMIC_Cluster,
            y = -0.05,
            label = paste0("N = ", Total)
          ),
          inherit.aes = FALSE,
          size = sample_size_text_size,
          fontface = "bold"
        )
      } else {
        NULL
      } } +
    scale_fill_manual(
      values = c(
        "Positive Induction" = "#3fb949ff",
        "Negative/No Induction" = "#f76661ff"
      ),
      drop = FALSE
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      limits = c(-0.10, 1.14),
      breaks = seq(0, 1, 0.2),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_x_discrete(
      expand = expansion(add = 0.5),
      drop = FALSE
    ) +
    labs(
      title = gene,
      y = "Percentage",
      x = NULL,
      fill = "Induction"
    ) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = plot_title_size_penetrance
      ),
      legend.position = "top",
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1.5
      ),
      axis.line.x = element_line(colour = "black", linewidth = 1.5),
      axis.line.y = element_line(colour = "black", linewidth = 1.5),
      axis.text.x = element_text(face = "bold", size = axis_text_size_x),
      axis.text.y = element_text(face = "bold", size = axis_text_size_y),
      axis.title.y = element_text(face = "bold", size = axis_title_size_y),
      legend.text = element_text(size = legend_text_size),
      legend.title = element_text(size = legend_title_size, face = "bold")
    )
}

# ============================================================================
# 9. Optional point‑biserial correlation plot
# ============================================================================
compare_pointbiserial_plot <- function(df, gene) {
  gene_column <- find_gene_column(df, gene)
  df2 <- df %>%
    dplyr::transmute(
      CIMIC_Cluster = CIMIC_Cluster,
      Cluster_bin = ifelse(CIMIC_Cluster == "Fun-CIM", 1, 0),
      Expression = as.numeric(.data[[gene_column]])
    ) %>%
    dplyr::filter(
      !is.na(CIMIC_Cluster),
      is.finite(Expression)
    )
  if (nrow(df2) < 3L || length(unique(df2$Cluster_bin)) != 2L || stats::sd(df2$Expression) == 0) {
    correlation <- NA_real_
    p_value <- NA_real_
  } else {
    correlation_result <- stats::cor.test(df2$Cluster_bin, df2$Expression, method = "pearson")
    correlation <- unname(correlation_result$estimate)
    p_value <- correlation_result$p.value
  }
  correlation_text <- paste0(
    "r = ",
    ifelse(is.na(correlation), "NA", sprintf("%.3f", correlation)),
    "\nP = ",
    format_p4(p_value)
  )
  y_position <- annotation_y(df2$Expression, fraction = 0.12)
  ggplot(
    df2,
    aes(
      x = Cluster_bin,
      y = Expression,
      colour = CIMIC_Cluster
    )
  ) +
    geom_jitter(width = 0.1, height = 0, size = 2, alpha = 0.6) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 1.2) +
    scale_x_continuous(
      breaks = c(0, 1),
      labels = c("Dys-CIM", "Fun-CIM")
    ) +
    scale_colour_manual(
      values = c(
        "Dys-CIM" = "#2980B9",
        "Fun-CIM" = "#C0392B"
      ),
      drop = FALSE
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    annotate(
      geom = "text",
      x = 0.5,
      y = y_position,
      label = correlation_text,
      size = annotation_text_size_pointbiserial,
      fontface = "bold"
    ) +
    labs(
      title = gene,
      x = NULL,
      y = expression(Delta * "(log"[2] * "[TPM+1])")
    ) +
    theme_classic(base_size = base_size_pointbiserial) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = plot_title_size_pointbiserial
      ),
      legend.position = "none",
      axis.text.x = element_text(face = "bold", size = axis_text_size_x),
      axis.text.y = element_text(face = "bold", size = axis_text_size_y),
      axis.title.y = element_text(face = "bold", size = axis_title_size_y),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
      axis.line = element_line(colour = "black", linewidth = 1.5)
    )
}

# ============================================================================
# 10. Panel combiners
# ============================================================================
make_three_panel <- function(gene) {
  induction_plot <- compare_induction_boxplot(clustered_plot_df_for_merge, gene)
  penetrance_plot <- compare_penetrance_plot(clustered_plot_df_for_merge, gene)
  pointbiserial_plot <- compare_pointbiserial_plot(clustered_plot_df_for_merge, gene)
  induction_plot <- induction_plot + theme(plot.title = element_blank())
  penetrance_plot <- penetrance_plot + theme(plot.title = element_blank(), axis.title.x = element_blank())
  pointbiserial_plot <- pointbiserial_plot + theme(plot.title = element_blank())
  combined <- ggpubr::ggarrange(
    induction_plot,
    penetrance_plot,
    pointbiserial_plot,
    ncol = 3,
    nrow = 1,
    common.legend = FALSE
  )
  ggpubr::annotate_figure(
    combined,
    top = ggpubr::text_grob(gene, face = "bold", size = main_title_size)
  )
}

make_two_panel <- function(gene) {
  induction_plot <- compare_induction_boxplot(clustered_plot_df_for_merge, gene)
  penetrance_plot <- compare_penetrance_plot(clustered_plot_df_for_merge, gene)
  induction_plot <- induction_plot + theme(plot.title = element_blank())
  penetrance_plot <- penetrance_plot + theme(plot.title = element_blank(), axis.title.x = element_blank(), legend.position = "none")
  combined <- ggpubr::ggarrange(
    induction_plot,
    penetrance_plot,
    ncol = 2,
    nrow = 1,
    common.legend = FALSE
  )
  ggpubr::annotate_figure(
    combined,
    top = ggpubr::text_grob(gene, face = "bold", size = main_title_size)
  )
}

# ============================================================================
# 11. Build and save figures
# ============================================================================


pyroptosis <- c(
 "AIM2", "CASP1", "CASP4", "CASP5", "GSDMD", "GSDME", "IL1B", "IL18", "NLRC4", "NLRP3", "PYCARD"
)

ferroptosis <- c(
 "ACSL4", "LPCAT3", "PTGS2", "SAT1", "TFRC"
)

panoptosis <- c(
 "AIM2", "CASP8", "FADD", "IRF1", "RIPK1", "RIPK3", "PYCARD", "ZBP1"
)

necroptosis <- c(
 "MLKL", "RIPK1", "RIPK3", "TICAM1", "ZBP1"
)

apoptosis <- c(
 "APAF1", "BAD", "BAK1", "BAX", "BBC3", "BCL2L11", "BID", "BMF", 
 "CASP3", "CASP7", "CASP8", "CASP9", "CYCS", "DIABLO", "FADD", "FAS", 
 "HRK", "PMAIP1", "TNFRSF10B", "TP53AIP1"
)

nki_smc_cim_sensors <- c(
  "IFI16", "TLR9", "TLR3", "TLR7", "TLR8", "CGAS",
  "P2RX7", "FPR1", "TLR4", "IFNAR1", "IFNAR2", "IFNGR1", "STING1"
)

nki_smc_cim_checkpoints_and_immuno <- c(
"NFKB1", "Avg_MHC1", "PDCD1", "CD274", "CTLA4", 
"LAG3", "BTLA", "HAVCR2", "TIGIT", "VSIR", "VTCN1", "HHLA2"
)

icd_genes <- c(
  "CALR", "HMGB1", "EIF2AK3", "ATF4", "PDIA3", "ANXA1", "EIF2A"
)


tumor_cell_checkpoints <- c(
  # B7 Family Ligands (Tumor-expressed)
  "CD274",     # PD-L1
  "PDCD1LG2",  # PD-L2
  "VTCN1",     # B7-H4
  "HHLA2",     # B7-H7
  "VSIR",      # VISTA
  "CD276",     # B7-H3

  # TIGIT / DNAM-1 Ligands
  "PVR",       # CD155
  "NECTIN2",   # CD112

  # LAG3 Ligand
  "FGL1",      # Fibrinogen-like protein 1

  # CD200 Axis & Other Tumor-expressed Checkpoints
  "CD200",     # OX-2
  "LGALS9",     # Galectin-9 (TIM-3 ligand)

  # autonomous tumor cell checkpoints (not immune checkpoints)
  "CTLA4", "LAG3", "TIGIT", "HAVCR2", "BTLA"
)

tnbc_cl_nucleic_acid_sensors <- c(
  "CGAS", "STING1", "AIM2", "IFIH1", "DHX58", "RIGI"
)

genes_of_interest <- tnbc_cl_nucleic_acid_sensors

for (gene in genes_of_interest) {
  gene_column <- find_gene_column(clustered_plot_df_for_merge, gene)
  limma_stat <- get_limma_statistics(gene_column)
  message(
    sprintf(
      "%s: limma logFC(Fun-Dys) = %.4f, moderated t = %.4f, P = %s, adjusted P = %s",
      gene,
      limma_stat$logFC,
      limma_stat$moderated_t,
      format_p4(limma_stat$p_value),
      format_p4(limma_stat$adjusted_p_value)
    )
  )
  figure <- make_two_panel(gene)
  output_file <- file.path(
    output_dir,
    sprintf("tnbc_cl_figure_gene_panel_%s.png", gene)
  )
  ggplot2::ggsave(
    filename = output_file,
    plot = figure,
    width = 7.5,
    height = 5,
    dpi = 300,
    bg = "white"
  )
  message(sprintf("Saved: %s", output_file))
}

message("Done.")
