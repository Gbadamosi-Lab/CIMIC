# Gbadamosi Lab - I.L on 8-21-26
# Complete NEO PAM50 Classification Pipeline
# Adapted from NKI pipeline
# NEO has one replicate per sample (_B_Rep1 = Baseline, _S_Rep1 = Surgery)
# Cluster IDs use pattern: "NEO11 _Patient_11_delta"
# Gene expression IDs use pattern: "NEO11 _Patient_11_B_Rep1"

# --- Packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("data.table","dplyr")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
for (p in c("genefu","org.Hs.eg.db","AnnotationDbi")) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("bhklab/genefu")

library(data.table)
library(dplyr)
library(genefu)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Load PAM50 data
data(pam50)
data(pam50.robust)

# ==============================================================================
# --- File Paths ---
# ==============================================================================

# CHANGED BLOCK: portable script-directory detection for both RStudio and Rscript.
# No new mandatory dependency is introduced: rstudioapi is used only when available.
get_script_dir <- function() {
  script_path <- NULL

  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    active_path <- tryCatch(
      rstudioapi::getActiveDocumentContext()$path,
      error = function(e) ""
    )

    if (nzchar(active_path)) {
      script_path <- active_path
    }
  }

  if (is.null(script_path)) {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)

    if (length(file_arg) > 0) {
      script_path <- sub("^--file=", "", file_arg[1])
    }
  }

  if (is.null(script_path) || !nzchar(script_path)) {
    stop(
      "Unable to determine the script location. ",
      "Run this file from RStudio with the script open, or via Rscript.",
      call. = FALSE
    )
  }

  dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
}

# CHANGED: anchor all repository paths to this script's own directory.
SCRIPT_DIR <- get_script_dir()

# CHANGED: figure_generation_scripts/ -> oncoimmunology_paper/
PAPER_DIR <- normalizePath(
  file.path(SCRIPT_DIR, ".."),
  winslash = "/",
  mustWork = TRUE
)

# CHANGED: repository input paths.
tpm_path <- normalizePath(
  file.path(PAPER_DIR, "Datasets", "NEO", "neo_tpm.csv"),
  winslash = "/",
  mustWork = TRUE
)

# CHANGED: cluster filename updated from neo_clusters.csv to neo_clustered_plot_df.csv.
cluster_path <- normalizePath(
  file.path(PAPER_DIR, "Datasets", "NEO", "neo_clustered_plot_df.csv"),
  winslash = "/",
  mustWork = TRUE
)

# CHANGED BLOCK: repository output root and existing output subfolder structure.
DIRS <- list(
  base = normalizePath(
    file.path(PAPER_DIR, "Results", "neo_PAM50"),
    winslash = "/",
    mustWork = FALSE
  )
)

DIRS$results <- file.path(DIRS$base, "results")
DIRS$plots   <- file.path(DIRS$base, "plots")
DIRS$logs    <- file.path(DIRS$base, "logs")

dir.create(DIRS$base,    recursive = TRUE, showWarnings = FALSE)
dir.create(DIRS$results, recursive = TRUE, showWarnings = FALSE)
dir.create(DIRS$plots,   recursive = TRUE, showWarnings = FALSE)
dir.create(DIRS$logs,    recursive = TRUE, showWarnings = FALSE)

# CHANGED: normalize created output directories consistently with winslash = "/".
DIRS$base    <- normalizePath(DIRS$base,    winslash = "/", mustWork = TRUE)
DIRS$results <- normalizePath(DIRS$results, winslash = "/", mustWork = TRUE)
DIRS$plots   <- normalizePath(DIRS$plots,   winslash = "/", mustWork = TRUE)
DIRS$logs    <- normalizePath(DIRS$logs,    winslash = "/", mustWork = TRUE)

# --- Read Data ---
cat("Loading data files...\n")
RNAseq_TPM       <- as.data.frame(fread(tpm_path, check.names = FALSE))
cluster_ICD_data <- read.csv(cluster_path, stringsAsFactors = FALSE, check.names = FALSE)

# Validate cluster data
stopifnot(all(c("sample_id","cluster_assignments") %in% names(cluster_ICD_data)))
cat("Data loaded:", nrow(RNAseq_TPM), "genes,", ncol(RNAseq_TPM)-1, "samples\n")
cat("Cluster data loaded:", nrow(cluster_ICD_data), "sample-level assignments\n")

# --- Extract Patient ID from cluster sample_id ---

# We extract the core patient code e.g. "NEO11" for matching
cluster_ICD_data <- cluster_ICD_data %>%
  dplyr::mutate(patient_id = trimws(sub("\\s*_Patient_.*$", "", sample_id)))

cat("Cluster patient IDs extracted:\n")
print(cluster_ICD_data[, c("sample_id", "patient_id", "cluster_assignments")])

# CHANGED BLOCK: current repository neo_tpm.csv contains an exported row-index
# column before gene_id. Remove only that index column and retain gene_id as Gene.
first_col_name <- names(RNAseq_TPM)[1]
if (first_col_name %in% c("", "V1", "...1", "X")) {
  first_col_values <- RNAseq_TPM[[1]]
  if (all(suppressWarnings(!is.na(as.numeric(first_col_values))))) {
    RNAseq_TPM[[1]] <- NULL
  }
}

# CHANGED: current repository TPM uses gene_id as the gene identifier column.
gene_col <- if ("gene_id" %in% names(RNAseq_TPM)) "gene_id" else names(RNAseq_TPM)[1]
cat("Gene column detected as:", gene_col, "\n")

# CHANGED BLOCK: normalize only sample-column names to the canonical NEO naming
# convention already expected by the original downstream regexes.
# Example:
#   NEO11 _Patient_11_T1_B_Rep1 -> NEO11_Patient_11_B_Rep1
#   NEO11 _Patient_11_T4_S_Rep1 -> NEO11_Patient_11_S_Rep1
sample_name_idx <- setdiff(seq_along(names(RNAseq_TPM)), match(gene_col, names(RNAseq_TPM)))
names(RNAseq_TPM)[sample_name_idx] <- names(RNAseq_TPM)[sample_name_idx] %>%
  trimws() %>%
  gsub("\\s+", "", .) %>%
  gsub("_T[0-9]+_(B|S)_Rep", "_\\1_Rep", .)

# --- Keep Only Baseline/Surgery Columns ---
sample_cols <- grep("_(B|S)_Rep1$", names(RNAseq_TPM), value = TRUE)
cat("Found", length(sample_cols), "B/S samples\n")

RNAseq_TPM <- RNAseq_TPM[, c(gene_col, sample_cols), drop = FALSE]
names(RNAseq_TPM)[1] <- "Gene"

# Make gene names unique before setting as rownames
RNAseq_TPM$Gene <- make.unique(RNAseq_TPM$Gene)
rownames(RNAseq_TPM) <- RNAseq_TPM$Gene
cat("After processing:", nrow(RNAseq_TPM), "genes (duplicates made unique)\n")

# --- Build Expression Matrix and Log Transform ---
tpm_mat <- as.matrix(RNAseq_TPM[, sample_cols, drop = FALSE])
mode(tpm_mat) <- "numeric"
tpm_log <- log2(tpm_mat + 1)

# --- Split by Condition ---
expr_B <- tpm_log[, grepl('_B_Rep1$', colnames(tpm_log)), drop = FALSE]
expr_S <- tpm_log[, grepl('_S_Rep1$', colnames(tpm_log)), drop = FALSE]
cat("Baseline samples:", ncol(expr_B), "| Surgery samples:", ncol(expr_S), "\n")

# --- PAM50 Classification Helper ---
run_pam50_neo <- function(expr_mat, cluster_data, label) {
  if (ncol(expr_mat) == 0) return(NULL)

  # Transpose: samples as rows, genes as columns
  texp_mat <- t(expr_mat)
  gene_info <- data.frame(Genes = colnames(texp_mat))
  rownames(gene_info) <- gene_info$Genes

  cat("Running PAM50 for", label, ":", nrow(texp_mat), "samples,", ncol(texp_mat), "genes\n")

  # PAM50 classification
  pam50_predictions <- molecular.subtyping(
    sbt.model = "pam50",
    data = texp_mat,
    annot = gene_info,
    do.mapping = FALSE
  )

  # Build output
  out <- data.frame(
    Sample = rownames(texp_mat),
    PAM50_Subtype = pam50_predictions$subtype,
    stringsAsFactors = FALSE
  )

  # Add probabilities and correlations
  if (!is.null(pam50_predictions$subtype.proba)) {
    out <- cbind(out, pam50_predictions$subtype.proba)
  }
  if (!is.null(pam50_predictions$correlations)) {
    colnames(pam50_predictions$correlations) <- paste0("cor_", colnames(pam50_predictions$correlations))
    out <- cbind(out, pam50_predictions$correlations)
  }

  # We extract "NEO11" (trimmed) to match cluster patient_id
  out <- out %>%
    dplyr::mutate(
      patient_id = trimws(sub("\\s*_Patient_.*$", "", Sample)),
      Condition  = ifelse(grepl('_B_Rep1$', Sample), 'Baseline', 'Surgery')
    ) %>%
    dplyr::left_join(
      cluster_data %>% dplyr::select(patient_id, sample_id, cluster_assignments),
      by = "patient_id"
    )

  # Check for unmatched samples
  unmatched <- sum(is.na(out$cluster_assignments))
  if (unmatched > 0) {
    cat("WARNING:", unmatched, "samples could not be matched to a cluster in", label, "\n")
    print(out[is.na(out$cluster_assignments), c("Sample", "patient_id")])
  }

  # Summary by cluster
  summary_by_cluster <- out %>%
    dplyr::count(Condition, cluster_assignments, PAM50_Subtype, name = "n") %>%
    dplyr::group_by(Condition, cluster_assignments) %>%
    dplyr::mutate(prop = round(n / sum(n), 3)) %>%
    dplyr::ungroup()

  # Save files
  # CHANGED: route result CSVs to Results/PAM50_neo/results/ without changing filenames.
  write.csv(
    out,
    file.path(DIRS$results, paste0("NEO_TPM_PAM50_", label, "_samples.csv")),
    row.names = FALSE
  )
  write.csv(
    summary_by_cluster,
    file.path(DIRS$results, paste0("NEO_TPM_PAM50_", label, "_by_cluster.csv")),
    row.names = FALSE
  )

  list(samples = out, summary = summary_by_cluster)
}

# --- Run PAM50 Classification ---
cat("\n=== Running PAM50 Classification ===\n")
pam50_B <- run_pam50_neo(expr_B, cluster_ICD_data, label = "Baseline")
pam50_S <- run_pam50_neo(expr_S, cluster_ICD_data, label = "Surgery")

# --- Combine Results ---
cat("\n=== Combining Results ===\n")
summary_combined <- dplyr::bind_rows(
  if (!is.null(pam50_B)) pam50_B$summary,
  if (!is.null(pam50_S)) pam50_S$summary
)

if (nrow(summary_combined) > 0) {
  # CHANGED: route combined result CSV to Results/PAM50_neo/results/.
  write.csv(
    summary_combined,
    file.path(DIRS$results, "NEO_TPM_PAM50_by_cluster_B_vs_S.csv"),
    row.names = FALSE
  )

  cat("PAM50 classification completed successfully!\n")
  cat("\n=== Summary by Cluster and Condition ===\n")
  print(summary_combined)

  cat("\n=== Overall PAM50 Distribution ===\n")
  overall_dist <- summary_combined %>%
    dplyr::group_by(Condition, PAM50_Subtype) %>%
    dplyr::summarise(total_n = sum(n), .groups = "drop") %>%
    dplyr::group_by(Condition) %>%
    dplyr::mutate(total_prop = round(total_n / sum(total_n), 3)) %>%
    dplyr::ungroup()
  print(overall_dist)

  cat("\n=== PAM50 Distribution by Cluster ===\n")
  cluster_dist <- summary_combined %>%
    dplyr::group_by(cluster_assignments, PAM50_Subtype) %>%
    dplyr::summarise(total_n = sum(n), .groups = "drop") %>%
    dplyr::group_by(cluster_assignments) %>%
    dplyr::mutate(total_prop = round(total_n / sum(total_n), 3)) %>%
    dplyr::ungroup()
  print(cluster_dist)

  cat("\n=== Files Saved ===\n")
  cat("1. NEO_TPM_PAM50_Baseline_samples.csv\n")
  cat("2. NEO_TPM_PAM50_Surgery_samples.csv\n")
  cat("3. NEO_TPM_PAM50_Baseline_by_cluster.csv\n")
  cat("4. NEO_TPM_PAM50_Surgery_by_cluster.csv\n")
  cat("5. NEO_TPM_PAM50_by_cluster_B_vs_S.csv\n")
} else {
  cat("ERROR: No results generated. Check data compatibility.\n")
}

##########################################################################################
# PAM50 Stacked Plot
##########################################################################################
library(ggplot2)
library(dplyr)
library(tidyr)

create_pam50_plot_neo <- function(pam50_data, condition_name) {

  plot_data <- pam50_data$samples %>%
    dplyr::filter(!is.na(cluster_assignments)) %>%  # remove unmatched samples
    dplyr::mutate(
      Cluster_Name = dplyr::case_when(
        cluster_assignments == 1 ~ "Dys-CIM",
        cluster_assignments == 2 ~ "Fun-CIM",
        TRUE ~ as.character(cluster_assignments)
      )
    ) %>%
    dplyr::group_by(Cluster_Name, PAM50_Subtype) %>%
    dplyr::summarise(total_n = n(), .groups = "drop") %>%
    dplyr::group_by(Cluster_Name) %>%
    dplyr::mutate(
      total_samples = sum(total_n),
      percentage = round((total_n / total_samples) * 100)
    ) %>%
    dplyr::ungroup()

  # Check if plot_data has rows
  if (nrow(plot_data) == 0) {
    cat("WARNING: No data to plot for", condition_name, "\n")
    return(NULL)
  }

  cat("Plot data for", condition_name, ":\n")
  print(plot_data)

  sample_sizes <- plot_data %>%
    dplyr::group_by(Cluster_Name) %>%
    dplyr::summarise(N = sum(total_n), .groups = "drop")

  subtype_colors <- c(
    "Basal"  = "#C44E52",
    "Her2"   = "#DD8DA6",
    "LumA"   = "#0072B2",
    "LumB"   = "#5C8C45",
    "Normal" = "#9467BD"
  )

  # Fisher's exact test
  p_value <- 1.0
  if (length(unique(plot_data$Cluster_Name)) > 1 && length(unique(plot_data$PAM50_Subtype)) > 1) {
    tryCatch({
      contingency_wide <- plot_data %>%
        dplyr::select(Cluster_Name, PAM50_Subtype, total_n) %>%
        tidyr::pivot_wider(names_from = PAM50_Subtype, values_from = total_n, values_fill = 0)
      contingency_table <- as.matrix(contingency_wide[,-1])
      rownames(contingency_table) <- contingency_wide$Cluster_Name
      fisher_test <- fisher.test(contingency_table)
      p_value <- fisher_test$p.value
    }, error = function(e) {
      cat("Statistical test failed:", conditionMessage(e), "\n")
      p_value <<- 1.0
    })
  }

  fish_stats_P <- paste0("P = ", format(round(p_value, 4), nsmall = 4))

  labels_df <- data.frame(
    x     = sample_sizes$Cluster_Name,
    y     = -8,
    label = paste0("N = ", sample_sizes$N)
  )

  outcomes_plot <- ggplot(plot_data,
                          aes(y = percentage,
                              x = factor(Cluster_Name, levels = c("Dys-CIM", "Fun-CIM")),
                              fill = PAM50_Subtype)) +
    geom_bar(stat = "identity", width = 0.85) +
    geom_text(aes(label = ifelse(percentage > 0, paste0(round(percentage, 1), "%"), "")),
              position = position_stack(vjust = 0.5),
              size = 6, color = "black", fontface = "bold") +
    labs(x = "", y = "Percentage", fill = "",
         title = paste0("PAM50 Subtype - ", condition_name)) +
    geom_text(data = labels_df, aes(x = x, y = y, label = label),
              size = 6, inherit.aes = FALSE, color = "black", fontface = "bold") +
    theme_classic() +
    theme(
      plot.title       = element_text(hjust = 0.5, size = 20, face = "bold"),
      axis.text        = element_text(colour = "black", size = 20, face = "bold"),
      axis.text.x      = element_text(colour = "black", size = 20, face = "bold"),
      axis.text.y      = element_text(colour = "black", size = 20, face = "bold"),
      axis.title       = element_text(colour = "black", size = 20, face = "bold"),
      axis.ticks       = element_line(linewidth = 1.5),
      legend.text      = element_text(colour = "black", size = 17, face = "bold"),
      legend.title     = element_text(colour = "black", size = 17, face = "bold"),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1),
      legend.position  = "right",
      legend.direction = "vertical"
    ) +
    guides(fill = guide_legend(
      title = NULL,
      keywidth  = unit(1, "cm"),
      keyheight = unit(0.8, "cm"),
      override.aes = list(color = "black", linewidth = 1.2)
    )) +
    annotate(geom = "text", x = 1.5, y = 110, label = fish_stats_P,
             size = 6, fontface = "bold") +
    scale_fill_manual(values = subtype_colors) +
    scale_y_continuous(limits = c(-15, 117), breaks = seq(0, 100, 20)) +
    scale_x_discrete(expand = expansion(mult = c(0.55, 0.55)))

  return(list(plot = outcomes_plot, data = plot_data, p_value = p_value))
}

# --- Generate and Save Plots ---
if (!is.null(pam50_B) && !is.null(pam50_S)) {

  baseline_result <- create_pam50_plot_neo(pam50_B, "Baseline")
  if (!is.null(baseline_result)) {
    print(baseline_result$plot)
    # CHANGED: route plot to Results/PAM50_neo/plots/ without changing filename or ggsave parameters.
    ggsave(
      file.path(DIRS$plots, "NEO_PAM50_Baseline_by_Cluster.png"),
      plot = baseline_result$plot, width = 5.5, height = 6, dpi = 300
    )
    cat("\n=== Baseline Sample Summary ===\n")
    print(baseline_result$data)
  }

  surgery_result <- create_pam50_plot_neo(pam50_S, "Surgery")
  if (!is.null(surgery_result)) {
    print(surgery_result$plot)
    # CHANGED: route plot to Results/PAM50_neo/plots/ without changing filename or ggsave parameters.
    ggsave(
      file.path(DIRS$plots, "NEO_PAM50_Surgery_by_Cluster.png"),
      plot = surgery_result$plot, width = 5.5, height = 6, dpi = 300
    )
    cat("\n=== Surgery Sample Summary ===\n")
    print(surgery_result$data)
  }

} else {
  cat("Error: PAM50 results not available\n")
  cat("pam50_B is NULL:", is.null(pam50_B), "\n")
  cat("pam50_S is NULL:", is.null(pam50_S), "\n")
}
