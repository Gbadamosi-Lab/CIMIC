# ==============================================================================
# NEO PAM50 CLASSIFICATION PIPELINE
# genefu PAM50 classification + CIM cluster distribution plots
#
# Adapted from the original NEO/NKI PAM50 workflow
# Gbadamosi Lab - I.L.
# Updated: 08/17/2026
#
# Expected repository structure:
#   NEO_PAM50/
#   ├── data/
#   │   ├── neo_tpm.csv
#   │   └── neo_clusters.csv
#   ├── results/      # created automatically
#   ├── plots/        # created automatically
#   └── logs/         # created automatically
#
# The new NEO TPM input uses sample IDs such as:
#   NEO11 _Patient_11_T1_B_Rep1
#   NEO11 _Patient_11_T4_S_Rep1
#
# These are harmonized at import to:
#   NEO11_Patient_11_B_Rep1
#   NEO11_Patient_11_S_Rep1
#
# PAM50 classification logic, cluster mapping, statistical test, plot design,
# and output filenames are retained from the original workflow.
# ==============================================================================

# ==============================================================================
# 0) REPRODUCIBILITY + GLOBAL OPTIONS
# ==============================================================================

set.seed(2025)
options(stringsAsFactors = FALSE, scipen = 999)

# ==============================================================================
# 1) PROJECT CONFIGURATION + INPUT PATHS
# ==============================================================================

PROJECT_DIR <- Sys.getenv("NEO_PAM50_DIR", unset = ".")
PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = FALSE)

DATA_DIR <- file.path(PROJECT_DIR, "data")

tpm_path     <- file.path(DATA_DIR, "neo_tpm.csv")
cluster_path <- file.path(DATA_DIR, "neo_clusters.csv")

DIRS <- list(
  base    = PROJECT_DIR,
  results = file.path(PROJECT_DIR, "results"),
  plots   = file.path(PROJECT_DIR, "plots"),
  logs    = file.path(PROJECT_DIR, "logs")
)

invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

required_input_files <- c(tpm_path, cluster_path)
missing_input_files <- required_input_files[!file.exists(required_input_files)]

if (length(missing_input_files) > 0) {
  stop(
    paste0(
      "Required input file(s) not found:\n  ",
      paste(missing_input_files, collapse = "\n  "),
      "\n\nRun the script from the NEO PAM50 project root, ",
      "or set NEO_PAM50_DIR to the project directory."
    ),
    call. = FALSE
  )
}

# ==============================================================================
# 2) LOGGING
# ==============================================================================

LOG_FILE <- file.path(
  DIRS$logs,
  paste0("run_log_NEO_PAM50_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

sink(LOG_FILE, split = TRUE)

cat("=== NEO PAM50 pipeline started ===\n")
cat("Timestamp:", format(Sys.time()), "\n")
cat("Project directory:", PROJECT_DIR, "\n")
cat("TPM input:", tpm_path, "\n")
cat("Cluster input:", cluster_path, "\n\n")

# ==============================================================================
# 3) PACKAGES
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_pkgs <- c("data.table", "dplyr", "ggplot2", "tidyr")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}

bioc_pkgs <- c("genefu", "org.Hs.eg.db", "AnnotationDbi")

for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(genefu)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(tidyr)
})

# Load PAM50 model data used by genefu.
data(pam50)
data(pam50.robust)

# Record package versions used for the run.
package_names <- c(
  "R",
  "data.table",
  "dplyr",
  "ggplot2",
  "tidyr",
  "genefu",
  "org.Hs.eg.db",
  "AnnotationDbi"
)

package_versions <- data.frame(
  package = package_names,
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    vapply(
      package_names[-1],
      function(p) as.character(utils::packageVersion(p)),
      character(1)
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  package_versions,
  file.path(DIRS$logs, "NEO_PAM50_package_versions.csv"),
  row.names = FALSE
)

# Input manifest.
input_manifest <- data.frame(
  input = c("TPM expression", "CIM cluster assignments"),
  file = normalizePath(c(tpm_path, cluster_path), winslash = "/", mustWork = TRUE),
  size_bytes = file.info(c(tpm_path, cluster_path))$size,
  modified = as.character(file.info(c(tpm_path, cluster_path))$mtime),
  stringsAsFactors = FALSE
)

write.csv(
  input_manifest,
  file.path(DIRS$logs, "NEO_PAM50_input_manifest.csv"),
  row.names = FALSE
)

# ==============================================================================
# 4) NEO SAMPLE-ID HARMONIZATION
# ==============================================================================

clean_sample_ids <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s+", "", x)

  # New TPM nomenclature:
  # *_T1_B_Rep1 -> *_B_Rep1
  # *_T4_S_Rep1 -> *_S_Rep1
  #
  # The general pattern also handles other T-number labels if introduced later.
  x <- sub("_T[0-9]+_(B|S)_Rep([0-9]+)$", "_\\1_Rep\\2", x)

  x
}

get_condition <- function(sample_id) {
  ifelse(
    grepl("_B_Rep[0-9]+$", sample_id),
    "Baseline",
    ifelse(grepl("_S_Rep[0-9]+$", sample_id), "Surgery", NA_character_)
  )
}

# Preserve the original NEO matching philosophy: use the leading NEO patient
# code (e.g., NEO11) to link expression samples to the cluster table.
get_neo_patient_id <- function(x) {
  x <- clean_sample_ids(x)
  sub("_Patient_.*$", "", x)
}

# ==============================================================================
# 5) READ + VALIDATE INPUT DATA
# ==============================================================================

cat("Loading data files...\n")

RNAseq_TPM <- as.data.frame(
  data.table::fread(tpm_path, data.table = FALSE, check.names = FALSE)
)

cluster_ICD_data <- read.csv(
  cluster_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Validate cluster data.
required_cluster_cols <- c("sample_id", "cluster_assignments")
missing_cluster_cols <- setdiff(required_cluster_cols, names(cluster_ICD_data))

if (length(missing_cluster_cols) > 0) {
  stop(
    "Cluster file is missing required column(s): ",
    paste(missing_cluster_cols, collapse = ", ")
  )
}

# Remove accidental CSV index columns such as "Unnamed: 0" or "X".
index_cols <- names(RNAseq_TPM)[
  grepl("^Unnamed:", names(RNAseq_TPM)) |
    names(RNAseq_TPM) %in% c("X")
]

if (length(index_cols) > 0) {
  cat(
    "Removing non-biological CSV index column(s):",
    paste(index_cols, collapse = ", "),
    "\n"
  )
  RNAseq_TPM[index_cols] <- NULL
}

# Explicitly identify the gene-symbol column.
gene_col_candidates <- c(
  "gene_id", "gene", "Gene", "symbol", "SYMBOL",
  "hgnc_symbol", "GeneSymbol"
)

gene_col <- intersect(gene_col_candidates, names(RNAseq_TPM))

if (length(gene_col) == 0) {
  stop(
    "Could not identify a gene-symbol column. Expected one of: ",
    paste(gene_col_candidates, collapse = ", ")
  )
}

gene_col <- gene_col[[1]]
cat("Gene column detected as:", gene_col, "\n")

# Harmonize sample-column names before any filtering.
sample_candidate_cols <- setdiff(names(RNAseq_TPM), gene_col)
canonical_sample_names <- clean_sample_ids(sample_candidate_cols)

if (anyDuplicated(canonical_sample_names)) {
  duplicated_ids <- unique(
    canonical_sample_names[
      duplicated(canonical_sample_names) |
        duplicated(canonical_sample_names, fromLast = TRUE)
    ]
  )

  stop(
    "Sample-ID harmonization created duplicate canonical sample IDs: ",
    paste(duplicated_ids, collapse = ", ")
  )
}

names(RNAseq_TPM)[match(sample_candidate_cols, names(RNAseq_TPM))] <-
  canonical_sample_names

# Normalize cluster IDs and derive patient IDs.
cluster_ICD_data <- cluster_ICD_data %>%
  dplyr::mutate(
    sample_id = clean_sample_ids(sample_id),
    patient_id = get_neo_patient_id(sample_id)
  )

if (anyDuplicated(cluster_ICD_data$patient_id)) {
  dup_cluster_patients <- unique(
    cluster_ICD_data$patient_id[
      duplicated(cluster_ICD_data$patient_id) |
        duplicated(cluster_ICD_data$patient_id, fromLast = TRUE)
    ]
  )

  stop(
    "Cluster file contains duplicate NEO patient IDs after harmonization: ",
    paste(dup_cluster_patients, collapse = ", ")
  )
}

# Keep only biological biopsy/surgery columns.
sample_cols <- grep("_(B|S)_Rep[0-9]+$", names(RNAseq_TPM), value = TRUE)

if (length(sample_cols) == 0) {
  stop(
    "No biopsy/surgery expression columns were detected after ID harmonization."
  )
}

cat("Data loaded:", nrow(RNAseq_TPM), "genes,", length(sample_cols), "B/S samples\n")
cat("Cluster data loaded:", nrow(cluster_ICD_data), "patient-level assignments\n")

# ==============================================================================
# 6) PREFLIGHT SAMPLE MATCHING
# ==============================================================================

sample_preflight <- data.frame(
  Sample = sample_cols,
  patient_id = get_neo_patient_id(sample_cols),
  Condition = get_condition(sample_cols),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    cluster_ICD_data %>%
      dplyr::select(patient_id, cluster_sample_id = sample_id, cluster_assignments),
    by = "patient_id"
  )

write.csv(
  sample_preflight,
  file.path(DIRS$logs, "NEO_PAM50_sample_preflight.csv"),
  row.names = FALSE
)

n_unmatched <- sum(is.na(sample_preflight$cluster_assignments))
n_bad_condition <- sum(is.na(sample_preflight$Condition))

cat("\n=== Preflight sample matching ===\n")
cat("Expression samples:", nrow(sample_preflight), "\n")
cat("Unique expression patients:", dplyr::n_distinct(sample_preflight$patient_id), "\n")
cat("Baseline samples:", sum(sample_preflight$Condition == "Baseline", na.rm = TRUE), "\n")
cat("Surgery samples:", sum(sample_preflight$Condition == "Surgery", na.rm = TRUE), "\n")
cat("Samples matched to cluster metadata:", nrow(sample_preflight) - n_unmatched,
    "/", nrow(sample_preflight), "\n")

if (n_bad_condition > 0) {
  stop(
    n_bad_condition,
    " sample(s) could not be assigned to Baseline/Surgery. ",
    "See logs/NEO_PAM50_sample_preflight.csv."
  )
}

if (n_unmatched > 0) {
  cat("\nUnmatched sample(s):\n")
  print(
    sample_preflight %>%
      dplyr::filter(is.na(cluster_assignments)) %>%
      dplyr::select(Sample, patient_id, Condition)
  )

  stop(
    "Expression and cluster inputs are not fully matched. ",
    "See logs/NEO_PAM50_sample_preflight.csv."
  )
}

# Detect unexpected patient-level pairing problems.
pair_check <- sample_preflight %>%
  dplyr::count(patient_id, Condition, name = "n") %>%
  tidyr::pivot_wider(
    names_from = Condition,
    values_from = n,
    values_fill = 0
  )

if (!all(c("Baseline", "Surgery") %in% names(pair_check))) {
  stop("Preflight failed: both Baseline and Surgery samples are required.")
}

bad_pairs <- pair_check %>%
  dplyr::filter(Baseline != 1 | Surgery != 1)

if (nrow(bad_pairs) > 0) {
  cat("\nWARNING: Some patients do not have exactly one Baseline and one Surgery sample:\n")
  print(bad_pairs)
}

cat("Preflight matching complete.\n\n")

cat("Cluster patient IDs extracted:\n")
print(
  cluster_ICD_data[, c("sample_id", "patient_id", "cluster_assignments")]
)

# ==============================================================================
# 7) PREPARE EXPRESSION MATRIX
# ==============================================================================

RNAseq_TPM <- RNAseq_TPM[, c(gene_col, sample_cols), drop = FALSE]
names(RNAseq_TPM)[1] <- "Gene"

RNAseq_TPM$Gene <- as.character(RNAseq_TPM$Gene)

if (anyNA(RNAseq_TPM$Gene) || any(RNAseq_TPM$Gene == "")) {
  stop("The gene-symbol column contains missing or empty values.")
}

# Retain the original behavior for duplicated gene symbols.
RNAseq_TPM$Gene <- make.unique(RNAseq_TPM$Gene)
rownames(RNAseq_TPM) <- RNAseq_TPM$Gene

cat(
  "After processing:",
  nrow(RNAseq_TPM),
  "genes (duplicate symbols, if present, made unique)\n"
)

tpm_mat <- as.matrix(RNAseq_TPM[, sample_cols, drop = FALSE])
storage.mode(tpm_mat) <- "numeric"

if (anyNA(tpm_mat)) {
  stop("Expression matrix contains NA values after numeric conversion.")
}

if (any(tpm_mat < 0, na.rm = TRUE)) {
  warning(
    "Negative TPM values detected. TPM input is expected to be non-negative.",
    call. = FALSE
  )
}

# Retain original transformation for PAM50.
tpm_log <- log2(tpm_mat + 1)

expr_B <- tpm_log[, grepl("_B_Rep[0-9]+$", colnames(tpm_log)), drop = FALSE]
expr_S <- tpm_log[, grepl("_S_Rep[0-9]+$", colnames(tpm_log)), drop = FALSE]

cat(
  "Baseline samples:", ncol(expr_B),
  "| Surgery samples:", ncol(expr_S), "\n"
)

# ==============================================================================
# 8) PAM50 CLASSIFICATION HELPER
# ==============================================================================

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
      patient_id = get_neo_patient_id(Sample),
      Condition  = ifelse(grepl('_B_Rep[0-9]+$', Sample), 'Baseline', 'Surgery')
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
  write.csv(out, file.path(DIRS$results, paste0("NEO_TPM_PAM50_", label, "_samples.csv")), row.names = FALSE)
  write.csv(summary_by_cluster, file.path(DIRS$results, paste0("NEO_TPM_PAM50_", label, "_by_cluster.csv")), row.names = FALSE)
  
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
  write.csv(summary_combined, file.path(DIRS$results, "NEO_TPM_PAM50_by_cluster_B_vs_S.csv"), row.names = FALSE)
  
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
  cat("1. results/NEO_TPM_PAM50_Baseline_samples.csv\n")
  cat("2. results/NEO_TPM_PAM50_Surgery_samples.csv\n")
  cat("3. results/NEO_TPM_PAM50_Baseline_by_cluster.csv\n")
  cat("4. results/NEO_TPM_PAM50_Surgery_by_cluster.csv\n")
  cat("5. results/NEO_TPM_PAM50_by_cluster_B_vs_S.csv\n")
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
      fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE, B = 10000)
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
    ggsave(file.path(DIRS$plots, "NEO_PAM50_Baseline_by_Cluster.png"),
           plot = baseline_result$plot, width = 5.5, height = 6, dpi = 300)
    cat("\n=== Baseline Sample Summary ===\n")
    print(baseline_result$data)
  }
  
  surgery_result <- create_pam50_plot_neo(pam50_S, "Surgery")
  if (!is.null(surgery_result)) {
    print(surgery_result$plot)
    ggsave(file.path(DIRS$plots, "NEO_PAM50_Surgery_by_Cluster.png"),
           plot = surgery_result$plot, width = 5.5, height = 6, dpi = 300)
    cat("\n=== Surgery Sample Summary ===\n")
    print(surgery_result$data)
  }
  
} else {
  cat("Error: PAM50 results not available\n")
  cat("pam50_B is NULL:", is.null(pam50_B), "\n")
  cat("pam50_S is NULL:", is.null(pam50_S), "\n")
}

# ==============================================================================
# 10) REPRODUCIBILITY RECORD
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(DIRS$logs, "NEO_PAM50_sessionInfo.txt")
)

cat("\n=== Reproducibility files ===\n")
cat("logs/NEO_PAM50_input_manifest.csv\n")
cat("logs/NEO_PAM50_sample_preflight.csv\n")
cat("logs/NEO_PAM50_package_versions.csv\n")
cat("logs/NEO_PAM50_sessionInfo.txt\n")
cat("Run log:", LOG_FILE, "\n")
cat("\n=== NEO PAM50 pipeline completed ===\n")

sink()
