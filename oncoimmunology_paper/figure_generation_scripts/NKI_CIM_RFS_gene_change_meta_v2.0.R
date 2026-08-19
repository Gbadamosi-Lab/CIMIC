# ==============================================================================
# NKI CIM RFS Master Analysis
# ==============================================================================
#
# Standalone R conversion of NKI_CIM_RFS_gene_change_and_ssGSEA_master_v1.1.Rmd
#
# The YAML header became the `params` list below. Run from the repository
# root so the repo-relative paths in params resolve:
#
#   setwd("<repo root>"); source("<this file>", echo = TRUE)
#   or:  Rscript <this file>
#
# Differences from the Rmd, confined to setup:
#   * no knitr::opts_chunk$set() and no knitr::kable() rendering
#   * no rstudioapi YAML reload -- `params` is a plain list, always current
#   * figures and CSVs are written to disk exactly as before; no HTML report
# ==============================================================================

# ------------------------------------------------------------------------------
# ANALYSIS PARAMETERS  (was the Rmd YAML header)
# ------------------------------------------------------------------------------
params <- list(
  # Written under Results/ alongside the other analysis outputs.
  output_dir = "oncoimmunology_paper/Results/NKI_CIM_RFS_meta_NKI_CL_kAllC_rho030",

  # ---------------------------------------------------------------------------
  # INPUTS -- repo-relative, resolved against data_dir (the repository root).
  # v1.1 hunted for a flat data/ folder containing all three files; in this
  # repository they live in three different places, so each carries its subpath.
  # ---------------------------------------------------------------------------
  data_dir = ".",
  survival_file =
    "oncoimmunology_paper/Datasets/NKI_SMC/Clinical_data/NKI_survival_data.xlsx",
  expression_file =
    "oncoimmunology_paper/Datasets/NKI_SMC/nki_smc_combine_TPM_lognorm.csv",
  # v2.0: the CIM gene sets now come from the cross-dataset meta-analysis
  # (cimic_overlap_meta_rankbiserial.R) instead of the single-cohort
  # nki_smc_correlation_analysis.csv thresholded at rho >= 0.60.
  signature_file =
    "oncoimmunology_paper/Results/cimic_meta_analysis/meta_results.csv",

  # ---------------------------------------------------------------------------
  # GENE-SET SELECTION FROM THE META-ANALYSIS
  # ---------------------------------------------------------------------------
  # Same toggle surface as CIM_ssGSEA_survival_external_cohorts_v3.0_meta.R, so
  # the discovery-cohort RFS analysis and the external-cohort survival analysis
  # can be run on exactly the same gene sets.
  #
  # meta_scope       ALL3_NKI.NEO.TNBCCL | NKI.NEO | NKI.TNBCCL | NEO.TNBCCL
  # meta_k_mode      k_min2 | k_all | k_all_concordant
  #                  k_datasets is MAGNITUDE ONLY (|r_rb| >= rho, either
  #                  direction). Only k_all_concordant also requires the
  #                  datasets to agree on direction.
  # meta_rank_metric stouffer_padj | fisher_padj | mean_abs_rrb
  meta_scope = "NKI.TNBCCL",
  meta_rho_threshold = 0.30,
  meta_k_mode = "k_all_concordant",
  meta_require_fisher = TRUE,
  meta_require_stouffer = TRUE,
  meta_padj_threshold = 0.05,
  meta_min_mean_abs_rrb = NULL,   # e.g. 0.5 for an extra effect-size floor
  meta_top_n_genes = NULL,        # e.g. 200 for the top N per signature
  meta_rank_metric = "stouffer_padj",
  # A gene can be meta-significant toward BOTH CIMs. "drop" excludes it from
  # both signatures; "assign_best" gives it to the smaller ranking metric.
  meta_ambiguous_gene_policy = "drop",

  # Gene-level RFS analysis.
  cox_fdr_threshold = 0.05,
  top_n_forest = 15,
  top_n_heatmap = 50,
  run_firth_sensitivity = TRUE,

  # ssGSEA analysis.
  run_ssgsea = TRUE,
  alpha = 0.25,
  normalize_scores = TRUE,
  minimum_genes = 5,
  minimum_coverage = 0.50,
  make_ssgsea_figure = TRUE
)

# ------------------------------------------------------------------------------
# SETUP
# ------------------------------------------------------------------------------
options(stringsAsFactors = FALSE, scipen = 999)
analysis_parameter_source <- "parameters defined in the .R script params list"
analysis_rmd_path <- ""   # no Rmd; the configuration section sets paths explicitly
# knitr::kable() only pretty-prints in an HTML report. In a script, print().
kable <- function(x, ...) print(as.data.frame(x))

# [setup chunk replaced above -- knitr/YAML-reload logic does not apply]


# ==============================================================================
# PURPOSE AND STATISTICAL SCOPE
# ==============================================================================

# This R Markdown file is the **RFS-focused master workflow** for the NKI paired
# baseline-to-surgery CIM analyses. It deliberately combines the two original
# relapse-focused analyses into one reproducible GitHub-ready entry point:

# 1. **Gene-level paired expression change:** each selected Fun-CIM or Dys-CIM
#    gene is evaluated separately using the standardized baseline-to-surgery
#    change in expression.
# 2. **Signature-level ssGSEA change:** Fun-CIM and Dys-CIM ssGSEA scores are
#    calculated independently for baseline and surgery samples, followed by
#    paired subtraction and RFS Cox modeling.

# Both modules use the **same selected CIM gene sets**, the same paired NKI
# patients, the same RFS definition, and the same input files. This prevents
# gene-selection drift between scripts and creates one audit trail for a
# repository release.

# For patient \(i\) and gene \(g\), gene-expression change is

# $$
# \Delta_{ig} = \log_2(\mathrm{TPM}_{S,ig}+1) -
# \log_2(\mathrm{TPM}_{B,ig}+1).
# $$

# For each CIM signature, ssGSEA is first computed from the original
# transcriptome-wide sample expression profiles and the paired score change is

# $$
# \Delta\mathrm{ssGSEA}_i =
# \mathrm{ssGSEA}_{S,i} - \mathrm{ssGSEA}_{B,i}.
# $$

# The survival endpoint in this master workflow is **relapse-free survival
# (RFS)**. RFS time is the recorded RFS time for patients with recurrence and
# follow-up time for patients without recurrence. Gene-level and ssGSEA Cox
# models are univariable exploratory analyses. Predictors are standardized so
# hazard ratios represent a one-standard-deviation greater paired change.

# Because the NKI paired cohort has relatively few events, inference should be
# treated as exploratory and hypothesis-generating. Firth-penalized Cox models
# are retained only as a sensitivity analysis for the gene-level module.


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

# The workflow checks its dependencies at startup. Missing CRAN packages are
# installed from the configured CRAN repository. Package installation therefore
# requires internet access and permission to write to the active R library.

# ------------------------------------------------------------------------------
# SECTION: packages
# ------------------------------------------------------------------------------

if (identical(getOption("repos")[["CRAN"]], "@CRAN@") ||
    is.null(getOption("repos")[["CRAN"]])) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

cran_packages <- c(
  "tidyverse", "readxl", "janitor", "survival", "broom", "pheatmap",
  "gridExtra", "RColorBrewer", "scales", "knitr", "patchwork",
  "BiocManager", "rmarkdown", "rstudioapi"
)

if (isTRUE(params$run_firth_sensitivity)) {
  cran_packages <- c(cran_packages, "coxphf")
}

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0L) {
  install.packages(missing_cran, dependencies = TRUE)
}

if (!requireNamespace("GSVA", quietly = TRUE)) {
  BiocManager::install("GSVA", ask = FALSE, update = FALSE)
}

required_packages <- unique(c(cran_packages, "GSVA"))
still_missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(still_missing) > 0L) {
  stop(
    "These required packages remain unavailable: ",
    paste(still_missing, collapse = ", "),
    ". Check internet access and write permission for .libPaths()."
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(survival)
  library(broom)
  library(pheatmap)
  library(gridExtra)
  library(RColorBrewer)
  library(scales)
  library(knitr)
  library(patchwork)
})

# Firth-penalized Cox sensitivity models use the CRAN package `coxphf` and run
# only when `params$run_firth_sensitivity` is set to `true`. Signature-level
# ssGSEA is maintained in the dedicated `NKI_CIM_ssGSEA_analysis.Rmd` workflow.


# ------------------------------------------------------------------------------
# GENE-SELECTION PARAMETER EXAMPLES
# ------------------------------------------------------------------------------

# The correlation and FDR filters are changed only in the YAML header. Their
# values are not repeated inside the filtering code. With `correlation_operator:
# ">="`, the selected Fun-CIM rule is `r_rb >= correlation_threshold` and the
# selected Dys-CIM rule is `r_rb <= -correlation_threshold`. To include genes
# equal to the chosen FDR boundary, change `fdr_operator` from `"<"` to `"<="`.
# Any finite cutoff from 0 to 1 is accepted and applied at run time. The independent
# `cox_fdr_threshold` parameter controls only the downstream Cox-results flag and
# does not determine which genes enter the analysis.


# ==============================================================================
# 2. CONFIGURATION AND OUTPUT DIRECTORIES
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION: configuration
# ------------------------------------------------------------------------------

# v2.0: explicit repo-relative paths. v1.1 searched for a flat data/ folder
# holding all three inputs; in this repository the survival table, the
# expression matrix and the meta-analysis results live in three different
# places, so each params entry carries its own subpath and data_dir is just the
# root they are resolved against.
is_absolute_path <- function(x) grepl("^[A-Za-z]:|^/|^\\\\\\\\", x)
data_dir <- normalizePath(as.character(params$data_dir), winslash = "/",
                          mustWork = FALSE)
if (!dir.exists(data_dir)) {
  stop("data_dir does not exist: ", data_dir,
       "\n  working directory = ", getwd(),
       "\n  data_dir is relative, so run this from the repository root (the ",
       "folder containing oncoimmunology_paper/).")
}
input_path <- function(x) if (is_absolute_path(x)) x else file.path(data_dir, x)

input_files <- c(
  survival   = input_path(as.character(params$survival_file)),
  expression = input_path(as.character(params$expression_file)),
  signature  = input_path(as.character(params$signature_file))
)
missing_inputs <- names(input_files)[!file.exists(input_files)]
if (length(missing_inputs)) {
  stop(
    "Missing input file(s):\n",
    paste0("  ", missing_inputs, " = ", input_files[missing_inputs],
           collapse = "\n"),
    if ("signature" %in% missing_inputs)
      paste0("\n  `signature` must be meta_results.csv from ",
             "cimic_overlap_meta_rankbiserial.R; run that script first.")
    else "",
    "\n  Working directory = ", getwd(),
    " (run from the repository root)."
  )
}

output_dir <- as.character(params$output_dir)
qc_dir <- file.path(output_dir, "01_QC")
change_dir <- file.path(output_dir, "02_gene_expression_change")
gene_result_dir <- file.path(output_dir, "03_gene_RFS_results")
gene_figure_dir <- file.path(output_dir, "04_gene_RFS_figures")
ssgsea_score_dir <- file.path(output_dir, "05_ssGSEA_scores")
ssgsea_result_dir <- file.path(output_dir, "06_ssGSEA_RFS_results")
ssgsea_figure_dir <- file.path(output_dir, "07_ssGSEA_RFS_figures")
sensitivity_dir <- file.path(output_dir, "08_gene_RFS_sensitivity")
reproducibility_dir <- file.path(output_dir, "09_reproducibility")

purrr::walk(
  c(
    output_dir, qc_dir, change_dir, gene_result_dir, gene_figure_dir,
    ssgsea_score_dir, ssgsea_result_dir, ssgsea_figure_dir,
    sensitivity_dir, reproducibility_dir
  ),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
)

# Aliases retained so the original gene-level plotting/export functions remain
# unchanged while the master output tree is explicit.
result_dir <- gene_result_dir
figure_dir <- gene_figure_dir

fun_color <- "#159D9C"
dys_color <- "#F07B72"
increase_color <- "#B85C38"
decrease_color <- "#3B6FB6"

input_provenance <- purrr::imap_dfr(
  input_files,
  function(path, role) {
    info <- file.info(path)
    tibble::tibble(
      input_role = role,
      configured_filename = basename(path),
      resolved_path = normalizePath(path, winslash = "/", mustWork = TRUE),
      file_size_bytes = as.numeric(info$size),
      modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
      md5 = unname(tools::md5sum(path))
    )
  }
)
readr::write_csv(input_provenance, file.path(qc_dir, "input_file_provenance.csv"))


# ==============================================================================
# 3. IMPORT AND CLEAN THE GENE-SIGNATURE TABLE
# ==============================================================================

# The correlation file is expected to contain `gene_id`, `r_rb`, `p_value`,
# `padj`, and `trajectory`. Gene selection is controlled in the YAML header by
# four independent parameters: `correlation_threshold`, `correlation_operator`,
# `fdr_threshold`, and `fdr_operator`. The values are parsed directly from the
# YAML header, validated, and exported to `01_QC/analysis_parameters.csv`.
# The input `trajectory` annotation is retained for auditing but is not used as a
# selection gate because it may have been generated upstream using a different
# correlation cutoff. Instead, direction is derived from the correlation sign:
# positive `r_rb` is Fun-CIM and negative `r_rb` is Dys-CIM. The user-selected
# threshold and operator are then applied directionally. Correlation statistics
# are retained with explicit prefixes so they cannot be confused with the new
# Cox-model p-values.

# ------------------------------------------------------------------------------
# SECTION: import-signature
# ------------------------------------------------------------------------------
# Defensive package loading allows this and subsequent chunks to run after the
# configuration chunk even when the user did not manually execute the earlier
# package chunk. Knitting the complete RMD still remains the recommended route.
runtime_packages <- c(
  "tidyverse", "readxl", "janitor", "survival", "broom", "pheatmap",
  "gridExtra", "RColorBrewer", "scales", "knitr"
)
missing_runtime_packages <- runtime_packages[
  !vapply(runtime_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_runtime_packages)) {
  stop(
    "Required package(s) are not installed: ",
    paste(missing_runtime_packages, collapse = ", "),
    ". Run the 'packages' chunk first or knit the RMD from the beginning."
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(survival)
  library(broom)
  library(pheatmap)
  library(gridExtra)
  library(RColorBrewer)
  library(scales)
  library(knitr)
})

signature_raw <- readr::read_csv(
  input_files[["signature"]],
  show_col_types = FALSE,
  name_repair = "minimal"
) %>%
  janitor::clean_names()

signature_file_info <- file.info(input_files[["signature"]])
signature_input_provenance <- tibble::tibble(
  input_role = "correlation signature",
  configured_filename = as.character(params$signature_file),
  resolved_path = normalizePath(
    input_files[["signature"]], winslash = "/", mustWork = TRUE
  ),
  file_size_bytes = as.numeric(signature_file_info$size),
  modified_time = format(signature_file_info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
  md5 = unname(tools::md5sum(input_files[["signature"]])),
  imported_rows = nrow(signature_raw)
)
readr::write_csv(
  signature_input_provenance,
  file.path(qc_dir, "signature_input_provenance.csv")
)

# =============================================================================
# v2.0 GENE SELECTION FROM THE CIMIC META-ANALYSIS
# -----------------------------------------------------------------------------
# v1.1 read nki_smc_correlation_analysis.csv and kept genes with
#   r_rb >= +0.60 (Fun-CIM) or r_rb <= -0.60 (Dys-CIM), padj < 0.05,
# i.e. a single-cohort rule. v2.0 instead reads meta_results.csv from
# cimic_overlap_meta_rankbiserial.R, so the discovery-cohort RFS analysis uses
# the same cross-dataset consensus signatures as the external-cohort survival
# analysis.
#
# Everything downstream is unchanged. This block reproduces the exact contract
# the rest of the script depends on:
#   signature                 gene, direction ("Fun-CIM"/"Dys-CIM"), r_rb,
#                             correlation_p_value, correlation_padj
#   signature_filter_audit    every candidate row with its pass/fail flags
#   selection_verification    per-signature re-check of the active rules
#   signature_summary         stage counts
#   correlation_stage_label / correlation_fdr_stage_label  (stage labels reused
#                             by the gene-availability summary further down)
#
# COLUMN MAPPING. The meta table has no per-cohort r_rb, so the v1.1 column
# names are filled with their meta analogues, keeping the same sign convention
# (positive = Fun-CIM, negative = Dys-CIM) that downstream abs()/sign() code
# assumes:
#   r_rb                <- mean_abs_rrb, signed by the CIM the gene belongs to
#   correlation_p_value <- stouffer_p    (or fisher_p if only Fisher required)
#   correlation_padj    <- stouffer_padj (or fisher_padj likewise)
# =============================================================================

meta_scope <- as.character(params$meta_scope)
meta_rho_threshold <- as.numeric(params$meta_rho_threshold)
meta_k_mode <- as.character(params$meta_k_mode)
meta_require_fisher <- isTRUE(params$meta_require_fisher)
meta_require_stouffer <- isTRUE(params$meta_require_stouffer)
meta_padj_threshold <- as.numeric(params$meta_padj_threshold)
meta_rank_metric <- as.character(params$meta_rank_metric)
meta_ambiguous_gene_policy <- as.character(params$meta_ambiguous_gene_policy)
meta_top_n_genes <- if (is.null(params$meta_top_n_genes)) NULL else
  as.integer(params$meta_top_n_genes)
meta_effect_floor <- if (is.null(params$meta_min_mean_abs_rrb)) -Inf else
  as.numeric(params$meta_min_mean_abs_rrb)

# Unrelated to gene selection, but validated in the v1.1 block this replaces,
# and used by the gene-level Cox section further down.
cox_fdr_threshold <- as.numeric(params$cox_fdr_threshold)
if (length(cox_fdr_threshold) != 1L || !is.finite(cox_fdr_threshold) ||
    cox_fdr_threshold <= 0 || cox_fdr_threshold > 1)
  stop("params$cox_fdr_threshold must be one number in (0, 1].")

if (!meta_k_mode %in% c("k_min2", "k_all", "k_all_concordant"))
  stop("params$meta_k_mode must be k_min2, k_all or k_all_concordant.")
if (!meta_rank_metric %in% c("stouffer_padj", "fisher_padj", "mean_abs_rrb"))
  stop("params$meta_rank_metric must be stouffer_padj, fisher_padj or mean_abs_rrb.")
if (!meta_ambiguous_gene_policy %in% c("drop", "assign_best"))
  stop("params$meta_ambiguous_gene_policy must be drop or assign_best.")
if (!meta_require_fisher && !meta_require_stouffer)
  stop("At least one of meta_require_fisher / meta_require_stouffer must be TRUE.")

required_signature_columns <- c(
  "scope", "cim", "rho_threshold", "gene", "k_datasets", "n_scope_datasets",
  "direction_concordant", "fisher_padj", "stouffer_padj", "mean_abs_rrb"
)
missing_signature_columns <- setdiff(required_signature_columns,
                                     names(signature_raw))
if (length(missing_signature_columns) > 0) {
  stop(
    "meta_results.csv is missing required column(s): ",
    paste(missing_signature_columns, collapse = ", "),
    "\n  Re-run cimic_overlap_meta_rankbiserial.R."
  )
}

if (!meta_scope %in% unique(signature_raw$scope))
  stop("meta_scope '", meta_scope, "' is not in meta_results.csv. Available: ",
       paste(unique(signature_raw$scope), collapse = ", "))
available_rho <- sort(unique(signature_raw$rho_threshold))
if (!any(abs(available_rho - meta_rho_threshold) < 1e-8))
  stop("meta_rho_threshold ", meta_rho_threshold,
       " is not in meta_results.csv. Available: ",
       paste(available_rho, collapse = ", "))

# Stage labels, reused by the gene-availability summary further down.
correlation_fdr_stage_label <- sprintf(
  "Meta significance (%s < %s)",
  paste(c(if (meta_require_fisher) "fisher_padj",
          if (meta_require_stouffer) "stouffer_padj"), collapse = " and "),
  format(meta_padj_threshold, scientific = FALSE, trim = TRUE))
correlation_stage_label <- sprintf(
  "Dataset support (%s at rho >= %s, scope %s)",
  meta_k_mode, format(meta_rho_threshold, scientific = FALSE, trim = TRUE),
  meta_scope)

# Full audit: every meta row for this scope+rho with each rule evaluated
# separately, so QC shows exactly which filter removed what.
signature_filter_audit <- signature_raw %>%
  dplyr::filter(scope == meta_scope,
                abs(rho_threshold - meta_rho_threshold) < 1e-8) %>%
  dplyr::transmute(
    gene = stringr::str_trim(as.character(gene)),
    direction = dplyr::case_when(cim == "Fun" ~ "Fun-CIM",
                                 cim == "Dys" ~ "Dys-CIM",
                                 TRUE ~ "Unassigned"),
    meta_scope_used = scope,
    rho_threshold = rho_threshold,
    k_datasets = k_datasets,
    n_scope_datasets = n_scope_datasets,
    k_concordant = if ("k_concordant" %in% names(signature_raw))
                     k_concordant else NA_integer_,
    direction_concordant = direction_concordant,
    fisher_padj = suppressWarnings(as.numeric(fisher_padj)),
    stouffer_padj = suppressWarnings(as.numeric(stouffer_padj)),
    mean_abs_rrb = suppressWarnings(as.numeric(mean_abs_rrb)),
    # Signed so the sign convention matches v1.1: + = Fun-CIM, - = Dys-CIM.
    r_rb = ifelse(direction == "Dys-CIM", -mean_abs_rrb, mean_abs_rrb),
    correlation_p_value = if (meta_require_stouffer)
      suppressWarnings(as.numeric(stouffer_p)) else
      suppressWarnings(as.numeric(fisher_p)),
    correlation_padj = if (meta_require_stouffer) stouffer_padj else fisher_padj,
    directional_correlation_strength = mean_abs_rrb,
    passes_fisher = !meta_require_fisher |
      (is.finite(fisher_padj) & fisher_padj < meta_padj_threshold),
    passes_stouffer = !meta_require_stouffer |
      (is.finite(stouffer_padj) & stouffer_padj < meta_padj_threshold),
    passes_k_mode = switch(
      meta_k_mode,
      k_min2           = k_datasets >= 2,
      k_all            = k_datasets == n_scope_datasets,
      k_all_concordant = k_datasets == n_scope_datasets &
                           direction_concordant == "YES"),
    passes_effect_floor = is.finite(mean_abs_rrb) &
      mean_abs_rrb >= meta_effect_floor,
    passes_correlation_threshold = passes_k_mode & passes_effect_floor,
    selected_for_analysis = !is.na(gene) & gene != "" &
      direction %in% c("Fun-CIM", "Dys-CIM") &
      passes_fisher & passes_stouffer & passes_k_mode & passes_effect_floor
  )

readr::write_csv(
  signature_filter_audit,
  file.path(qc_dir, "signature_gene_filtering_audit.csv")
)

selected_signature_rows <- signature_filter_audit %>%
  dplyr::filter(selected_for_analysis)

# A gene can be meta-significant toward BOTH CIMs (impossible with a single
# signed r_rb, but the meta table tests Fun and Dys separately). Resolve it
# explicitly rather than letting the gene land in both signatures.
ambiguous_genes <- selected_signature_rows %>%
  dplyr::distinct(gene, direction) %>%
  dplyr::count(gene, name = "n_directions") %>%
  dplyr::filter(n_directions > 1) %>%
  dplyr::pull(gene)
if (length(ambiguous_genes)) {
  message("Genes meta-significant toward BOTH CIMs: ", length(ambiguous_genes),
          "  (policy: ", meta_ambiguous_gene_policy, ")")
  readr::write_csv(
    selected_signature_rows %>%
      dplyr::filter(gene %in% ambiguous_genes) %>% dplyr::arrange(gene, direction),
    file.path(qc_dir, "signature_ambiguous_genes.csv"))
  if (identical(meta_ambiguous_gene_policy, "drop"))
    selected_signature_rows <- selected_signature_rows %>%
      dplyr::filter(!gene %in% ambiguous_genes)
}

signature <- selected_signature_rows %>%
  dplyr::arrange(gene, .data[[meta_rank_metric]],
                 dplyr::desc(mean_abs_rrb)) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

if (!is.null(meta_top_n_genes)) {
  signature <- signature %>%
    dplyr::group_by(direction) %>%
    dplyr::arrange(.data[[meta_rank_metric]], dplyr::desc(mean_abs_rrb),
                   .by_group = TRUE) %>%
    dplyr::slice_head(n = meta_top_n_genes) %>%
    dplyr::ungroup()
  message("Applied meta_top_n_genes = ", meta_top_n_genes, " per signature.")
}

if (!all(c("Fun-CIM", "Dys-CIM") %in% unique(signature$direction)))
  stop("Both Fun-CIM and Dys-CIM must be non-empty after selection. Active: ",
       "scope ", meta_scope, " | rho ", meta_rho_threshold, " | ", meta_k_mode,
       ". Loosen meta_k_mode or lower meta_rho_threshold.")
if (anyDuplicated(signature$gene))
  stop("A gene appears in both signatures; see signature_ambiguous_genes.csv.")

# Re-check every active rule against the FINAL selection, independently of the
# flags computed above. Any non-zero count means a filter did not apply.
selection_verification <- signature %>%
  dplyr::group_by(direction) %>%
  dplyr::summarise(
    selected_genes = dplyr::n(),
    minimum_selected_mean_abs_rrb = min(mean_abs_rrb),
    maximum_selected_fisher_padj = max(fisher_padj),
    maximum_selected_stouffer_padj = max(stouffer_padj),
    minimum_k_datasets = min(k_datasets),
    genes_failing_fisher_rule = sum(
      meta_require_fisher & !(fisher_padj < meta_padj_threshold)),
    genes_failing_stouffer_rule = sum(
      meta_require_stouffer & !(stouffer_padj < meta_padj_threshold)),
    genes_failing_k_mode_rule = sum(switch(
      meta_k_mode,
      k_min2           = !(k_datasets >= 2),
      k_all            = !(k_datasets == n_scope_datasets),
      k_all_concordant = !(k_datasets == n_scope_datasets &
                             direction_concordant == "YES"))),
    genes_failing_effect_floor_rule = sum(!(mean_abs_rrb >= meta_effect_floor)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(genes_failing_runtime_rule =
    genes_failing_fisher_rule + genes_failing_stouffer_rule +
    genes_failing_k_mode_rule + genes_failing_effect_floor_rule)

readr::write_csv(selection_verification,
                 file.path(qc_dir, "signature_selection_verification.csv"))
if (any(selection_verification$genes_failing_runtime_rule > 0L))
  stop("Internal gene-selection verification failed; inspect 01_QC outputs.")

readr::write_csv(signature, file.path(qc_dir, "selected_signature_genes.csv"))

# Stage counts, in the same shape v1.1 produced.
signature_summary <- dplyr::bind_rows(
  signature_filter_audit %>%
    dplyr::filter(direction %in% c("Fun-CIM", "Dys-CIM")) %>%
    dplyr::count(direction, name = "n_genes") %>%
    dplyr::mutate(stage = "Direction assigned from r_rb sign"),
  signature_filter_audit %>%
    dplyr::filter(direction %in% c("Fun-CIM", "Dys-CIM"),
                  passes_fisher, passes_stouffer) %>%
    dplyr::count(direction, name = "n_genes") %>%
    dplyr::mutate(stage = correlation_fdr_stage_label),
  signature %>%
    dplyr::count(direction, name = "n_genes") %>%
    dplyr::mutate(stage = correlation_stage_label)
) %>%
  dplyr::select(stage, direction, n_genes) %>%
  dplyr::arrange(factor(stage, levels = c(
    "Direction assigned from r_rb sign", correlation_fdr_stage_label,
    correlation_stage_label)), direction)

readr::write_csv(signature_summary,
                 file.path(qc_dir, "signature_filtering_stage_counts.csv"))

message("CIM gene sets from meta_results.csv | scope ", meta_scope,
        " | rho ", meta_rho_threshold, " | ", meta_k_mode)
kable(signature_summary, caption = "Genes supplied in each CIM signature")

# The selected genes above are used identically by the gene-change and ssGSEA
# modules. The following object is the single source of truth for both signatures.

# ------------------------------------------------------------------------------
# SECTION: shared-gene-sets
# ------------------------------------------------------------------------------
gene_sets_supplied <- signature %>%
  filter(direction %in% c("Fun-CIM", "Dys-CIM")) %>%
  distinct(direction, gene) %>%
  arrange(direction, gene) %>%
  split(.$direction) %>%
  purrr::map(~ .x$gene)

shared_signature_summary <- tibble::tibble(
  signature = names(gene_sets_supplied),
  n_genes = purrr::map_int(gene_sets_supplied, length)
)

readr::write_csv(
  shared_signature_summary,
  file.path(qc_dir, "shared_signature_gene_counts.csv")
)

knitr::kable(
  shared_signature_summary,
  caption = "Shared CIM gene sets used by both RFS analysis modules"
)


# ==============================================================================
# 4. IMPORT AND CONSTRUCT THE NKI SURVIVAL ENDPOINTS
# ==============================================================================

# For patients with an event, the recorded event time is used. For patients
# without an event, follow-up time is used as the censoring time.

# ------------------------------------------------------------------------------
# SECTION: import-survival
# ------------------------------------------------------------------------------

survival_raw <- readxl::read_excel(input_files[["survival"]]) %>%
  janitor::clean_names()

required_survival_columns <- c(
  "study_id", "fu_months", "recurrence", "rfs_months"
)

missing_survival_columns <- setdiff(
  required_survival_columns,
  names(survival_raw)
)

if (length(missing_survival_columns) > 0L) {
  stop(
    "Survival file is missing required RFS column(s) after clean_names(): ",
    paste(missing_survival_columns, collapse = ", "),
    "\nObserved columns: ", paste(names(survival_raw), collapse = ", ")
  )
}

survival_data <- survival_raw %>%
  transmute(
    patient_id = stringr::str_trim(as.character(study_id)),
    follow_up_months = suppressWarnings(as.numeric(fu_months)),
    recurrence = suppressWarnings(as.integer(recurrence)),
    recorded_rfs_months = suppressWarnings(as.numeric(rfs_months))
  ) %>%
  mutate(
    rfs_time = if_else(
      recurrence == 1L,
      recorded_rfs_months,
      follow_up_months
    ),
    rfs_event = recurrence
  )

if (anyDuplicated(survival_data$patient_id)) {
  stop("Duplicated StudyID values were found in the survival file.")
}

survival_qc <- survival_data %>%
  mutate(
    invalid_recurrence_code = !recurrence %in% c(0L, 1L),
    missing_rfs_time = is.na(rfs_time),
    nonpositive_rfs_time = !is.na(rfs_time) & rfs_time <= 0,
    recurrence_time_missing = recurrence == 1L & is.na(recorded_rfs_months),
    event_after_follow_up = recurrence == 1L &
      !is.na(recorded_rfs_months) &
      !is.na(follow_up_months) &
      recorded_rfs_months > follow_up_months
  )

readr::write_csv(
  survival_qc,
  file.path(qc_dir, "RFS_survival_data_qc.csv")
)

fatal_survival_issues <- survival_qc %>%
  filter(
    invalid_recurrence_code |
      missing_rfs_time |
      nonpositive_rfs_time |
      recurrence_time_missing
  )

if (nrow(fatal_survival_issues) > 0L) {
  stop(
    "Fatal RFS data issues were detected. Review ",
    file.path(qc_dir, "RFS_survival_data_qc.csv")
  )
}


# ==============================================================================
# 5. IMPORT EXPRESSION AND PARSE SAMPLE NAMES
# ==============================================================================

# The combined matrix contains NKI and SMC samples. NKI identifiers are numeric
# (for example, `2046_B`), whereas SMC identifiers begin with `OB_` (for example,
# `OB_0139_B`). Only NKI patients with survival information are retained.

# The parser removes only the final `_B` or `_S`, preserving underscores within
# patient identifiers.

# ------------------------------------------------------------------------------
# SECTION: import-expression
# ------------------------------------------------------------------------------
expression_raw <- readr::read_csv(
  input_files[["expression"]],
  show_col_types = FALSE,
  name_repair = "minimal"
)

if (ncol(expression_raw) < 3) {
  stop("Expression matrix must contain one gene column and multiple sample columns.")
}

# Detect the gene identifier column by position. Position-based renaming is
# intentional because expression CSVs exported with row names commonly have a
# blank first-column header; tidyselect::all_of("") cannot select that column.
expression_column_names <- names(expression_raw)
normalized_expression_names <- str_to_lower(str_trim(expression_column_names))
candidate_gene_indices <- which(
  normalized_expression_names %in%
    c("gene", "genes", "gene_id", "gene_symbol", "symbol", "hgnc_symbol")
)

gene_column_index <- if (length(candidate_gene_indices) > 0L) {
  candidate_gene_indices[[1]]
} else {
  1L
}
gene_column_original_name <- expression_column_names[[gene_column_index]]
gene_column_detection <- if (length(candidate_gene_indices) > 0L) {
  "recognized gene-column header"
} else if (is.na(gene_column_original_name) ||
           !nzchar(str_trim(gene_column_original_name))) {
  "blank first-column header; used column 1"
} else {
  "no recognized gene-column header; used column 1"
}

# Rename directly by position so an empty original header is handled safely.
names(expression_raw)[gene_column_index] <- "gene"

# dplyr cannot transform a data frame while any column name is NA or blank.
# A trailing delimiter in a CSV commonly creates a completely empty unnamed
# column. Audit all remaining unnamed columns, remove only those that are
# entirely empty or verified sequential row indices, and stop if any other
# unnamed column contains data (it may be a sample whose header was lost).
current_expression_names <- names(expression_raw)
blank_name_flags <- is.na(current_expression_names) |
  trimws(ifelse(is.na(current_expression_names), "", current_expression_names)) == ""
unnamed_column_indices <- which(blank_name_flags)

unnamed_column_is_empty <- if (length(unnamed_column_indices) > 0L) {
  vapply(
    expression_raw[unnamed_column_indices],
    function(column_values) {
      column_text <- as.character(column_values)
      all(is.na(column_text) |
            trimws(ifelse(is.na(column_text), "", column_text)) == "")
    },
    logical(1)
  )
} else {
  logical(0)
}

# Excel/CSV exports may also prepend an unnamed row-number column. Treat it as
# an index only when every value exactly matches 1, 2, ..., nrow(expression_raw).
# This narrow rule prevents a genuine unnamed sample column from being removed.
unnamed_column_is_row_index <- if (length(unnamed_column_indices) > 0L) {
  vapply(
    expression_raw[unnamed_column_indices],
    function(column_values) {
      column_text <- trimws(as.character(column_values))
      column_numeric <- suppressWarnings(as.numeric(column_text))
      length(column_numeric) == nrow(expression_raw) &&
        !anyNA(column_numeric) &&
        identical(column_numeric, as.numeric(seq_len(nrow(expression_raw))))
    },
    logical(1)
  )
} else {
  logical(0)
}

unnamed_expression_column_audit <- tibble::tibble(
  column_index = unnamed_column_indices,
  entirely_empty = unname(unnamed_column_is_empty),
  sequential_row_index = unname(unnamed_column_is_row_index),
  action = case_when(
    entirely_empty ~ "removed: unnamed and entirely empty",
    sequential_row_index ~ "removed: verified sequential row index (1:n)",
    TRUE ~ "analysis stopped: unnamed column contains non-index data"
  )
)
readr::write_csv(
  unnamed_expression_column_audit,
  file.path(qc_dir, "unnamed_expression_columns.csv")
)

invalid_unnamed_indices <- unnamed_column_indices[
  !unnamed_column_is_empty & !unnamed_column_is_row_index
]
if (length(invalid_unnamed_indices) > 0L) {
  stop(
    "The expression matrix contains unnamed column(s) with non-index data at ",
    "position(s): ", paste(invalid_unnamed_indices, collapse = ", "),
    ". Add valid sample header(s) to the source CSV and rerun. Review ",
    file.path(qc_dir, "unnamed_expression_columns.csv"), "."
  )
}

dropped_unnamed_indices <- unnamed_column_indices[
  unnamed_column_is_empty | unnamed_column_is_row_index
]
if (length(dropped_unnamed_indices) > 0L) {
  # Explicit row indices avoid an RStudio 2026.07 IDE bug triggered by x[, j].
  expression_raw <- expression_raw[
    seq_len(nrow(expression_raw)),
    -dropped_unnamed_indices,
    drop = FALSE
  ]
}

if (ncol(expression_raw) < 3) {
  stop(
    "After removing non-data unnamed columns, the expression matrix must contain ",
    "one gene column and at least two sample columns."
  )
}

if (anyDuplicated(names(expression_raw))) {
  stop(
    "Renaming the detected gene column created duplicated column names. ",
    "Original column position: ", gene_column_index,
    "; original header: '", gene_column_original_name, "'."
  )
}

expression_gene_column_audit <- tibble::tibble(
  gene_column_index = gene_column_index,
  original_header = if_else(
    is.na(gene_column_original_name) ||
      !nzchar(str_trim(gene_column_original_name)),
    "<blank>", gene_column_original_name
  ),
  assigned_header = "gene",
  detection_rule = gene_column_detection,
  removed_unnamed_columns = length(dropped_unnamed_indices),
  removed_column_positions = if (length(dropped_unnamed_indices) > 0L) {
    paste(dropped_unnamed_indices, collapse = ";")
  } else {
    NA_character_
  }
)
readr::write_csv(
  expression_gene_column_audit,
  file.path(qc_dir, "expression_gene_column_detection.csv")
)

expression_data <- expression_raw %>%
  mutate(gene = str_trim(as.character(gene)))

if (any(is.na(expression_data$gene) | expression_data$gene == "")) {
  stop(
    "The detected gene column contains missing or blank identifiers. Review ",
    file.path(qc_dir, "expression_gene_column_detection.csv"), "."
  )
}

if (anyDuplicated(expression_data$gene)) {
  duplicated_expression_genes <- expression_data %>%
    count(gene, name = "n_rows") %>%
    filter(n_rows > 1)
  readr::write_csv(
    duplicated_expression_genes,
    file.path(qc_dir, "duplicated_expression_genes.csv")
  )
  stop(
    "Duplicated gene identifiers were found in the expression matrix. ",
    "Review 01_QC/duplicated_expression_genes.csv before continuing."
  )
}

sample_columns <- setdiff(names(expression_data), "gene")

sample_metadata <- tibble(sample_id = sample_columns) %>%
  tidyr::extract(
    sample_id,
    into = c("patient_id", "timepoint_code"),
    regex = "^(.+)_([BS])$",
    remove = FALSE
  ) %>%
  mutate(
    cohort_from_name = if_else(str_starts(patient_id, "OB_"), "smc", "nki"),
    timepoint = recode(timepoint_code, B = "Baseline", S = "Surgery")
  )

unparsed_samples <- sample_metadata %>%
  filter(is.na(patient_id) | is.na(timepoint_code))

readr::write_csv(unparsed_samples, file.path(qc_dir, "unparsed_expression_samples.csv"))
if (nrow(unparsed_samples) > 0) {
  stop(
    "Some expression sample names do not end in _B or _S. Review ",
    file.path(qc_dir, "unparsed_expression_samples.csv")
  )
}

nki_sample_metadata <- sample_metadata %>%
  filter(
    cohort_from_name == "nki",
    patient_id %in% survival_data$patient_id
  )

pair_qc <- nki_sample_metadata %>%
  count(patient_id, timepoint, name = "n_samples") %>%
  tidyr::complete(
    patient_id,
    timepoint = c("Baseline", "Surgery"),
    fill = list(n_samples = 0L)
  ) %>%
  pivot_wider(names_from = timepoint, values_from = n_samples) %>%
  mutate(complete_unique_pair = Baseline == 1L & Surgery == 1L)

readr::write_csv(pair_qc, file.path(qc_dir, "nki_pair_qc.csv"))

invalid_pairs <- pair_qc %>% filter(!complete_unique_pair)
if (nrow(invalid_pairs) > 0) {
  stop(
    "Some NKI patients do not have exactly one baseline and one surgery sample. ",
    "Review ", file.path(qc_dir, "nki_pair_qc.csv")
  )
}

included_patient_ids <- pair_qc %>%
  filter(complete_unique_pair) %>%
  pull(patient_id)

unmatched_survival_patients <- survival_data %>%
  filter(!patient_id %in% included_patient_ids)

readr::write_csv(
  unmatched_survival_patients,
  file.path(qc_dir, "survival_patients_without_paired_expression.csv")
)


# ==============================================================================
# 6. RESTRICT TO SIGNATURE GENES AND CALCULATE PAIRED CHANGES
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION: calculate-change
# ------------------------------------------------------------------------------
genes_found <- intersect(signature$gene, expression_data$gene)
genes_missing <- setdiff(signature$gene, expression_data$gene)

gene_availability <- signature %>%
  mutate(in_expression_matrix = gene %in% genes_found)

readr::write_csv(gene_availability, file.path(qc_dir, "signature_gene_availability.csv"))

signature_analysis_counts <- bind_rows(
  signature_summary,
  gene_availability %>%
    filter(in_expression_matrix) %>%
    count(direction, name = "n_genes") %>%
    mutate(stage = "Found in expression matrix") %>%
    dplyr::select(stage, direction, n_genes)
) %>%
  arrange(
    factor(stage, levels = c(
      "Direction assigned from r_rb sign", correlation_fdr_stage_label,
      correlation_stage_label, "Found in expression matrix"
    )),
    direction
  )
readr::write_csv(signature_analysis_counts,
                 file.path(qc_dir, "signature_analysis_stage_counts.csv"))

if (length(genes_found) == 0) {
  stop("None of the signature genes matched the expression matrix.")
}

nki_expression_long <- expression_data %>%
  filter(gene %in% genes_found) %>%
  dplyr::select(gene, all_of(nki_sample_metadata$sample_id)) %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "expression"
  ) %>%
  left_join(
    nki_sample_metadata %>% dplyr::select(sample_id, patient_id, timepoint),
    by = "sample_id"
  ) %>%
  mutate(expression = suppressWarnings(as.numeric(expression)))

non_numeric_expression <- nki_expression_long %>%
  filter(is.na(expression))

readr::write_csv(
  non_numeric_expression,
  file.path(qc_dir, "missing_or_non_numeric_expression.csv")
)

if (nrow(non_numeric_expression) > 0) {
  stop(
    "Missing or nonnumeric expression values were detected. Review ",
    file.path(qc_dir, "missing_or_non_numeric_expression.csv")
  )
}

gene_changes <- nki_expression_long %>%
  dplyr::select(gene, patient_id, timepoint, expression) %>%
  pivot_wider(names_from = timepoint, values_from = expression) %>%
  mutate(delta_expression = Surgery - Baseline) %>%
  group_by(gene) %>%
  mutate(
    delta_mean = mean(delta_expression, na.rm = TRUE),
    delta_sd = sd(delta_expression, na.rm = TRUE),
    delta_z = if_else(
      is.finite(delta_sd) & delta_sd > 0,
      (delta_expression - delta_mean) / delta_sd,
      NA_real_
    )
  ) %>%
  ungroup() %>%
  left_join(signature, by = "gene") %>%
  left_join(
    survival_data %>% filter(patient_id %in% included_patient_ids),
    by = "patient_id"
  )

zero_variance_genes <- gene_changes %>%
  distinct(gene, delta_sd) %>%
  filter(is.na(delta_sd) | delta_sd == 0)

readr::write_csv(zero_variance_genes, file.path(qc_dir, "zero_variance_genes.csv"))
readr::write_csv(gene_changes, file.path(change_dir, "NKI_gene_changes_long.csv"))

gene_change_matrix <- gene_changes %>%
  dplyr::select(gene, patient_id, delta_expression) %>%
  pivot_wider(names_from = patient_id, values_from = delta_expression)

readr::write_csv(
  gene_change_matrix,
  file.path(change_dir, "NKI_gene_change_matrix_log2TPMplus1.csv")
)


# ==============================================================================
# 7. COHORT SUMMARY AND ANALYTIC CHECKS
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION: cohort-summary
# ------------------------------------------------------------------------------

included_patients <- survival_data %>%
  filter(patient_id %in% included_patient_ids) %>%
  arrange(suppressWarnings(as.numeric(patient_id)))

readr::write_csv(
  included_patients,
  file.path(qc_dir, "included_NKI_patients.csv")
)

cohort_summary <- tibble(
  metric = c(
    "NKI patients with paired expression",
    "NKI expression samples",
    "RFS events",
    "Signature genes selected",
    "Signature genes found in expression matrix",
    "Signature genes missing from expression matrix",
    "Zero-variance gene changes"
  ),
  value = c(
    length(included_patient_ids),
    nrow(nki_sample_metadata),
    sum(included_patients$rfs_event == 1L),
    nrow(signature),
    length(genes_found),
    length(genes_missing),
    nrow(zero_variance_genes)
  )
)

readr::write_csv(
  cohort_summary,
  file.path(qc_dir, "analysis_cohort_summary.csv")
)

knitr::kable(
  cohort_summary,
  caption = "Final paired NKI RFS analysis cohort"
)


# ==============================================================================
# 8. GENE-LEVEL COX PROPORTIONAL-HAZARDS MODELS FOR RFS
# ==============================================================================

# The helper below fits one gene at a time. It records model failures rather than
# stopping the complete analysis. The proportional-hazards test is based on
# scaled Schoenfeld residuals (`cox.zph`).

# ------------------------------------------------------------------------------
# SECTION: cox-functions
# ------------------------------------------------------------------------------

fit_gene_rfs_cox <- function(data) {
  model_data <- data %>%
    transmute(
      time = rfs_time,
      event = rfs_event,
      delta_z = delta_z
    ) %>%
    filter(complete.cases(.))

  n_patients <- nrow(model_data)
  n_events <- sum(model_data$event == 1L)

  empty_result <- tibble(
    endpoint = "RFS",
    n_patients = n_patients,
    n_events = n_events,
    coefficient = NA_real_,
    hazard_ratio = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    p_value = NA_real_,
    ph_p_value = NA_real_,
    model_status = NA_character_,
    warning_text = NA_character_
  )

  if (
    n_patients < 3L ||
      n_events < 2L ||
      !is.finite(sd(model_data$delta_z)) ||
      sd(model_data$delta_z) == 0
  ) {
    return(
      empty_result %>%
        mutate(
          model_status = "not_fitted",
          warning_text =
            "Insufficient patients/events or zero predictor variance"
        )
    )
  }

  captured_warnings <- character(0)

  fit <- tryCatch(
    withCallingHandlers(
      survival::coxph(
        survival::Surv(time, event) ~ delta_z,
        data = model_data,
        ties = "efron",
        x = TRUE,
        y = TRUE
      ),
      warning = function(w) {
        captured_warnings <<- c(
          captured_warnings,
          conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(
      empty_result %>%
        mutate(
          model_status = "error",
          warning_text = conditionMessage(fit)
        )
    )
  }

  model_tidy <- broom::tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE
  )

  ph_p <- tryCatch(
    unname(
      survival::cox.zph(fit)$table["delta_z", "p"]
    ),
    error = function(e) NA_real_
  )

  tibble(
    endpoint = "RFS",
    n_patients = n_patients,
    n_events = n_events,
    coefficient = unname(stats::coef(fit)[[1]]),
    hazard_ratio = model_tidy$estimate[[1]],
    conf_low = model_tidy$conf.low[[1]],
    conf_high = model_tidy$conf.high[[1]],
    p_value = model_tidy$p.value[[1]],
    ph_p_value = ph_p,
    model_status = if_else(
      length(captured_warnings) == 0L,
      "converged",
      "warning"
    ),
    warning_text = if_else(
      length(captured_warnings) == 0L,
      NA_character_,
      paste(unique(captured_warnings), collapse = " | ")
    )
  )
}

# ------------------------------------------------------------------------------
# SECTION: run-cox-models
# ------------------------------------------------------------------------------

gene_annotations <- signature %>%
  dplyr::select(
    gene,
    direction,
    r_rb,
    correlation_p_value,
    correlation_padj
  ) %>%
  rename_with(
    ~ paste0("signature_", .x),
    -gene
  )

cox_results <- gene_changes %>%
  group_by(gene) %>%
  group_modify(~ fit_gene_rfs_cox(.x)) %>%
  ungroup() %>%
  left_join(gene_annotations, by = "gene") %>%
  mutate(
    fdr_all_genes = p.adjust(p_value, method = "BH")
  ) %>%
  group_by(signature_direction) %>%
  mutate(
    fdr_within_signature = p.adjust(p_value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    association_direction = case_when(
      hazard_ratio > 1 ~
        "Higher change associated with higher relapse hazard",
      hazard_ratio < 1 ~
        "Higher change associated with lower relapse hazard",
      TRUE ~ NA_character_
    ),
    fdr_supported = !is.na(fdr_all_genes) &
      fdr_all_genes < cox_fdr_threshold,
    nominal_p_below_0_05 =
      !is.na(p_value) & p_value < 0.05,
    expected_survival_direction = case_when(
      signature_direction == "Dys-CIM" ~
        "Detrimental (HR > 1)",
      signature_direction == "Fun-CIM" ~
        "Beneficial (HR < 1)",
      TRUE ~ NA_character_
    ),
    observed_survival_direction = case_when(
      hazard_ratio > 1 ~ "Detrimental",
      hazard_ratio < 1 ~ "Beneficial",
      hazard_ratio == 1 ~ "Null",
      TRUE ~ NA_character_
    ),
    direction_concordant = case_when(
      signature_direction == "Dys-CIM" &
        hazard_ratio > 1 ~ TRUE,
      signature_direction == "Fun-CIM" &
        hazard_ratio < 1 ~ TRUE,
      is.finite(hazard_ratio) ~ FALSE,
      TRUE ~ NA
    )
  ) %>%
  arrange(p_value, gene)

readr::write_csv(
  cox_results,
  file.path(gene_result_dir, "RFS_gene_level_cox.csv")
)

complete_rfs_results <- cox_results %>%
  arrange(signature_direction, p_value, gene)

readr::write_csv(
  complete_rfs_results,
  file.path(
    gene_result_dir,
    "RFS_complete_correlation_survival_report.csv"
  )
)

readr::write_csv(
  cox_results %>%
    filter(
      model_status != "converged" |
        !is.na(warning_text)
    ),
  file.path(
    gene_result_dir,
    "RFS_gene_model_diagnostics.csv"
  )
)


# ==============================================================================
# 9. GENE-LEVEL RFS RESULTS AND DIRECTIONAL-CONCORDANCE SUMMARIES
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION: results-summary
# ------------------------------------------------------------------------------

model_summary <- cox_results %>%
  summarise(
    endpoint = "RFS",
    genes_attempted = n(),
    models_converged_without_warning =
      sum(model_status == "converged"),
    models_with_warning =
      sum(model_status == "warning"),
    models_not_fitted_or_error =
      sum(model_status %in% c("not_fitted", "error")),
    nominal_p_below_0_05 =
      sum(nominal_p_below_0_05, na.rm = TRUE),
    fdr_supported =
      sum(fdr_supported, na.rm = TRUE),
    ph_p_below_0_05 =
      sum(ph_p_value < 0.05, na.rm = TRUE)
  )

knitr::kable(
  model_summary,
  caption = "Gene-level RFS Cox model summary"
)

top_results_table <- cox_results %>%
  filter(
    is.finite(hazard_ratio),
    is.finite(conf_low),
    is.finite(conf_high)
  ) %>%
  slice_min(
    p_value,
    n = 10,
    with_ties = FALSE
  ) %>%
  dplyr::select(
    endpoint, gene, signature_direction,
    n_patients, n_events,
    hazard_ratio, conf_low, conf_high,
    p_value, fdr_all_genes,
    ph_p_value, model_status
  )

knitr::kable(
  top_results_table,
  digits = 3,
  caption = "Top exploratory gene-level RFS associations"
)

# Directional concordance is prespecified as Dys-CIM HR > 1 (detrimental) and
# Fun-CIM HR < 1 (beneficial). It is summarized descriptively among all genes and
# within three evidence subsets. These counts do not replace gene-level p-values
# or multiple-testing correction.

# ------------------------------------------------------------------------------
# SECTION: direction-summary
# ------------------------------------------------------------------------------
evidence_definitions <- tribble(
  ~evidence_subset, ~subset_order,
  "All genes", 1L,
  "Nominal p < 0.05", 2L,
  "Global FDR < 0.05", 3L,
  "Signature FDR < 0.05", 4L
)

summarise_direction_subset <- function(data, subset_name, subset_order) {
  subset_data <- switch(
    subset_name,
    "All genes" = data,
    "Nominal p < 0.05" = filter(data, p_value < 0.05),
    "Global FDR < 0.05" = filter(data, fdr_all_genes < 0.05),
    "Signature FDR < 0.05" = filter(data, fdr_within_signature < 0.05)
  ) %>%
    filter(is.finite(hazard_ratio), !is.na(direction_concordant)) %>%
    arrange(p_value, gene)

  concordant <- subset_data %>% filter(direction_concordant)

  tibble(
    evidence_subset = subset_name,
    subset_order = subset_order,
    genes_tested = nrow(subset_data),
    concordant_n = nrow(concordant),
    concordant_percent = if_else(
      nrow(subset_data) > 0,
      100 * nrow(concordant) / nrow(subset_data),
      NA_real_
    ),
    concordant_n_percent = if_else(
      nrow(subset_data) > 0,
      sprintf("%d (%.1f%%)", nrow(concordant), 100 * nrow(concordant) / nrow(subset_data)),
      "0 (NA)"
    ),
    concordant_genes = if_else(
      nrow(concordant) > 0,
      paste(concordant$gene, collapse = ", "),
      "None"
    )
  )
}

direction_summary <- complete_rfs_results %>%
  filter(signature_direction %in% c("Dys-CIM", "Fun-CIM")) %>%
  group_split(signature_direction) %>%
  map_dfr(function(signature_data) {
    signature_name <- unique(signature_data$signature_direction)
    evidence_definitions %>%
      pmap_dfr(~ summarise_direction_subset(signature_data, ..1, ..2)) %>%
      mutate(signature = signature_name, .before = 1)
  }) %>%
  mutate(
    signature = factor(signature, levels = c("Dys-CIM", "Fun-CIM"))
  ) %>%
  arrange(signature, subset_order) %>%
  mutate(signature = as.character(signature))

readr::write_csv(
  direction_summary,
  file.path(result_dir, "RFS_directional_concordance_summary.csv")
)

kable(
  direction_summary %>%
    dplyr::select(signature, evidence_subset, genes_tested,
           concordant_n_percent, concordant_genes),
  caption = paste(
    "RFS directional concordance. Gene lists are ordered by increasing",
    "standard-Cox p-value."
  )
)


# For the figure, genes are classified by their observed RFS association rather
# than by signature concordance. Genes with raw Cox p >= 0.05 are shown as NS.
# Among genes with raw Cox p < 0.05, HR < 1 is labeled Beneficial and HR > 1 is
# labeled Detrimental. Segment labels report both the number and percentage of
# tested genes; labels for segments below 10% are placed outside the bar with a
# leader line. No between-signature hypothesis test is displayed.

# ------------------------------------------------------------------------------
# SECTION: direction-plot
# ------------------------------------------------------------------------------
# ==============================================================================
# Directional concordance of gene-level RFS associations
# Publication version for multi-panel figure
# ==============================================================================

direction_plot_data <- complete_rfs_results %>%
  filter(
    signature_direction %in% c("Dys-CIM", "Fun-CIM"),
    is.finite(hazard_ratio),
    is.finite(p_value)
  ) %>%
  mutate(
    signature = factor(
      signature_direction,
      levels = c("Dys-CIM", "Fun-CIM")
    ),
    outcome_class = case_when(
      p_value >= 0.05 ~ "NS",
      hazard_ratio < 1 ~ "Beneficial",
      hazard_ratio > 1 ~ "Detrimental",
      TRUE ~ "NS"
    ),
    outcome_class = factor(
      outcome_class,
      levels = c("Beneficial", "Detrimental", "NS")
    )
  ) %>%
  count(signature, outcome_class, name = "n") %>%
  tidyr::complete(
    signature,
    outcome_class,
    fill = list(n = 0L)
  ) %>%
  group_by(signature) %>%
  arrange(outcome_class, .by_group = TRUE) %>%
  mutate(
    total_n = sum(n),
    percent = 100 * n / total_n,

    y_max = cumsum(percent),
    y_min = lag(y_max, default = 0),
    y_mid = (y_min + y_max) / 2,

    # Cleaner labels for a small multi-panel figure.
    inside_label = if_else(
      n > 0,
      sprintf("%d\n%.1f%%", n, percent),
      ""
    ),

    # Only very small segments are labeled externally.
    external_label = if_else(
      n > 0,
      sprintf("%d (%.1f%%)", n, percent),
      ""
    ),

    signature_x = as.numeric(signature),

    # The 9.8% Fun-CIM beneficial segment remains inside.
    label_inside = percent >= 7.5 & n > 0
  ) %>%
  ungroup()


# ------------------------------------------------------------------------------
# Total n incorporated directly into the x-axis labels
# ------------------------------------------------------------------------------

signature_totals <- direction_plot_data %>%
  distinct(signature, total_n) %>%
  arrange(signature)

signature_axis_labels <- setNames(
  paste0(
    signature_totals$signature,
    "\n(n = ",
    signature_totals$total_n,
    ")"
  ),
  as.character(signature_totals$signature)
)


# ------------------------------------------------------------------------------
# External annotation data
# Only tiny non-zero segments should appear here.
# ------------------------------------------------------------------------------

external_data <- direction_plot_data %>%
  filter(!label_inside, n > 0) %>%
  mutate(
    external_y = case_when(
      outcome_class == "Beneficial" ~ pmax(y_mid - 1.5, 2),
      outcome_class == "Detrimental" ~ pmin(y_mid + 1.5, 98),
      TRUE ~ y_mid
    )
  )


# ==============================================================================
# Plot
# ==============================================================================

direction_plot <- ggplot(direction_plot_data) +

  # --------------------------------------------------------------------------
  # Stacked bars
  # Slightly wider bars improve visibility after multi-panel reduction.
  # --------------------------------------------------------------------------
  geom_rect(
    aes(
      xmin = signature_x - 0.36,
      xmax = signature_x + 0.36,
      ymin = y_min,
      ymax = y_max,
      fill = outcome_class
    ),
    color = NA
  ) +

  # --------------------------------------------------------------------------
  # Internal labels
  # Larger than before for final assembled-figure readability.
  # --------------------------------------------------------------------------
  geom_text(
    data = direction_plot_data %>%
      filter(label_inside),
    aes(
      x = signature_x,
      y = y_mid,
      label = inside_label
    ),
    size = 2.75,
    lineheight = 0.90,
    family = "sans",
    fontface = "bold",
    color = "black"
  ) +

  # --------------------------------------------------------------------------
  # Leader line for tiny segments
  # --------------------------------------------------------------------------
  geom_segment(
    data = external_data,
    aes(
      x = signature_x + 0.36,
      xend = signature_x + 0.45,
      y = y_mid,
      yend = external_y
    ),
    linewidth = 0.30,
    color = "grey40"
  ) +

  geom_segment(
    data = external_data,
    aes(
      x = signature_x + 0.45,
      xend = signature_x + 0.50,
      y = external_y,
      yend = external_y
    ),
    linewidth = 0.30,
    color = "grey40"
  ) +

  # Tiny category label: readable but intentionally less visually dominant.
  geom_text(
    data = external_data,
    aes(
      x = signature_x + 0.52,
      y = external_y,
      label = external_label
    ),
    hjust = 0,
    size = 2.45,
    family = "sans",
    fontface = "plain",
    color = "grey15"
  ) +

  # --------------------------------------------------------------------------
  # Colors
  # --------------------------------------------------------------------------
  scale_fill_manual(
    values = c(
      "Beneficial" = fun_color,
      "Detrimental" = dys_color,
      "NS" = "#E5E5E5"
    ),
    breaks = c(
      "Beneficial",
      "Detrimental",
      "NS"
    ),
    labels = c(
      "Beneficial",
      "Detrimental",
      "NS"
    )
  ) +

  # --------------------------------------------------------------------------
  # X-axis: signature + total n
  # --------------------------------------------------------------------------
  scale_x_continuous(
    breaks = c(1, 2),
    labels = signature_axis_labels
  ) +

  # --------------------------------------------------------------------------
  # Percentage axis
  # --------------------------------------------------------------------------
  scale_y_continuous(
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +

  # Enough right margin for the tiny external label.
  coord_cartesian(
    xlim = c(0.45, 2.72),
    ylim = c(0, 101),
    clip = "off"
  ) +

  labs(
    x = NULL,
    y = "Percentage of tested genes",
    fill = NULL
  ) +

  # --------------------------------------------------------------------------
  # Publication styling
  # Typography deliberately enlarged for multi-panel assembly.
  # --------------------------------------------------------------------------
  theme_classic(
    base_size = 8.5,
    base_family = "sans"
  ) +

  theme(
    axis.line.x = element_line(
      linewidth = 0.40,
      color = "black"
    ),

    axis.line.y = element_line(
      linewidth = 0.40,
      color = "black"
    ),

    axis.ticks = element_line(
      linewidth = 0.35,
      color = "black"
    ),

    axis.text.x = element_text(
      color = "black",
      face = "bold",
      size = 8.8,
      lineheight = 0.92,
      margin = margin(t = 4)
    ),

    axis.text.y = element_text(
      color = "black",
      size = 8.0
    ),

    axis.title.y = element_text(
      color = "black",
      size = 8.8,
      face = "bold",
      margin = margin(r = 5)
    ),

    legend.position = "top",
    legend.justification = "center",

    legend.text = element_text(
      size = 7.8,
      color = "black"
    ),

    legend.key.size = grid::unit(
      0.095,
      "in"
    ),

    legend.spacing.x = grid::unit(
      0.05,
      "in"
    ),

    # Slightly more room above the bars so the legend does not feel crowded.
    plot.margin = margin(
      t = 4,
      r = 16,
      b = 4,
      l = 3
    )
  ) +

  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )


print(direction_plot)


# ==============================================================================
# Export
# Slightly larger than the previous standalone panel to preserve readability.
# ==============================================================================

ggsave(
  file.path(
    figure_dir,
    "RFS_directional_concordance.pdf"
  ),
  direction_plot,
  device = grDevices::cairo_pdf,
  width = 2.80,
  height = 2.90,
  units = "in"
)

ggsave(
  file.path(
    figure_dir,
    "RFS_directional_concordance.png"
  ),
  direction_plot,
  width = 2.80,
  height = 2.90,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    figure_dir,
    "RFS_directional_concordance.tiff"
  ),
  direction_plot,
  width = 2.80,
  height = 2.90,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ==============================================================================
# 10. MAIN RFS TOP-15 FOREST PLOTS
# ==============================================================================

# Separate Dys-CIM and Fun-CIM forest panels show the top
# `r params$top_n_forest` RFS associations within each signature, ranked
# exclusively by increasing raw standard-Cox p-value **after restricting to the
# direction expected for that signature**: HR > 1 for Dys-CIM and HR < 1 for
# Fun-CIM. Opposite-direction genes remain in the complete results and
# directional-concordance tables but are excluded from these forest panels.
# Ranking by raw p-value is a visualization rule and does not imply
# multiple-testing-adjusted significance.
# The forest axis shows the HR on a logarithmic scale, with HR = 1 as the null.
# To keep the compact panels readable, graphical confidence intervals are capped
# at prespecified display limits without arrowheads.
# The adjacent table always reports the complete, uncapped HR (95% CI) and raw
# p-value.

# ------------------------------------------------------------------------------
# SECTION: rfs-forest-functions
# ------------------------------------------------------------------------------
format_effect_p <- function(x) {
  case_when(
    is.na(x) ~ "NA",
    x < 0.001 ~ formatC(x, format = "e", digits = 1),
    TRUE ~ formatC(x, format = "f", digits = 3)
  )
}

format_effect_number <- function(x) {
  ifelse(
    is.na(x), "NA",
    ifelse(abs(x) >= 100, formatC(x, format = "e", digits = 1),
           formatC(x, format = "f", digits = 2))
  )
}

select_forest_genes <- function(results, signature_name, top_n = 15) {
  results %>%
    filter(
      endpoint == "RFS",
      signature_direction == signature_name,
      is.finite(coefficient),
      is.finite(hazard_ratio),
      is.finite(p_value),
      p_value > 0,
      hazard_ratio > 0
    ) %>%
    mutate(
      survival_direction = case_when(
        hazard_ratio < 1 ~ "Beneficial",
        hazard_ratio > 1 ~ "Detrimental",
        TRUE ~ "Neutral"
      )
    ) %>%
    filter(
      (signature_name == "Dys-CIM" & survival_direction == "Detrimental") |
        (signature_name == "Fun-CIM" & survival_direction == "Beneficial")
    ) %>%
    slice_min(p_value, n = top_n, with_ties = FALSE) %>%
    # Order rows by raw Cox p-value so the strongest association is shown
    # at the top of the forest plot (smallest to largest p-value).
    arrange(p_value, gene) %>%
    mutate(
      p_label = format_effect_p(p_value),
      hr_ci_label = paste0(
        format_effect_number(hazard_ratio), " (",
        format_effect_number(conf_low), "-",
        format_effect_number(conf_high), ")"
      ),
      gene_plot = factor(gene, levels = rev(unique(gene)))
    )
}

make_forest_plot <- function(results, signature_name, top_n = 15,
                             display_limits = c(0.05, 20)) {
  plot_data <- select_forest_genes(
    results, signature_name, top_n
  )

  if (nrow(plot_data) == 0) {
    return(
      ggplot() +
        annotate(
          "text", x = 0, y = 0,
          label = paste0(signature_name, ": no estimable RFS associations"),
          family = "sans", size = 3
        ) +
        xlim(-1, 1) + ylim(-1, 1) + theme_void(base_family = "sans")
    )
  }

  signature_color <- if (signature_name == "Fun-CIM") fun_color else dys_color
  plot_data <- plot_data %>%
    mutate(
      hr_display = pmin(pmax(hazard_ratio, display_limits[1]), display_limits[2]),
      low_display = pmax(conf_low, display_limits[1]),
      high_display = pmin(conf_high, display_limits[2])
    )

  forest_panel <- ggplot(plot_data, aes(y = gene_plot)) +
    geom_vline(
      xintercept = 1, linetype = "dashed", linewidth = 0.32,
      color = "grey45"
    ) +
    geom_segment(
      aes(x = low_display, xend = high_display, yend = gene_plot),
      linewidth = 0.48, color = "#303030"
    ) +
    geom_point(
      aes(x = hr_display), shape = 21, size = 2.15, stroke = 0.35,
      fill = signature_color, color = "black"
    ) +
    scale_x_log10(
      limits = display_limits,
      breaks = c(0.05, 0.2, 1, 5, 20),
      labels = c("0.05", "0.2", "1", "5", "20"),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      title = signature_name,
      subtitle = paste0("Top ", top_n, " concordant by raw Cox p-value"),
      x = "RFS hazard ratio (log scale)", y = NULL
    ) +
    theme_classic(base_size = 8, base_family = "sans") +
    theme(
      axis.text.x = element_text(color = "black", size = 6.4),
      axis.text.y = element_text(color = "black", face = "bold", size = 7.0),
      axis.title.x = element_text(size = 6.8, face = "bold", margin = margin(t = 2)),
      plot.title = element_text(face = "bold", size = 9.2, color = signature_color,
                                margin = margin(b = 1)),
      plot.subtitle = element_text(size = 6.5, color = "grey25",
                                   margin = margin(b = 3)),
      plot.margin = margin(t = 2, r = 1, b = 2, l = 1)
    )

  table_panel <- ggplot(plot_data, aes(y = gene_plot)) +
    geom_rect(
      aes(ymin = as.numeric(gene_plot) - 0.5,
          ymax = as.numeric(gene_plot) + 0.5),
      xmin = -Inf, xmax = Inf,
      fill = rep(c("#F4F4F4", "white"), length.out = nrow(plot_data)),
      color = NA
    ) +
    geom_text(
      aes(x = 0, label = hr_ci_label), hjust = 0,
      size = 2.05, family = "sans", fontface = "plain", color = "#111111"
    ) +
    geom_text(
      aes(x = 1, label = p_label), hjust = 1,
      size = 2.10, family = "sans", fontface = "bold", color = "#111111"
    ) +
    annotate(
      "text", x = 0, y = Inf, label = "HR (95% CI)", hjust = 0,
      vjust = 1.25, size = 2.15, family = "sans", fontface = "bold",
      color = signature_color
    ) +
    annotate(
      "text", x = 1, y = Inf, label = "p value", hjust = 1,
      vjust = 1.25, size = 2.15, family = "sans", fontface = "bold",
      color = signature_color
    ) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
    theme_void(base_family = "sans") +
    theme(plot.margin = margin(t = 27, r = 1, b = 18, l = 1))

  gridExtra::arrangeGrob(
    forest_panel, table_panel,
    ncol = 2, widths = c(1.18, 0.82)
  )
}

# ------------------------------------------------------------------------------
# SECTION: rfs-forest-plots
# ------------------------------------------------------------------------------
export_forest <- function(signature_name) {
  selected <- select_forest_genes(
    cox_results, signature_name, params$top_n_forest
  )
  plot_object <- make_forest_plot(
    cox_results, signature_name, params$top_n_forest
  )
  filename_stub <- str_replace_all(signature_name, "-", "_")

  ggsave(
    file.path(figure_dir, paste0("RFS_", filename_stub, "_top", params$top_n_forest, "_forest.pdf")),
    plot_object, device = grDevices::cairo_pdf,
    width = 2.85, height = 4.00, units = "in"
  )
  ggsave(
    file.path(figure_dir, paste0("RFS_", filename_stub, "_top", params$top_n_forest, "_forest.png")),
    plot_object, width = 2.85, height = 4.00, units = "in",
    dpi = 600, bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0("RFS_", filename_stub, "_top", params$top_n_forest, "_forest.tiff")),
    plot_object, width = 2.85, height = 4.00, units = "in",
    dpi = 600, compression = "lzw", bg = "white"
  )

  readr::write_csv(
    selected %>%
      dplyr::select(signature_direction, gene, survival_direction, coefficient,
             hazard_ratio, conf_low, conf_high, p_value,
             fdr_all_genes, fdr_within_signature),
    file.path(result_dir, paste0("RFS_", filename_stub, "_top", params$top_n_forest,
                                "_forest_gene_table.csv"))
  )
  plot_object
}

dys_rfs_forest <- export_forest("Dys-CIM")
fun_rfs_forest <- export_forest("Fun-CIM")

print(dys_rfs_forest)
print(fun_rfs_forest)


# ==============================================================================
# 11. PATIENT-LEVEL RFS GENE-CHANGE HEATMAP
# ==============================================================================

# The heatmap shows the gene-wise z-score of the paired expression change
# ($\Delta GE$). First, paired change is calculated as
# $\Delta GE = \mathrm{Surgery} - \mathrm{Baseline}$ using the
# $\log_2(\mathrm{TPM}+1)$ values. Then, for each gene, these patient-level
# changes are standardized across patients:

# $$
# Z(\Delta GE)_{ig} =
# \frac{\Delta GE_{ig} - \overline{\Delta GE}_{g}}
# {SD(\Delta GE_{g})}.
# $$

# Positive values therefore indicate that a patient's expression change is
# higher than the mean change for that gene across patients; negative values
# indicate a lower-than-average change; and zero indicates a change close to the
# gene-specific mean. Because the plotted values are gene-wise z-scores, color
# intensity should not be interpreted as the original expression-unit change.
# This color scale does not represent hazard or statistical significance. Raw
# $\log_2(\mathrm{TPM}+1)$ changes remain available in the exported tables.
# Patients are clustered by similarity in their gene-wise $Z(\Delta GE)$
# profiles. The only clinical annotations shown are recurrence and RFS months.

# ------------------------------------------------------------------------------
# SECTION: heatmap-functions
# ------------------------------------------------------------------------------
heatmap_patient_annotation <- included_patients %>%
  transmute(
    patient_id,
    Recurrence = factor(if_else(rfs_event == 1L, "Yes", "No"),
                        levels = c("No", "Yes")),
    `RFS (months)` = rfs_time
  ) %>%
  arrange(Recurrence, desc(`RFS (months)`)) %>%
  column_to_rownames("patient_id")

annotation_colors <- list(
  Recurrence = c("No" = "#D9D9D9", "Yes" = "#C94747"),
  `RFS (months)` = c("#F2F2F2", "#252525")
)

heatmap_breaks <- seq(-2.5, 2.5, length.out = 101)
heatmap_colors <- colorRampPalette(c(decrease_color, "white", increase_color))(100)

prepare_heatmap_matrix <- function(selected_genes) {
  gene_changes %>%
    filter(gene %in% selected_genes) %>%
    dplyr::select(gene, patient_id, delta_z) %>%
    pivot_wider(names_from = patient_id, values_from = delta_z) %>%
    column_to_rownames("gene") %>%
    as.matrix()
}

save_change_heatmap <- function(
    selected_genes,
    filename_stem,
    title,
    width = 8,
    height = 8,
    cluster_rows = TRUE) {

  heatmap_matrix <- prepare_heatmap_matrix(selected_genes)
  patient_order <- intersect(rownames(heatmap_patient_annotation),
                             colnames(heatmap_matrix))
  # Explicit row indices avoid an RStudio 2026.07 IDE bug triggered by x[, j].
  heatmap_matrix <- heatmap_matrix[
    seq_len(nrow(heatmap_matrix)),
    patient_order,
    drop = FALSE
  ]
  annotation <- heatmap_patient_annotation[patient_order, , drop = FALSE]

  row_annotation <- signature %>%
    filter(gene %in% rownames(heatmap_matrix)) %>%
    dplyr::select(gene, Signature = direction) %>%
    column_to_rownames("gene")

  row_annotation <- row_annotation[rownames(heatmap_matrix), , drop = FALSE]

  heatmap_plot <- pheatmap::pheatmap(
    heatmap_matrix,
    color = heatmap_colors,
    breaks = heatmap_breaks,
    border_color = NA,
    cluster_rows = cluster_rows,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = annotation,
    annotation_row = row_annotation,
    annotation_colors = c(
      annotation_colors,
      list(Signature = c("Fun-CIM" = fun_color, "Dys-CIM" = dys_color))
    ),
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize = 8.5,
    fontsize_row = 7.5,
    fontsize_col = 7,
    angle_col = 90,
    main = title,
    legend_breaks = c(-2, -1, 0, 1, 2),
    legend_labels = c(
      "≤ -2  lower Z(ΔGE)", "-1", "0  mean Z(ΔGE)", "1", "≥ 2  higher Z(ΔGE)"
    ),
    width = width,
    height = height,
    silent = TRUE
  )

  ggsave(
    file.path(figure_dir, paste0(filename_stem, ".pdf")),
    heatmap_plot$gtable, device = grDevices::cairo_pdf,
    width = width, height = height, units = "in", bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(filename_stem, ".png")),
    heatmap_plot$gtable,
    width = width, height = height, units = "in", dpi = 600, bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(filename_stem, ".tiff")),
    heatmap_plot$gtable,
    width = width, height = height, units = "in", dpi = 600,
    compression = "lzw", bg = "white"
  )

  invisible(heatmap_plot)
}

# ------------------------------------------------------------------------------
# SECTION: create-rfs-heatmap
# ------------------------------------------------------------------------------
rfs_top_heatmap_genes <- cox_results %>%
  filter(
    endpoint == "RFS",
    is.finite(hazard_ratio),
    is.finite(p_value)
  ) %>%
  slice_min(p_value, n = params$top_n_heatmap, with_ties = FALSE) %>%
  arrange(p_value, gene) %>%
  pull(gene)

save_change_heatmap(
  rfs_top_heatmap_genes,
  filename_stem = "RFS_top_genes_change_heatmap",
  title = paste0("Top ", params$top_n_heatmap, " RFS-associated genes: Z-score of ΔGE"),
  width = 7.2,
  height = 6.4
)


# ==============================================================================
# 12. OPTIONAL FIRTH-PENALIZED COX SENSITIVITY ANALYSIS
# ==============================================================================

# Firth-penalized Cox regression can reduce monotone-likelihood bias in sparse
# event settings. It is treated as a sensitivity analysis, not as a substitute
# for external validation.

# ------------------------------------------------------------------------------
# SECTION: firth-sensitivity
# ------------------------------------------------------------------------------

if (isTRUE(params$run_firth_sensitivity)) {
  if (!requireNamespace("coxphf", quietly = TRUE)) {
    stop(
      "run_firth_sensitivity is TRUE, but package 'coxphf' is not installed."
    )
  }

  fit_firth_gene_rfs <- function(data) {
    model_data <- data %>%
      transmute(
        time = rfs_time,
        event = rfs_event,
        delta_z = delta_z
      ) %>%
      filter(complete.cases(.))

    fit <- tryCatch(
      coxphf::coxphf(
        survival::Surv(time, event) ~ delta_z,
        data = model_data
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(
        tibble(
          endpoint = "RFS",
          firth_hr = NA_real_,
          firth_conf_low = NA_real_,
          firth_conf_high = NA_real_,
          firth_p_value = NA_real_,
          firth_status = conditionMessage(fit)
        )
      )
    }

    tibble(
      endpoint = "RFS",
      firth_hr = exp(fit$coefficients[[1]]),
      firth_conf_low = fit$ci.lower[[1]],
      firth_conf_high = fit$ci.upper[[1]],
      firth_p_value = fit$prob[[1]],
      firth_status = "converged"
    )
  }

  firth_results <- gene_changes %>%
    group_by(gene) %>%
    group_modify(~ fit_firth_gene_rfs(.x)) %>%
    ungroup() %>%
    mutate(
      firth_fdr = p.adjust(
        firth_p_value,
        method = "BH"
      )
    ) %>%
    left_join(gene_annotations, by = "gene")

  readr::write_csv(
    firth_results,
    file.path(
      sensitivity_dir,
      "RFS_firth_gene_level_cox_results.csv"
    )
  )

  firth_rfs <- firth_results %>%
    dplyr::select(
      gene,
      firth_hr,
      firth_conf_low,
      firth_conf_high,
      firth_p_value,
      firth_fdr,
      firth_status
    )

  complete_rfs_results <- complete_rfs_results %>%
    left_join(
      firth_rfs,
      by = "gene"
    ) %>%
    mutate(
      standard_firth_direction_agreement =
        case_when(
          !is.finite(hazard_ratio) |
            !is.finite(firth_hr) ~ NA,
          hazard_ratio > 1 &
            firth_hr > 1 ~ TRUE,
          hazard_ratio < 1 &
            firth_hr < 1 ~ TRUE,
          hazard_ratio == 1 &
            firth_hr == 1 ~ TRUE,
          TRUE ~ FALSE
        )
    ) %>%
    arrange(
      signature_direction,
      p_value,
      gene
    )

  readr::write_csv(
    complete_rfs_results,
    file.path(
      gene_result_dir,
      "RFS_complete_correlation_survival_report.csv"
    )
  )

  firth_comparison_summary <- complete_rfs_results %>%
    summarise(
      genes_with_both_estimates = sum(
        is.finite(hazard_ratio) &
          is.finite(firth_hr),
        na.rm = TRUE
      ),
      same_effect_direction = sum(
        standard_firth_direction_agreement %in% TRUE,
        na.rm = TRUE
      ),
      opposite_effect_direction = sum(
        standard_firth_direction_agreement %in% FALSE,
        na.rm = TRUE
      )
    )

  readr::write_csv(
    firth_comparison_summary,
    file.path(
      sensitivity_dir,
      "RFS_standard_vs_firth_summary.csv"
    )
  )
}


# ==============================================================================
# 13. SIGNATURE-LEVEL SSGSEA RFS ANALYSIS
# ==============================================================================

# The ssGSEA module uses the **same selected Fun-CIM and Dys-CIM genes** as the
# gene-level module. ssGSEA is run on each original sample expression profile;
# paired subtraction occurs only after baseline and surgery scores have been
# calculated. This preserves the intended single-sample rank-based scoring
# procedure.


# ------------------------------------------------------------------------------
# 13.1 SSGSEA GENE-COVERAGE QC
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# SECTION: ssgsea-coverage-qc
# ------------------------------------------------------------------------------
ssgsea_gene_availability <- purrr::imap_dfr(
  gene_sets_supplied,
  ~ tibble(
    signature = .y,
    gene = unique(.x),
    present = unique(.x) %in% expression_data$gene
  )
)

ssgsea_coverage_summary <- ssgsea_gene_availability %>%
  group_by(signature) %>%
  summarise(
    genes_supplied = n(),
    genes_retained = sum(present),
    genes_missing = sum(!present),
    coverage = genes_retained / genes_supplied,
    .groups = "drop"
  )

readr::write_csv(
  ssgsea_gene_availability,
  file.path(qc_dir, "ssGSEA_signature_gene_availability.csv")
)

readr::write_csv(
  ssgsea_coverage_summary,
  file.path(qc_dir, "ssGSEA_signature_coverage_summary.csv")
)

failed_ssgsea_coverage <- ssgsea_coverage_summary %>%
  filter(
    genes_retained < as.integer(params$minimum_genes) |
      coverage < as.numeric(params$minimum_coverage)
  )

if (nrow(failed_ssgsea_coverage) > 0L) {
  stop(
    "At least one CIM signature failed the prespecified ssGSEA gene-coverage ",
    "threshold. Review 01_QC/ssGSEA_signature_coverage_summary.csv."
  )
}

gene_sets_filtered <- split(
  filter(ssgsea_gene_availability, present)$gene,
  filter(ssgsea_gene_availability, present)$signature
)

knitr::kable(
  ssgsea_coverage_summary,
  digits = 3,
  caption = "CIM signature coverage for ssGSEA"
)


# ------------------------------------------------------------------------------
# 13.2 CALCULATE SAMPLE-LEVEL SSGSEA SCORES
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# SECTION: run-ssgsea
# ------------------------------------------------------------------------------
ssgsea_scores_long <- NULL
paired_ssgsea_scores <- NULL
ssgsea_rfs_results <- NULL

if (isTRUE(params$run_ssgsea)) {
  expr_df_ssgsea <- expression_data %>%
    dplyr::select(
      gene,
      all_of(nki_sample_metadata$sample_id)
    )

  expr_mat_ssgsea <- as.matrix(
    expr_df_ssgsea[, -1, drop = FALSE]
  )
  storage.mode(expr_mat_ssgsea) <- "numeric"
  rownames(expr_mat_ssgsea) <- expr_df_ssgsea$gene

  if (
    anyNA(expr_mat_ssgsea) ||
      any(!is.finite(expr_mat_ssgsea))
  ) {
    stop(
      "Expression contains missing or non-finite values; ssGSEA was not run."
    )
  }

  sample_sd <- apply(
    expr_mat_ssgsea,
    2,
    stats::sd
  )

  if (any(!is.finite(sample_sd) | sample_sd == 0)) {
    stop(
      "At least one NKI expression sample has zero or invalid variance; ",
      "ssGSEA was not run."
    )
  }

  ssgsea_param <- GSVA::ssgseaParam(
    exprData = expr_mat_ssgsea,
    geneSets = gene_sets_filtered,
    minSize = 1,
    maxSize = Inf,
    alpha = as.numeric(params$alpha),
    normalize = isTRUE(params$normalize_scores),
    checkNA = "auto",
    use = "everything"
  )

  ssgsea_matrix <- GSVA::gsva(
    ssgsea_param,
    verbose = TRUE
  )

  ssgsea_scores_long <- as.data.frame(
    ssgsea_matrix
  ) %>%
    tibble::rownames_to_column("signature") %>%
    pivot_longer(
      -signature,
      names_to = "sample_id",
      values_to = "ssgsea_score"
    ) %>%
    left_join(
      nki_sample_metadata %>%
        dplyr::select(
          sample_id,
          patient_id,
          timepoint
        ),
      by = "sample_id"
    ) %>%
    arrange(
      signature,
      patient_id,
      timepoint
    )

  readr::write_csv(
    ssgsea_scores_long,
    file.path(
      ssgsea_score_dir,
      "CIM_ssGSEA_scores_by_sample.csv"
    )
  )
}


# ------------------------------------------------------------------------------
# 13.3 CALCULATE PAIRED SSGSEA CHANGES
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# SECTION: paired-ssgsea-change
# ------------------------------------------------------------------------------
if (isTRUE(params$run_ssgsea)) {
  paired_ssgsea_scores <- ssgsea_scores_long %>%
    dplyr::select(
      signature,
      patient_id,
      timepoint,
      ssgsea_score
    ) %>%
    pivot_wider(
      names_from = timepoint,
      values_from = ssgsea_score
    ) %>%
    mutate(
      delta_ssgsea = Surgery - Baseline
    ) %>%
    group_by(signature) %>%
    mutate(
      delta_ssgsea_mean = mean(delta_ssgsea, na.rm = TRUE),
      delta_ssgsea_sd = stats::sd(delta_ssgsea, na.rm = TRUE),
      delta_ssgsea_z = if_else(
        is.finite(delta_ssgsea_sd) & delta_ssgsea_sd > 0,
        (delta_ssgsea - delta_ssgsea_mean) / delta_ssgsea_sd,
        NA_real_
      )
    ) %>%
    ungroup() %>%
    left_join(
      survival_data %>%
        dplyr::select(
          patient_id,
          recurrence,
          rfs_time,
          rfs_event
        ),
      by = "patient_id"
    ) %>%
    arrange(
      signature,
      patient_id
    )

  readr::write_csv(
    paired_ssgsea_scores,
    file.path(
      ssgsea_score_dir,
      "CIM_ssGSEA_paired_changes.csv"
    )
  )
}


# ------------------------------------------------------------------------------
# 13.4 CONTINUOUS SSGSEA-CHANGE COX MODELS FOR RFS
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# SECTION: ssgsea-rfs-models
# ------------------------------------------------------------------------------
if (isTRUE(params$run_ssgsea)) {
  fit_ssgsea_rfs_cox <- function(data) {
    analysis_data <- data %>%
      filter(
        is.finite(delta_ssgsea_z),
        is.finite(rfs_time),
        rfs_time > 0,
        rfs_event %in% c(0L, 1L)
      )

    if (
      nrow(analysis_data) < 3L ||
        sum(analysis_data$rfs_event) < 2L ||
        !is.finite(stats::sd(
          analysis_data$delta_ssgsea_z
        )) ||
        stats::sd(
          analysis_data$delta_ssgsea_z
        ) == 0
    ) {
      return(
        tibble(
          n = nrow(analysis_data),
          events = sum(
            analysis_data$rfs_event,
            na.rm = TRUE
          ),
          HR_per_1SD_delta = NA_real_,
          conf_low = NA_real_,
          conf_high = NA_real_,
          p_value = NA_real_,
          proportional_hazards_p = NA_real_,
          model_status = "not_fitted"
        )
      )
    }

    fit <- tryCatch(
      survival::coxph(
        survival::Surv(
          rfs_time,
          rfs_event
        ) ~ delta_ssgsea_z,
        data = analysis_data,
        ties = "efron",
        x = TRUE,
        y = TRUE
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(
        tibble(
          n = nrow(analysis_data),
          events = sum(
            analysis_data$rfs_event,
            na.rm = TRUE
          ),
          HR_per_1SD_delta = NA_real_,
          conf_low = NA_real_,
          conf_high = NA_real_,
          p_value = NA_real_,
          proportional_hazards_p = NA_real_,
          model_status = paste0(
            "error: ",
            conditionMessage(fit)
          )
        )
      )
    }

    tidy_fit <- broom::tidy(
      fit,
      exponentiate = TRUE,
      conf.int = TRUE
    )

    ph_p <- tryCatch(
      unname(
        survival::cox.zph(fit)$table[
          "delta_ssgsea_z",
          "p"
        ]
      ),
      error = function(e) NA_real_
    )

    tibble(
      n = nrow(analysis_data),
      events = sum(
        analysis_data$rfs_event
      ),
      HR_per_1SD_delta =
        tidy_fit$estimate[[1]],
      conf_low =
        tidy_fit$conf.low[[1]],
      conf_high =
        tidy_fit$conf.high[[1]],
      p_value =
        tidy_fit$p.value[[1]],
      proportional_hazards_p = ph_p,
      model_status = "converged"
    )
  }

  ssgsea_rfs_results <- paired_ssgsea_scores %>%
    group_by(signature) %>%
    group_modify(
      ~ fit_ssgsea_rfs_cox(.x)
    ) %>%
    ungroup() %>%
    mutate(
      p_adjust_bh = p.adjust(
        p_value,
        method = "BH"
      )
    ) %>%
    arrange(
      factor(
        signature,
        levels = c(
          "Dys-CIM",
          "Fun-CIM"
        )
      )
    )

  readr::write_csv(
    ssgsea_rfs_results,
    file.path(
      ssgsea_result_dir,
      "CIM_delta_ssGSEA_RFS_Cox.csv"
    )
  )

  knitr::kable(
    ssgsea_rfs_results,
    digits = 4,
    caption =
      "Continuous CIM ΔssGSEA Cox models for RFS"
  )
}


# ------------------------------------------------------------------------------
# 13.5 PUBLICATION FOREST PLOT FOR SSGSEA RFS
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# SECTION: ssgsea-publication-forest
# ------------------------------------------------------------------------------
if (
  isTRUE(params$run_ssgsea) &&
    isTRUE(params$make_ssgsea_figure) &&
    !is.null(ssgsea_rfs_results)
) {
  signature_colors <- c(
    "Fun-CIM" = fun_color,
    "Dys-CIM" = dys_color
  )

  forest_data_ssgsea <- ssgsea_rfs_results %>%
    filter(
      is.finite(HR_per_1SD_delta),
      is.finite(conf_low),
      is.finite(conf_high)
    ) %>%
    mutate(
      signature = factor(
        signature,
        levels = c(
          "Dys-CIM",
          "Fun-CIM"
        )
      ),
      result_label = sprintf(
        "%.2f (%.2f\u2013%.2f)",
        HR_per_1SD_delta,
        conf_low,
        conf_high
      ),
      adjusted_p_label = case_when(
        is.na(p_adjust_bh) ~ "NA",
        p_adjust_bh < 0.001 ~ "<0.001",
        TRUE ~ sprintf(
          "%.3f",
          p_adjust_bh
        )
      )
    )

  signature_axis_labels <- c(
    "Dys-CIM" = "\u0394 ssGSEA Dys-CIM",
    "Fun-CIM" = "\u0394 ssGSEA Fun-CIM"
  )

  forest_axis_ssgsea <- ggplot(
    forest_data_ssgsea,
    aes(y = signature)
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.35,
      color = "grey65"
    ) +
    geom_errorbarh(
      aes(
        xmin = conf_low,
        xmax = conf_high,
        color = signature
      ),
      height = 0,
      linewidth = 0.55
    ) +
    geom_point(
      aes(
        x = HR_per_1SD_delta,
        color = signature
      ),
      shape = 16,
      size = 2.25
    ) +
    scale_color_manual(
      values = signature_colors,
      guide = "none"
    ) +
    scale_x_log10(
      limits = c(0.15, 12),
      breaks = c(
        0.25, 0.50, 1.00,
        2.00, 4.00, 8.00
      ),
      labels =
        scales::label_number(
          accuracy = 0.01
        )
    ) +
    scale_y_discrete(
      drop = FALSE,
      labels = signature_axis_labels
    ) +
    labs(
      x =
        "RFS hazard ratio per 1-SD increase in \u0394ssGSEA score",
      y = NULL
    ) +
    theme_classic(base_size = 8) +
    theme(
      axis.text.x =
        element_text(
          color = "black",
          size = 7.5
        ),
      axis.text.y =
        element_text(
          color = "black",
          face = "bold",
          size = 7.8,
          margin = margin(r = 3)
        ),
      axis.title.x =
        element_text(
          color = "black",
          face = "bold",
          size = 8.0,
          margin = margin(t = 4)
        ),
      axis.line =
        element_line(
          linewidth = 0.40,
          color = "black"
        ),
      axis.ticks =
        element_line(
          linewidth = 0.35,
          color = "black"
        ),
      plot.margin =
        margin(3, 1, 2, 2)
    )

  hr_table_ssgsea <- ggplot(
    forest_data_ssgsea,
    aes(y = signature)
  ) +
    geom_text(
      aes(
        x = 0.5,
        label = result_label
      ),
      hjust = 0.5,
      size = 2.8,
      color = "grey10"
    ) +
    annotate(
      "text",
      x = 0.5,
      y = 2.58,
      label = "HR (95% CI)",
      hjust = 0.5,
      fontface = "bold",
      size = 2.9,
      color = "black"
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    scale_y_discrete(drop = FALSE) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      plot.margin =
        margin(3, 0, 2, 0)
    )

  p_table_ssgsea <- ggplot(
    forest_data_ssgsea,
    aes(y = signature)
  ) +
    geom_text(
      aes(
        x = 0.5,
        label = adjusted_p_label
      ),
      hjust = 0.5,
      size = 2.8,
      color = "grey10"
    ) +
    annotate(
      "text",
      x = 0.5,
      y = 2.58,
      label = "Adjusted p",
      hjust = 0.5,
      fontface = "bold",
      size = 2.9,
      color = "black"
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    scale_y_discrete(drop = FALSE) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      plot.margin =
        margin(3, 1, 2, 0)
    )

  ssgsea_forest_panel <-
    forest_axis_ssgsea +
    hr_table_ssgsea +
    p_table_ssgsea +
    patchwork::plot_layout(
      widths = c(
        1.55,
        0.82,
        0.42
      )
    )

  ggsave(
    file.path(
      ssgsea_figure_dir,
      "CIM_delta_ssGSEA_RFS_forest.pdf"
    ),
    ssgsea_forest_panel,
    device = grDevices::cairo_pdf,
    width = 4.25,
    height = 2.15,
    units = "in"
  )

  ggsave(
    file.path(
      ssgsea_figure_dir,
      "CIM_delta_ssGSEA_RFS_forest.png"
    ),
    ssgsea_forest_panel,
    width = 4.25,
    height = 2.15,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  ggsave(
    file.path(
      ssgsea_figure_dir,
      "CIM_delta_ssGSEA_RFS_forest.tiff"
    ),
    ssgsea_forest_panel,
    width = 4.25,
    height = 2.15,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )

  print(ssgsea_forest_panel)
}


# ==============================================================================
# 14. INTERPRETATION GUIDANCE
# ==============================================================================

# The two RFS modules answer related but distinct questions. The gene-level
# analysis asks whether the paired change of an individual CIM-associated gene is
# associated with relapse hazard. The ssGSEA analysis asks whether the paired
# change in the aggregate rank-based enrichment of the full Fun-CIM or Dys-CIM
# signature is associated with relapse hazard.

# A positive delta always denotes a higher surgery/post-treatment value than the
# paired baseline value. For the predefined trajectory interpretation, increasing
# Dys-CIM change is expected to be detrimental (HR > 1), whereas increasing
# Fun-CIM change is expected to be beneficial (HR < 1). These directional
# expectations should not be treated as constraints on the fitted models.

# The Cox analyses are exploratory. Raw p-values, BH-adjusted p-values,
# confidence intervals, event counts, model diagnostics, and proportional-hazards
# tests are exported so the full evidence can be evaluated rather than relying
# only on selected figures.


# ==============================================================================
# 15. REPRODUCIBILITY
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION: session-info
# ------------------------------------------------------------------------------

version_packages <- unique(c(
  "tidyverse", "dplyr", "tidyr", "purrr", "readr", "stringr",
  "tibble", "forcats", "ggplot2", "readxl", "janitor", "survival",
  "broom", "pheatmap", "gridExtra", "RColorBrewer", "scales", "knitr",
  "patchwork", "BiocManager", "GSVA", "rmarkdown"
))

if (isTRUE(params$run_firth_sensitivity)) {
  version_packages <- c(
    version_packages,
    "coxphf"
  )
}

package_versions <- tibble::tibble(
  software = c(
    "R",
    "Bioconductor",
    version_packages
  ),
  version = c(
    paste(
      R.version$major,
      R.version$minor,
      sep = "."
    ),
    as.character(
      BiocManager::version()
    ),
    vapply(
      version_packages,
      function(package_name) {
        if (
          requireNamespace(
            package_name,
            quietly = TRUE
          )
        ) {
          as.character(
            utils::packageVersion(
              package_name
            )
          )
        } else {
          NA_character_
        }
      },
      character(1)
    )
  ),
  role = c(
    "Statistical computing environment",
    "Bioconductor release associated with the R environment",
    rep(
      "Analysis dependency",
      length(version_packages)
    )
  )
)

readr::write_csv(
  package_versions,
  file.path(
    reproducibility_dir,
    "package_versions.csv"
  )
)

# v2.0: the correlation_threshold / correlation_operator / fdr_* parameters no
# longer exist -- gene selection is driven by the meta_* settings instead.
or_null <- function(x) if (is.null(x)) "null" else as.character(x)
master_parameters <- tibble::tibble(
  parameter = c(
    "gene_set_source",
    "survival_file",
    "expression_file",
    "signature_file",
    "meta_scope",
    "meta_rho_threshold",
    "meta_k_mode",
    "meta_require_fisher",
    "meta_require_stouffer",
    "meta_padj_threshold",
    "meta_min_mean_abs_rrb",
    "meta_top_n_genes",
    "meta_rank_metric",
    "meta_ambiguous_gene_policy",
    "effective_meta_significance_stage",
    "effective_dataset_support_stage",
    "ambiguous_genes_found",
    "cox_fdr_threshold",
    "top_n_forest",
    "top_n_heatmap",
    "run_firth_sensitivity",
    "run_ssgsea",
    "ssgsea_alpha",
    "ssgsea_normalize_scores",
    "ssgsea_minimum_genes",
    "ssgsea_minimum_coverage",
    "make_ssgsea_figure"
  ),
  value = as.character(c(
    "cimic meta_results.csv",
    params$survival_file,
    params$expression_file,
    params$signature_file,
    meta_scope,
    meta_rho_threshold,
    meta_k_mode,
    meta_require_fisher,
    meta_require_stouffer,
    meta_padj_threshold,
    or_null(params$meta_min_mean_abs_rrb),
    or_null(params$meta_top_n_genes),
    meta_rank_metric,
    meta_ambiguous_gene_policy,
    correlation_fdr_stage_label,
    correlation_stage_label,
    length(ambiguous_genes),
    params$cox_fdr_threshold,
    params$top_n_forest,
    params$top_n_heatmap,
    isTRUE(params$run_firth_sensitivity),
    isTRUE(params$run_ssgsea),
    params$alpha,
    isTRUE(params$normalize_scores),
    params$minimum_genes,
    params$minimum_coverage,
    isTRUE(params$make_ssgsea_figure)
  ))
)

readr::write_csv(
  master_parameters,
  file.path(
    reproducibility_dir,
    "master_analysis_parameters.csv"
  )
)

session_info <- capture.output(
  sessionInfo()
)

writeLines(
  session_info,
  file.path(
    reproducibility_dir,
    "sessionInfo.txt"
  )
)

knitr::kable(
  package_versions,
  caption =
    "Exact software versions recorded for this analysis run"
)

cat(
  paste(
    session_info,
    collapse = "\n"
  )
)

