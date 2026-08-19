# ==============================================================================
# CIM ssGSEA Adjusted Survival Analysis: Meta-Analysis Signatures
# ==============================================================================
#
# Standalone R conversion of
#   CIM_ssGSEA_survival_external_cohorts_v3.0_meta.Rmd
#
# Behaviour is identical to knitting the Rmd. The YAML header became the
# `params` list below -- edit it there, then run the script top to bottom:
#
#   setwd("<repo root>")   # paths in params are resolved from here
#   source("oncoimmunology_paper/figure_generation_scripts/
#           CIM_ssGSEA_survival_external_cohorts_v3.0_meta.R", echo = TRUE)
#
# or from a shell:  Rscript <this file>
#
# Differences from the Rmd, all confined to the setup section:
#   * no knitr::opts_chunk$set() -- there are no chunks
#   * no rstudioapi YAML reload  -- `params` is a plain list, always current
#   * figures are written to disk exactly as before, but nothing is
#     rendered inline and no HTML report is produced
#
# Author: Gbadamosi Lab - IL/MG
# ==============================================================================

# ------------------------------------------------------------------------------
# ANALYSIS PARAMETERS  (was the Rmd YAML header)
# ------------------------------------------------------------------------------
# Every value below is read through `params$<name>`, exactly as when
# knitting. Validation is unchanged and happens in the configuration
# section further down.

params <- list(
  # Root that every *_file path below is resolved against. Repo-relative, so the
  # working directory must be the repository root (the folder containing
  # oncoimmunology_paper/). An absolute path also works if the data lives
  # elsewhere. The four cohort files sit in two different Datasets subfolders and
  # meta_results.csv sits under Results, so each entry carries its own subpath
  # rather than all sharing one flat folder.
  data_dir = "oncoimmunology_paper",
  # Written under Results/ alongside the other analysis outputs. A trailing
  # rhoNNN must agree with meta_rho_threshold (checked at run time).
  output_dir = "oncoimmunology_paper/Results/CIM_ssGSEA_survival_meta_NKI_CL_kAllC_rho030",

  # ===========================================================================
  # GENE-SET SOURCE: cimic_overlap_meta_rankbiserial.R meta_results.csv
  # ===========================================================================
  # v3.0 replaces v2.0's single-cohort correlation file with the cross-dataset
  # meta-analysis. Everything below this banner is the toggle surface for
  # choosing WHICH genes enter the two signatures. Nothing downstream of the
  # gene-sets chunk changed from v2.0.
  #
  # Produced by cimic_overlap_meta_rankbiserial.R. Resolved against data_dir.
  signature_file = "Results/cimic_meta_analysis/meta_results.csv",

  # -- WHICH DATASET COMBINATION ---------------------------------------------
  # ALL3_NKI.NEO.TNBCCL | NKI.NEO | NKI.TNBCCL | NEO.TNBCCL
  meta_scope = "ALL3_NKI.NEO.TNBCCL",

  # -- EFFECT-SIZE CUTOFF ----------------------------------------------------
  # Must be one of the rho_threshold values present in meta_results.csv (0.3, 0.6).
  meta_rho_threshold = 0.30,

  # -- HOW MANY DATASETS MUST BACK A GENE ------------------------------------
  # k_min2            k_datasets >= 2 (for ALL3: backed by 2 OR 3 datasets)
  # k_all             cleared rho in EVERY dataset of the scope
  # k_all_concordant  every dataset AND all pointing the same way (strictest)
  #
  # NOTE: k_datasets is MAGNITUDE ONLY (|r_rb| >= rho, either direction). Only
  # k_all_concordant additionally requires direction agreement.
  #
  # Signature sizes actually available (Fun / Dys), both meta tests at 0.05.
  # Combinations marked (!) put a signature below minimum_signature_genes (10)
  # and will stop the run.
  #   scope       rho   k_min2      k_all       k_all_concordant
  #   ALL3        0.3   5015 / 221  1096 / 14   550 / 6   (!)
  #   ALL3        0.6    232 /   4  (!)  10 / 0 (!)  1 / 0 (!)
  #   NKI.NEO     0.3   3065 /  23  same        3063 / 22
  #   NKI.NEO     0.6    135 /   0  (!) same (!) same (!)
  #   NKI.TNBCCL  0.3   1273 / 516  same        1079 / 402
  #   NKI.TNBCCL  0.6     20 /   7  (!) same (!)   7 / 4 (!)
  #   NEO.TNBCCL  0.3   2495 /   0  (!) same (!) 2366 / 0 (!)
  #   NEO.TNBCCL  0.6     85 /   0  (!) same (!)   85 / 0 (!)
  # For a strict "same direction in all three" run use ALL3 / rho 0.3 /
  # k_all_concordant and lower minimum_signature_genes, or use NKI.TNBCCL.
  meta_k_mode = "k_all_concordant",

  # -- WHICH META TESTS MUST PASS --------------------------------------------
  # Both true reproduces the "both tests" filter used by the ORA and volcano.
  meta_require_fisher = TRUE,
  meta_require_stouffer = TRUE,
  meta_padj_threshold = 0.05,

  # -- OPTIONAL EXTRA RESTRICTIONS -------------------------------------------
  # null disables each of these.
  meta_min_mean_abs_rrb = NULL,   # e.g. 0.5 for an additional effect-size floor
  meta_top_n_genes = NULL,   # e.g. 200 for the top N per signature
  meta_rank_metric = "stouffer_padj",   # stouffer_padj | fisher_padj | mean_abs_rrb
  # A gene can be meta-significant toward BOTH CIMs (cim_class == "Both").
  # "drop"        exclude it from both signatures (conservative; recommended)
  # "assign_best" give it to whichever CIM has the smaller ranking metric
  meta_ambiguous_gene_policy = "drop",

  # -- FULL MANUAL OVERRIDE --------------------------------------------------
  # Supply explicit symbols to bypass meta selection entirely. When either list
  # is non-null BOTH must be provided, and no meta filtering is applied.
  custom_fun_cim_genes = NULL,
  custom_dys_cim_genes = NULL,

  # Cohort files, each relative to data_dir. METABRIC and SCAN-B live in
  # separate Datasets subfolders, so the subfolder is part of each path.
  metabric_clinical_file = "Datasets/METABRIC/METABRIC_clinical_data.xlsx",
  metabric_expression_file = "Datasets/METABRIC/METABRIC_expression_data.csv",
  scanb_clinical_file = "Datasets/SCAN-B/SCANB_clinical_harmonized_OS_DRFi_RFi.rds",
  scanb_expression_file = "Datasets/SCAN-B/SCANB_expression_harmonized_OS_DRFi_RFi.rds",
  # First transcriptome-wide gene-expression column in the METABRIC CSV.
  # Change this value if the source-file layout changes; null is not accepted
  # because ssGSEA requires the full expression background.
  metabric_gene_start_col = 39,
  metabric_chemo_only = TRUE,
  scanb_chemo_only = TRUE,
  metabric_chemo_column = NULL,
  scanb_chemo_column = "Chemo_Combined",
  chemo_positive_values = c("true", "1", "yes", "y", "treated", "chemotherapy"),
  # Optional manuscript-version safeguards. Set either value to null when
  # intentionally using another selection. These are checked against the gene
  # counts produced by the meta-selection settings above.
  expected_fun_cim_genes = NULL,
  expected_dys_cim_genes = NULL,
  expected_signature_file_md5 = NULL,
  # Separate threshold used only to flag Cox-model FDR-supported results.
  cox_fdr_threshold = 0.05,
  cox_fdr_operator = "<",
  minimum_signature_coverage = 0.70,
  minimum_signature_genes = 5,
  duplicate_gene_method = "max_variance",
  missing_expression_gene_policy = "drop_gene",
  # ssGSEA settings (GSVA). These values are passed explicitly to both the
  # current and legacy GSVA interfaces.
  ssgsea_alpha = 0.25,
  ssgsea_normalize = TRUE,
  ssgsea_min_size = 1,
  ssgsea_max_size = NULL,
  bioc_parallel_workers = 1,
  min_n_model = 15,
  min_events_model = 5,
  events_per_parameter = 10,
  enforce_events_per_parameter = FALSE,
  sparse_factor_level_n = 5,
  alpha = 0.05,
  significance_metric = "p_adjust_signature",
  run_common_sample_sensitivity = TRUE,
  # Adjusted Cox exposure parameterizations.
  score_parameterizations = c("continuous", "median_split"),
  m2_covariates = c("age_5y", "tumor_size", "grade", "clinical_subtype"),
  primary_adjusted_model = "M2",
  # Median cutoffs are fixed within cohort and signature before endpoint/model exclusions.
  median_tie_rule = "high_includes_cutoff",
  median_group_warning_fraction = 0.20,
  # Optional descriptive KM/log-rank analysis; disabled unless explicitly requested.
  # Descriptive Kaplan-Meier curves + log-rank, using the fixed cohort-signature
  # median groups. Was FALSE (inherited from v2.0, where KM was opt-in), which
  # is why 08_KM_Plots came out empty. Requires the survminer package and
  # "median_split" in score_parameterizations; both are satisfied above.
  run_km_median_split = TRUE,
  km_break_months = 48,
  metabric_km_ylim = c(0.00, 1.00),
  scanb_km_ylim = c(0.75, 1.00),
  metabric_km_pval_coord = c(0.15, 0.20),
  scanb_km_pval_coord = c(0.15, 0.80),
  km_width = 7,
  km_height = 6.6,
  km_dpi = 600,
  seed = 2026,
  # Publication forest plots. The default reproduces the continuous M2 panels.
  make_forest_plots = TRUE,
  forest_model = "M2",
  forest_parameterization = "continuous",
  hr_axis_mode = "manual",
  hr_axis_min = 0.62,
  hr_axis_max = 1.65,
  hr_axis_breaks = c(0.65, 0.80, 1.00, 1.25, 1.60),
  forest_width = 7.2,
  forest_height = 3.8,
  combined_forest_width = 9.0,
  score_qc_plot_width = 7,
  score_qc_plot_height = 4.5,
  figure_dpi = 600
)

# ------------------------------------------------------------------------------
# SETUP
# ------------------------------------------------------------------------------
options(stringsAsFactors = FALSE, scipen = 999)
# The Rmd reloaded its own YAML when chunks were run interactively. In a
# plain script `params` above is always the active configuration, so the
# reload logic is replaced by a fixed provenance string.
analysis_parameter_source <- "parameters defined in the .R script params list"
analysis_rmd_path <- "CIM_ssGSEA_survival_external_cohorts_v3.0_meta.R"
analysis_seed <- suppressWarnings(as.integer(params$seed))
if (length(analysis_seed) != 1L || is.na(analysis_seed)) {
  stop("params$seed must be one integer.")
}
set.seed(analysis_seed)

# [setup chunk replaced above -- knitr/YAML-reload logic does not apply to a script]

# This standalone analysis evaluates pretreatment Fun-CIM and Dys-CIM ssGSEA
# scores using two Cox parameterizations in METABRIC and SCAN-B. By default,
# both cohorts are restricted to chemotherapy-treated patients. The exposure is
# modeled as either a continuous within-cohort standardized score or fixed
# cohort-signature median groups (High versus Low). Continuous hazard ratios are reported per
# one-standard-deviation increase. Median-split hazard ratios compare High with
# Low, with Low as the reference. Age is reported per five-year increase.
# Fun-CIM and Dys-CIM are fitted separately.

# The model set is:

# | Model | Formula after the survival outcome | Status |
# |---|---|---|
# | M0: Unadjusted | CIM score | Primary univariable |
# | M2: Prespecified adjusted | CIM score + covariates listed in `m2_covariates` | Primary adjusted |

# M2 is the prespecified primary adjusted model and uses the covariates supplied
# through `m2_covariates`. M0 is the unadjusted univariable reference. No
# intermediate adjusted models are fitted in this analysis.
# Optional Kaplan-Meier curves are descriptive and are disabled by default with
# `run_km_median_split: false`.

# The raw-score median is calculated once for each cohort-signature pair using
# all finite scores after cohort restriction and clinical-expression matching,
# before endpoint or model complete-case exclusions. The same cutoff and group
# assignment are reused across every endpoint, adjusted model, and common-sample
# sensitivity analysis. Medians are never pooled across cohorts or recalculated
# inside endpoint-specific analysis populations.

# METABRIC endpoints are OS, DSS, and RFI. SCAN-B endpoints are OS, DRFi, and
# RFi. Direct endpoint replication is limited to OS versus OS and RFI versus RFi;
# DSS and DRFi remain cohort-specific.


# ------------------------------------------------------------------------------
# META-ANALYSIS SIGNATURE DEFINITION
# ------------------------------------------------------------------------------

# **This is the only analytical difference from v2.0.** Version 2.0 built the two
# signatures by thresholding `r_rb` in a single cohort's
# `*_correlation_analysis.csv`. Version 3.0 instead draws them from
# `meta_results.csv`, the cross-dataset Fisher + Stouffer meta-analysis produced by
# `cimic_overlap_meta_rankbiserial.R` over NKI/SMC, NEO, and the TNBC cell-line
# panel. Everything after the gene-sets chunk — ssGSEA, the Cox engine, the forest
# plots, the workbook — is unchanged.

# A gene enters a signature when, for the selected `meta_scope` and
# `meta_rho_threshold`, it satisfies all of:

# 1. `cim` matches the signature (Fun or Dys);
# 2. the required meta tests pass at `meta_padj_threshold` (`meta_require_fisher`
#    and/or `meta_require_stouffer`);
# 3. the `meta_k_mode` dataset-support rule;
# 4. any optional `meta_min_mean_abs_rrb` floor;
# 5. it survives the `meta_ambiguous_gene_policy` check.


# ------------------------------------------------------------------------------
# WHAT `META_K_MODE` DOES, PRECISELY
# ------------------------------------------------------------------------------

# `k_datasets` in `meta_results.csv` counts the datasets where
# `|r_rb| >= rho_threshold`. That is **magnitude only** — it does not require
# significance within any individual dataset, and it does **not** require the
# datasets to agree on direction. A gene at `r_rb = -0.7` in one dataset and
# `+0.7` in another clears `rho = 0.6` in both while pointing at opposite CIMs.
# Direction consistency is otherwise enforced only statistically, through the
# one-sided Fisher/Stouffer p.

# Consequently:

# | `meta_k_mode` | Requires |
# |---|---|
# | `k_min2` | `k_datasets >= 2` |
# | `k_all` | cleared rho in every dataset of the scope (magnitude) |
# | `k_all_concordant` | every dataset **and** all agreeing on direction |

# Only `k_all_concordant` means "same direction in all three". Expect it to be far
# more restrictive. The YAML block for `meta_k_mode` tabulates the signature size
# for every scope / rho / k-mode combination and flags those that fall below
# `minimum_signature_genes`. The gene-sets chunk stops immediately with that
# guidance if the active combination is undersized, rather than failing later
# during cohort QC.

# The shipped default is `ALL3` / `rho 0.30` / `k_all` (1096 Fun, 14 Dys), the
# strictest combination that still gives both signatures a workable size. For a
# true same-direction-in-all-three run, use `ALL3` / `rho 0.30` /
# `k_all_concordant` (550 Fun, 6 Dys) and lower `minimum_signature_genes`, or
# switch to `NKI.TNBCCL` (1079 Fun, 402 Dys).


# ------------------------------------------------------------------------------
# AMBIGUOUS GENES
# ------------------------------------------------------------------------------

# Unlike a signed `r_rb`, the meta table can call one gene significant toward
# **both** CIMs (`cim_class == "Both"`), because the Fun and Dys tests are run
# separately. `meta_ambiguous_gene_policy: "drop"` excludes such genes from both
# signatures; `"assign_best"` gives each to whichever CIM has the smaller ranking
# metric. The count is always reported and written to the QC folder.


# ------------------------------------------------------------------------------
# RUNNING VARIANTS
# ------------------------------------------------------------------------------

# Put each parameter combination in its own `output_dir`. If the folder name ends
# in `rhoNNN`, the RMD verifies the suffix matches `meta_rho_threshold` before
# writing, e.g. `..._meta_ALL3_kAllC_rho030` with `meta_rho_threshold: 0.30`.


# ==============================================================================
# PACKAGES, PATHS, AND ANALYSIS CONFIGURATION
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION: packages-paths-configuration
# ------------------------------------------------------------------------------
required_packages <- c(
  "data.table", "readr", "readxl", "dplyr", "tidyr", "purrr", "tibble",
  "stringr", "survival", "broom", "GSVA", "BiocParallel", "ggplot2",
  "patchwork", "openxlsx", "scales", "rmarkdown"
)
if (isTRUE(params$run_km_median_split)) {
  required_packages <- unique(c(required_packages, "survminer"))
}
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(stringr); library(ggplot2)
})

validate_numeric <- function(value, parameter_name, lower = -Inf, upper = Inf,
                             lower_inclusive = TRUE, upper_inclusive = TRUE) {
  value <- suppressWarnings(as.numeric(value))
  lower_ok <- if (lower_inclusive) value >= lower else value > lower
  upper_ok <- if (upper_inclusive) value <= upper else value < upper
  if (length(value) != 1L || !is.finite(value) || !lower_ok || !upper_ok) {
    stop(parameter_name, " must be one finite numeric value in the permitted range.")
  }
  value
}
validate_integer <- function(value, parameter_name, lower = 1L) {
  value <- validate_numeric(value, parameter_name, lower = lower)
  if (value != as.integer(value)) stop(parameter_name, " must be an integer.")
  as.integer(value)
}
validate_numeric_vector <- function(value, parameter_name, length_required,
                                    lower = -Inf, upper = Inf) {
  value <- suppressWarnings(as.numeric(unlist(
    value, recursive = TRUE, use.names = FALSE)))
  if (length(value) != length_required || any(!is.finite(value)) ||
      any(value < lower) || any(value > upper)) {
    stop(parameter_name, " must contain exactly ", length_required,
         " finite numeric values in the permitted range.")
  }
  value
}
validate_choice <- function(value, choices, parameter_name) {
  value <- trimws(as.character(value))
  if (length(value) != 1L || is.na(value) || !value %in% choices) {
    stop(parameter_name, " must be one of: ", paste(choices, collapse = ", "), ".")
  }
  value
}
validate_flag <- function(value, parameter_name) {
  if (length(value) != 1L || is.na(value) || !is.logical(value)) {
    stop(parameter_name, " must be YAML true or false (without quotation marks).")
  }
  isTRUE(value)
}
validate_optional_integer <- function(value, parameter_name, lower = 0L) {
  if (is.null(value) || length(value) == 0L || all(is.na(value))) return(NULL)
  validate_integer(value, parameter_name, lower = lower)
}
apply_comparator <- function(x, cutoff, operator) {
  switch(operator,
    "<" = x < cutoff, "<=" = x <= cutoff,
    ">" = x > cutoff, ">=" = x >= cutoff,
    stop("Unsupported comparison operator: ", operator))
}
format_cutoff <- function(x) {
  format(x, scientific = FALSE, trim = TRUE, digits = 15)
}

# ── v3.0: meta-analysis gene-selection settings ─────────────────────────────
validate_optional_numeric <- function(value, parameter_name,
                                      lower = -Inf, upper = Inf) {
  if (is.null(value) || length(value) == 0L || all(is.na(value))) return(NULL)
  validate_numeric(value, parameter_name, lower = lower, upper = upper)
}
validate_optional_character_vector <- function(value, parameter_name) {
  if (is.null(value) || length(value) == 0L || all(is.na(value))) return(NULL)
  out <- toupper(trimws(as.character(unlist(
    value, recursive = TRUE, use.names = FALSE))))
  out <- unique(out[!is.na(out) & nzchar(out)])
  if (!length(out)) return(NULL)
  out
}

meta_scope <- validate_choice(
  params$meta_scope,
  c("ALL3_NKI.NEO.TNBCCL", "NKI.NEO", "NKI.TNBCCL", "NEO.TNBCCL"),
  "params$meta_scope")
meta_rho_threshold <- validate_numeric(
  params$meta_rho_threshold, "params$meta_rho_threshold", 0, 1)
meta_k_mode <- validate_choice(
  params$meta_k_mode, c("k_min2", "k_all", "k_all_concordant"),
  "params$meta_k_mode")
meta_require_fisher <- validate_flag(
  params$meta_require_fisher, "params$meta_require_fisher")
meta_require_stouffer <- validate_flag(
  params$meta_require_stouffer, "params$meta_require_stouffer")
if (!meta_require_fisher && !meta_require_stouffer) {
  stop("At least one of meta_require_fisher / meta_require_stouffer must be ",
       "true, otherwise no significance filter is applied at all.")
}
meta_padj_threshold <- validate_numeric(
  params$meta_padj_threshold, "params$meta_padj_threshold", 0, 1)
meta_min_mean_abs_rrb <- validate_optional_numeric(
  params$meta_min_mean_abs_rrb, "params$meta_min_mean_abs_rrb", 0, 1)
meta_top_n_genes <- validate_optional_integer(
  params$meta_top_n_genes, "params$meta_top_n_genes", lower = 1L)
meta_rank_metric <- validate_choice(
  params$meta_rank_metric, c("stouffer_padj", "fisher_padj", "mean_abs_rrb"),
  "params$meta_rank_metric")
meta_ambiguous_gene_policy <- validate_choice(
  params$meta_ambiguous_gene_policy, c("drop", "assign_best"),
  "params$meta_ambiguous_gene_policy")
custom_fun_cim_genes <- validate_optional_character_vector(
  params$custom_fun_cim_genes, "params$custom_fun_cim_genes")
custom_dys_cim_genes <- validate_optional_character_vector(
  params$custom_dys_cim_genes, "params$custom_dys_cim_genes")
use_custom_gene_sets <- !is.null(custom_fun_cim_genes) ||
  !is.null(custom_dys_cim_genes)
if (use_custom_gene_sets &&
    (is.null(custom_fun_cim_genes) || is.null(custom_dys_cim_genes))) {
  stop("custom_fun_cim_genes and custom_dys_cim_genes must both be supplied, ",
       "or both be null.")
}
if (use_custom_gene_sets) {
  overlapping_custom <- base::intersect(custom_fun_cim_genes,
                                        custom_dys_cim_genes)
  if (length(overlapping_custom)) {
    stop("These genes appear in BOTH custom signatures: ",
         paste(utils::head(overlapping_custom, 20), collapse = ", "),
         if (length(overlapping_custom) > 20) " ..." else "")
  }
}

cox_fdr_threshold <- validate_numeric(
  params$cox_fdr_threshold, "params$cox_fdr_threshold", 0, 1)
cox_fdr_operator <- validate_choice(
  params$cox_fdr_operator, c("<", "<="), "params$cox_fdr_operator")

# Human-readable descriptions of the active selection, reused in the QC exports,
# the reproducibility manifest, and the console summary. These replace v2.0's
# fun_correlation_rule / dys_correlation_rule / correlation_fdr_rule strings.
meta_k_mode_description <- c(
  k_min2           = "k_datasets >= 2",
  k_all            = "k_datasets == n_scope_datasets (all datasets in scope)",
  k_all_concordant = paste("k_datasets == n_scope_datasets AND",
                           "direction_concordant == YES (same direction in all)")
)[[meta_k_mode]]
meta_test_description <- paste(
  c(if (meta_require_fisher) "fisher_padj",
    if (meta_require_stouffer) "stouffer_padj"),
  collapse = " and ")
meta_significance_rule <- paste(
  meta_test_description, "<", format_cutoff(meta_padj_threshold))
fun_selection_rule <- if (use_custom_gene_sets) {
  paste0("custom list supplied via params (", length(custom_fun_cim_genes),
         " symbols)")
} else {
  paste0("cim == Fun | scope ", meta_scope,
         " | rho_threshold == ", format_cutoff(meta_rho_threshold),
         " | ", meta_k_mode, " | ", meta_significance_rule)
}
dys_selection_rule <- if (use_custom_gene_sets) {
  paste0("custom list supplied via params (", length(custom_dys_cim_genes),
         " symbols)")
} else {
  paste0("cim == Dys | scope ", meta_scope,
         " | rho_threshold == ", format_cutoff(meta_rho_threshold),
         " | ", meta_k_mode, " | ", meta_significance_rule)
}
# v2.0 names kept as aliases so no downstream chunk needs editing.
fun_correlation_rule <- fun_selection_rule
dys_correlation_rule <- dys_selection_rule
correlation_fdr_rule <- meta_significance_rule

allowed_parameterizations <- c("continuous", "median_split")
score_parameterizations <- unique(as.character(unlist(
  params$score_parameterizations, recursive = TRUE, use.names = FALSE)))
if (!length(score_parameterizations) ||
    any(!score_parameterizations %in% allowed_parameterizations)) {
  stop("score_parameterizations must contain continuous and/or median_split.")
}
median_tie_rule <- validate_choice(
  params$median_tie_rule,
  c("high_includes_cutoff", "low_includes_cutoff"),
  "params$median_tie_rule")
median_group_warning_fraction <- validate_numeric(
  params$median_group_warning_fraction,
  "params$median_group_warning_fraction", 0, 0.5)
significance_metric <- validate_choice(
  params$significance_metric, c("p_value", "p_adjust_signature"),
  "params$significance_metric")

minimum_signature_coverage <- validate_numeric(
  params$minimum_signature_coverage,
  "params$minimum_signature_coverage", 0, 1, lower_inclusive = FALSE)
minimum_signature_genes <- validate_integer(
  params$minimum_signature_genes, "params$minimum_signature_genes")
duplicate_gene_method <- validate_choice(
  params$duplicate_gene_method, c("max_variance", "mean", "median"),
  "params$duplicate_gene_method")
missing_expression_gene_policy <- validate_choice(
  params$missing_expression_gene_policy, c("drop_gene", "stop"),
  "params$missing_expression_gene_policy")
ssgsea_alpha <- validate_numeric(
  params$ssgsea_alpha, "params$ssgsea_alpha", 0, Inf,
  lower_inclusive = FALSE)
ssgsea_normalize <- validate_flag(
  params$ssgsea_normalize, "params$ssgsea_normalize")
ssgsea_min_size <- validate_integer(
  params$ssgsea_min_size, "params$ssgsea_min_size")
ssgsea_max_size <- if (is.null(params$ssgsea_max_size) ||
                       length(params$ssgsea_max_size) == 0L ||
                       all(is.na(params$ssgsea_max_size))) {
  Inf
} else {
  validate_integer(params$ssgsea_max_size, "params$ssgsea_max_size")
}
if (ssgsea_max_size < ssgsea_min_size) {
  stop("ssgsea_max_size must be null/Inf or at least ssgsea_min_size.")
}
bioc_parallel_workers <- validate_integer(
  params$bioc_parallel_workers, "params$bioc_parallel_workers")
ssgsea_bpparam <- if (bioc_parallel_workers == 1L) {
  BiocParallel::SerialParam()
} else {
  BiocParallel::SnowParam(workers = bioc_parallel_workers, type = "SOCK")
}

expected_fun_cim_genes <- validate_optional_integer(
  params$expected_fun_cim_genes, "params$expected_fun_cim_genes")
expected_dys_cim_genes <- validate_optional_integer(
  params$expected_dys_cim_genes, "params$expected_dys_cim_genes")
metabric_chemo_only <- validate_flag(
  params$metabric_chemo_only, "params$metabric_chemo_only")
scanb_chemo_only <- validate_flag(
  params$scanb_chemo_only, "params$scanb_chemo_only")
enforce_events_per_parameter <- validate_flag(
  params$enforce_events_per_parameter,
  "params$enforce_events_per_parameter")
run_common_sample_sensitivity <- validate_flag(
  params$run_common_sample_sensitivity,
  "params$run_common_sample_sensitivity")
run_km_median_split <- validate_flag(
  params$run_km_median_split, "params$run_km_median_split")
make_forest_plots <- validate_flag(
  params$make_forest_plots, "params$make_forest_plots")

is_absolute_path <- function(x) {
  grepl("^[A-Za-z]:|^/|^\\\\\\\\", x)
}
data_dir <- normalizePath(params$data_dir, winslash = "/", mustWork = FALSE)
# A repo-relative data_dir is resolved against the working directory, so catch a
# wrong wd here with one clear message rather than five missing-file paths.
if (!dir.exists(data_dir)) {
  stop("data_dir does not exist: ", data_dir,
       "\n  params$data_dir = ", params$data_dir,
       "\n  working directory = ", getwd(),
       if (!is_absolute_path(params$data_dir))
         paste0("\n  data_dir is relative, so run this from the repository ",
                "root (the folder containing oncoimmunology_paper/), or set ",
                "params$data_dir to an absolute path.")
       else "")
}
input_path <- function(x) {
  if (is_absolute_path(x)) x else file.path(data_dir, x)
}
files <- c(
  # v3.0: this is meta_results.csv, not a single-cohort correlation table. The
  # name `correlation` is retained so downstream chunks need no edits.
  correlation = input_path(params$signature_file),
  metabric_clinical = input_path(params$metabric_clinical_file),
  metabric_expression = input_path(params$metabric_expression_file),
  scanb_clinical = input_path(params$scanb_clinical_file),
  scanb_expression = input_path(params$scanb_expression_file)
)
# meta_results.csv is only required when the gene sets are actually derived from
# it; a fully custom gene-set run does not need the file at all.
required_file_names <- if (use_custom_gene_sets) {
  base::setdiff(names(files), "correlation")
} else {
  names(files)
}
missing_files <- required_file_names[!file.exists(files[required_file_names])]
if (length(missing_files)) {
  stop("Missing input file(s): ",
       paste(missing_files, files[missing_files],
             sep = " = ", collapse = "; "),
       if (!use_custom_gene_sets && "correlation" %in% missing_files)
         paste0("\n  `correlation` must point at meta_results.csv from ",
                "cimic_overlap_meta_rankbiserial.R (repo location: ",
                "oncoimmunology_paper/Results/cimic_meta_analysis/).")
       else "")
}
expected_md5 <- params$expected_signature_file_md5
if (!use_custom_gene_sets &&
    !is.null(expected_md5) && length(expected_md5) &&
    !all(is.na(expected_md5)) && nzchar(trimws(as.character(expected_md5)))) {
  observed_md5 <- unname(tools::md5sum(files[["correlation"]]))
  if (!identical(tolower(observed_md5),
                 tolower(trimws(as.character(expected_md5))))) {
    stop("The signature-file MD5 does not match expected_signature_file_md5.")
  }
}

output_dir <- if (is_absolute_path(params$output_dir)) {
  params$output_dir
} else {
  file.path(getwd(), params$output_dir)
}
# A trailing rhoNNN (or the v2.0 corrNNN) in the folder name must agree with the
# active rho, so variant runs cannot silently overwrite each other's folders.
output_rho_match <- regexec("(?:rho|corr)([0-9]{3})$", basename(output_dir),
                            ignore.case = TRUE)
output_rho_parts <- regmatches(basename(output_dir), output_rho_match)[[1]]
if (length(output_rho_parts) == 2L) {
  output_rho_threshold <- as.numeric(output_rho_parts[2]) / 100
  if (!isTRUE(all.equal(output_rho_threshold, meta_rho_threshold,
                        tolerance = 1e-12))) {
    stop("output_dir suffix ", output_rho_parts[1], " represents ",
         format_cutoff(output_rho_threshold),
         ", but meta_rho_threshold is ",
         format_cutoff(meta_rho_threshold),
         ". Rename output_dir or change meta_rho_threshold.")
  }
}
subdirs <- setNames(file.path(output_dir, c(
  "01_QC", "02_Scores", "03_Models", "04_Forest_Plots",
  "05_Consolidated_Figures", "06_Sensitivity", "07_Reproducibility",
  "08_KM_Plots")),
  c("qc", "scores", "models", "forests", "figures", "sensitivity",
    "repro", "km"))
invisible(lapply(c(output_dir, unname(subdirs)), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

endpoint_config <- list(
  METABRIC = tribble(
    ~endpoint_original, ~endpoint_harmonized, ~time_col,     ~event_col,
    "OS",               "OS",                 "OS_MONTHS",   "OS_EVENT",
    "DSS",              "DSS",                "DSS_MONTHS",  "DSS_EVENT",
    "RFI",              "RFI",                "RFI_MONTHS",  "RFI_EVENT"
  ),
  SCANB = tribble(
    ~endpoint_original, ~endpoint_harmonized, ~time_col,      ~event_col,
    "OS",               "OS",                 "OS_MONTHS",    "OS_EVENT",
    "DRFi",             "DRFi",               "DRFI_MONTHS",  "DRFI_EVENT",
    "RFi",              "RFI",                "RFI_MONTHS",   "RFI_EVENT"
  )
)
# ── Model definitions: M0 (unadjusted) and M2 (primary adjusted) only ──────
allowed_model_covariates <- c(
  "age_5y", "tumor_size", "grade", "pam50", "clinical_subtype")
m2_covariates <- unique(as.character(unlist(
  params$m2_covariates, recursive = TRUE, use.names = FALSE)))
if (!length(m2_covariates) || any(!m2_covariates %in% allowed_model_covariates)) {
  stop("m2_covariates must contain one or more of: ",
       paste(allowed_model_covariates, collapse = ", "), ".")
}
model_definitions <- tribble(
  ~model_id, ~model_label,                          ~covariates,              ~model_role,             ~color,
  "M0",      "M0: Unadjusted",                      list(character(0)),        "Primary univariable",   "#4DAF4A",
  "M2",      "M2: Prespecified adjusted",            list(m2_covariates),       "Primary adjusted",       "#E41A1C"
)
model_levels <- model_definitions$model_id
model_colors <- setNames(model_definitions$color, model_definitions$model_id)
primary_adjusted_model <- as.character(params$primary_adjusted_model)
if (length(primary_adjusted_model) != 1L ||
    !primary_adjusted_model %in% model_levels) {
  stop("primary_adjusted_model must identify one configured adjusted model.")
}
primary_covariates <- model_definitions$covariates[
  match(primary_adjusted_model, model_definitions$model_id)][[1]]
if (length(primary_covariates) == 0L) {
  stop("primary_adjusted_model must identify an adjusted, not unadjusted, model.")
}
forest_model <- validate_choice(
  params$forest_model, model_levels, "params$forest_model")
forest_parameterization <- validate_choice(
  params$forest_parameterization, allowed_parameterizations,
  "params$forest_parameterization")
if (make_forest_plots &&
    !forest_parameterization %in% score_parameterizations) {
  stop("forest_parameterization must also be included in score_parameterizations.")
}
hr_axis_mode <- validate_choice(
  params$hr_axis_mode, c("manual", "automatic"), "params$hr_axis_mode")
hr_axis_min <- validate_numeric(
  params$hr_axis_min, "params$hr_axis_min", 0, Inf,
  lower_inclusive = FALSE)
hr_axis_max <- validate_numeric(
  params$hr_axis_max, "params$hr_axis_max", 0, Inf,
  lower_inclusive = FALSE)
if (hr_axis_max <= hr_axis_min) {
  stop("hr_axis_max must be greater than hr_axis_min.")
}
manual_hr_breaks <- sort(unique(suppressWarnings(as.numeric(unlist(
  params$hr_axis_breaks, recursive = TRUE, use.names = FALSE)))))
if (!length(manual_hr_breaks) || any(!is.finite(manual_hr_breaks)) ||
    any(manual_hr_breaks <= 0)) {
  stop("hr_axis_breaks must contain positive finite numbers.")
}
forest_width <- validate_numeric(
  params$forest_width, "params$forest_width", 0, Inf,
  lower_inclusive = FALSE)
forest_height <- validate_numeric(
  params$forest_height, "params$forest_height", 0, Inf,
  lower_inclusive = FALSE)
combined_forest_width <- validate_numeric(
  params$combined_forest_width, "params$combined_forest_width", 0, Inf,
  lower_inclusive = FALSE)
score_qc_plot_width <- validate_numeric(
  params$score_qc_plot_width, "params$score_qc_plot_width", 0, Inf,
  lower_inclusive = FALSE)
score_qc_plot_height <- validate_numeric(
  params$score_qc_plot_height, "params$score_qc_plot_height", 0, Inf,
  lower_inclusive = FALSE)
figure_dpi <- validate_integer(params$figure_dpi, "params$figure_dpi")
km_break_months <- validate_numeric(
  params$km_break_months, "params$km_break_months", 0, Inf,
  lower_inclusive = FALSE)
km_width <- validate_numeric(params$km_width, "params$km_width", 0, Inf,
                             lower_inclusive = FALSE)
km_height <- validate_numeric(params$km_height, "params$km_height", 0, Inf,
                              lower_inclusive = FALSE)
km_dpi <- validate_integer(params$km_dpi, "params$km_dpi")
metabric_km_ylim <- validate_numeric_vector(
  params$metabric_km_ylim, "params$metabric_km_ylim", 2L, 0, 1)
scanb_km_ylim <- validate_numeric_vector(
  params$scanb_km_ylim, "params$scanb_km_ylim", 2L, 0, 1)
if (metabric_km_ylim[1] >= metabric_km_ylim[2] ||
    scanb_km_ylim[1] >= scanb_km_ylim[2]) {
  stop("Each KM y-limit vector must be strictly increasing.")
}
metabric_km_pval_coord <- validate_numeric_vector(
  params$metabric_km_pval_coord, "params$metabric_km_pval_coord", 2L,
  lower = 0)
scanb_km_pval_coord <- validate_numeric_vector(
  params$scanb_km_pval_coord, "params$scanb_km_pval_coord", 2L,
  lower = 0)
if (metabric_km_pval_coord[2] < metabric_km_ylim[1] ||
    metabric_km_pval_coord[2] > metabric_km_ylim[2] ||
    scanb_km_pval_coord[2] < scanb_km_ylim[1] ||
    scanb_km_pval_coord[2] > scanb_km_ylim[2]) {
  stop("Each KM p-value y coordinate must fall within its cohort y limits.")
}
min_n_model <- validate_integer(params$min_n_model, "params$min_n_model")
min_events_model <- validate_integer(
  params$min_events_model, "params$min_events_model")
events_per_parameter_target <- validate_numeric(
  params$events_per_parameter, "params$events_per_parameter", 0, Inf,
  lower_inclusive = FALSE)
sparse_factor_level_n <- validate_integer(
  params$sparse_factor_level_n, "params$sparse_factor_level_n")
alpha_threshold <- validate_numeric(params$alpha, "params$alpha", 0, 1)
chemo_positive_values <- unique(str_to_lower(str_squish(as.character(unlist(
  params$chemo_positive_values, recursive = TRUE, use.names = FALSE)))))
if (!length(chemo_positive_values) || any(!nzchar(chemo_positive_values))) {
  stop("chemo_positive_values must contain at least one nonempty value.")
}
metabric_gene_start_col <- validate_optional_integer(
  params$metabric_gene_start_col, "params$metabric_gene_start_col", lower = 1L)
if (is.null(metabric_gene_start_col)) {
  stop("metabric_gene_start_col cannot be null: ssGSEA requires the full ",
       "transcriptome-wide expression background, which cannot be inferred ",
       "safely from signature overlap alone.")
}
covariate_dictionary <- tribble(
  ~variable,            ~display_base,             ~unit_or_reference,
  "cim_z",              "CIM score",               "per 1 SD",
  "median_group",       "CIM median group",        "High versus Low",
  "age_5y",             "Age",                     "per 5 years",
  "tumor_size",         "Tumor size",              "reference: ≤20 mm",
  "grade",              "Histologic grade",        "reference: 1-2",
  "pam50",              "PAM50 subtype",           "reference: LumA",
  "clinical_subtype",   "Clinical subtype",        "reference: HR+/HER2-"
)
input_manifest <- tibble(
  input = names(files), path = unname(files),
  size_bytes = unname(file.info(files)$size),
  modified = format(unname(file.info(files)$mtime), "%Y-%m-%d %H:%M:%S")
)
readr::write_csv(input_manifest,
  file.path(subdirs["qc"], "input_manifest.csv"))
readr::write_csv(
  model_definitions %>% mutate(covariates = map_chr(covariates,
    ~ paste(as.character(unlist(.x, recursive = TRUE)), collapse = " + "))),
  file.path(subdirs["qc"], "model_definitions.csv"))
readr::write_csv(covariate_dictionary,
  file.path(subdirs["qc"], "harmonized_covariate_dictionary.csv"))
or_null <- function(x) if (is.null(x)) "null" else as.character(x)
gene_selection_parameters <- tibble(
  parameter = c(
    "gene_set_source", "signature_file",
    "meta_scope", "meta_rho_threshold", "meta_k_mode",
    "meta_k_mode_description",
    "meta_require_fisher", "meta_require_stouffer", "meta_padj_threshold",
    "effective_meta_significance_rule",
    "meta_min_mean_abs_rrb", "meta_top_n_genes", "meta_rank_metric",
    "meta_ambiguous_gene_policy",
    "effective_fun_cim_rule", "effective_dys_cim_rule",
    "expected_fun_cim_genes", "expected_dys_cim_genes",
    "cox_fdr_threshold", "cox_fdr_operator", "output_dir", "parameter_source"
  ),
  value = c(
    if (use_custom_gene_sets) "custom gene lists from params"
      else "cimic meta_results.csv",
    as.character(params$signature_file),
    meta_scope, format_cutoff(meta_rho_threshold), meta_k_mode,
    meta_k_mode_description,
    as.character(meta_require_fisher), as.character(meta_require_stouffer),
    format_cutoff(meta_padj_threshold),
    meta_significance_rule,
    or_null(meta_min_mean_abs_rrb), or_null(meta_top_n_genes),
    meta_rank_metric, meta_ambiguous_gene_policy,
    fun_selection_rule, dys_selection_rule,
    or_null(expected_fun_cim_genes), or_null(expected_dys_cim_genes),
    format_cutoff(cox_fdr_threshold),
    cox_fdr_operator, output_dir, analysis_parameter_source
  )
)
readr::write_csv(
  gene_selection_parameters,
  file.path(subdirs["qc"], "gene_selection_parameters.csv")
)
cat(
  "Gene-selection settings parsed from YAML:\n",
  "  Source:      ",
  if (use_custom_gene_sets) "custom gene lists" else "meta_results.csv", "\n",
  "  Scope:       ", meta_scope, "\n",
  "  rho:         ", format_cutoff(meta_rho_threshold), "\n",
  "  k-mode:      ", meta_k_mode, "  [", meta_k_mode_description, "]\n",
  "  Significance:", meta_significance_rule, "\n",
  "  Fun-CIM:     ", fun_selection_rule, "\n",
  "  Dys-CIM:     ", dys_selection_rule, "\n",
  "  Output dir:  ", output_dir, "\n",
  sep = ""
)


# ==============================================================================
# GENERAL AND HARMONIZATION HELPERS
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION: helpers
# ------------------------------------------------------------------------------
# Core data-harmonization helpers
safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))
clean_id <- function(x) stringr::str_squish(as.character(x))
collapse_values <- function(x) paste(sort(unique(na.omit(as.character(x)))), collapse = "; ")
subset_rows <- function(x, rows) x[rows, seq_len(ncol(x)), drop = FALSE]
subset_columns <- function(x, columns) x[seq_len(nrow(x)), columns, drop = FALSE]
extract_column <- function(x, column) x[seq_len(nrow(x)), column, drop = TRUE]
get_first_available <- function(df, candidates, numeric = FALSE) {
  hits <- candidates[candidates %in% names(df)]
  if (!length(hits)) {
    value <- if (numeric) rep(NA_real_, nrow(df)) else rep(NA_character_, nrow(df))
    return(list(value = value, source = NA_character_))
  }
  usable <- vapply(hits, function(candidate) {
    value <- df[[candidate]]
    if (numeric) value <- safe_numeric(value)
    any(!is.na(value) & str_squish(as.character(value)) != "")
  }, logical(1))
  hit <- if (any(usable)) hits[which(usable)[1]] else hits[1]
  value <- df[[hit]]
  if (numeric) value <- safe_numeric(value)
  list(value = value, source = hit)
}
standardize_event <- function(x, variable, cohort) {
  if (is.logical(x)) return(as.integer(x))
  z <- str_to_lower(str_squish(as.character(x)))
  map <- c("0" = 0, "1" = 1, "false" = 0, "true" = 1, "no" = 0,
           "yes" = 1, "alive" = 0, "living" = 0, "censored" = 0,
           "dead" = 1, "deceased" = 1, "event" = 1,
           "0:living" = 0, "1:deceased" = 1,
           "0:not recurred" = 0, "1:recurred" = 1,
           "no recurrence" = 0, "recurrence" = 1,
           "died of disease" = 1, "died of other causes" = 0)
  bad <- setdiff(unique(z[!is.na(z) & z != ""]), names(map))
  if (length(bad)) stop(cohort, " ", variable, " has unmapped event values: ",
                        paste(bad, collapse = ", "))
  unname(map[z])
}
clean_tumor_size <- function(x) {
  z <- str_to_lower(str_squish(as.character(x))) %>%
    str_replace_all("[−–—]", "-")
  numeric_x <- safe_numeric(x)
  case_when(
    is.finite(numeric_x) & numeric_x <= 20 ~ "≤20 mm",
    is.finite(numeric_x) & numeric_x > 20 & numeric_x <= 50 ~ "21-50 mm",
    is.finite(numeric_x) & numeric_x > 50 ~ ">50 mm",
    z %in% c("≤20 mm", "≤ 20 mm", "<=20 mm", "<= 20 mm", "0-20", "0-20 mm") ~ "≤20 mm",
    z %in% c("21-50", "21-50 mm", "21 - 50 mm", ">20-50 mm") ~ "21-50 mm",
    z %in% c(">50", ">50 mm", "> 50 mm", "50+", "51+ mm") ~ ">50 mm",
    is.na(z) | z %in% c("", "na", "n/a", "unknown") ~ NA_character_,
    TRUE ~ paste0("UNMAPPED:", x)
  )
}
clean_grade <- function(x) {
  z <- str_to_lower(str_squish(as.character(x)))
  case_when(
    z %in% c("1", "2", "1-2", "1 & 2", "grade 1", "grade 2", "grade 1-2") ~ "1-2",
    z %in% c("3", "grade 3") ~ "3",
    is.na(z) | z %in% c("", "na", "n/a", "unknown") ~ NA_character_,
    TRUE ~ paste0("UNMAPPED:", x)
  )
}
clean_pam50 <- function(x) {
  z <- str_to_lower(str_squish(as.character(x)))
  case_when(
    z %in% c("luma", "luminal a", "luminal_a") ~ "LumA",
    z %in% c("lumb", "luminal b", "luminal_b") ~ "LumB",
    z %in% c("her2", "her2-enriched", "her2 enriched") ~ "Her2",
    z %in% c("basal", "basal-like", "basal like") ~ "Basal",
    z %in% c("normal", "normal-like", "normal like") ~ "Normal",
    is.na(z) | z %in% c("", "na", "n/a", "unknown") ~ NA_character_,
    TRUE ~ paste0("UNMAPPED:", x)
  )
}
clean_clinical_subtype <- function(x) {
  z <- str_to_lower(str_replace_all(str_squish(as.character(x)), "\\s+", ""))
  case_when(
    z %in% c("hr+/her2-", "her2-/hr+", "hrpositive/her2negative") ~ "HR+/HER2-",
    z %in% c("hr+/her2+", "her2+/hr+", "hrpositive/her2positive") ~ "HR+/HER2+",
    z %in% c("hr-/her2+", "her2+/hr-", "hrnegative/her2positive") ~ "HR-/HER2+",
    z %in% c("hr-/her2-", "her2-/hr-", "hrnegative/her2negative", "tnbc", "triple-negative") ~ "HR-/HER2-",
    is.na(z) | z %in% c("", "na", "n/a", "unknown") ~ NA_character_,
    TRUE ~ paste0("UNMAPPED:", x)
  )
}
strict_factor <- function(x, cleaner, levels, variable, cohort) {
  cleaned <- cleaner(x)
  bad <- unique(cleaned[!is.na(cleaned) & str_starts(cleaned, "UNMAPPED:")])
  if (length(bad)) stop(cohort, " ", variable, " has unmapped values: ",
                        paste(sub("^UNMAPPED:", "", bad), collapse = ", "))
  factor(cleaned, levels = levels)
}
harmonize_clinical <- function(df, cohort) {
  id_col <- if (cohort == "METABRIC") "PATIENT_ID" else "GEX.assay"
  if (!id_col %in% names(df)) stop(cohort, " clinical data must contain ", id_col, ".")
  audit <- tibble(cohort = character(), target = character(), source = character(),
                  derivation = character(), nonmissing_n = integer())
  record <- function(target, source, derivation, value) {
    audit <<- bind_rows(audit, tibble(cohort, target,
      source = ifelse(is.na(source), "not_available", source), derivation,
      nonmissing_n = sum(!is.na(value))))
  }
  aliases <- if (cohort == "METABRIC") list(
    OS_MONTHS  = c("OS_MONTHS",   "OS_months"),
    OS_EVENT   = c("OS_EVENT",    "OS_event"),
    DSS_MONTHS = c("DSS_MONTHS",  "DSS_months"),
    DSS_EVENT  = c("DSS_EVENT",   "DSS_event"),
    RFI_MONTHS = c("RFI_MONTHS",  "RFI_months"),
    RFI_EVENT  = c("RFI_EVENT",   "RFI_event")
  ) else list(
    OS_MONTHS   = c("OS_MONTHS",   "OS_months"),
    OS_EVENT    = c("OS_EVENT",    "OS_event"),
    DRFI_MONTHS = c("DRFI_MONTHS", "DRFi_months", "DRFI_months"),
    DRFI_EVENT  = c("DRFI_EVENT",  "DRFi_event",  "DRFI_event"),
    RFI_MONTHS  = c("RFI_MONTHS",  "RFi_months",  "RFI_months"),
    RFI_EVENT   = c("RFI_EVENT",   "RFi_event",   "RFI_event")
  )
  for (target in names(aliases)) {
    found <- get_first_available(df, aliases[[target]],
                                 numeric = str_ends(target, "_MONTHS"))
    df[[target]] <- found$value
    record(target, found$source, "first available endpoint alias", found$value)
  }
  age_candidates <- if (cohort == "METABRIC") {
    c("AGE_AT_DIAGNOSIS", "Age", "age")
  } else c("AGE_mid", "Age (5-year range)", "Age", "age")
  age <- get_first_available(df, age_candidates, numeric = TRUE)
  df$age_years <- age$value
  df$age_5y    <- df$age_years / 5
  record("age_5y", age$source, "numeric age divided by 5", df$age_5y)
  size <- get_first_available(df,
    c("Tumor_Size_Cat",
      if (cohort == "METABRIC") "TUMOR_SIZE" else "Size.mm",
      "TUMOR_SIZE", "Size.mm"))
  df$tumor_size <- strict_factor(size$value, clean_tumor_size,
    c("≤20 mm", "21-50 mm", ">50 mm"), "tumor_size", cohort)
  record("tumor_size", size$source, "harmonized or derived size category", df$tumor_size)
  grade <- get_first_available(df,
    c("NHG_collapsed", "NHG_Cat", "GRADE", "NHG", "Grade"))
  df$grade <- strict_factor(grade$value, clean_grade,
    c("1-2", "3"), "grade", cohort)
  record("grade", grade$source, "grade 1/2 versus grade 3", df$grade)
  pam <- get_first_available(df,
    if (cohort == "METABRIC")
      c("PAM50_metabric", "PAM50", "SSP.PAM50")
    else
      c("PAM50", "SSP.PAM50", "PAM50_metabric"))
  df$pam50 <- strict_factor(pam$value, clean_pam50,
    c("LumA", "LumB", "Her2", "Basal", "Normal"), "pam50", cohort)
  record("pam50", pam$source, "harmonized PAM50", df$pam50)
  subtype <- get_first_available(df,
    c("Clinical_Subtype", "Clinical_Subtype_full", "clinical_subtype"))
  df$clinical_subtype <- strict_factor(subtype$value, clean_clinical_subtype,
    c("HR+/HER2-", "HR+/HER2+", "HR-/HER2+", "HR-/HER2-"),
    "clinical_subtype", cohort)
  record("clinical_subtype", subtype$source,
         "harmonized receptor subtype", df$clinical_subtype)
  df[[id_col]] <- clean_id(df[[id_col]])
  attr(df, "derivation_audit") <- audit
  df
}
resolve_duplicates <- function(mat, method = "max_variance") {
  if (!method %in% c("max_variance", "mean", "median"))
    stop("Unknown duplicate method.")
  symbols <- toupper(trimws(colnames(mat)))
  groups  <- split(seq_along(symbols), symbols)
  resolved <- lapply(groups, function(idx) {
    x <- subset_columns(mat, idx)
    if (ncol(x) == 1L) return(extract_column(x, 1L))
    if (method == "max_variance") {
      variance <- apply(x, 2, var, na.rm = TRUE)
      variance[!is.finite(variance)] <- -Inf
      return(extract_column(x, which.max(variance)))
    }
    apply(x, 1, if (method == "mean") mean else median, na.rm = TRUE)
  })
  out <- do.call(cbind, resolved)
  rownames(out) <- rownames(mat)
  colnames(out) <- names(groups)
  storage.mode(out) <- "double"
  list(matrix = out,
       audit  = tibble(gene = names(groups),
                       n_input_columns = lengths(groups),
                       method = method))
}
term_display <- function(term) {
  case_when(
    term == "cim_z"                          ~ "CIM score (per 1 SD)",
    term == "median_groupHigh"               ~ "CIM score: High vs Low",
    term == "age_5y"                         ~ "Age (per 5 years)",
    term == "tumor_size21-50 mm"             ~ "Tumor size: 21-50 vs ≤20 mm",
    term == "tumor_size>50 mm"               ~ "Tumor size: >50 vs ≤20 mm",
    term == "grade3"                         ~ "Grade: 3 vs 1-2",
    term == "pam50LumB"                      ~ "PAM50: LumB vs LumA",
    term == "pam50Her2"                      ~ "PAM50: HER2-enriched vs LumA",
    term == "pam50Basal"                     ~ "PAM50: Basal-like vs LumA",
    term == "pam50Normal"                    ~ "PAM50: Normal-like vs LumA",
    term == "clinical_subtypeHR+/HER2+"      ~ "HR+/HER2+ vs HR+/HER2-",
    term == "clinical_subtypeHR-/HER2+"      ~ "HR-/HER2+ vs HR+/HER2-",
    term == "clinical_subtypeHR-/HER2-"      ~ "HR-/HER2- vs HR+/HER2-",
    TRUE ~ term
  )
}
term_variable <- function(term) case_when(
  term == "cim_z"                    ~ "cim_z",
  term == "median_groupHigh"         ~ "median_group",
  term == "age_5y"                   ~ "age_5y",
  str_starts(term, "tumor_size")     ~ "tumor_size",
  str_starts(term, "grade")          ~ "grade",
  str_starts(term, "pam50")          ~ "pam50",
  str_starts(term, "clinical_subtype") ~ "clinical_subtype",
  TRUE ~ NA_character_
)
term_reference <- function(term) case_when(
  term == "cim_z"                    ~ "per 1-SD increase",
  term == "median_groupHigh"         ~ "Low",
  term == "age_5y"                   ~ "per 5-year increase",
  str_starts(term, "tumor_size")     ~ "≤20 mm",
  str_starts(term, "grade")          ~ "1-2",
  str_starts(term, "pam50")          ~ "LumA",
  str_starts(term, "clinical_subtype") ~ "HR+/HER2-",
  TRUE ~ NA_character_
)
term_comparison <- function(term) case_when(
  term == "cim_z"                    ~ "+1 SD",
  term == "median_groupHigh"         ~ "High",
  term == "age_5y"                   ~ "+5 years",
  str_starts(term, "tumor_size")     ~ str_remove(term, "^tumor_size"),
  str_starts(term, "grade")          ~ str_remove(term, "^grade"),
  str_starts(term, "pam50")          ~ str_remove(term, "^pam50"),
  str_starts(term, "clinical_subtype") ~ str_remove(term, "^clinical_subtype"),
  TRUE ~ NA_character_
)


# ==============================================================================
# CONSTRUCT FUN-CIM AND DYS-CIM GENE SETS
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION: gene-sets
# ------------------------------------------------------------------------------
# ── v3.0 signature construction from the CIMIC meta-analysis ────────────────
# Produces exactly the three objects v2.0 produced, so nothing downstream
# changes:
#   gene_selection_audit  every candidate row with its pass/fail flags
#   gene_definitions      the selected genes, columns `gene` and `signature`
#   gene_sets             list("Fun-CIM" = <chr>, "Dys-CIM" = <chr>)
if (use_custom_gene_sets) {

  # ---- Manual override: skip meta selection entirely -----------------------
  gene_selection_audit <- dplyr::bind_rows(
    tibble(gene = custom_fun_cim_genes, signature = "Fun-CIM"),
    tibble(gene = custom_dys_cim_genes, signature = "Dys-CIM")
  ) %>%
    dplyr::mutate(
      selection_source = "custom params list",
      # The columns the shared code path below references must exist even in the
      # custom branch, otherwise the ranking arrange() and the workbook sheet
      # fail on a missing variable.
      meta_scope_used = NA_character_, rho_threshold = NA_real_,
      k_datasets = NA_integer_, n_scope_datasets = NA_integer_,
      k_concordant = NA_integer_, direction_concordant = NA_character_,
      fisher_padj = NA_real_, stouffer_padj = NA_real_,
      mean_abs_rrb = NA_real_,
      r_rb = NA_real_, source_p_adjust = NA_real_,
      selected = !is.na(gene) & gene != ""
    )

} else {

  # ---- Selection from meta_results.csv ------------------------------------
  meta_raw <- readr::read_csv(files[["correlation"]], show_col_types = FALSE,
                              name_repair = "minimal")
  required_meta <- c("scope", "cim", "rho_threshold", "gene", "k_datasets",
                     "n_scope_datasets", "direction_concordant",
                     "fisher_padj", "stouffer_padj", "mean_abs_rrb")
  missing_meta <- base::setdiff(required_meta, names(meta_raw))
  if (length(missing_meta))
    stop("meta_results.csv is missing column(s): ",
         paste(missing_meta, collapse = ", "),
         "\n  Re-run cimic_overlap_meta_rankbiserial.R; `n_scope_datasets`, ",
         "`k_concordant` and `direction_concordant` were added when the ",
         "k-mode toggle was introduced.")

  if (!meta_scope %in% unique(meta_raw$scope))
    stop("meta_scope '", meta_scope, "' is not present in meta_results.csv. ",
         "Available: ", paste(unique(meta_raw$scope), collapse = ", "))

  available_rho <- sort(unique(meta_raw$rho_threshold))
  if (!any(abs(available_rho - meta_rho_threshold) < 1e-8))
    stop("meta_rho_threshold ", format_cutoff(meta_rho_threshold),
         " is not present in meta_results.csv. Available: ",
         paste(available_rho, collapse = ", "),
         ". Add the value to RHO_CUTS in cimic_overlap_meta_rankbiserial.R ",
         "and re-run it, or choose an available cutoff.")

  # Scalar floor precomputed: `mean_abs_rrb >= NULL` would collapse to
  # logical(0) inside transmute rather than recycling.
  meta_effect_floor <- if (is.null(meta_min_mean_abs_rrb)) -Inf else
    meta_min_mean_abs_rrb

  # Audit table: every row for this scope+rho, with each rule evaluated
  # separately so a QC reader can see exactly which filter removed what.
  gene_selection_audit <- meta_raw %>%
    dplyr::filter(scope == meta_scope,
                  abs(rho_threshold - meta_rho_threshold) < 1e-8) %>%
    dplyr::transmute(
      gene = toupper(trimws(as.character(gene))),
      signature = dplyr::case_when(cim == "Fun" ~ "Fun-CIM",
                                   cim == "Dys" ~ "Dys-CIM",
                                   TRUE ~ "NS"),
      meta_scope_used   = scope,
      rho_threshold     = rho_threshold,
      k_datasets        = k_datasets,
      n_scope_datasets  = n_scope_datasets,
      k_concordant      = if ("k_concordant" %in% names(meta_raw))
                            k_concordant else NA_integer_,
      direction_concordant = direction_concordant,
      fisher_padj       = safe_numeric(fisher_padj),
      stouffer_padj     = safe_numeric(stouffer_padj),
      mean_abs_rrb      = safe_numeric(mean_abs_rrb),
      # `r_rb` / `source_p_adjust` retain v2.0's column names so the QC
      # verification block and the workbook sheet keep working unchanged.
      r_rb              = mean_abs_rrb,
      # Scalar `if`, not dplyr::if_else: the condition is a single setting, and
      # if_else() requires its condition to be the same length as the branches.
      source_p_adjust   = if (meta_require_stouffer) stouffer_padj
                          else fisher_padj,
      passes_fisher     = !meta_require_fisher |
                            (is.finite(fisher_padj) &
                               fisher_padj < meta_padj_threshold),
      passes_stouffer   = !meta_require_stouffer |
                            (is.finite(stouffer_padj) &
                               stouffer_padj < meta_padj_threshold),
      passes_k_mode = switch(
        meta_k_mode,
        k_min2           = k_datasets >= 2,
        k_all            = k_datasets == n_scope_datasets,
        k_all_concordant = k_datasets == n_scope_datasets &
                             direction_concordant == "YES"
      ),
      passes_effect_floor = is.finite(mean_abs_rrb) &
        mean_abs_rrb >= meta_effect_floor,
      selected = !is.na(gene) & gene != "" &
        signature %in% c("Fun-CIM", "Dys-CIM") &
        passes_fisher & passes_stouffer & passes_k_mode & passes_effect_floor,
      selection_source = "cimic meta_results.csv"
    )
}

# ---- Resolve genes that qualify for BOTH signatures -----------------------
# Impossible in v2.0 (a single signed r_rb), but the meta table tests Fun and
# Dys separately, so a gene can pass both. Handle it explicitly rather than
# letting split() silently duplicate the gene across signatures.
selected_rows <- gene_selection_audit %>% dplyr::filter(selected)
rank_col <- if (use_custom_gene_sets) "gene" else meta_rank_metric
ambiguous_genes <- selected_rows %>%
  dplyr::distinct(gene, signature) %>%
  dplyr::count(gene, name = "n_signatures") %>%
  dplyr::filter(n_signatures > 1) %>%
  dplyr::pull(gene)

if (length(ambiguous_genes)) {
  message("Genes meta-significant toward BOTH CIMs: ", length(ambiguous_genes),
          "  (policy: ", meta_ambiguous_gene_policy, ")")
  readr::write_csv(
    selected_rows %>% dplyr::filter(gene %in% ambiguous_genes) %>%
      dplyr::arrange(gene, signature),
    file.path(subdirs["qc"], "gene_set_ambiguous_genes.csv"))
  if (identical(meta_ambiguous_gene_policy, "drop")) {
    selected_rows <- selected_rows %>%
      dplyr::filter(!gene %in% ambiguous_genes)
  }
  # "assign_best" needs no filtering here: the arrange + distinct below keeps
  # only the row with the smallest ranking metric for each gene.
}

gene_definitions <- selected_rows %>%
  dplyr::arrange(gene, .data[[rank_col]], dplyr::desc(mean_abs_rrb)) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

# Optional top-N per signature, applied after de-duplication so N counts
# unique genes.
if (!is.null(meta_top_n_genes)) {
  gene_definitions <- gene_definitions %>%
    dplyr::group_by(signature) %>%
    dplyr::arrange(.data[[rank_col]], dplyr::desc(mean_abs_rrb),
                   .by_group = TRUE) %>%
    dplyr::slice_head(n = meta_top_n_genes) %>%
    dplyr::ungroup()
  message("Applied meta_top_n_genes = ", meta_top_n_genes,
          " per signature.")
}

gene_sets <- split(gene_definitions$gene, gene_definitions$signature)
# A signature can be emptied entirely by the selection (e.g. Dys-CIM is empty for
# NEO.TNBCCL at either rho). split() then drops the name, so report which one
# vanished together with the settings that caused it.
emptied <- base::setdiff(c("Fun-CIM", "Dys-CIM"), names(gene_sets))
if (length(emptied)) {
  stop(
    "No genes selected for: ", paste(emptied, collapse = " and "),
    ". Both signatures are required.\n",
    "  Active selection: ",
    if (use_custom_gene_sets) "custom gene lists" else
      paste0("scope ", meta_scope, " | rho ",
             format_cutoff(meta_rho_threshold), " | ", meta_k_mode, " | ",
             meta_significance_rule), "\n",
    "  Loosen meta_k_mode (k_all_concordant -> k_all -> k_min2), lower ",
    "meta_rho_threshold, or choose a meta_scope where both CIMs have signal ",
    "(Dys-CIM is empty for NEO.TNBCCL at every setting)."
  )
}
observed_signature_counts <- lengths(gene_sets[c("Fun-CIM", "Dys-CIM")])
if (!is.null(expected_fun_cim_genes) &&
    observed_signature_counts[["Fun-CIM"]] != expected_fun_cim_genes) {
  stop("Fun-CIM selected-gene count is ", observed_signature_counts[["Fun-CIM"]],
       ", but expected_fun_cim_genes is ", expected_fun_cim_genes,
       ". Confirm the meta_results.csv version and the active meta_* ",
       "selection settings, or set ",
       "expected_fun_cim_genes to null for an intentional alternative analysis.")
}
if (!is.null(expected_dys_cim_genes) &&
    observed_signature_counts[["Dys-CIM"]] != expected_dys_cim_genes) {
  stop("Dys-CIM selected-gene count is ", observed_signature_counts[["Dys-CIM"]],
       ", but expected_dys_cim_genes is ", expected_dys_cim_genes,
       ". Confirm the meta_results.csv version and the active meta_* ",
       "selection settings, or set ",
       "expected_dys_cim_genes to null for an intentional alternative analysis.")
}
readr::write_csv(gene_selection_audit,
  file.path(subdirs["qc"], "gene_set_selection_audit.csv"))
readr::write_csv(gene_definitions,
  file.path(subdirs["qc"], "gene_set_definitions.csv"))
# Re-check every rule against the FINAL selection, independently of the flags
# computed above. Any non-zero count means a filter did not actually apply.
gene_selection_verification <- if (use_custom_gene_sets) {
  gene_definitions %>%
    dplyr::group_by(signature) %>%
    dplyr::summarise(selected_genes = dplyr::n(),
                     genes_failing_runtime_rule = 0L, .groups = "drop")
} else {
  gene_definitions %>%
    dplyr::group_by(signature) %>%
    dplyr::summarise(
      selected_genes            = dplyr::n(),
      minimum_selected_mean_abs_rrb = min(mean_abs_rrb),
      maximum_selected_fisher_padj  = max(fisher_padj),
      maximum_selected_stouffer_padj = max(stouffer_padj),
      minimum_k_datasets        = min(k_datasets),
      minimum_k_concordant      = suppressWarnings(min(k_concordant)),
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
                    genes_failing_k_mode_rule +
                    genes_failing_effect_floor_rule)
}
readr::write_csv(gene_selection_verification,
  file.path(subdirs["qc"], "gene_set_selection_verification.csv"))
if (any(gene_selection_verification$genes_failing_runtime_rule > 0L))
  stop("Internal gene-selection verification failed; inspect 01_QC outputs.")
# Genes that qualified for both CIMs must never survive into both signatures.
if (anyDuplicated(gene_definitions$gene))
  stop("A gene appears in more than one signature after selection; ",
       "inspect gene_set_ambiguous_genes.csv in 01_QC.")
print(gene_selection_verification)
message("Signature sizes: ",
        paste(names(gene_sets), lengths(gene_sets), sep = " = ",
              collapse = " | "))

# Fail fast, with actionable guidance. The cohort-expression QC further down
# already enforces minimum_signature_genes on MATCHED genes, but that happens
# after both cohorts have been loaded and harmonized. A signature that starts
# below the minimum can never satisfy it, so stop here instead.
undersized <- lengths(gene_sets)[lengths(gene_sets) < minimum_signature_genes]
if (length(undersized)) {
  stop(
    "Signature too small before cohort matching: ",
    paste(names(undersized), undersized, sep = " = ", collapse = "; "),
    " (minimum_signature_genes = ", minimum_signature_genes, ").\n",
    "  Active selection: scope ", meta_scope,
    " | rho ", format_cutoff(meta_rho_threshold),
    " | ", meta_k_mode, " | ", meta_significance_rule, "\n",
    "  Loosen the selection (meta_k_mode k_all_concordant -> k_all -> k_min2, ",
    "or meta_rho_threshold 0.6 -> 0.3), choose a different meta_scope, ",
    "or lower minimum_signature_genes if a small signature is intended.\n",
    "  Note that matched-gene counts can only be lower than these, so this ",
    "run would fail the cohort-expression QC regardless."
  )
}


# ==============================================================================
# READ AND ALIGN COHORTS
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION: read-cohorts
# ------------------------------------------------------------------------------
# Cohort import, restriction, matching, and expression QC
met_expr_raw <- data.table::fread(files[["metabric_expression"]],
                                  data.table = FALSE, check.names = FALSE)
met_clin <- readxl::read_xlsx(files[["metabric_clinical"]]) %>% as.data.frame()
if (!"PATIENT_ID" %in% names(met_expr_raw))
  stop("METABRIC expression lacks PATIENT_ID.")
met_clin  <- harmonize_clinical(met_clin, "METABRIC")
met_audit <- attr(met_clin, "derivation_audit")
if (metabric_gene_start_col > ncol(met_expr_raw)) {
  stop("metabric_gene_start_col exceeds the number of expression columns.")
}
gene_columns <- names(met_expr_raw)[
  seq.int(metabric_gene_start_col, ncol(met_expr_raw))]
metabric_gene_detection_method <- paste0(
  "user-specified first gene-column index: ", metabric_gene_start_col)
if (!length(gene_columns)) {
  stop("No METABRIC gene-expression columns were detected.")
}
if (!any(toupper(trimws(gene_columns)) %in% gene_definitions$gene)) {
  stop("The detected METABRIC gene columns do not overlap the selected signatures.")
}
readr::write_csv(
  tibble(column = gene_columns,
         normalized_symbol = toupper(trimws(gene_columns)),
         selected_signature_gene = normalized_symbol %in% gene_definitions$gene,
         detection_method = metabric_gene_detection_method),
  file.path(subdirs["qc"], "METABRIC_expression_gene_columns.csv"))
met_mat <- as.matrix(data.frame(
  lapply(met_expr_raw[gene_columns], safe_numeric), check.names = FALSE))
rownames(met_mat) <- clean_id(met_expr_raw$PATIENT_ID)
scan_clin <- readRDS(files[["scanb_clinical"]]) %>% as.data.frame()
scan_raw  <- as.matrix(readRDS(files[["scanb_expression"]]))
if (!"GEX.assay" %in% names(scan_clin))
  stop("SCAN-B clinical data lack GEX.assay.")
scan_clin$GEX.assay  <- clean_id(scan_clin$GEX.assay)
row_id_overlap  <- sum(clean_id(rownames(scan_raw)) %in% scan_clin$GEX.assay)
col_id_overlap  <- sum(clean_id(colnames(scan_raw)) %in% scan_clin$GEX.assay)
row_gene_overlap <- sum(toupper(trimws(rownames(scan_raw))) %in% gene_definitions$gene)
col_gene_overlap <- sum(toupper(trimws(colnames(scan_raw))) %in% gene_definitions$gene)
if (row_id_overlap > col_id_overlap && col_gene_overlap >= row_gene_overlap) {
  scan_mat <- scan_raw; scan_orientation <- "samples_by_genes"
} else if (col_id_overlap > row_id_overlap && row_gene_overlap >= col_gene_overlap) {
  scan_mat <- t(scan_raw); scan_orientation <- "genes_by_samples_transposed"
} else stop("SCAN-B expression orientation is ambiguous.")
rownames(scan_mat) <- clean_id(rownames(scan_mat))
storage.mode(scan_mat) <- "double"
filter_chemotherapy <- function(df, cohort, requested_column = NULL) {
  candidate_columns <- if (!is.null(requested_column) &&
                           length(requested_column) &&
                           !all(is.na(requested_column)) &&
                           nzchar(trimws(as.character(requested_column)))) {
    trimws(as.character(requested_column))
  } else {
    c("Chemo_Combined", "CHEMOTHERAPY", "Chemotherapy", "chemotherapy")
  }
  chemo_column <- candidate_columns[candidate_columns %in% names(df)][1]
  if (is.na(chemo_column)) {
    stop(cohort, " chemotherapy restriction was requested, but none of these ",
         "columns was found: ", paste(candidate_columns, collapse = ", "), ".")
  }
  normalized <- str_to_lower(str_squish(as.character(df[[chemo_column]])))
  keep <- !is.na(normalized) & normalized %in% chemo_positive_values
  audit <- tibble(
    cohort, chemo_column,
    input_n = nrow(df), retained_n = sum(keep), excluded_n = sum(!keep),
    positive_values = paste(chemo_positive_values, collapse = "; "))
  if (!any(keep)) {
    stop(cohort, " chemotherapy filtering retained zero patients. Values found: ",
         paste(sort(unique(normalized[!is.na(normalized)])), collapse = ", "), ".")
  }
  list(data = df[keep, , drop = FALSE], audit = audit)
}
chemo_audits <- list()
if (scanb_chemo_only) {
  scan_chemo <- filter_chemotherapy(
    scan_clin, "SCANB", params$scanb_chemo_column)
  scan_clin <- scan_chemo$data
  chemo_audits[["SCANB"]] <- scan_chemo$audit
}
scan_clin  <- harmonize_clinical(scan_clin, "SCANB")
scan_audit <- attr(scan_clin, "derivation_audit")
if (metabric_chemo_only) {
  met_chemo <- filter_chemotherapy(
    met_clin, "METABRIC", params$metabric_chemo_column)
  met_clin <- met_chemo$data
  chemo_audits[["METABRIC"]] <- met_chemo$audit
}
if (length(chemo_audits)) {
  readr::write_csv(bind_rows(chemo_audits),
    file.path(subdirs["qc"], "chemotherapy_filter_audit.csv"))
}
if (anyDuplicated(rownames(met_mat))  || anyDuplicated(met_clin$PATIENT_ID) ||
    anyDuplicated(rownames(scan_mat)) || anyDuplicated(scan_clin$GEX.assay))
  stop("Duplicate sample identifiers detected; resolve them before analysis.")
met_common  <- intersect(met_clin$PATIENT_ID, rownames(met_mat))
scan_common <- intersect(scan_clin$GEX.assay,  rownames(scan_mat))
if (!length(met_common)) stop("No METABRIC clinical-expression sample matches were found.")
if (!length(scan_common)) stop("No SCAN-B clinical-expression sample matches were found.")
met_clin  <- subset_rows(met_clin,  match(met_common,  met_clin$PATIENT_ID))
met_mat   <- subset_rows(met_mat,   met_common)
scan_clin <- subset_rows(scan_clin, match(scan_common, scan_clin$GEX.assay))
scan_mat  <- subset_rows(scan_mat,  scan_common)
met_resolved  <- resolve_duplicates(met_mat,  duplicate_gene_method)
scan_resolved <- resolve_duplicates(scan_mat, duplicate_gene_method)
met_mat  <- met_resolved$matrix
scan_mat <- scan_resolved$matrix
expression_gene_qc <- function(mat, cohort) {
  qc <- tibble(
    cohort,
    gene = colnames(mat),
    nonfinite_n = vapply(seq_len(ncol(mat)), function(j)
      sum(!is.finite(mat[, j])), integer(1)),
    variance = vapply(seq_len(ncol(mat)), function(j)
      var(mat[, j], na.rm = TRUE), numeric(1))) %>%
    mutate(
      finite_complete = nonfinite_n == 0L,
      variable = is.finite(variance) & variance > 0,
      retained = finite_complete & variable,
      exclusion_reason = case_when(
        retained ~ "retained",
        !finite_complete ~ "nonfinite_expression",
        !variable ~ "zero_or_undefined_variance",
        TRUE ~ "excluded"))
  if (missing_expression_gene_policy == "stop" && any(!qc$finite_complete)) {
    stop(cohort, " expression contains nonfinite values. Inspect expression_gene_QC.csv ",
         "or use missing_expression_gene_policy: drop_gene.")
  }
  list(matrix = subset_columns(mat, qc$retained), audit = qc)
}
met_expression_qc <- expression_gene_qc(met_mat, "METABRIC")
scan_expression_qc <- expression_gene_qc(scan_mat, "SCANB")
met_mat <- met_expression_qc$matrix
scan_mat <- scan_expression_qc$matrix
readr::write_csv(bind_rows(met_expression_qc$audit, scan_expression_qc$audit),
  file.path(subdirs["qc"], "expression_gene_QC.csv"))
cohorts <- list(
  METABRIC = list(clinical = met_clin, expression = met_mat,
                  id_col = "PATIENT_ID", endpoints = endpoint_config$METABRIC),
  SCANB    = list(clinical = scan_clin, expression = scan_mat,
                  id_col = "GEX.assay",  endpoints = endpoint_config$SCANB)
)
endpoint_input_audit <- imap_dfr(cohorts, function(obj, cohort) {
  map_dfr(seq_len(nrow(obj$endpoints)), function(i) {
    endpoint_row <- obj$endpoints[i, , drop = FALSE]
    time <- safe_numeric(obj$clinical[[endpoint_row$time_col]])
    event <- standardize_event(obj$clinical[[endpoint_row$event_col]],
                               endpoint_row$event_col, cohort)
    valid <- is.finite(time) & time > 0 & event %in% c(0, 1)
    tibble(
      cohort,
      endpoint_original = endpoint_row$endpoint_original,
      time_column = endpoint_row$time_col,
      event_column = endpoint_row$event_col,
      cohort_n = nrow(obj$clinical),
      valid_endpoint_n = sum(valid),
      events = sum(event[valid] == 1),
      censored = sum(event[valid] == 0))
  })
})
readr::write_csv(endpoint_input_audit,
  file.path(subdirs["qc"], "endpoint_input_audit.csv"))
if (any(endpoint_input_audit$valid_endpoint_n < min_n_model |
        endpoint_input_audit$events < min_events_model)) {
  stop("At least one configured endpoint has insufficient valid observations or events; ",
       "inspect endpoint_input_audit.csv.")
}
required_covariate_input_audit <- imap_dfr(cohorts, function(obj, cohort) {
  map_dfr(m2_covariates, function(variable) {
    x <- obj$clinical[[variable]]
    tibble(
      cohort, variable,
      cohort_n = nrow(obj$clinical),
      nonmissing_n = sum(!is.na(x)),
      distinct_nonmissing = dplyr::n_distinct(x, na.rm = TRUE),
      missing_n = sum(is.na(x)))
  })
})
readr::write_csv(required_covariate_input_audit,
  file.path(subdirs["qc"], "required_covariate_input_audit.csv"))
if (any(required_covariate_input_audit$nonmissing_n < min_n_model |
        required_covariate_input_audit$distinct_nonmissing < 2L)) {
  stop("At least one M2 covariate is unavailable or nonvariable; inspect ",
       "required_covariate_input_audit.csv.")
}
readr::write_csv(bind_rows(met_audit, scan_audit),
  file.path(subdirs["qc"], "clinical_variable_derivation_audit.csv"))
readr::write_csv(tibble(cohort = "SCANB", row_id_overlap, col_id_overlap,
  row_gene_overlap, col_gene_overlap, selected_orientation = scan_orientation),
  file.path(subdirs["qc"], "scanb_orientation_audit.csv"))
readr::write_csv(
  bind_rows(mutate(met_resolved$audit,  cohort = "METABRIC"),
            mutate(scan_resolved$audit, cohort = "SCANB")),
  file.path(subdirs["qc"], "duplicate_gene_resolution.csv"))
coverage <- imap_dfr(cohorts, function(obj, cohort) {
  imap_dfr(gene_sets, function(genes, signature) tibble(
    cohort, signature,
    requested = length(genes),
    matched   = sum(genes %in% colnames(obj$expression)),
    coverage  = matched / requested,
    missing_genes = collapse_values(setdiff(genes, colnames(obj$expression)))))
})
readr::write_csv(coverage,
  file.path(subdirs["qc"], "signature_coverage.csv"))
if (any(coverage$coverage < minimum_signature_coverage |
        coverage$matched  < minimum_signature_genes))
  stop("A signature failed the prespecified gene-coverage threshold.")
if (any(coverage$matched < ssgsea_min_size |
        coverage$matched > ssgsea_max_size)) {
  stop("At least one matched signature falls outside ssgsea_min_size / ",
       "ssgsea_max_size and would be removed by GSVA.")
}


# ==============================================================================
# CALCULATE SSGSEA SCORES
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION: ssgsea-scores
# ------------------------------------------------------------------------------
# ssGSEA scoring and cohort-signature median assignment
run_ssgsea <- function(mat, cohort) {
  expression_genes_by_samples <- t(mat)
  sets <- lapply(gene_sets, intersect, y = rownames(expression_genes_by_samples))
  score_result <- tryCatch({
    if (exists("ssgseaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
      parameter <- GSVA::ssgseaParam(
        exprData  = expression_genes_by_samples,
        geneSets  = sets,
        alpha = ssgsea_alpha,
        minSize = ssgsea_min_size,
        maxSize = ssgsea_max_size,
        normalize = ssgsea_normalize)
      list(
        score = GSVA::gsva(parameter, verbose = FALSE,
                           BPPARAM = ssgsea_bpparam),
        interface = "ssgseaParam")
    } else {
      list(
        score = GSVA::gsva(expression_genes_by_samples, sets,
          method = "ssgsea", tau = ssgsea_alpha,
          ssgsea.norm = ssgsea_normalize,
          min.sz = ssgsea_min_size, max.sz = ssgsea_max_size,
          verbose = FALSE, BPPARAM = ssgsea_bpparam),
        interface = "legacy")
    }
  }, error = function(e) {
    stop(cohort, " ssGSEA calculation failed: ", conditionMessage(e),
         call. = FALSE)
  })
  score <- as.matrix(score_result$score)
  interface <- score_result$interface
  if (!all(names(sets) %in% rownames(score))) {
    stop(cohort, " ssGSEA output is missing one or more configured signatures.")
  }
  if (any(!is.finite(score))) {
    stop(cohort, " ssGSEA returned nonfinite scores.")
  }
  as.data.frame(score) %>%
    rownames_to_column("signature") %>%
    pivot_longer(-signature, names_to = "sample_id", values_to = "ssgsea_score") %>%
    group_by(signature) %>%
    mutate(cim_z = as.numeric(scale(ssgsea_score))) %>%
    ungroup() %>%
    mutate(cohort, gsva_version = as.character(packageVersion("GSVA")),
           interface, .before = 1)
}
ssgsea_scores <- imap_dfr(cohorts, ~ run_ssgsea(.x$expression, .y))
median_cutoffs <- ssgsea_scores %>%
  filter(is.finite(ssgsea_score), is.finite(cim_z)) %>%
  group_by(cohort, signature) %>%
  summarise(
    n_scored   = n(),
    raw_median = median(ssgsea_score),
    z_median   = median(cim_z),
    raw_min    = min(ssgsea_score),
    raw_max    = max(ssgsea_score),
    .groups = "drop"
  )
expected_cutoffs <- tidyr::crossing(
  cohort = names(cohorts), signature = names(gene_sets))
missing_cutoffs <- expected_cutoffs %>%
  anti_join(median_cutoffs, by = c("cohort", "signature"))
if (nrow(missing_cutoffs))
  stop("At least one cohort-signature combination has no finite ssGSEA cutoff.")
ssgsea_scores <- ssgsea_scores %>%
  left_join(median_cutoffs, by = c("cohort", "signature")) %>%
  mutate(
    median_group = case_when(
      !is.finite(ssgsea_score) ~ NA_character_,
      median_tie_rule == "high_includes_cutoff" &
        ssgsea_score >= raw_median ~ "High",
      median_tie_rule == "high_includes_cutoff" &
        ssgsea_score < raw_median ~ "Low",
      median_tie_rule == "low_includes_cutoff" &
        ssgsea_score > raw_median ~ "High",
      median_tie_rule == "low_includes_cutoff" &
        ssgsea_score <= raw_median ~ "Low",
      TRUE                                  ~ NA_character_
    ),
    median_group  = factor(median_group, levels = c("Low", "High")),
    at_raw_median = is.finite(ssgsea_score) & dplyr::near(ssgsea_score, raw_median)
  )
median_group_balance <- ssgsea_scores %>%
  group_by(cohort, signature, raw_median, z_median) %>%
  summarise(
    n_scored              = sum(!is.na(median_group)),
    low_n                 = sum(median_group == "Low",  na.rm = TRUE),
    high_n                = sum(median_group == "High", na.rm = TRUE),
    at_median_n           = sum(at_raw_median,           na.rm = TRUE),
    low_fraction          = low_n  / n_scored,
    high_fraction         = high_n / n_scored,
    minimum_group_fraction = min(low_fraction, high_fraction),
    balance_warning       = minimum_group_fraction < median_group_warning_fraction,
    .groups = "drop"
  )
median_assignment_audit <- ssgsea_scores %>%
  transmute(
    cohort, signature, sample_id, ssgsea_score, cim_z,
    raw_median, z_median, median_group, at_raw_median,
    assignment_rule = if_else(
      median_tie_rule == "high_includes_cutoff",
      "High if raw score >= cohort-signature median; Low otherwise",
      "High if raw score > cohort-signature median; Low otherwise"),
    assignment_valid = case_when(
      !is.finite(ssgsea_score) ~ is.na(median_group),
      median_tie_rule == "high_includes_cutoff" &
        ssgsea_score >= raw_median ~ as.character(median_group) == "High",
      median_tie_rule == "high_includes_cutoff" &
        ssgsea_score < raw_median ~ as.character(median_group) == "Low",
      median_tie_rule == "low_includes_cutoff" &
        ssgsea_score > raw_median ~ as.character(median_group) == "High",
      median_tie_rule == "low_includes_cutoff" &
        ssgsea_score <= raw_median ~ as.character(median_group) == "Low",
      TRUE ~ FALSE
    )
  )
if (any(!median_assignment_audit$assignment_valid))
  stop("Median-group assignment verification failed.")
if (any(median_group_balance$balance_warning))
  warning("At least one median split has a group smaller than the configured ",
          "warning fraction; inspect median_group_balance.csv.")
readr::write_csv(median_cutoffs,
  file.path(subdirs["scores"], "ssGSEA_median_cutoffs.csv"))
readr::write_csv(median_group_balance,
  file.path(subdirs["scores"], "median_group_balance.csv"))
readr::write_csv(median_assignment_audit,
  file.path(subdirs["scores"], "median_group_assignment_audit.csv"))
readr::write_csv(ssgsea_scores,
  file.path(subdirs["scores"], "sample_level_ssGSEA_scores_with_groups.csv"))
readr::write_csv(ssgsea_scores,
  file.path(subdirs["scores"], "sample_level_ssGSEA_scores.csv"))
score_qc <- ssgsea_scores %>%
  group_by(cohort, signature) %>%
  summarise(
    n          = n(),
    raw_mean   = mean(ssgsea_score),
    raw_sd     = sd(ssgsea_score),
    raw_median = median(ssgsea_score),
    z_mean     = mean(cim_z),
    z_sd       = sd(cim_z),
    .groups = "drop")
score_correlation <- ssgsea_scores %>%
  select(cohort, sample_id, signature, cim_z) %>%
  pivot_wider(names_from = signature, values_from = cim_z) %>%
  group_by(cohort) %>%
  summarise(
    n            = sum(complete.cases(`Fun-CIM`, `Dys-CIM`)),
    pearson_r    = cor(`Fun-CIM`, `Dys-CIM`, use = "complete.obs", method = "pearson"),
    spearman_rho = cor(`Fun-CIM`, `Dys-CIM`, use = "complete.obs", method = "spearman"),
    .groups = "drop")
readr::write_csv(score_qc,
  file.path(subdirs["scores"], "ssGSEA_score_QC.csv"))
readr::write_csv(score_correlation,
  file.path(subdirs["scores"], "Fun_Dys_score_correlation.csv"))
p_score <- ggplot(ssgsea_scores,
                  aes(x = signature, y = cim_z, fill = signature)) +
  geom_violin(trim = FALSE, alpha = .65, color = "grey30") +
  geom_boxplot(width = .16, outlier.shape = NA, fill = "white") +
  facet_wrap(~ cohort, scales = "free_y") +
  scale_fill_manual(values = c("Fun-CIM" = "#1B9E77", "Dys-CIM" = "#D95F02")) +
  labs(x = NULL, y = "Within-cohort standardized ssGSEA score", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
ggsave(file.path(subdirs["figures"], "ssGSEA_score_distributions.png"),
       p_score, width = score_qc_plot_width, height = score_qc_plot_height,
       dpi = figure_dpi, bg = "white")
score_pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
ggsave(file.path(subdirs["figures"], "ssGSEA_score_distributions.pdf"),
       p_score, width = score_qc_plot_width, height = score_qc_plot_height,
       device = score_pdf_device, bg = "white")


# ==============================================================================
# FULL-COEFFICIENT COX ENGINE
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION: cox-engine
# ------------------------------------------------------------------------------
# Cox model engine
model_covariates <- function(model_id) {
  idx <- match(model_id, model_definitions$model_id)
  if (is.na(idx)) stop("Unknown model: ", model_id)
  as.character(unlist(model_definitions$covariates[[idx]],
                      recursive = TRUE, use.names = FALSE))
}
exposure_specification <- function(parameterization) {
  parameterization <- match.arg(parameterization, allowed_parameterizations)
  if (parameterization == "continuous") {
    list(variable = "cim_z", coefficient_term = "cim_z",
         label = "CIM score (per 1 SD)", reference = "per 1-SD increase",
         comparison = "+1 SD")
  } else {
    list(variable = "median_group", coefficient_term = "median_groupHigh",
         label = "CIM score: High vs Low", reference = "Low",
         comparison = "High")
  }
}
prepare_model_data <- function(obj, cohort, endpoint_row, signature, model_id,
                               parameterization,
                               sample_mode = "available", common_ids = NULL) {
  score <- ssgsea_scores %>%
    filter(.data$cohort == .env$cohort, .data$signature == .env$signature)
  clinical <- obj$clinical
  clinical$sample_id <- clinical[[obj$id_col]]
  dat <- clinical %>%
    left_join(
      select(score, sample_id, ssgsea_score, cim_z, median_group,
             raw_median, z_median),
      by = "sample_id"
    ) %>%
    transmute(
      sample_id,
      time  = safe_numeric(.data[[endpoint_row$time_col]]),
      event = standardize_event(.data[[endpoint_row$event_col]],
                                endpoint_row$event_col, cohort),
      cim_z, ssgsea_score,
      median_group = factor(median_group, levels = c("Low", "High")),
      raw_median, z_median,
      age_5y, tumor_size, grade, pam50, clinical_subtype)
  covariates <- model_covariates(model_id)
  exposure   <- exposure_specification(parameterization)
  covariates_complete <- if (!length(covariates)) {
    rep(TRUE, nrow(dat))
  } else {
    complete.cases(dat[, covariates, drop = FALSE])
  }
  dat <- dat %>% mutate(
    endpoint_valid  = is.finite(time) & time > 0 & event %in% c(0, 1),
    score_available = if (exposure$variable == "cim_z") {
      is.finite(cim_z)
    } else {
      !is.na(median_group)
    },
    covariates_complete = .env$covariates_complete,
    included_available  = endpoint_valid & score_available & covariates_complete,
    exclusion_reason = case_when(
      !endpoint_valid      ~ "invalid_or_missing_endpoint",
      !score_available     ~ "missing_signature_score",
      !covariates_complete ~ "missing_model_covariate",
      TRUE                 ~ "included"))
  if (sample_mode == "common") {
    dat <- dat %>% mutate(
      included = included_available & sample_id %in% common_ids,
      exclusion_reason = if_else(included, "included",
        if_else(included_available, "not_in_common_sample", exclusion_reason)))
  } else {
    dat <- dat %>% mutate(included = included_available)
  }
  dat
}
fit_one_cox <- function(obj, cohort, endpoint_row, signature, model_id,
                        parameterization,
                        sample_mode = "available", common_ids = NULL) {
  dat_audit <- prepare_model_data(obj, cohort, endpoint_row, signature,
                                  model_id, parameterization,
                                  sample_mode, common_ids)
  dat <- dat_audit %>%
    filter(included) %>%
    mutate(across(where(is.factor), droplevels))
  covariates <- model_covariates(model_id)
  exposure   <- exposure_specification(parameterization)
  rhs          <- paste(c(exposure$variable, covariates), collapse = " + ")
  formula_text <- paste("survival::Surv(time, event) ~", rhs)
  metadata <- tibble(
    cohort,
    endpoint_original  = endpoint_row$endpoint_original,
    endpoint_harmonized = endpoint_row$endpoint_harmonized,
    signature, model_id,
    model_label = model_definitions$model_label[
      match(model_id, model_definitions$model_id)],
    model_role = model_definitions$model_role[
      match(model_id, model_definitions$model_id)],
    parameterization,
    exposure_variable    = exposure$variable,
    sample_mode, formula_text,
    configured_covariates = paste(covariates, collapse = "; "),
    raw_median_cutoff    = dplyr::first(dat_audit$raw_median),
    z_median_cutoff      = dplyr::first(dat_audit$z_median),
    n        = nrow(dat),
    events   = sum(dat$event == 1),
    censored = sum(dat$event == 0),
    low_n    = sum(dat$median_group == "Low",  na.rm = TRUE),
    high_n   = sum(dat$median_group == "High", na.rm = TRUE),
    low_events  = sum(dat$event == 1 & dat$median_group == "Low",  na.rm = TRUE),
    high_events = sum(dat$event == 1 & dat$median_group == "High", na.rm = TRUE),
    n_endpoint_score = sum(dat_audit$endpoint_valid & dat_audit$score_available),
    excluded_covariate_missing = sum(dat_audit$endpoint_valid &
      dat_audit$score_available & !dat_audit$covariates_complete))
  audit_out <- dat_audit %>% transmute(
    cohort,
    endpoint_original   = endpoint_row$endpoint_original,
    endpoint_harmonized = endpoint_row$endpoint_harmonized,
    signature, model_id, parameterization, sample_mode,
    sample_id, ssgsea_score, cim_z, raw_median, z_median, median_group,
    endpoint_valid, score_available, covariates_complete,
    included, exclusion_reason)
  fail <- function(status, warning_text, parameters_n = NA_integer_)
    list(
      coefficients = metadata %>% transmute(
        across(everything()),
        term = exposure$coefficient_term,
        harmonized_variable = exposure$variable,
        display_label       = exposure$label,
        reference_level     = exposure$reference,
        comparison_level    = exposure$comparison,
        estimate = NA_real_, hazard_ratio = NA_real_,
        conf_low = NA_real_, conf_high = NA_real_,
        p_value = NA_real_, term_ph_p = NA_real_,
        model_status = status, warning_text),
      summary = metadata %>% mutate(
        parameters_n,
        events_per_parameter = ifelse(parameters_n > 0,
                                      events / parameters_n, NA_real_),
        concordance    = NA_real_, concordance_se = NA_real_,
        aic            = NA_real_, log_likelihood = NA_real_,
        lr_chisq       = NA_real_, lr_df  = NA_real_, lr_p   = NA_real_,
        wald_chisq     = NA_real_, wald_df = NA_real_, wald_p = NA_real_,
        score_chisq    = NA_real_, score_df = NA_real_, score_p = NA_real_,
        global_ph_p    = NA_real_, max_vif = NA_real_,
        model_status   = status, warning_text),
      audit = audit_out, ph = tibble(), fit = NULL)
  if (nrow(dat) < min_n_model ||
      sum(dat$event) < min_events_model)
    return(fail("skipped", "insufficient n or events"))
  model_variables <- c(exposure$variable, covariates)
  nonestimable <- model_variables[vapply(model_variables, function(v)
    length(unique(dat[[v]][!is.na(dat[[v]])])) < 2L, logical(1))]
  if (length(nonestimable))
    return(fail("skipped",
                paste("fewer than two values/levels:",
                      paste(nonestimable, collapse = ", "))))
  mm <- model.matrix(as.formula(paste("~", rhs)), dat)
  parameters_n <- ncol(mm) - 1L
  if (qr(mm)$rank < ncol(mm))
    return(fail("skipped", "rank-deficient model matrix", parameters_n))
  low_epv <- sum(dat$event) < events_per_parameter_target * parameters_n
  if (low_epv && enforce_events_per_parameter)
    return(fail("skipped", "events-per-parameter safeguard not met", parameters_n))
  sparse_factors <- covariates[vapply(covariates, function(v) {
    x <- dat[[v]]
    (is.factor(x) || is.character(x)) &&
      any(table(x) < sparse_factor_level_n)
  }, logical(1))]
  warnings <- c(
    if (low_epv) paste0("events per parameter below configured target of ",
                        events_per_parameter_target),
    if (length(sparse_factors)) paste0("sparse factor level(s) in: ",
                                       paste(sparse_factors, collapse = ", "))
  )
  fit <- tryCatch(
    withCallingHandlers(
      survival::coxph(as.formula(formula_text), data = dat,
                      ties = "efron", x = TRUE, y = TRUE),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    error = function(e) e)
  if (inherits(fit, "error"))
    return(fail("failed", conditionMessage(fit), parameters_n))
  sm <- summary(fit)
  ph <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
  ph_table <- if (is.null(ph)) tibble() else
    as.data.frame(ph$table) %>%
    rownames_to_column("ph_term") %>% as_tibble() %>%
    rename(ph_chisq = chisq, ph_df = df, ph_p = p)
  term_ph <- if (nrow(ph_table))
    ph_table %>% filter(ph_term != "GLOBAL") else tibble()
  coefficient_table <- broom::tidy(fit, exponentiate = FALSE, conf.int = TRUE) %>%
    transmute(
      term,
      harmonized_variable = term_variable(term),
      display_label       = term_display(term),
      reference_level     = term_reference(term),
      comparison_level    = term_comparison(term),
      estimate,
      hazard_ratio = exp(estimate),
      conf_low     = exp(conf.low),
      conf_high    = exp(conf.high),
      p_value      = p.value) %>%
    left_join(
      term_ph %>% transmute(term = ph_term, term_ph_p = ph_p),
      by = "term") %>%
    mutate(
      model_status = "converged",
      warning_text = if (length(warnings)) collapse_values(warnings) else NA_character_) %>%
    bind_cols(metadata[rep(1, nrow(.)), ]) %>%
    relocate(all_of(names(metadata)))
  max_vif <- NA_real_
  if (ncol(mm) > 2L) {
    x <- mm[, -1, drop = FALSE]
    vif_values <- vapply(seq_len(ncol(x)), function(j) {
      others <- setdiff(seq_len(ncol(x)), j)
      if (!length(others)) return(1)
      r2 <- summary(lm(x[, j] ~ x[, others, drop = FALSE]))$r.squared
      ifelse(is.finite(r2) && r2 < 1, 1 / (1 - r2), Inf)
    }, numeric(1))
    max_vif <- max(vif_values, na.rm = TRUE)
  }
  ll <- logLik(fit)
  summary_table <- metadata %>% mutate(
    parameters_n,
    events_per_parameter = events / parameters_n,
    concordance    = unname(sm$concordance[1]),
    concordance_se = unname(sm$concordance[2]),
    aic            = AIC(fit),
    log_likelihood = as.numeric(ll),
    lr_chisq       = unname(sm$logtest["test"]),
    lr_df          = unname(sm$logtest["df"]),
    lr_p           = unname(sm$logtest["pvalue"]),
    wald_chisq     = unname(sm$waldtest["test"]),
    wald_df        = unname(sm$waldtest["df"]),
    wald_p         = unname(sm$waldtest["pvalue"]),
    score_chisq    = unname(sm$sctest["test"]),
    score_df       = unname(sm$sctest["df"]),
    score_p        = unname(sm$sctest["pvalue"]),
    global_ph_p    = if (nrow(ph_table))
      ph_table$ph_p[ph_table$ph_term == "GLOBAL"] else NA_real_,
    max_vif,
    model_status = "converged",
    warning_text = if (length(warnings)) collapse_values(warnings) else NA_character_)
  list(coefficients = coefficient_table, summary = summary_table,
       audit = audit_out,
       ph = ph_table %>% mutate(
         cohort,
         endpoint_original = endpoint_row$endpoint_original,
         signature, model_id, parameterization, sample_mode, .before = 1),
       fit = fit)
}


# ==============================================================================
# FIT AVAILABLE-CASE AND COMMON-SAMPLE MODELS
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION: fit-models
# ------------------------------------------------------------------------------
available_fits        <- list()
common_fits           <- list()
common_sample_manifest <- tibble()

for (cohort in names(cohorts)) {
  obj <- cohorts[[cohort]]
  for (i in seq_len(nrow(obj$endpoints))) {
    endpoint_row <- obj$endpoints[i, , drop = FALSE]
    for (signature in names(gene_sets)) {
      key_base <- paste(cohort, endpoint_row$endpoint_original,
                        signature, sep = "__")
      for (parameterization in score_parameterizations) {
        for (model_id in model_levels) {
          fit_key <- paste(key_base, parameterization, model_id,
                           "available", sep = "__")
          available_fits[[fit_key]] <- fit_one_cox(
            obj, cohort, endpoint_row, signature, model_id,
            parameterization, "available"
          )
        }
      }
      if (run_common_sample_sensitivity) {
        # Anchor every model to the complete-case population of the configured
        # primary adjusted model.
        full_audit <- prepare_model_data(
          obj, cohort, endpoint_row, signature,
          primary_adjusted_model, "continuous"
        )
        # ────────────────────────────────────────────────────────────────────
        common_ids <- full_audit$sample_id[full_audit$included_available]
        common_sample_manifest <- bind_rows(
          common_sample_manifest,
          tibble(
            cohort,
            endpoint_original = endpoint_row$endpoint_original,
            signature,
            raw_median_cutoff = first(full_audit$raw_median),
            n_common          = length(common_ids),
            common_low_n      = sum(full_audit$sample_id %in% common_ids &
                                      full_audit$median_group == "Low",  na.rm = TRUE),
            common_high_n     = sum(full_audit$sample_id %in% common_ids &
                                      full_audit$median_group == "High", na.rm = TRUE),
            common_ids        = paste(common_ids, collapse = "; ")
          )
        )
        for (parameterization in score_parameterizations) {
          for (model_id in model_levels) {
            fit_key <- paste(key_base, parameterization, model_id,
                             "common", sep = "__")
            common_fits[[fit_key]] <- fit_one_cox(
              obj, cohort, endpoint_row, signature, model_id,
              parameterization, "common", common_ids
            )
          }
        }
      }
    }
  }
}

all_fits <- c(available_fits, common_fits)

coefficient_results <- map_dfr(all_fits, "coefficients")
model_summary       <- map_dfr(all_fits, "summary")
patient_audit       <- map_dfr(all_fits, "audit")
ph_results          <- map_dfr(all_fits, "ph")

coefficient_results <- coefficient_results %>%
  group_by(cohort, endpoint_original, signature, parameterization,
           sample_mode, model_id) %>%
  mutate(p_adjust_model_terms = ifelse(
    model_status == "converged",
    p.adjust(p_value, method = "BH"), NA_real_)) %>%
  ungroup() %>%
  group_by(cohort, endpoint_original, parameterization,
           sample_mode, model_id, term) %>%
  mutate(p_adjust_signature = ifelse(
    term %in% c("cim_z", "median_groupHigh") & model_status == "converged",
    p.adjust(p_value, method = "BH"), NA_real_)) %>%
  ungroup() %>%
  mutate(
    exposure_coefficient = term %in% c("cim_z", "median_groupHigh"),
    cox_fdr_supported    = exposure_coefficient &
      is.finite(p_adjust_signature) &
      apply_comparator(p_adjust_signature, cox_fdr_threshold,
                       cox_fdr_operator),
    significance_value   = if (significance_metric == "p_value")
      p_value else p_adjust_signature,
    significance_cutoff  = if (significance_metric == "p_value")
      alpha_threshold else cox_fdr_threshold,
    significant = exposure_coefficient & is.finite(significance_value) &
      if (significance_metric == "p_value") {
        significance_value < significance_cutoff
      } else {
        apply_comparator(significance_value, significance_cutoff,
                         cox_fdr_operator)
      },
    model_id    = factor(model_id, levels = model_levels),
    model_color = unname(model_colors[as.character(model_id)])
  )

coefficient_source_map <- bind_rows(met_audit, scan_audit) %>%
  transmute(cohort, harmonized_variable = target,
            source_column    = source,
            source_derivation = derivation) %>%
  bind_rows(tibble(
    cohort = c("METABRIC", "SCANB"),
    harmonized_variable = "cim_z",
    source_column    = "ssGSEA score",
    source_derivation = "within-cohort z standardization")) %>%
  bind_rows(tibble(
    cohort = c("METABRIC", "SCANB"),
    harmonized_variable = "median_group",
    source_column    = "ssGSEA score",
    source_derivation = "fixed raw-score median within cohort and signature")) %>%
  distinct(cohort, harmonized_variable, .keep_all = TRUE)

coefficient_results <- coefficient_results %>%
  left_join(coefficient_source_map,
            by = c("cohort", "harmonized_variable"))

signature_results    <- coefficient_results %>% filter(exposure_coefficient)
continuous_results   <- signature_results   %>% filter(parameterization == "continuous")
median_split_results <- signature_results   %>% filter(parameterization == "median_split")

readr::write_csv(coefficient_results,
  file.path(subdirs["models"], "complete_coefficient_results.csv"))
readr::write_csv(signature_results,
  file.path(subdirs["models"], "signature_model_comparison.csv"))
readr::write_csv(continuous_results,
  file.path(subdirs["models"], "continuous_adjusted_cox_results.csv"))
readr::write_csv(median_split_results,
  file.path(subdirs["models"], "median_split_adjusted_cox_results.csv"))
readr::write_csv(model_summary,
  file.path(subdirs["models"], "model_level_summary.csv"))
readr::write_csv(patient_audit,
  file.path(subdirs["models"], "patient_inclusion_audit.csv"))
readr::write_csv(ph_results,
  file.path(subdirs["models"], "proportional_hazards_tests.csv"))
readr::write_csv(common_sample_manifest,
  file.path(subdirs["sensitivity"], "common_sample_manifest.csv"))

ph_plot_index <- imap_dfr(all_fits, function(result, fit_key) {
  if (is.null(result$fit) || !nrow(result$ph) ||
      !any(is.finite(result$ph$ph_p) & result$ph$ph_p < alpha_threshold))
    return(tibble())
  ph_object <- tryCatch(survival::cox.zph(result$fit), error = function(e) NULL)
  if (is.null(ph_object)) return(tibble())
  safe_key  <- str_replace_all(fit_key, "[^A-Za-z0-9_-]", "_")
  png_file  <- file.path(subdirs["sensitivity"],
                         paste0("PH_diagnostic__", safe_key, ".png"))
  pdf_file  <- file.path(subdirs["sensitivity"],
                         paste0("PH_diagnostic__", safe_key, ".pdf"))
  n_ph_terms  <- max(1L, nrow(ph_object$table) - 1L)
  panel_rows  <- ceiling(sqrt(n_ph_terms))
  panel_cols  <- ceiling(n_ph_terms / panel_rows)
  grDevices::png(png_file, width = 2400, height = 1800, res = 300)
  par(mfrow = c(panel_rows, panel_cols), mar = c(4, 4, 2, 1))
  for (j in seq_len(n_ph_terms)) plot(ph_object[j])
  grDevices::dev.off()
  grDevices::pdf(pdf_file, width = 8, height = 6)
  par(mfrow = c(panel_rows, panel_cols), mar = c(4, 4, 2, 1))
  for (j in seq_len(n_ph_terms)) plot(ph_object[j])
  grDevices::dev.off()
  tibble(fit_key, png_file, pdf_file,
         reason = "at least one term/global PH p < alpha")
})
readr::write_csv(ph_plot_index,
  file.path(subdirs["sensitivity"], "PH_diagnostic_plot_index.csv"))

covariate_level_audit <- imap_dfr(cohorts, function(obj, cohort) {
  map_dfr(c("age_5y", "tumor_size", "grade", "pam50", "clinical_subtype"),
    function(variable) {
      x <- obj$clinical[[variable]]
      if (is.factor(x) || is.character(x)) {
        as_tibble(table(level = as.character(x), useNA = "ifany"),
                  .name_repair = "minimal") %>%
          setNames(c("level", "n")) %>%
          mutate(cohort, variable,
                 missing_n = sum(is.na(x)), .before = 1)
      } else {
        tibble(cohort, variable, level = "continuous",
               n = sum(!is.na(x)), missing_n = sum(is.na(x)))
      }
    })
})
readr::write_csv(covariate_level_audit,
  file.path(subdirs["qc"], "covariate_level_and_missingness_audit.csv"))

model_covariate_audit <- imap_dfr(cohorts, function(obj, cohort) {
  map_dfr(seq_len(nrow(obj$endpoints)), function(i) {
    endpoint_row <- obj$endpoints[i, , drop = FALSE]
    map_dfr(names(gene_sets), function(signature) {
      map_dfr(model_levels, function(model_id) {
        dat <- prepare_model_data(
          obj, cohort, endpoint_row, signature, model_id, "continuous"
        ) %>% filter(included)
        covariates <- model_covariates(model_id)
        if (!length(covariates)) return(tibble())
        map_dfr(covariates, function(variable) {
          x <- dat[[variable]]
          if (is.factor(x) || is.character(x)) {
            observed_levels <- levels(droplevels(factor(x)))
            map_dfr(observed_levels, function(level) {
              in_level <- !is.na(x) & as.character(x) == level
              tibble(
                cohort,
                endpoint_original   = endpoint_row$endpoint_original,
                endpoint_harmonized = endpoint_row$endpoint_harmonized,
                signature, model_id, variable,
                variable_type   = "categorical",
                level,
                reference_level = levels(x)[1],
                n_model         = nrow(dat),
                events_model    = sum(dat$event == 1),
                level_n         = sum(in_level),
                level_events    = sum(dat$event[in_level] == 1),
                mean = NA_real_, sd = NA_real_,
                minimum = NA_real_, maximum = NA_real_
              )
            })
          } else {
            tibble(
              cohort,
              endpoint_original   = endpoint_row$endpoint_original,
              endpoint_harmonized = endpoint_row$endpoint_harmonized,
              signature, model_id, variable,
              variable_type   = "continuous",
              level           = NA_character_,
              reference_level = "per 5-year increase",
              n_model         = nrow(dat),
              events_model    = sum(dat$event == 1),
              level_n         = sum(is.finite(x)),
              level_events    = NA_integer_,
              mean    = mean(x, na.rm = TRUE),
              sd      = sd(x,   na.rm = TRUE),
              minimum = min(x,  na.rm = TRUE),
              maximum = max(x,  na.rm = TRUE)
            )
          }
        })
      })
    })
  })
}) %>% mutate(model_id = factor(model_id, levels = model_levels))
readr::write_csv(model_covariate_audit,
  file.path(subdirs["models"], "model_specific_covariate_audit.csv"))

model_covariate_inclusion <- model_definitions %>%
  transmute(
    model_id, model_label, model_role,
    included_covariates = map_chr(covariates,
      ~ paste(as.character(unlist(.x, recursive = TRUE)), collapse = "; ")),
    includes_age = map_lgl(covariates,
      ~ "age_5y"           %in% unlist(.x, recursive = TRUE)),
    includes_tumor_size = map_lgl(covariates,
      ~ "tumor_size"       %in% unlist(.x, recursive = TRUE)),
    includes_grade = map_lgl(covariates,
      ~ "grade"            %in% unlist(.x, recursive = TRUE)),
    includes_pam50 = map_lgl(covariates,
      ~ "pam50"            %in% unlist(.x, recursive = TRUE)),
    includes_clinical_subtype = map_lgl(covariates,
      ~ "clinical_subtype" %in% unlist(.x, recursive = TRUE)),
    color)
readr::write_csv(model_covariate_inclusion,
  file.path(subdirs["models"], "model_covariate_inclusion_map.csv"))


# ==============================================================================
# PUBLICATION FOREST PLOTS
# ==============================================================================

# This section produces publication-quality forest plots for the model and score
# parameterization selected by `forest_model` and `forest_parameterization`.
# Only the ssGSEA exposure coefficient is displayed; adjustment-variable
# coefficients are intentionally excluded. Three figures are generated:

# 1. **ssGSEA Fun-CIM score** — individual forest
# 2. **ssGSEA Dys-CIM score** — individual forest
# 3. **Combined Fun-CIM and Dys-CIM** — both signatures on a shared linear HR axis

# Individual figures include N and Events columns. The combined figure omits
# these columns to preserve readability; both fields remain in the exported CSV.
# The combined figure has no title or subtitle so it can be labeled at the
# multipanel level. All three figures share the same linear HR axis. The axis can
# be manually specified or calculated automatically from the plotted CIs.


# ------------------------------------------------------------------------------
# SHARED HR AXIS
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# SECTION: forest-axis
# ------------------------------------------------------------------------------
plot_data_forest <- tibble()
if (make_forest_plots) {
# ── Build plot_skeleton: the fixed six-row endpoint layout ──────────────────
plot_skeleton <- tibble::tribble(
  ~cohort,    ~endpoint_original, ~endpoint_group,        ~group_order, ~row_order, ~y,
  "METABRIC", "OS",               "Overall survival",              1,          1,  6,
  "SCANB",    "OS",               "Overall survival",              1,          2,  5,
  "METABRIC", "RFI",              "Recurrence-focused",            2,          3,  4,
  "SCANB",    "RFi",              "Recurrence-focused",            2,          4,  3,
  "METABRIC", "DSS",              "Cohort-specific",               3,          5,  2,
  "SCANB",    "DRFi",             "Cohort-specific",               3,          6,  1
) %>% mutate(
  cohort_display = if_else(cohort == "SCANB", "SCAN-B", cohort),
  endpoint_display = case_when(
    endpoint_original %in% c("RFI", "RFi")   ~ "RFI",
    endpoint_original %in% c("DRFI", "DRFi") ~ "DRFI",
    TRUE ~ endpoint_original
  ),
  row_label = paste(cohort_display, endpoint_display, sep = " | ")
)

# ── Build forest data from the configured model and parameterization ────────
forest_exposure <- exposure_specification(forest_parameterization)
plot_data_forest <- signature_results %>%
  filter(
    as.character(model_id) == forest_model,
    sample_mode   == "available",
    term          == forest_exposure$coefficient_term,
    parameterization == forest_parameterization,
    signature     %in% c("Fun-CIM", "Dys-CIM")
  ) %>%
  left_join(plot_skeleton,
            by = c("cohort", "endpoint_original"),
            relationship = "many-to-one") %>%
  arrange(signature, row_order) %>%
  mutate(
    estimable = model_status == "converged" &
      is.finite(hazard_ratio) & is.finite(conf_low) & is.finite(conf_high) &
      hazard_ratio > 0 & conf_low > 0 & conf_high > 0,
    valid_ci_order = !estimable |
      (conf_low <= hazard_ratio & hazard_ratio <= conf_high),
    hr_ci_text = if_else(
      estimable,
      sprintf("%.2f (%.2f\u2013%.2f)", hazard_ratio, conf_low, conf_high),
      "NE"),
    p_adjust_text = case_when(
      !estimable | !is.finite(p_adjust_signature) ~ "NE",
      p_adjust_signature < 0.001                  ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_adjust_signature)),
    significant_bh = estimable & is.finite(p_adjust_signature) &
      apply_comparator(p_adjust_signature, cox_fdr_threshold,
                       cox_fdr_operator),
    n_events_text  = if_else(!is.na(n) & !is.na(events),
                             paste0(n, " / ", events), "NE")
  )
expected_forest_rows <- tidyr::crossing(
  plot_skeleton %>% select(cohort, endpoint_original),
  signature = c("Fun-CIM", "Dys-CIM"))
missing_forest_rows <- expected_forest_rows %>%
  anti_join(plot_data_forest,
            by = c("cohort", "endpoint_original", "signature"))
if (nrow(missing_forest_rows)) {
  stop("The configured forest is missing cohort-endpoint-signature rows; ",
       "inspect the requested model and parameterization.")
}

invalid_ci_forest <- plot_data_forest %>% filter(!valid_ci_order)
if (nrow(invalid_ci_forest))
  stop("At least one forest estimate is not contained within its 95% CI.")

readr::write_csv(plot_data_forest,
  file.path(subdirs["forests"], paste0(
    forest_model, "_", forest_parameterization,
    "_forest_plot_source_data.csv")))
purrr::walk(c("Fun-CIM", "Dys-CIM"), function(sig) {
  stem <- paste0(stringr::str_replace_all(sig, "-", "_"),
                 "_", forest_model, "_", forest_parameterization,
                 "_consolidated_forest_source.csv")
  readr::write_csv(dplyr::filter(plot_data_forest, signature == sig),
                   file.path(subdirs["forests"], stem))
})

# ── Shared linear HR axis ────────────────────────────────────────────────────
finite_ci_vals <- plot_data_forest %>%
  filter(estimable) %>%
  select(conf_low, conf_high) %>%
  unlist(use.names = FALSE)
finite_ci_vals <- finite_ci_vals[is.finite(finite_ci_vals) & finite_ci_vals > 0]
if (!length(finite_ci_vals))
  stop("No finite positive confidence intervals are available for the configured forest.")

if (hr_axis_mode == "manual") {
  hr_limits <- c(hr_axis_min, hr_axis_max)
  outside_axis <- plot_data_forest %>%
    filter(estimable & (conf_low < hr_limits[1] | conf_high > hr_limits[2]))
  if (nrow(outside_axis)) {
    stop("Manual HR-axis limits would truncate one or more 95% CIs; ",
         "widen hr_axis_min / hr_axis_max or use hr_axis_mode: automatic.")
  }
  hr_breaks <- manual_hr_breaks[
    manual_hr_breaks >= hr_limits[1] & manual_hr_breaks <= hr_limits[2]]
} else {
  raw_range <- range(c(finite_ci_vals, 1), finite = TRUE)
  padding <- max(diff(raw_range) * 0.08, 0.03)
  hr_limits <- c(max(.001, raw_range[1] - padding), raw_range[2] + padding)
  hr_breaks <- scales::breaks_pretty(n = 6)(hr_limits)
  hr_breaks <- hr_breaks[hr_breaks >= hr_limits[1] & hr_breaks <= hr_limits[2]]
}
if (length(hr_breaks) < 2L) stop("Too few HR-axis breaks remain after validation.")
if (!1 %in% hr_breaks && 1 >= hr_limits[1] && 1 <= hr_limits[2])
  hr_breaks <- sort(unique(c(hr_breaks, 1)))

readr::write_csv(
  tibble(axis_min = hr_limits[1], axis_max = hr_limits[2],
         axis_breaks = paste(hr_breaks, collapse = "; "),
         axis_scale  = "linear",
         axis_source = hr_axis_mode,
         forest_model, forest_parameterization),
  file.path(subdirs["forests"], "shared_HR_axis_audit.csv"))
}


# ------------------------------------------------------------------------------
# FOREST PLOT FUNCTIONS
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# SECTION: forest-plot-functions
# ------------------------------------------------------------------------------
# ── Visual constants ─────────────────────────────────────────────────────────
signature_colors <- c("Fun-CIM" = "#1B9E77", "Dys-CIM" = "#D95F02")
signature_titles <- c(
  "Fun-CIM" = "ssGSEA Fun-CIM score",
  "Dys-CIM" = "ssGSEA Dys-CIM score"
)
combined_signature_headers <- c(
  "Fun-CIM" = "Fun-CIM score",
  "Dys-CIM" = "Dys-CIM score"
)
pt_to_mm <- function(pt) pt / 2.845276
forest_covariates <- model_covariates(forest_model)
forest_covariate_labels <- c(
  age_5y = "age", tumor_size = "tumor size",
  grade = "histologic grade", pam50 = "PAM50 subtype",
  clinical_subtype = "clinical subtype")
if (!length(forest_covariates)) {
  figure_subtitle <- "Unadjusted"
} else {
  covariate_text <- unname(forest_covariate_labels[forest_covariates])
  covariate_text[is.na(covariate_text)] <- forest_covariates[
    is.na(covariate_text)]
  if (length(covariate_text) > 1L) {
    covariate_text <- paste0(
      paste(covariate_text[-length(covariate_text)], collapse = ", "),
      ", and ", covariate_text[length(covariate_text)])
  }
  figure_subtitle <- stringr::str_wrap(
    paste("Adjusted for", covariate_text), width = 58)
}
if (forest_parameterization == "median_split") {
  figure_subtitle <- stringr::str_wrap(
    paste("High versus Low (Low reference);", figure_subtitle), width = 58)
}
forest_x_label <- if (forest_parameterization == "continuous") {
  "Hazard ratio per 1-SD increase in ssGSEA score"
} else {
  "Hazard ratio: High vs Low ssGSEA score"
}
base_theme <- theme_void(base_size = 8) +
  theme(
    plot.margin = margin(0, 1.5, 0, 1.5),
    text = element_text(family = "sans", color = "black")
  )

# ── Individual signature forest (with N and Events columns) ─────────────────
make_signature_forest <- function(signature) {
  signature <- match.arg(signature, c("Fun-CIM", "Dys-CIM"))
  d <- plot_data_forest %>% filter(.data$signature == .env$signature)
  sig_color <- unname(signature_colors[signature])

  label_panel <- ggplot() +
    geom_text(data = tibble(x = 0, y = 7, label = "Cohort | Endpoint"),
              aes(x, y, label = label), hjust = 0,
              size = pt_to_mm(11.5), fontface = "bold", color = "black") +
    geom_text(data = d, aes(x = 0.03, y = y, label = row_label),
              hjust = 0, size = pt_to_mm(11), color = "black") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.45, 7.25), clip = "off") +
    base_theme + theme(plot.margin = margin(0, 0, 0, 1.5))

  n_panel <- ggplot() +
    geom_text(data = tibble(x = 0.5, y = 7, label = "N"),
              aes(x, y, label = label), hjust = 0.5,
              size = pt_to_mm(10.8), fontface = "bold", color = "black") +
    geom_text(data = d %>% mutate(n_text = if_else(!is.na(n), as.character(n), "NE")),
              aes(x = 0.5, y = y, label = n_text),
              hjust = 0.5, size = pt_to_mm(10.5), color = "black") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.45, 7.25), clip = "off") +
    base_theme + theme(plot.margin = margin(0, 0.8, 0, 0.8))

  events_panel <- ggplot() +
    geom_text(data = tibble(x = 0.5, y = 7, label = "Events"),
              aes(x, y, label = label), hjust = 0.5,
              size = pt_to_mm(10.8), fontface = "bold", color = "black") +
    geom_text(
      data = d %>% mutate(
        events_text = if_else(!is.na(events), as.character(events), "NE")),
      aes(x = 0.5, y = y, label = events_text),
      hjust = 0.5, size = pt_to_mm(10.5), color = "black") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.45, 7.25), clip = "off") +
    base_theme + theme(plot.margin = margin(0, 1.0, 0, 0.8))

  forest_panel <- ggplot(d, aes(y = y)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               linewidth = 0.28, color = "grey72") +
    geom_segment(data = filter(d, estimable),
                 aes(x = conf_low, xend = conf_high, yend = y),
                 linewidth = 0.58, color = sig_color, lineend = "round") +
    geom_segment(data = filter(d, estimable),
                 aes(x = conf_low,  xend = conf_low,
                     y = y - 0.11,  yend = y + 0.11),
                 linewidth = 0.42, color = sig_color) +
    geom_segment(data = filter(d, estimable),
                 aes(x = conf_high, xend = conf_high,
                     y = y - 0.11,  yend = y + 0.11),
                 linewidth = 0.42, color = sig_color) +
    geom_point(data = filter(d, estimable),
               aes(x = hazard_ratio),
               shape = 16, size = 2.05, color = sig_color) +
    geom_text(data = filter(d, !estimable),
              aes(x = 1, label = "NE"),
              size = 2.5, color = "#555555") +
    scale_x_continuous(
      limits = hr_limits, breaks = hr_breaks,
      labels = scales::label_number(accuracy = 0.01, trim = TRUE),
      expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0.45, 7.25), breaks = NULL,
                       expand = expansion(mult = c(0, 0))) +
    labs(x = forest_x_label, y = NULL) +
    theme_classic(base_size = 8) +
    theme(
      axis.line.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.x  = element_line(color = "black", linewidth = 0.4),
      axis.ticks.x = element_line(color = "black", linewidth = 0.35),
      axis.text.x  = element_text(size = 10,   face = "plain", color = "black"),
      axis.title.x = element_text(size = 11.5, face = "bold",  color = "black",
                                  margin = margin(t = 4)),
      plot.margin  = margin(0, 1, 0, 0))

  hr_panel <- ggplot() +
    geom_text(data = tibble(x = 0, y = 7, label = "HR (95% CI)"),
              aes(x, y, label = label), hjust = 0,
              size = pt_to_mm(10.8), fontface = "bold", color = "black") +
    geom_text(data = d, aes(x = 0, y = y, label = hr_ci_text),
              hjust = 0, size = pt_to_mm(10.5),
              color = ifelse(d$estimable, "black", "#555555")) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.45, 7.25), clip = "off") +
    base_theme + theme(plot.margin = margin(0, 1.5, 0, 0))

  adjusted_p_panel <- ggplot() +
    geom_text(data = tibble(x = 1, y = 7, label = "Adjusted p"),
              aes(x, y, label = label), hjust = 1,
              size = pt_to_mm(10.8), fontface = "bold", color = "black") +
    geom_text(data = d, aes(x = 1, y = y, label = p_adjust_text),
              hjust = 1, size = pt_to_mm(10.5), fontface = "plain",
              color = "black") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.45, 7.25), clip = "off") +
    base_theme + theme(plot.margin = margin(0, 0, 0, 1.5))

  combined <- label_panel + n_panel + events_panel +
    forest_panel + hr_panel + adjusted_p_panel +
    patchwork::plot_layout(widths = c(1.20, 0.30, 0.42, 1.85, 1.22, 0.60)) +
    patchwork::plot_annotation(
      title    = unname(signature_titles[signature]),
      subtitle = figure_subtitle,
      theme = theme(
        plot.title    = element_text(family = "sans", size = 13.5,
                                     face = "bold", color = "black",
                                     hjust = 0.5, margin = margin(b = 2)),
        plot.subtitle = element_text(family = "sans", size = 9,
                                     face = "plain", color = "black",
                                     hjust = 0.5, lineheight = 1.0,
                                     margin = margin(t = 1, b = 5)),
        plot.margin = margin(5, 5, 4, 5)))
  list(plot = combined, data = d)
}

# ── Combined two-signature forest (no title/subtitle; no N/Events columns) ──
make_combined_forest <- function() {
  d <- plot_data_forest %>%
    mutate(
      signature = factor(signature, levels = c("Dys-CIM", "Fun-CIM")),
      forest_y  = y + if_else(signature == "Dys-CIM", 0.11, -0.11)
    )
  row_labels <- plot_skeleton %>% select(y, row_label)
  axis_span <- diff(hr_limits)
  legend_fun_point_x <- hr_limits[1] + 0.12 * axis_span
  legend_fun_text_x  <- hr_limits[1] + 0.16 * axis_span
  legend_dys_point_x <- hr_limits[1] + 0.58 * axis_span
  legend_dys_text_x  <- hr_limits[1] + 0.62 * axis_span

  label_panel <- ggplot() +
    geom_text(
      data = tibble(x = 0, y = 6.82, label = "Cohort | Endpoint"),
      aes(x, y, label = label), hjust = 0,
      size = pt_to_mm(11.8), fontface = "bold", color = "black") +
    geom_text(data = row_labels,
              aes(x = 0.03, y = y, label = row_label),
              hjust = 0, size = pt_to_mm(11.2), color = "black") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.42, 7.35), clip = "off") +
    base_theme + theme(plot.margin = margin(0, 0, 0, 1.5))

  forest_panel <- ggplot(d, aes(y = forest_y, color = signature)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               linewidth = 0.28, color = "grey72") +
    geom_segment(data = filter(d, estimable),
                 aes(x = conf_low, xend = conf_high, yend = forest_y),
                 linewidth = 0.56, lineend = "round", show.legend = FALSE) +
    geom_segment(data = filter(d, estimable),
                 aes(x = conf_low,  xend = conf_low,
                     y = forest_y - 0.075, yend = forest_y + 0.075),
                 linewidth = 0.40, show.legend = FALSE) +
    geom_segment(data = filter(d, estimable),
                 aes(x = conf_high, xend = conf_high,
                     y = forest_y - 0.075, yend = forest_y + 0.075),
                 linewidth = 0.40, show.legend = FALSE) +
    geom_point(data = filter(d, estimable),
               aes(x = hazard_ratio),
               shape = 16, size = 1.90, show.legend = FALSE) +
    geom_text(data = filter(d, !estimable),
              aes(x = 1, label = "NE"),
              size = 2.2, show.legend = FALSE) +
    annotate("point", x = legend_fun_point_x, y = 7.18,
             shape = 16, size = 1.85,
             color = signature_colors[["Fun-CIM"]]) +
    annotate("text", x = legend_fun_text_x, y = 7.18, label = "Fun-CIM",
             hjust = 0, size = pt_to_mm(10.5), color = "black") +
    annotate("point", x = legend_dys_point_x, y = 7.18,
             shape = 16, size = 1.85,
             color = signature_colors[["Dys-CIM"]]) +
    annotate("text", x = legend_dys_text_x, y = 7.18, label = "Dys-CIM",
             hjust = 0, size = pt_to_mm(10.5), color = "black") +
    scale_color_manual(values = signature_colors, guide = "none") +
    scale_x_continuous(
      limits = hr_limits, breaks = hr_breaks,
      labels = scales::label_number(accuracy = 0.01, trim = TRUE),
      expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0.42, 7.35), breaks = NULL,
                       expand = expansion(mult = c(0, 0))) +
    labs(x = forest_x_label, y = NULL) +
    theme_classic(base_size = 8) +
    theme(
      axis.line.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.x  = element_line(color = "black", linewidth = 0.4),
      axis.ticks.x = element_line(color = "black", linewidth = 0.35),
      axis.text.x  = element_text(size = 10.5, face = "plain", color = "black"),
      axis.title.x = element_text(size = 11.5, face = "bold",  color = "black",
                                  margin = margin(t = 4)),
      plot.margin  = margin(0, 2, 0, 0),
      legend.position = "none")

  make_hr_panel <- function(selected_signature) {
    stats_d     <- d %>% filter(signature == selected_signature)
    header_text <- unname(combined_signature_headers[selected_signature])
    ggplot() +
      annotate("text", x = 0.72, y = 7.20, label = header_text,
               size = pt_to_mm(11.8), fontface = "bold", color = "black") +
      geom_text(data = tibble(x = 0, y = 6.82, label = "HR (95% CI)"),
                aes(x, y, label = label), hjust = 0,
                size = pt_to_mm(10.5), fontface = "bold", color = "black") +
      geom_text(data = stats_d,
                aes(x = 0, y = y, label = hr_ci_text), hjust = 0,
                size = pt_to_mm(10.5),
                color = ifelse(stats_d$estimable, "black", "#555555")) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0.42, 7.35), clip = "off") +
      base_theme + theme(plot.margin = margin(0, 1.5, 0, 0))
  }

  make_adjusted_p_panel <- function(selected_signature) {
    stats_d <- d %>% filter(signature == selected_signature)
    ggplot() +
      geom_text(data = tibble(x = 1, y = 6.82, label = "Adjusted p"),
                aes(x, y, label = label), hjust = 1,
                size = pt_to_mm(10.5), fontface = "bold", color = "black") +
      geom_text(data = stats_d,
                aes(x = 1, y = y, label = p_adjust_text), hjust = 1,
                size = pt_to_mm(10.7), fontface = "plain", color = "black") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0.42, 7.35), clip = "off") +
      base_theme + theme(plot.margin = margin(0, 0, 0, 1.5))
  }

  dys_hr <- make_hr_panel("Dys-CIM")
  dys_p  <- make_adjusted_p_panel("Dys-CIM")
  fun_hr <- make_hr_panel("Fun-CIM")
  fun_p  <- make_adjusted_p_panel("Fun-CIM")

  combined <- label_panel + forest_panel + dys_hr + dys_p + fun_hr + fun_p +
    patchwork::plot_layout(widths = c(1.08, 1.95, 1.28, 0.58, 1.28, 0.58))

  list(plot = combined, data = d)
}


# ------------------------------------------------------------------------------
# GENERATE AND SAVE PUBLICATION FOREST FIGURES
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# SECTION: generate-forest-figures
# ------------------------------------------------------------------------------
forest_output_manifest <- tibble()
if (make_forest_plots) {
fun_result      <- make_signature_forest("Fun-CIM")
dys_result      <- make_signature_forest("Dys-CIM")
combined_result <- make_combined_forest()

save_forest <- function(plot, stem, width = forest_width,
                        height = forest_height) {
  png_file <- file.path(subdirs["figures"], paste0(stem, ".png"))
  pdf_file <- file.path(subdirs["figures"], paste0(stem, ".pdf"))
  ggsave(png_file, plot = plot, width = width, height = height,
         units = "in", dpi = figure_dpi, bg = "white")
  if (capabilities("cairo")) {
    ggsave(pdf_file, plot = plot, width = width, height = height,
           units = "in", device = grDevices::cairo_pdf, bg = "white")
  } else {
    ggsave(pdf_file, plot = plot, width = width, height = height,
           units = "in", device = "pdf", bg = "white")
  }
  c(png = png_file, pdf = pdf_file)
}

forest_tag <- paste(forest_model, forest_parameterization, sep = "_")
fun_files <- save_forest(fun_result$plot,
  paste0("Fun_CIM_", forest_tag, "_consolidated_forest"))
dys_files <- save_forest(dys_result$plot,
  paste0("Dys_CIM_", forest_tag, "_consolidated_forest"))
combined_files <- save_forest(combined_result$plot,
  paste0("Combined_CIM_", forest_tag, "_consolidated_forest"),
  width = combined_forest_width, height = forest_height)

readr::write_csv(
  combined_result$data %>% select(
    signature, cohort, cohort_display, endpoint_original, endpoint_display,
    row_order, y, forest_y, model_id, sample_mode, term, model_status,
    n, events, hazard_ratio, conf_low, conf_high,
    p_value, p_adjust_signature, significant_bh, estimable, formula_text,
    warning_text),
  file.path(subdirs["forests"], paste0(
    "Combined_CIM_", forest_tag, "_consolidated_forest.csv")))

forest_output_manifest <- tibble(
  signature = c("Fun-CIM", "Fun-CIM", "Dys-CIM", "Dys-CIM",
                "Combined", "Combined"),
  format    = c("PNG", "PDF", "PNG", "PDF", "PNG", "PDF"),
  path      = c(fun_files[["png"]], fun_files[["pdf"]],
                dys_files[["png"]], dys_files[["pdf"]],
                combined_files[["png"]], combined_files[["pdf"]]),
  exists     = file.exists(path),
  size_bytes = file.info(path)$size
)
readr::write_csv(forest_output_manifest,
  file.path(subdirs["forests"], "forest_output_manifest.csv"))
if (any(!forest_output_manifest$exists))
  warning("One or more forest figure files were not written; ",
          "inspect forest_output_manifest.csv.")

# Render figures inline in the HTML report
fun_result$plot
dys_result$plot
combined_result$plot
}


# ==============================================================================
# OPTIONAL MEDIAN-SPLIT KAPLAN-MEIER ANALYSIS
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION: optional-km-median-split
# ------------------------------------------------------------------------------
# Optional descriptive KM and log-rank analysis
km_logrank_results <- tibble()
km_plot_manifest   <- tibble()
km_analysis_data   <- tibble()

if (run_km_median_split) {
  km_palette <- c("Low" = "#0072B2", "High" = "#CC0000")
  for (cohort in names(cohorts)) {
    obj      <- cohorts[[cohort]]
    clinical <- obj$clinical
    clinical$sample_id <- clinical[[obj$id_col]]
    for (i in seq_len(nrow(obj$endpoints))) {
      endpoint_row <- obj$endpoints[i, , drop = FALSE]
      for (signature in names(gene_sets)) {
        score <- ssgsea_scores %>%
          filter(.data$cohort == .env$cohort,
                 .data$signature == .env$signature) %>%
          select(sample_id, ssgsea_score, raw_median, median_group)
        raw_cutoff <- dplyr::first(score$raw_median)
        km_data <- clinical %>%
          left_join(score, by = "sample_id") %>%
          transmute(
            sample_id,
            time  = safe_numeric(.data[[endpoint_row$time_col]]),
            event = standardize_event(.data[[endpoint_row$event_col]],
                                      endpoint_row$event_col, cohort),
            ssgsea_score, raw_median,
            median_group = factor(median_group, levels = c("Low", "High"))
          ) %>%
          filter(is.finite(time), time > 0, event %in% c(0, 1),
                 !is.na(median_group))
        km_analysis_data <- bind_rows(
          km_analysis_data,
          km_data %>% transmute(
            cohort,
            endpoint_original   = endpoint_row$endpoint_original,
            endpoint_harmonized = endpoint_row$endpoint_harmonized,
            signature, raw_median_cutoff = raw_cutoff,
            sample_id, time, event, ssgsea_score, median_group))
        if (nrow(km_data) < min_n_model ||
            sum(km_data$event == 1) < min_events_model ||
            nlevels(droplevels(km_data$median_group)) < 2L) {
          km_logrank_results <- bind_rows(km_logrank_results, tibble(
            cohort,
            endpoint_original   = endpoint_row$endpoint_original,
            endpoint_harmonized = endpoint_row$endpoint_harmonized,
            signature, raw_median_cutoff = raw_cutoff,
            n = nrow(km_data), events = sum(km_data$event == 1),
            low_n  = sum(km_data$median_group == "Low"),
            high_n = sum(km_data$median_group == "High"),
            logrank_chisq = NA_real_, logrank_df = NA_integer_,
            logrank_p = NA_real_, status = "skipped",
            warning_text = "insufficient n, events, or median-group levels"))
          next
        }
        km_fit <- survival::survfit(
          survival::Surv(time, event) ~ median_group, data = km_data)
        logrank    <- survival::survdiff(
          survival::Surv(time, event) ~ median_group, data = km_data)
        logrank_df <- length(logrank$n) - 1L
        logrank_p  <- stats::pchisq(logrank$chisq, df = logrank_df,
                                    lower.tail = FALSE)
        endpoint_title <- ifelse(cohort == "SCANB", "SCAN-B", cohort)
        plot_title <- paste(signature, endpoint_title,
                            endpoint_row$endpoint_original, sep = " \u2014 ")
        max_time <- max(km_data$time, na.rm = TRUE)
        km_ylim <- if (cohort == "METABRIC") {
          metabric_km_ylim
        } else {
          scanb_km_ylim
        }
        km_pval_coord <- if (cohort == "METABRIC") {
          metabric_km_pval_coord
        } else {
          scanb_km_pval_coord
        }
        if (km_pval_coord[1] > max_time) {
          stop(cohort, " KM p-value x coordinate exceeds the observed follow-up; ",
               "change the corresponding YAML parameter.")
        }
        km_plot <- survminer::ggsurvplot(
          km_fit, data = km_data, palette = km_palette,
          conf.int = FALSE, censor.shape = "|", censor.size = 2,
          xlab = "Time (months)", ylab = "Survival probability",
          xlim = c(0, max_time), ylim = km_ylim,
          break.x.by = km_break_months,
          pval = TRUE, pval.coord = km_pval_coord,
          risk.table = TRUE,
          risk.table.title  = "Number at risk",
          risk.table.height = 0.25,
          legend = "top", legend.title = "",
          legend.labs = c("Low", "High"),
          title = plot_title,
          ggtheme      = theme_classic(base_size = 13),
          tables.theme = theme_classic(base_size = 11))
        file_tag <- paste(cohort, endpoint_row$endpoint_original,
                          str_replace_all(signature, "-", "_"),
                          "median_KM", sep = "__")
        png_file <- file.path(subdirs["km"], paste0(file_tag, ".png"))
        pdf_file <- file.path(subdirs["km"], paste0(file_tag, ".pdf"))
        grDevices::png(png_file, width = km_width, height = km_height,
                       units = "in", res = km_dpi)
        print(km_plot); grDevices::dev.off()
        if (capabilities("cairo")) {
          grDevices::cairo_pdf(pdf_file, width = km_width, height = km_height)
        } else {
          grDevices::pdf(pdf_file, width = km_width, height = km_height)
        }
        print(km_plot); grDevices::dev.off()
        km_logrank_results <- bind_rows(km_logrank_results, tibble(
          cohort,
          endpoint_original   = endpoint_row$endpoint_original,
          endpoint_harmonized = endpoint_row$endpoint_harmonized,
          signature, raw_median_cutoff = raw_cutoff,
          n = nrow(km_data), events = sum(km_data$event == 1),
          low_n       = sum(km_data$median_group == "Low"),
          high_n      = sum(km_data$median_group == "High"),
          low_events  = sum(km_data$event[km_data$median_group == "Low"]  == 1),
          high_events = sum(km_data$event[km_data$median_group == "High"] == 1),
          logrank_chisq = unname(logrank$chisq), logrank_df,
          logrank_p, status = "completed", warning_text = NA_character_))
        km_plot_manifest <- bind_rows(km_plot_manifest, tibble(
          cohort, endpoint_original = endpoint_row$endpoint_original,
          signature, png_file, pdf_file,
          png_exists = file.exists(png_file),
          pdf_exists = file.exists(pdf_file)))
      }
    }
  }
  km_logrank_results <- km_logrank_results %>%
    group_by(cohort, endpoint_original) %>%
    mutate(logrank_p_adjust_signature = ifelse(
      status == "completed", p.adjust(logrank_p, method = "BH"), NA_real_)) %>%
    ungroup()
  readr::write_csv(km_logrank_results,
    file.path(subdirs["km"], "KM_logrank_results.csv"))
  readr::write_csv(km_plot_manifest,
    file.path(subdirs["km"], "KM_plot_manifest.csv"))
  readr::write_csv(km_analysis_data,
    file.path(subdirs["km"], "KM_analysis_data.csv"))
}


# ==============================================================================
# WORKBOOK AND REPRODUCIBILITY OUTPUTS
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION: workbook-reproducibility
# ------------------------------------------------------------------------------
# Workbook and reproducibility exports
workbook <- openxlsx::createWorkbook()
add_sheet <- function(name, data) {
  name <- substr(gsub("[\\[\\]*?/\\\\:]", "_", name), 1, 31)
  if (!ncol(data)) data <- tibble(note = "No records generated")
  openxlsx::addWorksheet(workbook, name)
  openxlsx::writeDataTable(workbook, name, as.data.frame(data),
                            tableStyle = "TableStyleMedium2")
  openxlsx::freezePane(workbook, name, firstRow = TRUE)
  openxlsx::setColWidths(workbook, name,
                         cols = seq_len(ncol(data)), widths = "auto")
}
add_sheet("Signature results",      signature_results)
add_sheet("Continuous Cox",         continuous_results)
add_sheet("Median split Cox",       median_split_results)
add_sheet("All coefficients",       coefficient_results)
add_sheet("Model summary",          model_summary)
add_sheet("PH tests",               ph_results)
add_sheet("Covariate audit",        covariate_level_audit)
add_sheet("Model covariates",       model_covariate_audit)
add_sheet("Covariate inclusion",    model_covariate_inclusion)
add_sheet("Gene sets",              gene_definitions)
add_sheet("Coverage",               coverage)
add_sheet("Score QC",               score_qc)
add_sheet("Median cutoffs",         median_cutoffs)
add_sheet("Median group balance",   median_group_balance)
add_sheet("Score group assignments", median_assignment_audit)
if (run_km_median_split) {
  add_sheet("KM logrank",       km_logrank_results)
  add_sheet("KM manifest",      km_plot_manifest)
  add_sheet("KM analysis data", km_analysis_data)
}
add_sheet("Model definitions",
  model_definitions %>% mutate(covariates = map_chr(covariates,
    ~ paste(as.character(unlist(.x, recursive = TRUE)), collapse = "; "))))
parameter_text <- function(x) {
  if (is.null(x) || !length(x) || all(is.na(x))) return("null")
  paste(as.character(unlist(x, recursive = TRUE, use.names = FALSE)),
        collapse = "; ")
}
parameter_values <- list(
  gene_set_source = if (use_custom_gene_sets) "custom gene lists from params"
    else "cimic meta_results.csv",
  signature_file = params$signature_file,
  signature_file_md5 = if (use_custom_gene_sets) NA_character_ else
    unname(tools::md5sum(files[["correlation"]])),
  expected_signature_file_md5 = params$expected_signature_file_md5,
  meta_scope = meta_scope,
  meta_rho_threshold = meta_rho_threshold,
  meta_k_mode = meta_k_mode,
  meta_k_mode_description = meta_k_mode_description,
  meta_require_fisher = meta_require_fisher,
  meta_require_stouffer = meta_require_stouffer,
  meta_padj_threshold = meta_padj_threshold,
  effective_meta_significance_rule = meta_significance_rule,
  meta_min_mean_abs_rrb = meta_min_mean_abs_rrb,
  meta_top_n_genes = meta_top_n_genes,
  meta_rank_metric = meta_rank_metric,
  meta_ambiguous_gene_policy = meta_ambiguous_gene_policy,
  ambiguous_genes_found = length(ambiguous_genes),
  effective_fun_cim_rule = fun_selection_rule,
  effective_dys_cim_rule = dys_selection_rule,
  expected_fun_cim_genes = expected_fun_cim_genes,
  expected_dys_cim_genes = expected_dys_cim_genes,
  selected_fun_cim_genes = observed_signature_counts[["Fun-CIM"]],
  selected_dys_cim_genes = observed_signature_counts[["Dys-CIM"]],
  cox_fdr_threshold = cox_fdr_threshold,
  cox_fdr_operator = cox_fdr_operator,
  minimum_signature_coverage = minimum_signature_coverage,
  minimum_signature_genes = minimum_signature_genes,
  duplicate_gene_method = duplicate_gene_method,
  missing_expression_gene_policy = missing_expression_gene_policy,
  ssgsea_alpha = ssgsea_alpha,
  ssgsea_normalize = ssgsea_normalize,
  ssgsea_min_size = ssgsea_min_size,
  ssgsea_max_size = ssgsea_max_size,
  bioc_parallel_workers = bioc_parallel_workers,
  score_parameterizations = score_parameterizations,
  median_scope = "cohort_signature",
  median_tie_rule = median_tie_rule,
  median_group_warning_fraction = median_group_warning_fraction,
  age_unit = "5 years",
  score_unit = "1 within-cohort SD",
  median_split_unit = "High versus Low",
  min_n_model = min_n_model,
  min_events_model = min_events_model,
  events_per_parameter_target = events_per_parameter_target,
  enforce_events_per_parameter = enforce_events_per_parameter,
  sparse_factor_level_n = sparse_factor_level_n,
  alpha = alpha_threshold,
  significance_metric = significance_metric,
  m2_covariates = m2_covariates,
  primary_adjusted_model = primary_adjusted_model,
  run_common_sample_sensitivity = run_common_sample_sensitivity,
  metabric_chemo_only = metabric_chemo_only,
  scanb_chemo_only = scanb_chemo_only,
  metabric_chemo_column = params$metabric_chemo_column,
  scanb_chemo_column = params$scanb_chemo_column,
  chemo_positive_values = chemo_positive_values,
  metabric_gene_start_col = metabric_gene_start_col,
  run_km_median_split = run_km_median_split,
  km_break_months = km_break_months,
  metabric_km_ylim = metabric_km_ylim,
  scanb_km_ylim = scanb_km_ylim,
  metabric_km_pval_coord = metabric_km_pval_coord,
  scanb_km_pval_coord = scanb_km_pval_coord,
  km_width = km_width,
  km_height = km_height,
  km_dpi = km_dpi,
  make_forest_plots = make_forest_plots,
  forest_model = forest_model,
  forest_parameterization = forest_parameterization,
  hr_axis_mode = hr_axis_mode,
  hr_axis_min = hr_axis_min,
  hr_axis_max = hr_axis_max,
  hr_axis_breaks = manual_hr_breaks,
  forest_width = forest_width,
  forest_height = forest_height,
  combined_forest_width = combined_forest_width,
  score_qc_plot_width = score_qc_plot_width,
  score_qc_plot_height = score_qc_plot_height,
  figure_dpi = figure_dpi,
  seed = analysis_seed,
  output_dir = output_dir,
  parameter_source = analysis_parameter_source
)
analysis_parameters <- tibble(
  parameter = names(parameter_values),
  value = map_chr(parameter_values, parameter_text))
add_sheet("Analysis parameters", analysis_parameters)
openxlsx::saveWorkbook(workbook,
  file.path(output_dir, "CIM_ssGSEA_survival_results.xlsx"), overwrite = TRUE)
readr::write_csv(analysis_parameters,
  file.path(subdirs["repro"], "analysis_parameters.csv"))
writeLines(capture.output(sessionInfo()),
  file.path(subdirs["repro"], "sessionInfo.txt"))
writeLines(capture.output(tools::md5sum(files)),
  file.path(subdirs["repro"], "input_md5_checksums.txt"))

completion_check <- tibble(
  check = c(
    "both_signatures", "all_models_attempted",
    "both_parameterizations_attempted", "all_endpoints_attempted",
    "one_cutoff_per_cohort_signature", "median_assignments_valid",
    "coefficient_table_written", "model_summary_written", "workbook_written",
    "forest_figures_written", "km_analysis"
  ),
  passed = c(
    all(c("Fun-CIM", "Dys-CIM") %in% unique(signature_results$signature)),
    all(model_levels %in% unique(as.character(signature_results$model_id))),
    all(score_parameterizations %in% unique(signature_results$parameterization)),
    all(c("OS", "DSS", "RFI", "DRFi", "RFi") %in%
          unique(signature_results$endpoint_original)),
    nrow(median_cutoffs) == length(names(cohorts)) * length(names(gene_sets)) &&
      !anyDuplicated(median_cutoffs[c("cohort", "signature")]),
    all(median_assignment_audit$assignment_valid),
    file.exists(file.path(subdirs["models"], "complete_coefficient_results.csv")),
    file.exists(file.path(subdirs["models"], "model_level_summary.csv")),
    file.exists(file.path(output_dir, "CIM_ssGSEA_survival_results.xlsx")),
    if (make_forest_plots) {
      nrow(forest_output_manifest) == 6L && all(forest_output_manifest$exists)
    } else TRUE,
    if (run_km_median_split) {
      nrow(km_logrank_results) > 0 &&
        file.exists(file.path(subdirs["km"], "KM_logrank_results.csv"))
    } else TRUE
  ),
  detail = c(
    "Fun-CIM and Dys-CIM",
    paste(model_levels, collapse = "; "),
    paste(score_parameterizations, collapse = "; "),
    "OS, DSS, RFI, DRFi, and RFi",
    "Cutoffs fixed within cohort and signature",
    "All patient assignments satisfy the fixed median rule",
    "Required", "Required", "Required",
    if (make_forest_plots) {
      paste("Fun-CIM, Dys-CIM, and Combined", forest_model,
            forest_parameterization, "forest PNGs and PDFs")
    } else "not requested",
    if (run_km_median_split) "requested" else "not requested"
  )
)
readr::write_csv(completion_check,
  file.path(subdirs["repro"], "completion_check.csv"))
if (any(!completion_check$passed))
  warning("One or more completion checks failed; inspect completion_check.csv.")


# ==============================================================================
# INTERPRETATION SAFEGUARDS
# ==============================================================================

# - Continuous hazard ratios are per one within-cohort standard deviation.
#   Median-split hazard ratios compare High with Low using fixed raw-score medians
#   calculated separately within each cohort and signature.
# - Median cutoffs are not pooled across cohorts and are not recalculated across
#   endpoints, adjusted models, available-case populations, or common samples.
# - M0 is the unadjusted univariable reference. M2 is the primary adjusted model.
#   No intermediate models are fitted in this analysis.
# - The common-sample sensitivity analysis uses the configured primary-adjusted
#   model complete-case population as its anchor, fitting all configured models
#   to the same set of samples for direct
#   comparison of unadjusted and adjusted estimates.
# - Available-case estimates can change because both covariate adjustment and
#   sample composition change. The common-sample sensitivity analysis separates
#   these effects by fitting all models to the primary-model-complete population.
# - Raw p-values, signature-family BH-adjusted p-values, and within-model
#   term-adjusted p-values are retained as separate columns. Continuous and
#   median-split signature p-values are adjusted in separate families.
# - The publication forest plots display the model and parameterization selected
#   by `forest_model` and `forest_parameterization`. All fitted results are
#   retained in the CSV and Excel outputs.
# - Kaplan-Meier and log-rank analyses are optional descriptive analyses and do
#   not run unless `run_km_median_split` is explicitly changed to `true`.
# - RFI/RFi provide the direct recurrence-focused cross-cohort comparison. DSS
#   and DRFi are reported as cohort-specific endpoints rather than equivalents.

