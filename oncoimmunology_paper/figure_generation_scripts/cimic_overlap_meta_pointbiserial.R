# =============================================================================
# cimic_overlap_meta_pointbiserial.R
#
# Cross-dataset analysis of the rank-biserial / point-biserial (r_rb) gene-vs-
# cluster correlations for the Fun-CIM and Dys-CIM signatures, using the three
# per-dataset correlation tables in oncoimmunology_paper/Results/ produced by
# the *_correlation_analysis.R scripts.
#
# Adapted from overlap_meta_rankbiserial.R (9_Projects/combined_cimic_analysis_bc
# 9-dataset version). The statistical logic is carried over unchanged; only the
# dataset registry, the CIM-direction derivation, and the analysis scopes differ.
#
# THREE complementary analyses, all per CIM (Fun / Dys) and at |r_rb| >= 0.3/0.6:
#
#   (1) META-ANALYSIS on RAW p  -- combine each gene's one-sided raw p across
#       datasets into a meta p-value, to surface aggregate signal even where no
#       single dataset clears threshold.
#         * Fisher's method  (PRIMARY)      : X2 = -2 * sum ln(p_one_sided), df=2k
#         * Stouffer's Z     (SECONDARY)    : signed, sqrt(n)-weighted
#       Directional: p_one_sided is taken TOWARD the CIM of interest, so a gene
#       that flips direction between datasets cancels out (unlike a naive Fisher
#       on two-sided p). A dataset contributes to a gene at a given rho level
#       only where |r_rb| >= cutoff; >= 2 contributing datasets required.
#       Meta p-values are BH-adjusted ACROSS GENES (the single, correct point of
#       multiplicity control -- inputs are RAW p, never BH-adjusted).
#
#   (2) CROSS-DATASET OVERLAP on BH-adjusted significance -- a gene counts as
#       significant in a dataset if padj < 0.05 there; then frequency (shared by
#       >= k datasets) + group intersections.
#
#   (3) CROSS-DATASET OVERLAP on RAW p -- significant if raw p < 0.05, and
#       separately < 0.01; same frequency + intersections.
#
# ----------------------------------------------------------------------------
# CIM DIRECTION -- IMPORTANT
# ----------------------------------------------------------------------------
# The per-dataset tables here carry no `direction` column (unlike the 9-dataset
# version), so direction is recovered from the SIGN of r_rb. All three upstream
# scripts fit `model.matrix(~ clusterf)` and read `coef = 2`, i.e. Cluster 2 vs
# Cluster 1, so:
#
#     r_rb > 0  ->  higher in Cluster 2
#     r_rb < 0  ->  higher in Cluster 1
#
# Combined with the per-dataset cluster coding:
#
#     nki_smc   C1 = Dys-CIM, C2 = Fun-CIM   ->  r_rb > 0 = Fun,  r_rb < 0 = Dys
#     neo       C1 = Dys-CIM, C2 = Fun-CIM   ->  r_rb > 0 = Fun,  r_rb < 0 = Dys
#     tnbc_cl   C1 = Fun-CIM, C2 = Dys-CIM   ->  r_rb > 0 = Dys,  r_rb < 0 = Fun
#
# This reproduces the `trajectory` labels already written by each upstream
# script, which is the cross-check that the mapping below is right. The
# `trajectory` column itself is NOT used: it bakes in a per-script rho cutoff
# (0.6 / 0.5 / 0.5) and a padj filter, whereas this script re-derives CIM
# membership at the requested RHO_CUTS.
#
# ----------------------------------------------------------------------------
# meta_results.csv COLUMNS
# ----------------------------------------------------------------------------
# A machine-readable copy of this table is written to
# meta_results_data_dictionary.csv next to the data.
#
#   scope                     Which dataset combination the row was computed over.
#   cim                       CIM the test was directed TOWARD (Fun or Dys).
#                             Every gene appears twice per scope+rho, once per
#                             CIM; the p-values differ because the one-sided p
#                             is taken toward this CIM.
#   cim_class                 Verdict across BOTH cim rows (Fisher, primary):
#                             Fun_only / Dys_only / Both / none.
#   rho_threshold             The |r_rb| cutoff a dataset had to clear.
#   gene                      Gene symbol.
#   n_scope_datasets          Datasets in this scope (3 for all-three, 2 pairs).
#
#   k_datasets                Number of datasets where |r_rb| >= rho_threshold.
#                             *** MAGNITUDE ONLY. *** Does NOT require
#                             significance in the individual dataset, and does
#                             NOT require the datasets to agree on direction.
#                             Identical between a gene's Fun and Dys rows.
#   all_datasets_contributed  YES when k_datasets == n_scope_datasets, i.e. the
#                             gene cleared rho in EVERY dataset of the scope.
#                             This is the "require all three" filter.
#   k_concordant              Of those contributors, how many point toward `cim`.
#                             This is the DIRECTION-agreement count.
#   direction_concordant      YES when k_concordant == k_datasets.
#                             Combine with all_datasets_contributed to get
#                             "same direction in all three".
#
#   <dataset>                 YES/NO -- did that dataset contribute?
#   <dataset>_dir             Fun/Dys -- which CIM its r_rb points to;
#                             NA where it did not contribute.
#
#   fisher_stat               -2 * sum(log(p_one_sided)) over contributors.
#   fisher_p                  PRIMARY meta p. Chi-square tail, 2*k df. Raw.
#   fisher_padj               fisher_p, BH-adjusted ACROSS GENES within
#                             scope+cim+rho. The single multiplicity control.
#   fisher_sig                YES when fisher_padj < META_ALPHA.
#   stouffer_Z                SECONDARY. sqrt(n)-weighted one-sided z toward
#                             `cim`. Signed: positive favours `cim`.
#   stouffer_p                One-sided normal tail of stouffer_Z. Raw.
#   stouffer_padj             BH-adjusted across genes within scope+cim+rho.
#   stouffer_sig              YES when stouffer_padj < META_ALPHA.
#   mean_abs_rrb              Mean |r_rb| over contributors (3 dp).
#
# Outputs -> oncoimmunology_paper/Results/cimic_meta_analysis/ :
#   meta_results.csv                 (per gene: fisher/stouffer p + BH padj)
#   meta_results_data_dictionary.csv (column meanings for the above)
#   meta_cim_class_summary.csv       (Fun_only / Dys_only / Both gene counts)
#   meta_significant_summary.csv     (counts of meta-significant genes)
#   overlap_summary_counts.csv       (intersection/frequency gene counts)
#   overlap_group_intersections.csv  (intersection gene lists)
#   overlap_frequency_per_gene.csv   (per-gene dataset frequency)
#   overlap_cim_class_summary.csv    (Fun_only / Dys_only / Both gene counts)
#   overlap_frequency_summary.csv    (shared-by->=k breakdown)
#
# Run from the repository root:
#   Rscript oncoimmunology_paper/cimic_overlap_meta_rankbiserial.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(purrr)
  library(tibble)
})

# ---- 0. Config --------------------------------------------------------------
# Repository-relative paths; run from the repo root.
base_dir <- "oncoimmunology_paper"
out_dir  <- file.path(base_dir, "Results", "cimic_meta_analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

RHO_CUTS   <- c(0.3)
META_ALPHA <- 0.05          # BH-padj cutoff for calling a meta gene significant

# ---- 1. Dataset registry ----------------------------------------------------
# corr = point/rank-biserial correlation table from *_correlation_analysis.R
# ind  = master induction df (used only to recover n via sd/sem for Stouffer)
# map  = which CIM the Cluster 1 / Cluster 2 labels mean in that dataset
datasets <- list(
  nki_smc = list(
    corr = file.path(base_dir, "Results", "nki_smc_correlation_analysis.csv"),
    ind  = file.path(base_dir, "Datasets", "NKI_SMC",
                     "nki_smc_master_induction_df.csv"),
    map  = c("Cluster 1" = "Dys", "Cluster 2" = "Fun")),
  neo = list(
    corr = file.path(base_dir, "Results", "neo_correlation_analysis.csv"),
    ind  = file.path(base_dir, "Datasets", "NEO",
                     "neo_master_induction_df.csv"),
    map  = c("Cluster 1" = "Dys", "Cluster 2" = "Fun")),
  tnbc_cl = list(
    corr = file.path(base_dir, "Results", "tnbc_cl_correlation_analysis.csv"),
    ind  = file.path(base_dir, "Datasets", "TNBC_CL_Epirubicin",
                     "tnbc_cl_epirubicin_master_induction_df.csv"),
    map  = c("Cluster 1" = "Fun", "Cluster 2" = "Dys"))
)
ALL_DATASETS <- names(datasets)

# ---- 2. Analysis scopes -----------------------------------------------------
# The four contexts available with these three datasets. Names spell out members.
analysis_scopes <- list(
  ALL3_NKI.NEO.TNBCCL = c("nki_smc", "neo", "tnbc_cl"),
  NKI.NEO             = c("nki_smc", "neo"),
  NKI.TNBCCL          = c("nki_smc", "tnbc_cl"),
  NEO.TNBCCL          = c("neo", "tnbc_cl")
)
# Meta-analysis scopes, overlap intersection groups, and frequency universes are
# all the same four contexts here (unlike the 9-dataset version, where natural
# and custom groupings differed).
meta_scopes    <- analysis_scopes
overlap_groups <- analysis_scopes
freq_universes <- analysis_scopes

# ---- 3. Recover sample sizes (for Stouffer weights) from master induction ---
# n is constant across genes; recover from sem = sd/sqrt(n) (mode over finite).
recover_n <- function(ind_path) {
  if (!file.exists(ind_path)) return(c(n1 = NA, n2 = NA, n = NA))
  m <- data.table::fread(ind_path, showProgress = FALSE, data.table = FALSE)
  mode_n <- function(sd, sem) {
    v <- round((sd / sem)^2); v <- v[is.finite(v) & v > 0]
    if (!length(v)) return(NA_integer_)
    as.integer(names(which.max(table(v))))
  }
  if (!all(c("cluster1_sd", "cluster1_sem", "cluster2_sd", "cluster2_sem")
           %in% names(m)))
    return(c(n1 = NA, n2 = NA, n = NA))
  n1 <- mode_n(m$cluster1_sd, m$cluster1_sem)
  n2 <- mode_n(m$cluster2_sd, m$cluster2_sem)
  c(n1 = n1, n2 = n2, n = ifelse(is.na(n1) | is.na(n2), NA, n1 + n2))
}
n_map <- vapply(ALL_DATASETS,
                function(d) unname(recover_n(datasets[[d]]$ind)["n"]),
                numeric(1))

# ---- Effective-n override for cell lines (Stouffer weighting) ---------------
# The dupcor LIMMA analysis in tnbc_cl_correlation_analysis.R treats the CELL
# LINE (not the 27 replicate samples) as the unit of replication. recover_n()
# reads sd/sem computed over the 27 samples (-> 27), which would over-weight the
# cell-line panel in the Stouffer meta. Override to the panel size (n = 9 cell
# lines x 3 replicates) so the weight reflects the effective evidence. Patient
# cohorts keep their sample n.
CELL_LINE_DATASETS <- c("tnbc_cl")
CELL_LINE_EFF_N    <- 9
n_map[base::intersect(names(n_map), CELL_LINE_DATASETS)] <- CELL_LINE_EFF_N

message("Sample sizes (n) per dataset (cell lines overridden to effective n=",
        CELL_LINE_EFF_N, "):")
print(n_map)

# ---- 4. Load + annotate all correlation tables ------------------------------
# Gene column is `gene_id` here (was `gene` in the 9-dataset tables).
# direction is derived from sign(r_rb) -- see the header note.
load_one <- function(name) {
  fp <- datasets[[name]]$corr
  if (!file.exists(fp))
    stop("Missing correlation table (run the *_correlation_analysis.R script ",
         "first): ", fp)
  d <- data.table::fread(fp, showProgress = FALSE, data.table = FALSE)

  gene_col <- base::intersect(c("gene_id", "gene"), names(d))[1]
  if (is.na(gene_col)) stop("No gene_id/gene column in ", fp)
  names(d)[names(d) == gene_col] <- "gene"

  needed <- c("gene", "r_rb", "p_value", "padj")
  if (!all(needed %in% names(d)))
    stop("Missing column(s) in ", fp, ": ",
         paste(base::setdiff(needed, names(d)), collapse = ", "))

  d %>%
    dplyr::mutate(
      dataset   = name,
      # sign(r_rb) -> cluster -> CIM, via the per-dataset mapping
      direction = ifelse(r_rb > 0, "Cluster 2", "Cluster 1"),
      cim       = unname(datasets[[name]]$map[direction]),
      abs_rrb   = abs(r_rb)
    ) %>%
    # r_rb == 0 has no direction; it also never clears any rho cutoff.
    dplyr::filter(!is.na(cim), !is.na(p_value), is.finite(r_rb), r_rb != 0) %>%
    dplyr::select(dataset, gene, cim, r_rb, abs_rrb, p_value, padj)
}
message("\nLoading ", length(ALL_DATASETS), " correlation tables ...")
all_data <- purrr::map_dfr(ALL_DATASETS, load_one)
message("Loaded ", nrow(all_data), " gene-rows.")

# Sanity check: the derived CIM assignment must agree with the sign convention.
qc <- all_data %>%
  dplyr::count(dataset, cim, sign_rrb = sign(r_rb)) %>%
  dplyr::arrange(dataset, cim)
message("\nCIM assignment by sign(r_rb) [expect one sign per dataset+cim]:")
print(as.data.frame(qc))

# Gene universe per dataset (used by the ORA script's universe options).
gene_universe <- all_data %>%
  dplyr::distinct(dataset, gene) %>%
  dplyr::count(dataset, name = "n_genes")
message("\nGenes per dataset:")
print(as.data.frame(gene_universe))

# =============================================================================
# ANALYSIS 1 -- META-ANALYSIS (Fisher primary + Stouffer secondary) on RAW p
# =============================================================================
CLAMP <- 1e-300
meta_one <- function(members, cim_target, rho) {
  w <- all_data %>%
    dplyr::filter(dataset %in% members, cim %in% c("Fun", "Dys"), abs_rrb >= rho)
  if (nrow(w) == 0) return(NULL)
  # one-sided raw p toward the CIM of interest; z toward the CIM (signed).
  w <- w %>%
    dplyr::mutate(
      p_os = ifelse(cim == cim_target, p_value / 2, 1 - p_value / 2),
      p_os = pmin(pmax(p_os, CLAMP), 1 - 1e-15),
      z_os = stats::qnorm(1 - pmin(pmax(p_os, 1e-15), 1 - 1e-15)),
      w_i  = { nn <- n_map[dataset]; ifelse(is.finite(nn), sqrt(nn), 1) }
    )
  w %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      k_datasets   = n(),
      datasets     = paste(sort(dataset), collapse = ";"),
      # k_concordant / concordant_datasets: of the contributing datasets, how
      # many actually point TOWARD cim_target. k_datasets counts magnitude only
      # (|r_rb| >= rho, either direction); these two count direction agreement.
      # A gene with k_datasets = 3 but k_concordant = 2 clears the rho cutoff in
      # all three datasets while disagreeing on direction in one of them.
      k_concordant = sum(cim == cim_target),
      concordant_datasets = paste(sort(dataset[cim == cim_target]),
                                  collapse = ";"),
      fisher_stat  = -2 * sum(log(p_os)),
      fisher_p     = stats::pchisq(fisher_stat, df = 2 * n(), lower.tail = FALSE),
      stouffer_Z   = sum(w_i * z_os) / sqrt(sum(w_i^2)),
      stouffer_p   = stats::pnorm(stouffer_Z, lower.tail = FALSE),
      mean_abs_rrb = round(mean(abs_rrb), 3),
      .groups = "drop"
    ) %>%
    dplyr::filter(k_datasets >= 2) %>%
    dplyr::mutate(fisher_padj   = stats::p.adjust(fisher_p,   method = "BH"),
           stouffer_padj = stats::p.adjust(stouffer_p, method = "BH"),
           scope = NA_character_, cim = cim_target, rho_threshold = rho,
           # Number of datasets in this scope -- lets a downstream script express
           # "contributed in ALL datasets of the scope" without hard-coding 2 vs 3.
           n_scope_datasets = length(members))
}

message("\n== Analysis 1: meta-analysis (Fisher + Stouffer) ==")
meta_rows <- list()
for (sc_name in names(meta_scopes)) {
  for (cim_t in c("Fun", "Dys")) {
    for (rc in RHO_CUTS) {
      res <- meta_one(meta_scopes[[sc_name]], cim_t, rc)
      if (!is.null(res) && nrow(res) > 0) {
        res$scope <- sc_name
        meta_rows[[length(meta_rows) + 1]] <- res
      }
    }
  }
}
if (!length(meta_rows))
  stop("Meta-analysis produced no rows -- check RHO_CUTS and the input tables.")

meta_results <- dplyr::bind_rows(meta_rows) %>%
  dplyr::arrange(scope, cim, rho_threshold, fisher_p)
# Per-dataset columns, two per dataset:
#   <ds>      YES/NO -- did that dataset CONTRIBUTE to this gene's meta signal,
#                       i.e. |r_rb| >= cutoff within the scope? (magnitude only)
#   <ds>_dir  Fun/Dys/NA -- which CIM that dataset's r_rb actually points to.
#                       NA where the dataset did not contribute.
# Exact token match on the ";"-padded lists avoids substring false positives.
ds_pad    <- paste0(";", meta_results$datasets, ";")
cd_pad    <- paste0(";", meta_results$concordant_datasets, ";")
other_cim <- ifelse(meta_results$cim == "Fun", "Dys", "Fun")
for (d in ALL_DATASETS) {
  contributed <- grepl(paste0(";", d, ";"), ds_pad, fixed = TRUE)
  concordant  <- grepl(paste0(";", d, ";"), cd_pad, fixed = TRUE)
  meta_results[[d]] <- ifelse(contributed, "YES", "NO")
  # A contributing dataset either agrees with cim (it is in concordant_datasets)
  # or, by construction, points to the other CIM.
  meta_results[[paste0(d, "_dir")]] <- ifelse(
    !contributed, NA_character_,
    ifelse(concordant, meta_results$cim, other_cim))
}
meta_results <- meta_results %>%
  dplyr::mutate(fisher_sig   = ifelse(fisher_padj   < META_ALPHA, "YES", "NO"),
         stouffer_sig = ifelse(stouffer_padj < META_ALPHA, "YES", "NO"),
         # Did every contributing dataset agree on direction?
         direction_concordant = ifelse(k_concordant == k_datasets, "YES", "NO"),
         # Did the gene clear the rho cutoff in EVERY dataset of the scope?
         all_datasets_contributed =
           ifelse(k_datasets == n_scope_datasets, "YES", "NO"))
# CIM class within each scope+rho: is the gene meta-significant (Fisher, primary)
# toward Fun ONLY, Dys ONLY, BOTH, or neither ("none")?
meta_cls <- meta_results %>%
  dplyr::group_by(scope, rho_threshold, gene) %>%
  dplyr::summarise(has_fun = any(cim == "Fun" & fisher_sig == "YES"),
            has_dys = any(cim == "Dys" & fisher_sig == "YES"),
            .groups = "drop") %>%
  dplyr::mutate(cim_class = dplyr::case_when(has_fun & has_dys ~ "Both",
                                      has_fun ~ "Fun_only",
                                      has_dys ~ "Dys_only",
                                      TRUE ~ "none")) %>%
  dplyr::select(scope, rho_threshold, gene, cim_class)
# Interleave each dataset's contribution flag with its direction, so the pair
# reads together in the CSV: nki_smc, nki_smc_dir, neo, neo_dir, ...
ds_cols <- as.vector(rbind(ALL_DATASETS, paste0(ALL_DATASETS, "_dir")))

meta_results <- meta_results %>%
  dplyr::left_join(meta_cls, by = c("scope", "rho_threshold", "gene")) %>%
  dplyr::select(scope, cim, cim_class, rho_threshold, gene,
         n_scope_datasets, k_datasets, all_datasets_contributed,
         k_concordant, direction_concordant,
         dplyr::all_of(ds_cols),
         fisher_stat, fisher_p, fisher_padj, fisher_sig,
         stouffer_Z, stouffer_p, stouffer_padj, stouffer_sig, mean_abs_rrb)
data.table::fwrite(meta_results, file.path(out_dir, "meta_results.csv"),
                   na = "NA")

# ---- Data dictionary for meta_results.csv -----------------------------------
# Written alongside the data so the column semantics travel with the file.
# Kept in sync with the dplyr::select() above; see also the header block of this script.
meta_dict <- tibble::tribble(
  ~column, ~type, ~meaning,

  "scope", "character",
  paste0("Which combination of datasets this row was computed over. One of: ",
         paste(names(meta_scopes), collapse = ", "), "."),

  "cim", "character",
  paste0("The CIM the meta-test was directed TOWARD for this row (Fun or Dys). ",
         "Every gene appears twice per scope+rho, once per CIM; the p-values ",
         "differ because the one-sided p is taken toward this CIM."),

  "cim_class", "character",
  paste0("Per scope+rho verdict for the gene across BOTH cim rows, using the ",
         "primary (Fisher) test: Fun_only / Dys_only / Both / none."),

  "rho_threshold", "numeric",
  paste0("The |r_rb| cutoff a dataset had to clear to contribute to this gene. ",
         "One of ", paste(RHO_CUTS, collapse = ", "), "."),

  "gene", "character", "Gene symbol.",

  "n_scope_datasets", "integer",
  "How many datasets are in this scope (3 for the all-three scope, 2 for pairs).",

  "k_datasets", "integer",
  paste0("Number of datasets in the scope where |r_rb| >= rho_threshold. ",
         "MAGNITUDE ONLY -- this does NOT require significance within the ",
         "dataset and does NOT require the datasets to agree on direction. ",
         "Identical between the Fun and Dys rows of the same gene."),

  "all_datasets_contributed", "character",
  paste0("YES when k_datasets == n_scope_datasets, i.e. the gene cleared the ",
         "rho cutoff in EVERY dataset of the scope. This is the ",
         "'require all three' filter."),

  "k_concordant", "integer",
  paste0("Of the k_datasets contributors, how many actually point toward `cim`. ",
         "This is the DIRECTION-agreement count, as distinct from k_datasets."),

  "direction_concordant", "character",
  paste0("YES when k_concordant == k_datasets, i.e. every contributing dataset ",
         "agrees on direction. Combine with all_datasets_contributed for ",
         "'same direction in all three'."),

  "<dataset>", "character",
  paste0("One column per dataset (", paste(ALL_DATASETS, collapse = ", "),
         "). YES if that dataset contributed (|r_rb| >= rho_threshold)."),

  "<dataset>_dir", "character",
  paste0("One column per dataset. Which CIM that dataset's r_rb points to ",
         "(Fun or Dys), or NA where the dataset did not contribute."),

  "fisher_stat", "numeric",
  "Fisher combined statistic: -2 * sum(log(p_one_sided)) over contributors.",

  "fisher_p", "numeric",
  paste0("PRIMARY meta p-value. Chi-square tail of fisher_stat on 2*k_datasets ",
         "df. Raw, unadjusted."),

  "fisher_padj", "numeric",
  paste0("fisher_p, BH-adjusted ACROSS GENES within this scope+cim+rho. This is ",
         "the single point of multiplicity control."),

  "fisher_sig", "character",
  paste0("YES when fisher_padj < ", META_ALPHA, "."),

  "stouffer_Z", "numeric",
  paste0("SECONDARY meta statistic. sqrt(n)-weighted sum of one-sided z toward ",
         "`cim`, normalised by sqrt(sum of squared weights). Signed: positive ",
         "favours `cim`, negative favours the other CIM. Cell-line weights use ",
         "an effective n of ", CELL_LINE_EFF_N, " (cell lines, not samples)."),

  "stouffer_p", "numeric", "One-sided normal tail of stouffer_Z. Raw.",

  "stouffer_padj", "numeric",
  "stouffer_p, BH-adjusted across genes within this scope+cim+rho.",

  "stouffer_sig", "character",
  paste0("YES when stouffer_padj < ", META_ALPHA, "."),

  "mean_abs_rrb", "numeric",
  paste0("Mean |r_rb| over the contributing datasets, rounded to 3 dp. Used as ",
         "the effect-size axis in the volcano and as a ranking tie-break.")
)
data.table::fwrite(meta_dict,
                   file.path(out_dir, "meta_results_data_dictionary.csv"))

# Small summary: unique genes meta-significant (Fisher) toward Fun_only /
# Dys_only / Both per scope + rho (excludes 'none').
meta_cim_class_summary <- meta_results %>%
  dplyr::filter(cim_class != "none") %>%
  dplyr::distinct(scope, rho_threshold, gene, cim_class) %>%
  dplyr::count(scope, rho_threshold, cim_class, name = "n_genes") %>%
  dplyr::arrange(scope, rho_threshold, cim_class)
data.table::fwrite(meta_cim_class_summary,
                   file.path(out_dir, "meta_cim_class_summary.csv"), na = "NA")

meta_significant_summary <- meta_results %>%
  dplyr::group_by(scope, cim, rho_threshold) %>%
  dplyr::summarise(n_tested       = n(),
            n_fisher_sig   = sum(fisher_padj   < META_ALPHA, na.rm = TRUE),
            n_stouffer_sig = sum(stouffer_padj < META_ALPHA, na.rm = TRUE),
            n_both_sig     = sum(fisher_padj < META_ALPHA &
                                 stouffer_padj < META_ALPHA, na.rm = TRUE),
            .groups = "drop") %>%
  dplyr::arrange(scope, cim, rho_threshold)
data.table::fwrite(meta_significant_summary,
                   file.path(out_dir, "meta_significant_summary.csv"), na = "NA")

# =============================================================================
# ANALYSIS 2 & 3 -- CROSS-DATASET OVERLAP (BH  and  raw-p 0.05 / 0.01)
# =============================================================================
# Significance modes applied per dataset.
sig_modes <- list(
  padj   = list(col = "padj",    cut = 0.05),   # analysis 2 (BH)
  rawp05 = list(col = "p_value", cut = 0.05),   # analysis 3
  rawp01 = list(col = "p_value", cut = 0.01)    # analysis 3
)

sig_genes <- function(ds, cim_t, mode, rho) {
  m <- sig_modes[[mode]]
  all_data %>%
    dplyr::filter(dataset == ds, cim == cim_t, abs_rrb >= rho,
           .data[[m$col]] < m$cut) %>%
    dplyr::pull(gene) %>% unique()
}

message("\n== Analysis 2 & 3: cross-dataset overlap ==")
freq_rows <- list(); group_rows <- list(); ov_summary <- list()

add_ov_summary <- function(analysis, mode, group, members, cim, rho, n) {
  ov_summary[[length(ov_summary) + 1]] <<- tibble::tibble(
    analysis = analysis, sig_mode = mode, group = group,
    group_members = paste(members, collapse = ";"), n_members = length(members),
    cim = cim, rho_threshold = rho, n_genes = n)
}

for (mode in names(sig_modes)) {
  for (rc in RHO_CUTS) {
    # per (dataset, cim) significant gene sets for this mode+rho
    sets <- setNames(lapply(ALL_DATASETS, function(ds)
      list(Fun = sig_genes(ds, "Fun", mode, rc),
           Dys = sig_genes(ds, "Dys", mode, rc))), ALL_DATASETS)

    for (cim_t in c("Fun", "Dys")) {
      # ---- frequency across each universe ----
      for (uname in names(freq_universes)) {
        ud <- freq_universes[[uname]]
        per <- dplyr::bind_rows(lapply(ud, function(ds) {
          g <- sets[[ds]][[cim_t]]
          if (!length(g)) return(NULL)
          all_data %>% dplyr::filter(dataset == ds, cim == cim_t, gene %in% g) %>%
            dplyr::select(dataset, gene, abs_rrb)
        }))
        if (nrow(per) > 0) {
          fr <- per %>% dplyr::group_by(gene) %>%
            dplyr::summarise(n_datasets = dplyr::n_distinct(dataset),
                      datasets = paste(sort(unique(dataset)), collapse = ";"),
                      mean_abs_rrb = round(mean(abs_rrb), 3),
                      .groups = "drop") %>%
            dplyr::arrange(dplyr::desc(n_datasets), dplyr::desc(mean_abs_rrb)) %>%
            dplyr::mutate(universe = uname, sig_mode = mode, cim = cim_t,
                   rho_threshold = rc)
          freq_rows[[length(freq_rows) + 1]] <- fr
        }
      }
      # ---- group intersections ----
      for (gname in names(overlap_groups)) {
        members <- overlap_groups[[gname]]
        ms <- lapply(members, function(ds) sets[[ds]][[cim_t]])
        inter <- if (length(ms) == 1) ms[[1]] else Reduce(intersect, ms)
        analysis <- if (mode == "padj") "overlap_BH" else "overlap_rawp"
        add_ov_summary(analysis, mode, gname, members, cim_t, rc, length(inter))
        if (length(inter) > 0)
          group_rows[[length(group_rows) + 1]] <- tibble::tibble(
            analysis = analysis, sig_mode = mode, group = gname,
            group_members = paste(members, collapse = ";"),
            cim = cim_t, rho_threshold = rc, gene = sort(inter))
      }
    }
  }
}

overlap_summary_counts <- dplyr::bind_rows(ov_summary) %>%
  dplyr::arrange(analysis, sig_mode, group, cim, rho_threshold)
data.table::fwrite(overlap_summary_counts,
                   file.path(out_dir, "overlap_summary_counts.csv"), na = "NA")

overlap_group_intersections <- dplyr::bind_rows(group_rows)
if (nrow(overlap_group_intersections) > 0)
  overlap_group_intersections <- overlap_group_intersections %>%
    dplyr::select(analysis, sig_mode, group, group_members, cim, rho_threshold, gene)
data.table::fwrite(overlap_group_intersections,
                   file.path(out_dir, "overlap_group_intersections.csv"),
                   na = "NA")

overlap_frequency_per_gene <- dplyr::bind_rows(freq_rows)
# Per-dataset YES/NO columns: is the gene a significant correlate of THIS cim in
# that dataset? (replaces the ";"-joined datasets string.)
fp <- paste0(";", overlap_frequency_per_gene$datasets, ";")
for (d in ALL_DATASETS)
  overlap_frequency_per_gene[[d]] <-
    ifelse(grepl(paste0(";", d, ";"), fp, fixed = TRUE), "YES", "NO")
# CIM class within each (universe, sig_mode, rho): is the gene a significant
# correlate toward Fun ONLY, Dys ONLY, or BOTH (direction flips across datasets)?
cim_pres <- overlap_frequency_per_gene %>%
  dplyr::group_by(universe, sig_mode, rho_threshold, gene) %>%
  dplyr::summarise(has_fun = any(cim == "Fun"), has_dys = any(cim == "Dys"),
            .groups = "drop") %>%
  dplyr::mutate(cim_class = dplyr::case_when(has_fun & has_dys ~ "Both",
                                      has_fun ~ "Fun_only",
                                      TRUE ~ "Dys_only"))
overlap_frequency_per_gene <- overlap_frequency_per_gene %>%
  dplyr::left_join(dplyr::select(cim_pres, universe, sig_mode, rho_threshold, gene,
                          cim_class),
            by = c("universe", "sig_mode", "rho_threshold", "gene")) %>%
  dplyr::select(universe, sig_mode, cim, cim_class, rho_threshold, gene, n_datasets,
         dplyr::all_of(ALL_DATASETS), mean_abs_rrb)
data.table::fwrite(overlap_frequency_per_gene,
                   file.path(out_dir, "overlap_frequency_per_gene.csv"),
                   na = "NA")

# Small summary: how many unique genes are Fun_only / Dys_only / Both.
overlap_cim_class_summary <- overlap_frequency_per_gene %>%
  dplyr::distinct(universe, sig_mode, rho_threshold, gene, cim_class) %>%
  dplyr::count(universe, sig_mode, rho_threshold, cim_class, name = "n_genes") %>%
  dplyr::arrange(universe, sig_mode, rho_threshold, cim_class)
data.table::fwrite(overlap_cim_class_summary,
                   file.path(out_dir, "overlap_cim_class_summary.csv"),
                   na = "NA")

overlap_frequency_summary <- overlap_frequency_per_gene %>%
  dplyr::count(universe, sig_mode, cim, rho_threshold, n_datasets,
        name = "n_genes_exactly") %>%
  dplyr::arrange(universe, sig_mode, cim, rho_threshold, dplyr::desc(n_datasets)) %>%
  dplyr::group_by(universe, sig_mode, cim, rho_threshold) %>%
  dplyr::mutate(n_genes_shared_by_at_least = cumsum(n_genes_exactly)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(min_n_datasets = n_datasets) %>%
  dplyr::arrange(universe, sig_mode, cim, rho_threshold, min_n_datasets)
data.table::fwrite(overlap_frequency_summary,
                   file.path(out_dir, "overlap_frequency_summary.csv"),
                   na = "NA")

# ---- Console summary --------------------------------------------------------
message("\n==================== DONE ====================")
message("Output folder: ", out_dir)
message("\nMeta-significant gene counts (Fisher padj<", META_ALPHA, "):")
print(as.data.frame(meta_significant_summary))
message("\nFiles written:")
message("  meta_results.csv                (", nrow(meta_results), " rows, ",
        ncol(meta_results), " cols)")
message("  meta_results_data_dictionary.csv(", nrow(meta_dict), " rows)")
message("  meta_cim_class_summary.csv      (", nrow(meta_cim_class_summary), " rows)")
message("  meta_significant_summary.csv    (", nrow(meta_significant_summary), " rows)")
message("  overlap_summary_counts.csv      (", nrow(overlap_summary_counts), " rows)")
message("  overlap_group_intersections.csv (", nrow(overlap_group_intersections), " rows)")
message("  overlap_frequency_per_gene.csv  (", nrow(overlap_frequency_per_gene), " rows)")
message("  overlap_cim_class_summary.csv   (", nrow(overlap_cim_class_summary), " rows)")
message("  overlap_frequency_summary.csv   (", nrow(overlap_frequency_summary), " rows)")
