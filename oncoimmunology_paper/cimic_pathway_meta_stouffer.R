# =============================================================================
# cimic_pathway_meta_stouffer.R
#
# Stouffer meta-analysis of the PATHWAY-LEVEL fry / camera statistics across
# NKI-SMC, NEO, and the TNBC cell-line panel.
#
# This is the pathway analogue of cimic_overlap_meta_rankbiserial.R (which
# combines gene-level r_rb). Here the per-dataset inputs are the fry + camera
# tables already written by the heatmap and ICD figure scripts, so no
# expression data is re-read and no gene-set test is re-run.
#
# TWO FAMILIES, each tested independently:
#   zscore : the 19 a priori immune/stress pathways behind the Z-score heatmaps
#            (fig2e / neo / fig5f). The cell-line panel measures only 13 of the
#            19, so k_datasets varies within the ALL3 scope -- see below.
#   icd    : the five a priori ICD programs behind the box-plot figures
#            (fig2f_2k / neo / fig5i_5l), plus the standalone
#            GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS pathway.
#
# TWO TESTS, run separately end to end:
#   Fry    : self-contained, from Fry_PValue + Fry_Direction
#   Camera : competitive, from Camera_PValue + Camera_Direction
#            (cameraPR with block adjustment in the cell-line cohort)
#
# FOUR SCOPES:
#   ALL3_NKI.NEO.TNBCCL, NKI.TNBCCL, NKI.NEO, NEO.TNBCCL
#
# ----------------------------------------------------------------------------
# DIRECTION -- why no per-cohort flipping is needed here
# ----------------------------------------------------------------------------
# Every input table was written with the SAME contrast, "Dys-CIM - Fun-CIM",
# and its Direction column is expressed against that contrast:
#     Up   = enriched in Dys-CIM
#     Down = enriched in Fun-CIM
# The upstream scripts already handled each cohort's cluster coding when they
# built the contrast (TNBC_CL codes clusters the other way round), so the
# Direction columns are directly comparable across all six files. This is unlike
# the gene-level meta, where the sign of r_rb had to be re-mapped per dataset.
#
# ----------------------------------------------------------------------------
# METHOD
# ----------------------------------------------------------------------------
# For each dataset the reported p-value is two-sided. It is converted to a
# one-sided p toward Dys-CIM,
#     p_os = p/2      if Direction == "Up"    (toward Dys-CIM)
#            1 - p/2  if Direction == "Down"  (toward Fun-CIM)
# and then to a signed z, z = qnorm(1 - p_os). Positive z favours Dys-CIM.
#
#   STOUFFER (PRIMARY, as requested)
#     Z = sum(w_i * z_i) / sqrt(sum(w_i^2)),  w_i = sqrt(n_i)
#     Two-sided p = 2 * pnorm(-|Z|). The sign of Z gives the consensus
#     direction. Because the inputs are signed, a pathway that flips direction
#     between datasets cancels rather than accumulating evidence.
#
#   FISHER (SECONDARY, direction-agnostic)
#     X2 = -2 * sum(log(p_two_sided)), df = 2k. Reported for completeness and
#     deliberately NOT directional: it answers "is there evidence of a shift in
#     any dataset, in any direction", which is a weaker and different question
#     from Stouffer's "is there a consistent shift toward one CIM". Where the
#     two disagree, Stouffer is the one to quote.
#
# Weights use the same effective sample sizes as the gene-level meta: patient
# cohorts contribute their sample n, the cell-line panel contributes n = 9
# (cell lines, the unit of replication under duplicateCorrelation) rather than
# its 27 sequenced samples.
#
# Fry_PValue.Mixed is deliberately excluded: the mixed fry test is
# non-directional by construction, so a signed Stouffer Z is undefined for it.
#
# ----------------------------------------------------------------------------
# MULTIPLICITY
# ----------------------------------------------------------------------------
# Meta p-values are BH-adjusted across pathways within each
# (family x test x scope). The standalone GOBP_INFLAMMATORY_CELL_APOPTOTIC_
# PROCESS row is excluded from the ICD adjustment family and its FDR reported
# as NA, matching how the upstream scripts treat it.
#
# Outputs -> oncoimmunology_paper/Results/cimic_pathway_meta/ :
#   pathway_meta_results.csv           per pathway x family x test x scope
#   pathway_meta_results_dictionary.csv
#   pathway_meta_summary.csv           significant counts per family/test/scope
#   pathway_meta_input_long.csv        the harmonised per-dataset inputs
#   heatmap_<family>_<test>.png/.pdf   signed -log10(p) per pathway x scope
#
# Run from the repository root:
#   Rscript oncoimmunology_paper/cimic_pathway_meta_stouffer.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(purrr)
  library(tibble); library(ggplot2); library(stringr)
})

# ---- 0. Config --------------------------------------------------------------
base_dir <- "oncoimmunology_paper"
res_dir  <- file.path(base_dir, "Results")
out_dir  <- file.path(res_dir, "cimic_pathway_meta")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

META_ALPHA <- 0.05
CLAMP      <- 1e-300

# ---- 1. Input registry ------------------------------------------------------
# family -> dataset -> the fry/camera statistics CSV written by the figure script
inputs <- list(
  zscore = c(
    nki_smc = file.path(res_dir, "fig2e_zscore_heatmap_fry_camera_statistics.csv"),
    neo     = file.path(res_dir, "neo_zscore_heatmap_fry_camera_statistics.csv"),
    tnbc_cl = file.path(res_dir, "fig5f_tnbc_zscore_heatmap_fry_camera_statistics.csv")
  ),
  icd = c(
    nki_smc = file.path(res_dir, "fig2f_2k_icd_programs_fry_camera_statistics.csv"),
    neo     = file.path(res_dir, "neo_icd_programs_fry_camera_statistics.csv"),
    tnbc_cl = file.path(res_dir, "fig5i_5l_tnbc_cl_icd_programs_fry_camera_statistics.csv")
  )
)
ALL_DATASETS <- c("nki_smc", "neo", "tnbc_cl")
TESTS        <- c("Fry") 
# camera

# Master induction tables, used only to recover n for the Stouffer weights.
induction_files <- c(
  nki_smc = file.path(base_dir, "Datasets", "NKI_SMC",
                      "nki_smc_master_induction_df.csv"),
  neo     = file.path(base_dir, "Datasets", "NEO",
                      "neo_master_induction_df.csv"),
  tnbc_cl = file.path(base_dir, "Datasets", "TNBC_CL_Epirubicin",
                      "tnbc_cl_epirubicin_master_induction_df.csv")
)

# ---- 2. Scopes --------------------------------------------------------------
analysis_scopes <- list(
  ALL3_NKI.NEO.TNBCCL = c("nki_smc", "neo", "tnbc_cl"),
  NKI.TNBCCL          = c("nki_smc", "tnbc_cl"),
  NKI.NEO             = c("nki_smc", "neo"),
  NEO.TNBCCL          = c("neo", "tnbc_cl")
)

# ---- 3. Effective sample sizes (Stouffer weights) ---------------------------
# n is constant across genes; recover it from sem = sd/sqrt(n) (mode over
# finite values), exactly as cimic_overlap_meta_rankbiserial.R does.
recover_n <- function(ind_path) {
  if (!file.exists(ind_path)) return(NA_real_)
  m <- data.table::fread(ind_path, showProgress = FALSE, data.table = FALSE)
  need <- c("cluster1_sd", "cluster1_sem", "cluster2_sd", "cluster2_sem")
  if (!all(need %in% names(m))) return(NA_real_)
  mode_n <- function(sd, sem) {
    v <- round((sd / sem)^2); v <- v[is.finite(v) & v > 0]
    if (!length(v)) return(NA_integer_)
    as.integer(names(which.max(table(v))))
  }
  n1 <- mode_n(m$cluster1_sd, m$cluster1_sem)
  n2 <- mode_n(m$cluster2_sd, m$cluster2_sem)
  if (is.na(n1) || is.na(n2)) NA_real_ else n1 + n2
}
n_map <- vapply(ALL_DATASETS, function(d) recover_n(induction_files[[d]]),
                numeric(1))

# The cell-line fry/camera tests block on CELL LINE, so the panel size (9), not
# the 27 sequenced samples, is the effective evidence. Patient cohorts keep n.
CELL_LINE_DATASETS <- c("tnbc_cl")
CELL_LINE_EFF_N    <- 9
n_map[base::intersect(names(n_map), CELL_LINE_DATASETS)] <- CELL_LINE_EFF_N
message("Effective n per dataset (cell lines overridden to ",
        CELL_LINE_EFF_N, "):")
print(n_map)

# ---- 4. Load + harmonise the six input tables -------------------------------
# Column names differ slightly between files (fig5f uses Analysis / Fry_Test,
# the others AnalysisPool / FryTest). Only GeneSet, the two p-value columns and
# the two Direction columns are needed, and those are consistently named.
first_present <- function(d, candidates, default = NA_character_) {
  hit <- base::intersect(candidates, names(d))
  if (!length(hit)) rep(default, nrow(d)) else d[[hit[1]]]
}

load_one <- function(family, dataset) {
  fp <- inputs[[family]][[dataset]]
  if (!file.exists(fp))
    stop("Missing input table for ", family, "/", dataset, ": ", fp,
         "\n  Run the corresponding figure script first.")
  d <- data.table::fread(fp, showProgress = FALSE, data.table = FALSE)

  need <- c("GeneSet", "Fry_PValue", "Fry_Direction",
            "Camera_PValue", "Camera_Direction")
  miss <- base::setdiff(need, names(d))
  if (length(miss))
    stop("Missing column(s) in ", fp, ": ", paste(miss, collapse = ", "))

  pool <- first_present(d, c("AnalysisPool", "Analysis"))
  base <- data.frame(
    family        = family,
    dataset       = dataset,
    GeneSet       = as.character(d$GeneSet),
    analysis_pool = as.character(pool),
    # Standalone pathways are excluded from the pooled BH family, matching the
    # upstream convention that also blanks their FDR.
    standalone    = grepl("standalone", as.character(pool), ignore.case = TRUE),
    n_genes       = suppressWarnings(as.numeric(
                      first_present(d, c("Fry_NGenes", "Camera_NGenes")))),
    stringsAsFactors = FALSE
  )

  dplyr::bind_rows(lapply(TESTS, function(tst) {
    p   <- suppressWarnings(as.numeric(d[[paste0(tst, "_PValue")]]))
    dir <- as.character(d[[paste0(tst, "_Direction")]])
    cbind(base, data.frame(test = tst, p_value = p, direction = dir,
                           stringsAsFactors = FALSE))
  }))
}

message("\nLoading ", length(unlist(inputs)), " input tables ...")
long <- dplyr::bind_rows(
  purrr::map_dfr(names(inputs), function(fam)
    purrr::map_dfr(ALL_DATASETS, function(ds) load_one(fam, ds)))
)

# Direction must be exactly Up/Down for the one-sided conversion to be valid.
bad_dir <- long %>%
  dplyr::filter(!is.na(p_value), !direction %in% c("Up", "Down"))
if (nrow(bad_dir))
  stop("Unexpected Direction value(s): ",
       paste(unique(bad_dir$direction), collapse = ", "),
       ". Expected only 'Up' or 'Down'.")

long <- long %>% dplyr::filter(!is.na(p_value), is.finite(p_value))
message("Loaded ", nrow(long), " dataset-pathway-test rows.")

# Coverage report: which pathways each dataset actually measured.
coverage <- long %>%
  dplyr::filter(test == "Fry") %>%
  dplyr::count(family, dataset, name = "n_pathways") %>%
  tidyr::pivot_wider(names_from = dataset, values_from = n_pathways)
message("\nPathways per dataset:")
print(as.data.frame(coverage))

data.table::fwrite(long, file.path(out_dir, "pathway_meta_input_long.csv"),
                   na = "NA")

# ---- 5. One-sided p and signed z toward Dys-CIM -----------------------------
long <- long %>%
  dplyr::mutate(
    p_os = ifelse(direction == "Up", p_value / 2, 1 - p_value / 2),
    p_os = pmin(pmax(p_os, CLAMP), 1 - 1e-15),
    # Positive z favours Dys-CIM, negative favours Fun-CIM.
    z    = stats::qnorm(1 - p_os),
    w    = { nn <- n_map[dataset]; ifelse(is.finite(nn), sqrt(nn), 1) }
  )

# ---- 6. Stouffer (primary) + Fisher (secondary) -----------------------------
meta_one <- function(fam, tst, scope_name, members) {
  w <- long %>%
    dplyr::filter(family == fam, test == tst, dataset %in% members)
  if (!nrow(w)) return(NULL)

  w %>%
    dplyr::group_by(GeneSet) %>%
    dplyr::summarise(
      analysis_pool = dplyr::first(analysis_pool),
      standalone    = any(standalone),
      k_datasets    = dplyr::n(),
      datasets      = paste(sort(dataset), collapse = ";"),
      # Signed Stouffer toward Dys-CIM.
      stouffer_Z    = sum(w * z) / sqrt(sum(w^2)),
      # Direction-agnostic Fisher on the ORIGINAL two-sided p-values.
      fisher_stat   = -2 * sum(log(pmax(p_value, CLAMP))),
      fisher_p      = stats::pchisq(fisher_stat, df = 2 * dplyr::n(),
                                    lower.tail = FALSE),
      # How many datasets point the same way as the meta consensus?
      n_up          = sum(direction == "Up"),
      n_down        = sum(direction == "Down"),
      min_p_dataset = min(p_value),
      max_p_dataset = max(p_value),
      mean_n_genes  = round(mean(n_genes, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    dplyr::filter(k_datasets >= 2) %>%
    dplyr::mutate(
      stouffer_p         = 2 * stats::pnorm(-abs(stouffer_Z)),
      meta_direction     = ifelse(stouffer_Z > 0, "Dys-CIM", "Fun-CIM"),
      k_concordant       = ifelse(stouffer_Z > 0, n_up, n_down),
      direction_concordant = ifelse(k_concordant == k_datasets, "YES", "NO"),
      n_scope_datasets   = length(members),
      all_datasets_present = ifelse(k_datasets == length(members), "YES", "NO"),
      family = fam, test = tst, scope = scope_name
    )
}

message("\n== Stouffer meta-analysis ==")
rows <- list()
for (fam in names(inputs))
  for (tst in TESTS)
    for (sc in names(analysis_scopes)) {
      r <- meta_one(fam, tst, sc, analysis_scopes[[sc]])
      if (!is.null(r) && nrow(r)) rows[[length(rows) + 1]] <- r
    }
if (!length(rows)) stop("Meta-analysis produced no rows.")
meta <- dplyr::bind_rows(rows)

# ---- 7. BH within (family, test, scope), standalone excluded ----------------
meta <- meta %>%
  dplyr::group_by(family, test, scope) %>%
  dplyr::mutate(
    stouffer_padj = ifelse(standalone, NA_real_,
      stats::p.adjust(ifelse(standalone, NA_real_, stouffer_p), method = "BH")),
    fisher_padj   = ifelse(standalone, NA_real_,
      stats::p.adjust(ifelse(standalone, NA_real_, fisher_p), method = "BH"))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    stouffer_sig = ifelse(!is.na(stouffer_padj) & stouffer_padj < META_ALPHA,
                          "YES", "NO"),
    fisher_sig   = ifelse(!is.na(fisher_padj) & fisher_padj < META_ALPHA,
                          "YES", "NO"),
    stouffer_sig_nominal = ifelse(stouffer_p < META_ALPHA, "YES", "NO")
  )

# Per-dataset direction and p, so each meta row can be audited at a glance.
wide_dir <- long %>%
  dplyr::select(family, test, GeneSet, dataset, direction, p_value) %>%
  tidyr::pivot_wider(names_from = dataset,
                     values_from = c(direction, p_value),
                     names_glue = "{dataset}_{.value}")
meta <- meta %>% dplyr::left_join(wide_dir, by = c("family", "test", "GeneSet"))

ds_cols <- as.vector(rbind(paste0(ALL_DATASETS, "_direction"),
                           paste0(ALL_DATASETS, "_p_value")))
ds_cols <- base::intersect(ds_cols, names(meta))

meta <- meta %>%
  dplyr::arrange(family, test, scope, stouffer_p) %>%
  dplyr::select(family, test, scope, GeneSet, analysis_pool, standalone,
                n_scope_datasets, k_datasets, all_datasets_present,
                meta_direction, k_concordant, direction_concordant,
                n_up, n_down,
                stouffer_Z, stouffer_p, stouffer_padj, stouffer_sig,
                stouffer_sig_nominal,
                fisher_stat, fisher_p, fisher_padj, fisher_sig,
                min_p_dataset, max_p_dataset, mean_n_genes, datasets,
                dplyr::all_of(ds_cols))

data.table::fwrite(meta, file.path(out_dir, "pathway_meta_results.csv"),
                   na = "NA")

# ---- 8. Data dictionary -----------------------------------------------------
dict <- tibble::tribble(
  ~column, ~meaning,
  "family", "zscore (19 immune/stress pathways) or icd (5 ICD programs + standalone).",
  "test", "Fry (self-contained) or Camera (competitive). Run independently.",
  "scope", "Which dataset combination was combined.",
  "GeneSet", "Pathway / programme name as written by the upstream figure script.",
  "analysis_pool", "AnalysisPool label carried from the input table.",
  "standalone", paste("TRUE for GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS in the",
                      "icd family, which is tested outside the pooled family;",
                      "its FDR columns are NA by design."),
  "n_scope_datasets", "Datasets in the scope (3 or 2).",
  "k_datasets", paste("Datasets that actually reported this pathway. Less than",
                      "n_scope_datasets where a cohort did not measure it (the",
                      "cell-line panel covers 13 of the 19 zscore pathways)."),
  "all_datasets_present", "YES when k_datasets == n_scope_datasets.",
  "meta_direction", "Sign of stouffer_Z: Dys-CIM if positive, Fun-CIM if negative.",
  "k_concordant", "Datasets whose own Direction agrees with meta_direction.",
  "direction_concordant", "YES when every contributing dataset agrees.",
  "n_up / n_down", "Datasets reporting Up (Dys-CIM) / Down (Fun-CIM).",
  "stouffer_Z", paste("PRIMARY. sqrt(n)-weighted sum of signed one-sided z",
                      "toward Dys-CIM, normalised by sqrt(sum of squared",
                      "weights). Cell-line weight uses n = 9 cell lines."),
  "stouffer_p", "Two-sided p for stouffer_Z: 2 * pnorm(-|Z|). Raw.",
  "stouffer_padj", "stouffer_p, BH-adjusted within family x test x scope.",
  "stouffer_sig", "YES when stouffer_padj < 0.05.",
  "stouffer_sig_nominal", "YES when the RAW stouffer_p < 0.05 (no FDR).",
  "fisher_stat / fisher_p", paste("SECONDARY, direction-agnostic. Fisher on the",
                      "original TWO-SIDED p-values, so it does not require the",
                      "datasets to agree on direction. Quote Stouffer when the",
                      "two disagree."),
  "fisher_padj / fisher_sig", "As above, BH-adjusted within family x test x scope.",
  "min_p_dataset / max_p_dataset", "Range of the contributing per-dataset p-values.",
  "mean_n_genes", "Mean measured gene-set size across contributing datasets.",
  "datasets", "Semicolon-separated contributing datasets.",
  "<dataset>_direction / <dataset>_p_value",
    "That dataset's own Direction and two-sided p for this pathway and test."
)
data.table::fwrite(dict,
                   file.path(out_dir, "pathway_meta_results_dictionary.csv"))

# ---- 9. Summary -------------------------------------------------------------
summary_tbl <- meta %>%
  dplyr::filter(!standalone) %>%
  dplyr::group_by(family, test, scope) %>%
  dplyr::summarise(
    n_pathways          = dplyr::n(),
    n_all_present       = sum(all_datasets_present == "YES"),
    n_stouffer_sig      = sum(stouffer_sig == "YES"),
    n_stouffer_nominal  = sum(stouffer_sig_nominal == "YES"),
    n_fisher_sig        = sum(fisher_sig == "YES"),
    n_concordant        = sum(direction_concordant == "YES"),
    n_toward_Dys        = sum(meta_direction == "Dys-CIM"),
    n_toward_Fun        = sum(meta_direction == "Fun-CIM"),
    .groups = "drop"
  ) %>%
  dplyr::arrange(family, test, scope)
data.table::fwrite(summary_tbl, file.path(out_dir, "pathway_meta_summary.csv"),
                   na = "NA")

# ---- 10. Heatmaps: signed -log10(p) per pathway x scope ---------------------
star <- function(p) dplyr::case_when(
  is.na(p)   ~ "",
  p <= 1e-4  ~ "****",
  p <= 1e-3  ~ "***",
  p <= 1e-2  ~ "**",
  p <= 0.05  ~ "*",
  TRUE       ~ ""
)

for (fam in names(inputs)) {
  for (tst in TESTS) {
    d <- meta %>% dplyr::filter(family == fam, test == tst)
    if (!nrow(d)) next

    d <- d %>%
      dplyr::mutate(
        # Signed so the colour encodes direction as well as strength.
        signed_score = -log10(pmax(stouffer_p, 1e-300)) * sign(stouffer_Z),
        label = paste0(star(stouffer_padj),
                       ifelse(direction_concordant == "NO", "†", "")),
        Program = GeneSet %>%
          stringr::str_remove_all("GOBP_|REACTOME_|HALLMARK_|WP_") %>%
          stringr::str_replace_all("_", " ") %>%
          stringr::str_wrap(width = 40),
        scope = factor(scope, levels = names(analysis_scopes))
      )
    ord <- d %>%
      dplyr::group_by(Program) %>%
      dplyr::summarise(m = mean(signed_score), .groups = "drop") %>%
      dplyr::arrange(m) %>% dplyr::pull(Program)
    d$Program <- factor(d$Program, levels = ord)

    lim <- max(abs(d$signed_score), na.rm = TRUE)
    p <- ggplot(d, aes(x = scope, y = Program, fill = signed_score)) +
      geom_tile(colour = "black", linewidth = 0.5) +
      geom_text(aes(label = label), size = 6, fontface = "bold") +
      scale_fill_gradient2(
        low = "#313695", mid = "white", high = "#A50026", midpoint = 0,
        limits = c(-lim, lim),
        name = expression(paste("signed  ", -log[10](p[Stouffer])))) +
      labs(
        x = NULL, y = NULL,
        title = paste0("Stouffer meta-analysis - ", fam, " / ", tst),
        subtitle = paste0("Red = toward Dys-CIM, blue = toward Fun-CIM. ",
                          "Stars = BH-adjusted p. † = datasets disagree ",
                          "on direction.")) +
      theme_classic(base_size = 14) +
      theme(
        axis.text.x = element_text(face = "bold", size = 12, colour = "black",
                                   angle = 30, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 11, colour = "black"),
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 10, colour = "grey30"),
        legend.position = "right",
        axis.line = element_blank()
      )

    h <- max(4, 0.42 * length(unique(d$Program)) + 2.5)
    stem <- file.path(out_dir, paste0("heatmap_", fam, "_", tolower(tst)))
    ggsave(paste0(stem, ".png"), p, width = 10, height = h, dpi = 300,
           bg = "white")
    ggsave(paste0(stem, ".pdf"), p, width = 10, height = h, bg = "white")
    message("  saved ", basename(stem), ".png / .pdf")
  }
}

# ---- 11. Console summary ----------------------------------------------------
message("\n==================== DONE ====================")
message("Output folder: ", out_dir)
message("\nSignificant pathway counts (BH < ", META_ALPHA,
        ", standalone excluded):")
print(as.data.frame(summary_tbl))
message("\nFiles written:")
message("  pathway_meta_results.csv            (", nrow(meta), " rows)")
message("  pathway_meta_results_dictionary.csv (", nrow(dict), " rows)")
message("  pathway_meta_summary.csv            (", nrow(summary_tbl), " rows)")
message("  pathway_meta_input_long.csv         (", nrow(long), " rows)")
message("  heatmap_<family>_<test>.png/.pdf    (",
        length(inputs) * length(TESTS), " figures)")
