# =============================================================================
# cimic_ora_signatures.R
# -----------------------------------------------------------------------------
# Over-representation analysis (ORA) on the Fun-CIM and Dys-CIM consensus gene
# signatures produced by cimic_overlap_meta_rankbiserial.R.
#
# Adapted from ora_cimic_signatures.R (9_Projects/combined_cimic_analysis_bc).
# The ORA logic, databases, thresholds, and plotting are carried over unchanged.
# The only structural change: gene lists come from meta_results.csv rather than
# the cimic_gene_analysis_FDR.xlsx workbook (which does not exist for this
# 3-dataset analysis), so "tier" is replaced by "scope".
#
# Scopes tested (the four contexts available with three datasets):
#   ALL3_NKI.NEO.TNBCCL, NKI.NEO, NKI.TNBCCL, NEO.TNBCCL
#
# Significance filter: fisher_sig == "YES" & stouffer_sig == "YES" (both)
# Rho thresholds: 0.3, 0.6
# Top 200 ranked by: stouffer_padj ASC, mean_abs_rrb DESC, k_datasets DESC
#
# K-modes (see section 2a): how many datasets must back a gene.
#   k_min2            k_datasets >= 2  (2 OR 3 datasets for the ALL3 scope)
#   k_all             k_datasets == n_scope_datasets  (all three)
#   k_all_concordant  all three AND all agreeing on direction
# Each mode writes to its own subdirectory, so modes never overwrite each other.
#
# Output layout (short names -- Windows MAX_PATH is 260 chars; see section 2b):
#   <OUT_ROOT>/<scope_short>/<kmode_short>/rho<rho>_<cim>_<all|top200>/
# e.g.
#   .../cimic_ora_signatures/ALL3/kAllC/rho0.3_Fun_all/ora_GO_BP.csv
# run_index.csv at the output root maps every short directory back to its full
# scope / k-mode / cim / rho description.
#
# Databases:
#   * GO Biological Process  (enrichGO + simplified)
#   * Reactome               (enrichPathway / ReactomePA)
#   * MSigDB Hallmark (H)    (enricher)
#   * Proteostasis branches  (custom enricher from proteostatic_network_genes.csv,
#                             skipped automatically if the file is unavailable)
#
# Universe: protein-coding genes (helper_code/protein_coding_genes.rds)
#           intersected with the genes testable in that scope, i.e. the genes
#           present in meta_results.csv for that scope. See section 4.
#
# Outputs per analysis directory:
#   CSV  genes_used.csv              the exact ranked gene list fed to ORA, with
#                                    its meta_results statistics and an
#                                    in_ora_universe / used_in_ora flag. Written
#                                    even when the list is too small to test.
#   CSV  ora_GO_BP.csv               enrichGO full result
#   CSV  ora_GO_BP_simplified.csv    after clusterProfiler::simplify()
#   CSV  ora_reactome.csv            ReactomePA::enrichPathway
#   CSV  ora_hallmark.csv            MSigDB Hallmark enricher
#   CSV  ora_proteostasis.csv        proteostasis-branch enricher (optional)
#   FIG  <same stem>.png / .pdf      barplot per database
#
# Output root: oncoimmunology_paper/Results/cimic_ora_signatures/
#
# Run from the repository root, AFTER cimic_overlap_meta_rankbiserial.R:
#   Rscript oncoimmunology_paper/cimic_ora_signatures.R
# =============================================================================

suppressMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(msigdbr)
})

# -- 1. PATHS -----------------------------------------------------------------
BASE      <- "oncoimmunology_paper"
META_CSV  <- file.path(BASE, "Results", "cimic_meta_analysis", "meta_results.csv")
PC_RDS    <- file.path(BASE, "helper_code", "protein_coding_genes.rds")
# Proteostasis branch annotation lives outside this repository. Optional: if the
# file is not found the proteostasis database is skipped and everything else runs.
PROT_CSV  <- paste0("D:/OneDrive - University of Florida/Gbadamosi Lab/",
                    "Mohammed Gbadamosi/9_Projects/tonic_cimic_analysis/",
                    "proteostatic_network_genes.csv")
OUT_ROOT  <- file.path(BASE, "Results", "cimic_ora_signatures")
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(META_CSV))
  stop("meta_results.csv not found. Run cimic_overlap_meta_rankbiserial.R ",
       "first.\n  Expected at: ", META_CSV)

# -- 2. PARAMETERS ------------------------------------------------------------
RHO_LEVELS  <- c(0.3) 
#c(0.3, 0.6)
PADJ_CUTOFF <- 0.05
N_TOP_TERMS <- 15      # terms per barplot
TOP_N_GENES <- 200     # top genes per scope
MIN_GENES   <- 3       # skip gene lists smaller than this

SCOPES      <- c("ALL3_NKI.NEO.TNBCCL", "NKI.NEO", "NKI.TNBCCL", "NEO.TNBCCL")
#c("ALL3_NKI.NEO.TNBCCL", "NKI.NEO", "NKI.TNBCCL", "NEO.TNBCCL")
CIM_CLASSES <- c("Fun", "Dys")

# -- 2b. SHORT NAMES FOR PATHS ------------------------------------------------
# Windows caps a full path at 260 characters (MAX_PATH). The previous layout
# repeated the whole scope/k-mode/rho/cim label in BOTH the directory tree and
# every filename, which pushed the deepest files past that limit and made the
# writes fail silently. Directories and filenames are now short and
# non-redundant; the full descriptive label is kept for plot titles and is
# recorded in run_index.csv at the output root, so nothing is lost.
SCOPE_SHORT <- c(
  "ALL3_NKI.NEO.TNBCCL" = "ALL3",
  "NKI.NEO"             = "NKI_NEO",
  "NKI.TNBCCL"          = "NKI_CL",
  "NEO.TNBCCL"          = "NEO_CL"
)
KMODE_SHORT <- c(
  k_min2           = "k2",
  k_all            = "kAll",
  k_all_concordant = "kAllC"
)
# Fall back to the full name if an unmapped scope appears, rather than NA.
short_scope <- function(sc) if (sc %in% names(SCOPE_SHORT)) SCOPE_SHORT[[sc]] else sc
short_kmode <- function(km) if (km %in% names(KMODE_SHORT)) KMODE_SHORT[[km]] else km

# -- 2a. K-MODE: how many datasets must back a gene ---------------------------
# NOTE ON WHAT k MEANS. `k_datasets` counts the datasets in which the gene's
# |r_rb| >= rho_threshold. It is MAGNITUDE ONLY: it does not require the gene to
# be significant within any individual dataset, and it does not require the
# datasets to agree on direction. Direction consistency is enforced separately,
# by the one-sided Fisher/Stouffer p (a gene pointing the wrong way in one
# dataset still counts toward k but is penalised in the meta statistic).
#
# `k_concordant` counts only the contributors that point toward the CIM, so
# direction agreement is expressible too. Available modes:
#
#   "k_min2"          k_datasets >= 2. The permissive default: for the all-three
#                     scope this admits genes backed by 2 OR 3 datasets.
#   "k_all"           k_datasets == n_scope_datasets. "Require all three" for the
#                     ALL3 scope (and all two for a pair scope). Magnitude in
#                     every dataset; direction still only enforced statistically.
#   "k_all_concordant" k_datasets == n_scope_datasets AND every contributor
#                     points the same way (direction_concordant == "YES"). This
#                     is the strict "same direction in all three" reading.
#
# Edit K_MODES to choose which run(s) to produce. Outputs are written to a
# separate subdirectory per mode, so modes never overwrite each other.
K_MODES <- c("k_all_concordant")

# For a 2-dataset scope, "k_all" is the same set as "k_min2" (k can only be 2),
# so running both is wasted work. TRUE skips the duplicate and logs it; set to
# FALSE to force both to be written anyway.
SKIP_REDUNDANT_K_MODE <- TRUE

k_mode_filter <- function(dt, mode) {
  switch(mode,
    k_min2 = dt[k_datasets >= 2],
    k_all  = dt[k_datasets == n_scope_datasets],
    k_all_concordant = dt[k_datasets == n_scope_datasets &
                            direction_concordant == "YES"],
    stop("Unknown K_MODE: ", mode)
  )
}

k_mode_label <- c(
  k_min2           = "k>=2 datasets",
  k_all            = "all datasets in scope",
  k_all_concordant = "all datasets, same direction"
)

bad_modes <- setdiff(K_MODES, names(k_mode_label))
if (length(bad_modes))
  stop("Unknown entry in K_MODES: ", paste(bad_modes, collapse = ", "),
       "\n  Valid: ", paste(names(k_mode_label), collapse = ", "))

# -- 3. LOAD META RESULTS -----------------------------------------------------
message("Loading meta results: ", META_CSV)
meta <- fread(META_CSV, showProgress = FALSE)

req_cols <- c("scope", "cim", "rho_threshold", "gene", "k_datasets",
              "n_scope_datasets", "k_concordant", "direction_concordant",
              "fisher_sig", "stouffer_sig", "stouffer_padj", "mean_abs_rrb")
missing_cols <- setdiff(req_cols, names(meta))
if (length(missing_cols))
  stop("meta_results.csv is missing column(s): ",
       paste(missing_cols, collapse = ", "),
       "\n  Re-run cimic_overlap_meta_rankbiserial.R -- these columns were ",
       "added when the k-mode toggle was introduced.")

scopes_present <- intersect(SCOPES, unique(meta$scope))
if (!length(scopes_present))
  stop("None of the expected scopes are present in meta_results.csv. Found: ",
       paste(unique(meta$scope), collapse = ", "))
if (length(scopes_present) < length(SCOPES))
  message("  [note] scopes absent from meta_results.csv: ",
          paste(setdiff(SCOPES, scopes_present), collapse = ", "))
message("  Scopes: ", paste(scopes_present, collapse = ", "))
message("  Rows: ", nrow(meta))

# -- 4. UNIVERSE --------------------------------------------------------------
# Two-part universe, mirroring the intent of the 9-dataset script:
#   (a) protein-coding restriction from the repo's curated symbol list, and
#   (b) the genes actually testable in that scope, i.e. those appearing in
#       meta_results.csv for the scope (union over cim and rho).
# The scope-specific part matters because the three datasets do not share a gene
# universe (21k / 39k / 79k measured symbols), so a single pooled background
# would over-state enrichment for the pairs that exclude the widest dataset.
protein_coding <- readRDS(PC_RDS)
message("  Protein-coding symbols: ", length(protein_coding))

universe_for_scope <- function(sc) {
  g <- unique(meta$gene[meta$scope == sc])
  intersect(g, protein_coding)
}

# Entrez mapping is needed for Reactome. Mapped once per scope in the main loop
# and passed down, rather than re-mapped for every (cim x rho x gene-set-size).
universe_entrez_for_scope <- function(universe_genes) {
  suppressWarnings(
    clusterProfiler::bitr(universe_genes, fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
}

# -- 5. HALLMARK TERM2GENE ----------------------------------------------------
# msigdbr renamed `category` -> `collection`; try the new name and fall back.
message("Loading Hallmark gene sets ...")
hallmark_t2g <- tryCatch(
  msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H")
) |>
  dplyr::select(gs_name, gene_symbol) |>
  dplyr::distinct() |>
  as.data.frame()
message("  Hallmark sets: ", length(unique(hallmark_t2g$gs_name)))

# -- 6. PROTEOSTASIS TERM2GENE (optional) -------------------------------------
# All branch rows per gene -- multi-branch genes contribute to all relevant sets.
prot_t2g <- NULL
if (file.exists(PROT_CSV)) {
  message("Building proteostasis gene sets ...")
  prot_raw <- fread(PROT_CSV, showProgress = FALSE)
  setnames(prot_raw,
           old = intersect(names(prot_raw),
                           c("Gene Symbol", "gene_symbol", "Gene.Symbol")),
           new = "gene_symbol", skip_absent = TRUE)
  setnames(prot_raw,
           old = intersect(names(prot_raw), c("Branch", "branch")),
           new = "branch", skip_absent = TRUE)

  if (all(c("gene_symbol", "branch") %in% names(prot_raw))) {
    prot_t2g <- unique(prot_raw[
      !is.na(branch) & branch != "" & !is.na(gene_symbol) & gene_symbol != "",
      .(gs_name = branch, gene_symbol)
    ])
    message("  Proteostasis branches: ", length(unique(prot_t2g$gs_name)),
            " | gene-branch rows: ", nrow(prot_t2g))
    print(prot_t2g[, .N, by = gs_name][order(-N)])
  } else {
    message("  [skip] proteostasis file lacks gene_symbol/branch columns")
  }
} else {
  message("[skip] proteostasis annotation not found -- that database will be ",
          "omitted.\n  Looked for: ", PROT_CSV)
}

# -- 7. HELPER: barplot -------------------------------------------------------
ora_barplot <- function(res_df, title, wrap = 45) {
  if (is.null(res_df) || nrow(res_df) == 0) return(NULL)
  top <- res_df |>
    dplyr::filter(Count >= 2, p.adjust < PADJ_CUTOFF) |>
    dplyr::slice_min(order_by = p.adjust, n = N_TOP_TERMS, with_ties = FALSE) |>
    dplyr::mutate(
      Description = str_wrap(Description, width = wrap),
      Description = factor(Description, levels = rev(unique(Description)))
    )
  if (nrow(top) == 0) return(NULL)
  ggplot(top, aes(x = Description, y = Count, fill = p.adjust)) +
    geom_col(color = "black") +
    coord_flip() +
    labs(title = title, x = NULL, y = "Gene count", fill = "Adj. p") +
    scale_fill_gradient(low = "red", high = "blue", trans = "reverse",
                        labels = scales::label_scientific(digits = 2)) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid   = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      plot.title   = element_text(face = "bold", size = 13),
      axis.text.y  = element_text(size = 20, colour = "black", face = "bold"),
      axis.text.x  = element_text(size = 20, colour = "black",  face = "bold"),
      axis.title.x = element_text(size = 20, colour = "black", face = "bold"),
      legend.title = element_text(size = 16, colour = "black", face = "bold"),
      legend.text = element_text(size = 16, colour = "black", face = "bold")
    )
}

# -- 8. HELPER: save plot -----------------------------------------------------
save_ora_plot <- function(p, path_no_ext) {
  if (is.null(p)) {
    message("    [no enriched terms] ", basename(path_no_ext))
    return(invisible())
  }
  n_bars <- length(unique(p$data$Description))
  w <- 10; h <- max(4, 0.45 * n_bars + 2)
  ggsave(paste0(path_no_ext, ".png"), p,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(paste0(path_no_ext, ".pdf"), p,
         width = w, height = h, device = cairo_pdf, bg = "white")
}

# -- 9. CORE ORA FUNCTION -----------------------------------------------------
# `label` is the long human-readable description used in PLOT TITLES only.
# Output filenames are short and fixed (ora_GO_BP.csv, ora_hallmark.csv, ...)
# because the directory path already encodes scope / k-mode / rho / cim / set.
# `gene_table` is the ranked meta_results subset the gene list came from; it is
# written to genes_used.csv so every analysis records its own input.
run_ora <- function(genes, label, out_dir, universe_genes, universe_entrez,
                    gene_table = NULL) {

  genes_requested <- unique(genes[!is.na(genes) & genes != ""])
  genes <- intersect(genes_requested, universe_genes)  # restrict to universe
  message("  ", label, ": ", length(genes), " genes (after universe filter)")

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # -- Record the exact gene list used, before any early return ---------------
  # Written even when the list is too small to test, so a SKIP is auditable.
  if (!is.null(gene_table) && nrow(gene_table) > 0) {
    genes_used <- as.data.table(gene_table)
    genes_used[, rank_in_list := seq_len(.N)]
    genes_used[, in_ora_universe := gene %in% universe_genes]
    genes_used[, used_in_ora := in_ora_universe]
    genes_used[, analysis_label := label]
    fwrite(genes_used, file.path(out_dir, "genes_used.csv"), na = "NA")
  } else {
    fwrite(data.table(gene = genes_requested,
                      rank_in_list = seq_along(genes_requested),
                      in_ora_universe = genes_requested %in% universe_genes,
                      used_in_ora = genes_requested %in% universe_genes,
                      analysis_label = label),
           file.path(out_dir, "genes_used.csv"), na = "NA")
  }
  n_dropped <- length(genes_requested) - length(genes)
  if (n_dropped > 0)
    message("    [note] ", n_dropped, " gene(s) not in the ORA universe; ",
            "see genes_used.csv (in_ora_universe = FALSE)")

  if (length(genes) < MIN_GENES) {
    message("    [skip] fewer than ", MIN_GENES, " protein-coding genes")
    return(invisible(length(genes)))
  }

  # -- GO BP (+ simplify) -----------------------------------------------------
  message("    Running GO BP ...")
  ego <- tryCatch(
    clusterProfiler::enrichGO(
      gene          = genes,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = PADJ_CUTOFF,
      qvalueCutoff  = PADJ_CUTOFF,
      universe      = universe_genes,
      keyType       = "SYMBOL",
      minGSSize     = 10,
      maxGSSize     = 500
    ),
    error = function(e) { message("    [GO error] ", conditionMessage(e)); NULL }
  )

  if (!is.null(ego) && nrow(ego@result) > 0) {
    out_stem <- file.path(out_dir, "ora_GO_BP")
    fwrite(ego@result, paste0(out_stem, ".csv"))

    ego_s <- tryCatch(
      clusterProfiler::simplify(ego, cutoff = 0.7, by = "p.adjust",
                                select_fun = min, measure = "Wang"),
      error = function(e) NULL
    )
    if (!is.null(ego_s) && nrow(ego_s@result) > 0) {
      fwrite(ego_s@result, paste0(out_stem, "_simplified.csv"))
      save_ora_plot(
        ora_barplot(ego_s@result, paste0("GO BP - ", label)),
        paste0(out_stem, "_simplified")
      )
    }
    save_ora_plot(
      ora_barplot(ego@result, paste0("GO BP (full) - ", label)),
      out_stem
    )
  } else {
    message("    [GO] no enriched terms")
  }

  # -- Reactome ---------------------------------------------------------------
  message("    Running Reactome ...")
  g_ent <- suppressWarnings(
    clusterProfiler::bitr(genes, fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )

  if (nrow(g_ent) >= MIN_GENES) {
    react <- tryCatch(
      ReactomePA::enrichPathway(
        gene          = g_ent$ENTREZID,
        organism      = "human",
        pvalueCutoff  = PADJ_CUTOFF,
        pAdjustMethod = "BH",
        qvalueCutoff  = 0.2,
        universe      = universe_entrez$ENTREZID,
        minGSSize     = 10,
        maxGSSize     = 500,
        readable      = TRUE
      ),
      error = function(e) {
        message("    [Reactome error] ", conditionMessage(e)); NULL }
    )

    if (!is.null(react) && nrow(react@result) > 0) {
      out_stem <- file.path(out_dir, "ora_reactome")
      fwrite(react@result, paste0(out_stem, ".csv"))
      save_ora_plot(
        ora_barplot(react@result, paste0("Reactome - ", label)),
        out_stem
      )
    } else {
      message("    [Reactome] no enriched terms")
    }
  }

  # -- Hallmark ---------------------------------------------------------------
  message("    Running Hallmark ...")
  hm <- tryCatch(
    clusterProfiler::enricher(
      gene          = genes,
      TERM2GENE     = hallmark_t2g,
      universe      = universe_genes,
      pvalueCutoff  = PADJ_CUTOFF,
      qvalueCutoff  = PADJ_CUTOFF,
      pAdjustMethod = "BH",
      minGSSize     = 5,
      maxGSSize     = 500
    ),
    error = function(e) {
      message("    [Hallmark error] ", conditionMessage(e)); NULL }
  )

  if (!is.null(hm) && nrow(hm@result) > 0) {
    out_stem <- file.path(out_dir, "ora_hallmark")
    fwrite(hm@result, paste0(out_stem, ".csv"))
    save_ora_plot(
      ora_barplot(hm@result, paste0("Hallmark - ", label), wrap = 35),
      out_stem
    )
  } else {
    message("    [Hallmark] no enriched terms")
  }

  # -- Proteostasis branches (optional) ---------------------------------------
  # Custom enricher: branches as gene sets, full protein-coding scope universe.
  # Multi-branch genes contribute to all relevant branches.
  if (!is.null(prot_t2g)) {
    message("    Running Proteostasis branches ...")
    prot <- tryCatch(
      clusterProfiler::enricher(
        gene          = genes,
        TERM2GENE     = as.data.frame(prot_t2g),
        universe      = universe_genes,
        pvalueCutoff  = 1,        # keep all for CSV; filter in plot
        qvalueCutoff  = 1,
        pAdjustMethod = "BH",
        minGSSize     = 5,        # small branches (EX n=24, PN n=41)
        maxGSSize     = 2000
      ),
      error = function(e) {
        message("    [Proteostasis error] ", conditionMessage(e)); NULL }
    )

    if (!is.null(prot) && nrow(prot@result) > 0) {
      out_stem <- file.path(out_dir, "ora_proteostasis")
      fwrite(prot@result, paste0(out_stem, ".csv"))
      # Only plot branches that pass padj cutoff
      sig_prot <- prot@result[prot@result$p.adjust < PADJ_CUTOFF, ]
      save_ora_plot(
        ora_barplot(sig_prot, paste0("Proteostasis - ", label), wrap = 40),
        out_stem
      )
    } else {
      message("    [Proteostasis] no enriched terms")
    }
  }

  invisible(length(genes))
}

# -- 10. EXTRACT GENE LIST FROM meta_results ----------------------------------
# Filter: fisher_sig == "YES" & stouffer_sig == "YES" (both tests), matching the
# "both" filter used by the 9-dataset workbook version.
# Top N: ranked by stouffer_padj ASC, mean_abs_rrb DESC, k_datasets DESC, then
# de-duplicated by gene BEFORE the head() so top_n counts UNIQUE genes.
# abs()<tol guards floating-point equality on rho_threshold.
# Returns the ranked meta_results subset (a data.table), not just gene symbols,
# so the caller can both take `$gene` and write it out as genes_used.csv.
get_scope_gene_table <- function(sc, cim_sel, rho_sel, k_mode, top_n = NULL) {
  empty <- meta[0]
  dt <- meta[
    scope == sc &
      cim == cim_sel &
      fisher_sig == "YES" &
      stouffer_sig == "YES" &
      abs(rho_threshold - rho_sel) < 1e-8
  ]
  if (nrow(dt) == 0) return(empty)

  dt <- k_mode_filter(dt, k_mode)
  if (nrow(dt) == 0) return(empty)

  setorder(dt, stouffer_padj, -mean_abs_rrb, -k_datasets)
  dt <- unique(dt, by = "gene")

  if (!is.null(top_n)) dt <- head(dt, top_n)
  dt
}

# -- 11. MAIN LOOP ------------------------------------------------------------
# Iterate: scope x k_mode x cim x rho x gene_set_size (all / top200)
log_lines <- c()
run_index <- list()   # short dir -> full description, written to run_index.csv

for (sc in scopes_present) {

  universe_genes  <- universe_for_scope(sc)
  message("\n########## SCOPE: ", sc, " ##########")
  message("  Universe (protein-coding & testable in scope): ",
          length(universe_genes), " genes")

  if (length(universe_genes) < MIN_GENES) {
    message("  [skip scope] universe too small")
    next
  }
  universe_entrez <- universe_entrez_for_scope(universe_genes)
  message("  Universe Entrez mapped: ", nrow(universe_entrez))

  # How many datasets this scope has, for the redundant-mode check below.
  n_sc <- unique(meta$n_scope_datasets[meta$scope == sc])[1]

  for (k_mode in K_MODES) {

    # For a 2-dataset scope, k_all selects exactly the same genes as k_min2.
    if (SKIP_REDUNDANT_K_MODE && k_mode == "k_all" && isTRUE(n_sc <= 2)) {
      message("\n  [skip k-mode '", k_mode, "'] scope has ", n_sc,
              " datasets, so it is identical to 'k_min2'")
      log_lines <- c(log_lines,
                     paste(sc, k_mode, "-", "-", "-", "-",
                           "SKIP_redundant_equals_k_min2", sep = "\t"))
      next
    }

    message("\n  ===== k-mode: ", k_mode, " (", k_mode_label[[k_mode]], ") =====")

    for (cim in CIM_CLASSES) {
      message("\n=== ", sc, " / ", cim, "-CIM / ", k_mode, " ===")

      for (rho in RHO_LEVELS) {
        message("\n  -- rho = ", rho, " --")

        # Short, non-redundant directory: <ALL3>/<kAllC>/rho0.3_Fun_all/
        # The full label lives in plot titles, genes_used.csv, and run_index.csv.
        set_dir <- function(set_tag)
          file.path(OUT_ROOT, short_scope(sc), short_kmode(k_mode),
                    paste0("rho", rho, "_", cim, "_", set_tag))

        # All genes
        tab_all   <- get_scope_gene_table(sc, cim, rho, k_mode, top_n = NULL)
        genes_all <- tab_all$gene
        n_all     <- length(genes_all)
        message("  All genes: ", n_all)

        # run_ora() is called whenever at least one gene was selected, even if
        # that is below MIN_GENES: it writes genes_used.csv first and then skips
        # the enrichment itself, so a too-few analysis is still auditable.
        if (n_all > 0) {
          label_all <- paste0(sc, " | ", k_mode, " | ", cim, "-CIM | rho>=", rho,
                              " | all genes")
          out_dir   <- set_dir("all")
          run_ora(genes_all, label_all, out_dir, universe_genes,
                  universe_entrez, gene_table = tab_all)
          run_index[[length(run_index) + 1]] <- data.table(
            dir = out_dir, scope = sc, k_mode = k_mode, cim = cim, rho = rho,
            gene_set = "all", n_genes = n_all,
            enrichment_run = n_all >= MIN_GENES, label = label_all)
        }
        log_lines <- c(log_lines,
                       paste(sc, k_mode, cim, rho, "all", n_all,
                             if (n_all >= MIN_GENES) "OK" else "SKIP_too_few",
                             sep = "\t"))

        # Top N
        tab_top   <- get_scope_gene_table(sc, cim, rho, k_mode,
                                          top_n = TOP_N_GENES)
        genes_top <- tab_top$gene
        n_top     <- length(genes_top)
        message("  Top ", TOP_N_GENES, " genes: ", n_top)

        # Only run topN if it differs meaningfully from all (n_all > TOP_N_GENES)
        if (n_all > TOP_N_GENES && n_top >= MIN_GENES) {
          label_top <- paste0(sc, " | ", k_mode, " | ", cim, "-CIM | rho>=", rho,
                              " | top ", TOP_N_GENES)
          out_dir   <- set_dir(paste0("top", TOP_N_GENES))
          run_ora(genes_top, label_top, out_dir, universe_genes,
                  universe_entrez, gene_table = tab_top)
          run_index[[length(run_index) + 1]] <- data.table(
            dir = out_dir, scope = sc, k_mode = k_mode, cim = cim, rho = rho,
            gene_set = paste0("top", TOP_N_GENES), n_genes = n_top,
            enrichment_run = TRUE, label = label_top)
          log_lines <- c(log_lines,
                         paste(sc, k_mode, cim, rho, paste0("top", TOP_N_GENES),
                               n_top, "OK", sep = "\t"))
        } else if (n_all <= TOP_N_GENES) {
          message("  [top", TOP_N_GENES, " skipped - all genes (", n_all,
                  ") <= ", TOP_N_GENES, "; all-genes result is equivalent]")
          log_lines <- c(log_lines,
                         paste(sc, k_mode, cim, rho, paste0("top", TOP_N_GENES),
                               n_top, "SKIP_same_as_all", sep = "\t"))
        }
      }
    }
  }
}

# -- 12. WRITE LOG ------------------------------------------------------------
log_header <- paste("scope", "k_mode", "cim", "rho", "gene_set", "n_genes",
                    "status", sep = "\t")
writeLines(
  c(log_header, log_lines),
  file.path(OUT_ROOT, "cimic_ora_run.log")
)

# Map each short output directory back to its full description, so the
# abbreviated path names remain unambiguous.
if (length(run_index)) {
  run_index_dt <- rbindlist(run_index)
  fwrite(run_index_dt, file.path(OUT_ROOT, "run_index.csv"), na = "NA")
  message("Run index: ", file.path(OUT_ROOT, "run_index.csv"),
          " (", nrow(run_index_dt), " analyses)")
  longest <- max(nchar(normalizePath(run_index_dt$dir, winslash = "/",
                                     mustWork = FALSE)))
  message("Longest output directory path: ", longest,
          " chars (Windows limit 260; filenames add ~30)")
}

message("\n== ORA complete ==")
message("Outputs written to: ", OUT_ROOT)
message("Log: ", file.path(OUT_ROOT, "cimic_ora_run.log"))
