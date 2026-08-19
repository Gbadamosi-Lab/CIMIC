# =============================================================================
# cimic_consensus_volcano.R
# -----------------------------------------------------------------------------
# Volcano plots of the Fun-CIM and Dys-CIM consensus genes produced by
# cimic_overlap_meta_rankbiserial.R.
#
# Adapted from figure2_consensus_volcano.R (9_Projects/combined_cimic_analysis_bc).
# Same axes, encodings, highlight scheme, and styling. Structural changes:
#   * source is meta_results.csv, not cimic_gene_analysis_FDR.xlsx
#   * "tier" -> "scope"; the script loops over every scope x rho rather than
#     plotting one hard-coded tier
#
# Filter : fisher_sig == "YES" & stouffer_sig == "YES" (both tests)
#
# X-axis : mean_abs_rrb           (effect size magnitude)
# Y-axis : -log10(stouffer_padj)  (direction-consistent significance)
# Size   : k_datasets
# Colour : Fun-CIM (red) / Dys-CIM (blue)
#
# Highlighted genes:
#   Fun-CIM  : immunostimulatory ISGs + immune genes (curated list below)
#   Dys-CIM  : mitochondrial proteostasis, cytonuclear proteostasis and
#              translation branch genes from proteostatic_network_genes.csv
#              (highlighting is skipped if that file is unavailable)
#   Top N per highlight category labelled
#
# K-modes (see section 2a): how many datasets must back a gene.
#   k_min2            k_datasets >= 2  (2 OR 3 datasets for the ALL3 scope)
#   k_all             k_datasets == n_scope_datasets  (all three)
#   k_all_concordant  all three AND all agreeing on direction
#
# Outputs -> oncoimmunology_paper/Results/cimic_overlap_figures/ :
#   volcano_<scope>_rho<rho>_<k_mode>.png / .pdf
#   volcano_<scope>_rho<rho>_<k_mode>_plotted_genes.csv
#   volcano_<scope>_rho<rho>_<k_mode>_labelled_genes.csv
#
# Run from the repository root, AFTER cimic_overlap_meta_rankbiserial.R:
#   Rscript oncoimmunology_paper/cimic_consensus_volcano.R
# =============================================================================

suppressMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
  library(patchwork)
})

# -- 1. PATHS -----------------------------------------------------------------
BASE     <- "oncoimmunology_paper"
META_CSV <- file.path(BASE, "Results", "cimic_meta_analysis", "meta_results.csv")
# Proteostasis branch annotation lives outside this repository. Optional: without
# it the Dys-CIM highlight groups collapse to "background" and the figure still
# renders (top genes are still labelled).
PROT_CSV <- paste0("D:/OneDrive - University of Florida/Gbadamosi Lab/",
                   "Mohammed Gbadamosi/9_Projects/tonic_cimic_analysis/",
                   "proteostatic_network_genes.csv")
OUT_DIR  <- file.path(BASE, "Results", "cimic_overlap_figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(META_CSV))
  stop("meta_results.csv not found. Run cimic_overlap_meta_rankbiserial.R ",
       "first.\n  Expected at: ", META_CSV)

# -- 2. PARAMETERS - all editable ---------------------------------------------
RHO_LEVELS <- c(0.3) # can change to 0.3 or 0.6
SCOPES     <- c("ALL3_NKI.NEO.TNBCCL", "NKI.NEO", "NKI.TNBCCL", "NEO.TNBCCL") 
# "ALL3_NKI.NEO.TNBCCL", "NKI.NEO", "NKI.TNBCCL", "NEO.TNBCCL" for scopes
N_LABEL    <- 20      # max genes to label per highlight group

# -- 2a. K-MODE: how many datasets must back a gene ---------------------------
# NOTE ON WHAT k MEANS. `k_datasets` counts the datasets in which the gene's
# |r_rb| >= rho_threshold. It is MAGNITUDE ONLY: it does not require the gene to
# be significant within any individual dataset, and it does not require the
# datasets to agree on direction -- a gene at r_rb = -0.7 in one dataset and
# +0.7 in another clears rho = 0.6 in both while pointing at opposite CIMs.
# Direction consistency is enforced separately, by the one-sided Fisher/Stouffer
# p (the wrong-way dataset gets p_os near 1, penalising the meta statistic).
#
# `k_concordant` counts only the contributors pointing toward the CIM, so
# direction agreement is expressible too. Available modes:
#
#   "k_min2"           k_datasets >= 2. Permissive: for the all-three scope this
#                      admits genes backed by 2 OR 3 datasets.
#   "k_all"            k_datasets == n_scope_datasets. "Require all three" for
#                      the ALL3 scope. Magnitude in every dataset; direction
#                      still only enforced statistically.
#   "k_all_concordant" all datasets AND direction_concordant == "YES". The strict
#                      "same direction in all three" reading.
#
# Every entry in K_MODES produces its own figure and tables, tagged with the mode.
K_MODES <- c("k_all_concordant")

# For a 2-dataset scope, "k_all" selects the same genes as "k_min2" (k can only
# be 2). TRUE skips the duplicate figure; FALSE writes it anyway.
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

# Font sizes - increased for publication quality
BASE_SIZE        <- 22
TITLE_SIZE       <- 26
AXIS_TITLE_SIZE  <- 24
AXIS_TEXT_SIZE   <- 20
LEGEND_TITLE_SZ  <- 20
LEGEND_TEXT_SZ   <- 18
LABEL_SIZE       <- 7

# Point sizes
PT_SIZE_RANGE    <- c(5, 8)   # min/max for k_datasets scaling
PT_ALPHA_BG      <- 0.25      # background (non-highlighted) points
PT_ALPHA_HL      <- 0.90      # highlighted points

# Colours
COL_FUN_BG      <- "#99CCFF"   # Fun-CIM
COL_DYS_BG      <- "#FF9999"   # Dys-CIM

# Portrait orientation dimensions (width < height)
FIG_WIDTH       <- 12
FIG_HEIGHT      <- 16

# -- 3. CURATED IMMUNOSTIMULATORY GENE LIST (Fun-CIM) -------------------------
# Well-established ISGs and immune activation genes; intersected with the data.
ISG_GENES <- c(
  # Type I/II interferon stimulated genes
  "STAT1", "STAT2", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL",
  "ISG15", "ISG20", "IFIT1", "IFIT2", "IFIT3", "IFIT5",
  "IFIH1", "DDX58", "DHX58", "HERC5", "HERC6",
  "IRF1", "IRF3", "IRF7", "IRF9",
  # Antigen presentation
  "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1",
  "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9",
  # Chemokines / cytokines
  "CXCL10", "CXCL9", "CXCL11", "CCL5", "CXCL8",
  "IL6", "IL1B", "TNF", "IFNG", "IL2",
  # Innate immune / NF-kB
  "TLR3", "TLR4", "STING1", "CGAS", "MAVS",
  "NFKB1", "NFKB2", "RELA", "RELB",
  # T/NK cell activation markers
  "CD8A", "CD8B", "GZMB", "PRF1", "GNLY",
  "CD274", "PDCD1LG2", "CD80", "CD86",
  # Additional immune
  "PTPRC", "ITGAM", "FCGR3A", "NCAM1",
  "IL18", "IL12A", "IL12B", "TNFSF10",
  "CXCL13", "CXCL12", "IFNB1", "IFNA1", "IFNA2", "LAG3", "CTLA4"
)

# Immediate-early genes (IEGs), labelled for Fun-CIM regardless of ranking
IEG_GENES <- c(
  "FOS", "JUN", "EGR1", "EGR2", "NR4A1", "NR4A2", "IER2", "IER3",
  "ATF3", "DUSP1", "DUSP2"
)

# -- 4. LOAD PROTEOSTASIS GENE SETS (Dys-CIM highlight) -----------------------
PROT_BRANCHES <- c("Mitochondrial proteostasis",
                   "Cytonuclear proteostasis",
                   "Translation")
prot_genes <- NULL
if (file.exists(PROT_CSV)) {
  message("Loading proteostasis gene sets ...")
  prot_raw <- fread(PROT_CSV, showProgress = FALSE)
  setnames(prot_raw,
           old = intersect(names(prot_raw),
                           c("Gene Symbol", "gene_symbol", "Gene.Symbol")),
           new = "gene_symbol", skip_absent = TRUE)
  setnames(prot_raw,
           old = intersect(names(prot_raw), c("Branch", "branch")),
           new = "branch", skip_absent = TRUE)

  if (all(c("gene_symbol", "branch") %in% names(prot_raw))) {
    prot_genes <- unique(prot_raw[
      branch %in% PROT_BRANCHES & !is.na(gene_symbol) & gene_symbol != "",
      .(gene_symbol, branch)
    ])
    message("  Proteostasis genes loaded: ", nrow(prot_genes),
            " across ", length(unique(prot_genes$branch)), " branches")
    print(prot_genes[, .N, by = branch])
  } else {
    message("  [skip] proteostasis file lacks gene_symbol/branch columns")
  }
} else {
  message("[note] proteostasis annotation not found - Dys-CIM branch ",
          "highlighting disabled.\n  Looked for: ", PROT_CSV)
}

branch_genes <- function(b) {
  if (is.null(prot_genes)) return(character(0))
  prot_genes[branch == b, gene_symbol]
}

# -- 5. LOAD META RESULTS -----------------------------------------------------
message("\nLoading meta results: ", META_CSV)
meta <- fread(META_CSV, showProgress = FALSE)

req_cols <- c("scope", "cim", "rho_threshold", "gene", "k_datasets",
              "n_scope_datasets", "k_concordant", "direction_concordant",
              "fisher_sig", "stouffer_sig", "fisher_padj", "stouffer_padj",
              "stouffer_Z", "mean_abs_rrb")
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
message("  Scopes: ", paste(scopes_present, collapse = ", "))

# -- 6. PANEL BUILDER ---------------------------------------------------------
# One panel for a single CIM class. `df` is the full (both-CIM) frame so the two
# panels share an x-axis range.
build_panel <- function(df, label_df, cim_sel, scope_name, rho_val, x_range,
                        k_mode) {

  cim_name  <- if (cim_sel == "Fun") "Fun-CIM" else "Dys-CIM"
  fill_col  <- if (cim_sel == "Fun") COL_FUN_BG else COL_DYS_BG
  caption_i <- if (cim_sel == "Fun") {
    paste0("Fun-CIM highlights: interferon-stimulated & immunostimulatory ",
           "genes\nTop ", N_LABEL, " genes labelled by Stouffer p")
  } else {
    paste0("Dys-CIM highlights: mitochondrial proteostasis, cytonuclear ",
           "proteostasis,\nand translation genes (proteostatic_network_genes",
           ".csv)\nTop ", N_LABEL, " genes labelled by Stouffer p")
  }

  d_cim <- df       %>% dplyr::filter(cim == cim_sel)
  d_bg  <- d_cim    %>% dplyr::filter(highlight == "background")
  d_hl  <- d_cim    %>% dplyr::filter(highlight != "background")
  d_lab <- label_df %>% dplyr::filter(cim == cim_sel)

  # k_datasets can only be 2 or 3 with three datasets; keep breaks integral so
  # the legend does not invent fractional dataset counts.
  k_breaks <- sort(unique(df$k_datasets))

  base_title    <- paste0("Consensus CIMIC Gene Signatures - ", cim_name)
  base_subtitle <- bquote(paste(.(scope_name), " | ", rho, " >= ", .(rho_val),
                                " | ", .(k_mode_label[[k_mode]]),
                                " | Fisher + Stouffer BH-padj < 0.05"))

  # A CIM can legitimately have zero genes at a given scope+rho (e.g. Dys-CIM in
  # NEO+TNBC_CL). Return an explicit placeholder rather than an empty panel full
  # of "no shared levels" scale warnings.
  if (nrow(d_cim) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, size = 8, fontface = "bold",
                 colour = "grey30",
                 label = paste0("No ", cim_name, " genes pass the\n",
                                "Fisher + Stouffer filter at rho >= ",
                                rho_val)) +
        labs(title = base_title, subtitle = base_subtitle) +
        theme_void(base_size = BASE_SIZE) +
        theme(
          plot.title = element_text(face = "bold", size = TITLE_SIZE,
                                    colour = "black", hjust = 0.5),
          plot.subtitle = element_text(size = BASE_SIZE - 2, colour = "black",
                                       hjust = 0.5),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(t = 10, r = 10, b = 15, l = 10, unit = "pt")
        )
    )
  }

  p <- ggplot(d_cim, aes(x = mean_abs_rrb, y = neg_log10_stouffer))
  # Layers are added conditionally: an empty data frame carries no fill levels,
  # which would otherwise trip a spurious scale_fill_manual warning.
  if (nrow(d_bg) > 0)
    p <- p + geom_point(data = d_bg, aes(fill = cim_label, size = k_datasets),
                        colour = "black", alpha = PT_ALPHA_BG,
                        shape = 21, stroke = 0.3)
  if (nrow(d_hl) > 0)
    p <- p + geom_point(data = d_hl, aes(fill = cim_label, size = k_datasets),
                        colour = "black", alpha = PT_ALPHA_HL,
                        shape = 21, stroke = 0.3)
  if (nrow(d_lab) > 0)
    p <- p + ggrepel::geom_text_repel(
      data = d_lab, aes(label = gene), colour = "black", size = LABEL_SIZE,
      fontface = "bold", max.overlaps = 300, box.padding = 1,
      point.padding = 0.5, segment.size = 0.3, segment.color = "grey50",
      min.segment.length = 0.3, show.legend = FALSE,  force = 20, force_pull = 0.5, 
      max.time = 2, max.iter = 100000
      )

  p +
    scale_fill_manual(values = setNames(fill_col, cim_name),
                      name = "Gene category",
                      guide = guide_legend(
                        override.aes = list(shape = 21, colour = "black",
                                            size = 5, alpha = 1),
                        title.position = "top")) +
    scale_size_continuous(range = PT_SIZE_RANGE, name = "Datasets\n(k)",
                          breaks = k_breaks,
                          guide = guide_legend(
                            title.position = "top",
                            override.aes = list(fill = fill_col,
                                                colour = "black",
                                                alpha = 0.8))) +
    scale_x_continuous(limits = x_range) +
    labs(title = base_title,
         subtitle = base_subtitle,
         x = expression(paste("Mean absolute rank-biserial correlation (|",
                              rho, "|)")),
         y = expression(-log[10](Stouffer~p[adj])),
         caption = caption_i) +
    theme_classic(base_size = BASE_SIZE) +
    theme(
      plot.title = element_text(face = "bold", size = TITLE_SIZE,
                                colour = "black", hjust = 0.5),
      plot.subtitle = element_text(size = BASE_SIZE - 2, colour = "black",
                                   hjust = 0.5),
      plot.caption = element_text(size = 10, colour = "grey30", hjust = 0),
      axis.title = element_text(face = "bold", size = AXIS_TITLE_SIZE,
                                colour = "black"),
      axis.text = element_text(face = "bold", size = AXIS_TEXT_SIZE,
                               colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.6),
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = LEGEND_TITLE_SZ,
                                  colour = "black"),
      legend.text = element_text(size = LEGEND_TEXT_SZ, colour = "black"),
      legend.background = element_rect(fill = "white", colour = "grey80",
                                       linewidth = 0.3),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(t = 10, r = 10, b = 15, l = 10, unit = "pt")
    )
}

# -- 7. MAIN LOOP: scope x k_mode x rho ---------------------------------------
for (sc in scopes_present) {

  # How many datasets this scope has, for the redundant-mode check below.
  n_sc <- unique(meta$n_scope_datasets[meta$scope == sc])[1]

  for (k_mode in K_MODES) {

  # For a 2-dataset scope, k_all selects exactly the same genes as k_min2.
  if (SKIP_REDUNDANT_K_MODE && k_mode == "k_all" && isTRUE(n_sc <= 2)) {
    message("\n[skip k-mode '", k_mode, "'] ", sc, " has ", n_sc,
            " datasets, so it is identical to 'k_min2'")
    next
  }

  for (rho_val in RHO_LEVELS) {

    tag <- paste0(sc, "_rho", rho_val, "_", k_mode)
    message("\n########## ", tag, " ##########")

    df <- meta[
      scope == sc &
        fisher_sig == "YES" &
        stouffer_sig == "YES" &
        abs(rho_threshold - rho_val) < 1e-8
    ] |>
      k_mode_filter(k_mode) |>
      as.data.frame() |>
      dplyr::filter(!is.na(stouffer_padj), !is.na(mean_abs_rrb)) |>
      dplyr::mutate(
        # +1e-300 keeps -log10 finite when padj underflows to 0
        neg_log10_stouffer = -log10(stouffer_padj + 1e-300),
        cim_label          = ifelse(cim == "Fun", "Fun-CIM", "Dys-CIM")
      )

    if (nrow(df) == 0) {
      message("  [skip] no genes pass the Fisher + Stouffer + ", k_mode,
              " filter")
      next
    }
    n_fun <- sum(df$cim == "Fun"); n_dys <- sum(df$cim == "Dys")
    message("  Genes for plot: ", nrow(df), " (Fun=", n_fun,
            ", Dys=", n_dys, ")")

    # -- highlight categories + proteostasis annotation ----------------------
    mito <- branch_genes("Mitochondrial proteostasis")
    cyto <- branch_genes("Cytonuclear proteostasis")
    tran <- branch_genes("Translation")

    df <- df |>
      dplyr::mutate(
        highlight = dplyr::case_when(
          cim == "Fun" & gene %in% ISG_GENES ~ "Fun-CIM: Immunostimulatory",
          cim == "Dys" & gene %in% mito ~ "Dys-CIM: Mitochondrial Proteostasis",
          cim == "Dys" & gene %in% cyto ~ "Dys-CIM: Cytonuclear Proteostasis",
          cim == "Dys" & gene %in% tran ~ "Dys-CIM: Translation",
          TRUE ~ "background"
        ),
        # Explicit proteostasis annotation for the exported table. (The upstream
        # script referenced these columns without ever creating them.)
        proteostasis_branch = dplyr::case_when(
          gene %in% mito ~ "Mitochondrial proteostasis",
          gene %in% cyto ~ "Cytonuclear proteostasis",
          gene %in% tran ~ "Translation",
          TRUE ~ NA_character_
        ),
        is_proteostatic = !is.na(proteostasis_branch)
      )

    message("  Highlight group sizes:")
    print(df |> dplyr::count(highlight) |> dplyr::arrange(dplyr::desc(n)))

    # -- save the full plotted-gene table -----------------------------------
    fwrite(df, file.path(OUT_DIR, paste0("volcano_", tag,
                                         "_plotted_genes.csv")), na = "NA")

    # -- label selection ----------------------------------------------------
    # (a) top N per highlight group
    label_df <- df |>
      dplyr::filter(highlight != "background") |>
      dplyr::group_by(highlight) |>
      dplyr::slice_min(order_by = stouffer_padj, n = N_LABEL,
                       with_ties = FALSE) |>
      dplyr::ungroup()

    # (b) all Fun-CIM IEGs regardless of ranking
    ieg_df <- df |> dplyr::filter(cim == "Fun", gene %in% IEG_GENES)

    # (c) top N Dys-CIM overall, so the most significant Dys markers appear
    #     even when they sit outside a highlight branch
    dys_top_df <- df |>
      dplyr::filter(cim == "Dys") |>
      dplyr::slice_min(order_by = stouffer_padj, n = N_LABEL,
                       with_ties = FALSE)

    # (d) top N from the background set, for a broader view
    top_overall_df <- df |>
      dplyr::filter(highlight == "background") |>
      dplyr::slice_min(order_by = stouffer_padj, n = N_LABEL,
                       with_ties = FALSE)

    label_df <- dplyr::bind_rows(label_df, ieg_df, dys_top_df,
                                 top_overall_df) |>
      dplyr::distinct()
    message("  Genes selected for labelling: ", nrow(label_df))

    # -- build + save figure -------------------------------------------------
    x_range <- range(df$mean_abs_rrb, na.rm = TRUE)
    # A single unique x value would make scale_x_continuous(limits=) degenerate.
    if (diff(x_range) == 0) x_range <- x_range + c(-0.01, 0.01)

    p_fun <- build_panel(df, label_df, "Fun", sc, rho_val, x_range, k_mode)
    p_dys <- build_panel(df, label_df, "Dys", sc, rho_val, x_range, k_mode)
    p <- p_fun / p_dys + patchwork::plot_layout(guides = "keep")

    out_stem <- file.path(OUT_DIR, paste0("volcano_", tag))
    ggsave(paste0(out_stem, ".png"), p, width = FIG_WIDTH,
           height = FIG_HEIGHT, dpi = 300, bg = "white")
    ggsave(paste0(out_stem, ".pdf"), p, width = FIG_WIDTH,
           height = FIG_HEIGHT, device = cairo_pdf, bg = "white")
    message("  Saved: ", out_stem, ".png / .pdf")

    # -- export labelled gene table -----------------------------------------
    out_table <- label_df |>
      dplyr::select(gene, cim, cim_label, highlight,
                    k_datasets, n_scope_datasets, k_concordant,
                    direction_concordant,
                    mean_abs_rrb, stouffer_padj, fisher_padj, stouffer_Z,
                    is_proteostatic, proteostasis_branch) |>
      dplyr::arrange(highlight, stouffer_padj)
    fwrite(out_table,
           file.path(OUT_DIR, paste0("volcano_", tag,
                                     "_labelled_genes.csv")), na = "NA")
  }
  }
}

message("\n== Done ==")
message("Figures and tables written to: ", OUT_DIR)
