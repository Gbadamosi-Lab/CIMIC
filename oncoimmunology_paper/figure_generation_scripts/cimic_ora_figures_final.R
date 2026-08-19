# =============================================================================
# cimic_ora_figures_final.R
# -----------------------------------------------------------------------------
# One publication ORA bar-plot from a pre-computed CSV under
# Results/cimic_ora_signatures/. Adapted from fig2c2d_ora_figure_final.R:
# same visual language, but every input and threshold is a toggle below.
#
# Edit the TOGGLES block, run, get one figure. To make another, change SCOPE /
# CIM / DATABASE (or point IN_CSV somewhere else) and run again.
#
# Defaults are set for the Dys-CIM signature, which is small (ALL3 = 6 genes,
# NKI_NEO = 22), so almost nothing clears FDR 0.05:
#
#   scope     database   terms   FDR<.05   FDR<.15
#   ALL3      GO BP        113         0         8
#   NKI_NEO   GO BP        307         0         4
#   NKI_CL    GO BP       2178        87       154
#
# Hence FDR_CUTOFF = 0.15. MIN_COUNT is 1, not 2: 105 of the 113 ALL3 GO BP
# terms have Count == 1, and a Count >= 2 filter at FDR 0.05 is exactly why
# cimic_ora_signatures.R produced no GO BP figure for that run.
#
# BH ties are heavy at this size (the top eight ALL3 GO BP terms all have
# p.adjust = 0.116), so sorting by p.adjust would order them arbitrarily.
# SORT_BY is therefore "pvalue".
#
# Colour matches the original: the MOST significant term is blue, the least is
# red. Swap FILL_SIG / FILL_NONSIG for the conventional "significant = red".
# (The original wrote low="red", high="blue", trans="reverse"; `trans=` is
# defunct in ggplot2 4.0.0, so the same mapping is expressed directly.)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(stringr); library(scales)
})

# =============================================================================
# TOGGLES
# =============================================================================
ORA_ROOT <- "oncoimmunology_paper/Results/cimic_ora_signatures"
OUT_DIR  <- "oncoimmunology_paper/Results/cimic_ora_figures"

# ---- Which result table --------------------------------------------------
SCOPE    <- "ALL3"      # ALL3 | NKI_NEO | NKI_CL | NEO_CL
K_MODE   <- "kAllC"     # k2 | kAll | kAllC
RHO      <- "0.3"
CIM      <- "Dys"       # Dys | Fun
GENE_SET <- "all"       # all | top200
DATABASE <- "ora_GO_BP" # ora_GO_BP | ora_GO_BP_simplified | ora_reactome |
                        # ora_hallmark | ora_proteostasis
DB_LABEL <- "GO Biological Process"   # used in the title only

# Built from the pieces above. Override directly to read any other CSV.
IN_CSV <- file.path(ORA_ROOT, SCOPE, K_MODE,
                    paste0("rho", RHO, "_", CIM, "_", GENE_SET),
                    paste0(DATABASE, ".csv"))

# ---- Term selection ------------------------------------------------------
FDR_CUTOFF    <- 0.15     # p.adjust < this. NULL to disable.
PVALUE_CUTOFF <- NULL     # pvalue   < this. NULL to disable.
MIN_COUNT     <- 1        # minimum genes per term
TOP_N_TERMS   <- 15       # bars kept after sorting
SORT_BY       <- "pvalue" # pvalue | p.adjust | Count | FoldEnrichment
# Optional exact-match term list, like the original script. NULL = use top N.
KEEP_TERMS    <- NULL

# ---- Encoding and layout -------------------------------------------------
X_METRIC    <- "Count"     # Count | FoldEnrichment | RichFactor
FILL_BY     <- "p.adjust"  # p.adjust | pvalue
FILL_SIG    <- "blue"      # colour of the most significant term
FILL_NONSIG <- "red"       # colour of the least significant term
WRAP_WIDTH  <- 40
FIG_WIDTH   <- 14
# NA scales height with the number of bars, which matters because wrapped
# labels take two lines at size 25 and collide in a fixed short canvas.
FIG_HEIGHT  <- NA
HEIGHT_PER_BAR <- 0.85     # only used when FIG_HEIGHT is NA
FIG_DPI     <- 300
SAVE_FORMATS <- c("png", "pdf")
SHOW_TITLE  <- TRUE
OUT_NAME    <- paste0("ora_fig_", SCOPE, "_", CIM, "_", DATABASE)

# =============================================================================
# BUILD THE FIGURE
# =============================================================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(IN_CSV))
  stop("Input CSV not found:\n  ", IN_CSV,
       "\n  Check SCOPE / K_MODE / RHO / CIM / GENE_SET / DATABASE, or see ",
       file.path(ORA_ROOT, "run_index.csv"), " for what exists.")

df <- as.data.frame(fread(IN_CSV, showProgress = FALSE))
need <- c("Description", "pvalue", "p.adjust", "Count")
if (!all(need %in% names(df)))
  stop("Missing column(s) in ", basename(IN_CSV), ": ",
       paste(setdiff(need, names(df)), collapse = ", "))
if (!X_METRIC %in% names(df))
  stop("X_METRIC '", X_METRIC, "' is not a column. Available: ",
       paste(names(df), collapse = ", "))
n_total <- nrow(df)

df <- df[!is.na(df$Description) & nzchar(df$Description), , drop = FALSE]
df <- df[df$Count >= MIN_COUNT, , drop = FALSE]
n_count <- nrow(df)
if (!is.null(FDR_CUTOFF))
  df <- df[!is.na(df$p.adjust) & df$p.adjust < FDR_CUTOFF, , drop = FALSE]
if (!is.null(PVALUE_CUTOFF))
  df <- df[!is.na(df$pvalue) & df$pvalue < PVALUE_CUTOFF, , drop = FALSE]
if (!is.null(KEEP_TERMS))
  df <- df[df$Description %in% KEEP_TERMS, , drop = FALSE]

if (nrow(df) == 0)
  stop("No terms pass the filters (", n_total, " rows -> ", n_count,
       " after Count >= ", MIN_COUNT, " -> 0 after significance).",
       "\n  Raise FDR_CUTOFF, lower MIN_COUNT, or pick another database.")

# Descending for effect-size keys, ascending for p-values.
df <- df[order(df[[SORT_BY]],
               decreasing = SORT_BY %in% c("Count", "FoldEnrichment",
                                           "RichFactor")), , drop = FALSE]
if (!is.null(TOP_N_TERMS) && nrow(df) > TOP_N_TERMS)
  df <- df[seq_len(TOP_N_TERMS), , drop = FALSE]

# unique() guards against two terms wrapping to the same string, which would
# otherwise silently drop a bar.
df$Description <- str_wrap(df$Description, width = WRAP_WIDTH)
df$Description <- factor(df$Description, levels = rev(unique(df$Description)))

x_label <- c(Count = "Gene count", FoldEnrichment = "Fold enrichment",
             RichFactor = "Rich factor")[[X_METRIC]]
fill_label <- c(p.adjust = "FDR-adjusted p", pvalue = "Nominal p")[[FILL_BY]]

# A single-valued fill range (common with heavy BH ties) gives an empty
# gradient, so fall back to one flat colour and say so in the legend.
fill_range <- range(df[[FILL_BY]], na.rm = TRUE)
fill_scale <- if (!all(is.finite(fill_range)) ||
                  isTRUE(all.equal(fill_range[1], fill_range[2]))) {
  scale_fill_gradient(low = FILL_SIG, high = FILL_SIG,
                      name = paste0(fill_label, "\n(all tied)"),
                      labels = label_scientific(digits = 2))
} else {
  scale_fill_gradient(low = FILL_SIG, high = FILL_NONSIG, name = fill_label,
                      labels = label_scientific(digits = 2))
}

p <- ggplot(df, aes(x = .data$Description, y = .data[[X_METRIC]])) +
  geom_col(aes(fill = .data[[FILL_BY]]), colour = "black") +
  coord_flip() +
  labs(
    x = NULL, y = x_label,
    title = if (SHOW_TITLE) paste0(DB_LABEL, " — ", CIM, "-CIM, ", SCOPE)
            else NULL,
    subtitle = if (SHOW_TITLE) paste0(
      "FDR<", if (is.null(FDR_CUTOFF)) "none" else FDR_CUTOFF,
      ", Count>=", MIN_COUNT, ", top ", TOP_N_TERMS, " by ", SORT_BY) else NULL
  ) +
  fill_scale +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text.y  = element_text(size = 25, colour = "black", face = "bold"),
    axis.text.x  = element_text(size = 25, colour = "black", face = "bold"),
    axis.title.x = element_text(size = 25, colour = "black", face = "bold"),
    legend.text  = element_text(size = 20),
    legend.title = element_text(size = 20),
    plot.title    = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 13, colour = "grey30")
  )

fig_h <- if (is.na(FIG_HEIGHT)) max(4, HEIGHT_PER_BAR * nrow(df) + 2.5) else
  FIG_HEIGHT
for (fmt in SAVE_FORMATS) {
  f <- file.path(OUT_DIR, paste0(OUT_NAME, ".", fmt))
  ggsave(f, plot = p, width = FIG_WIDTH, height = fig_h, dpi = FIG_DPI,
         bg = "white")
  message("Saved: ", f)
}

# The exact terms drawn, so the figure can be curated via KEEP_TERMS later.
fwrite(df[, intersect(c("ID", "Description", "Count", "GeneRatio", "BgRatio",
                        "FoldEnrichment", "pvalue", "p.adjust", "geneID"),
                      names(df))],
       file.path(OUT_DIR, paste0(OUT_NAME, "_terms.csv")), na = "NA")

message(nrow(df), " terms plotted (of ", n_total, " in the table). ",
        "FDR<", if (is.null(FDR_CUTOFF)) "none" else FDR_CUTOFF,
        ", Count>=", MIN_COUNT, ", sorted by ", SORT_BY, ".")
