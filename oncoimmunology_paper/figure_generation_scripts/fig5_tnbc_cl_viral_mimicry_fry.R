# ============================================================================
# fig5_tnbc_cl_viral_mimicry_fry.R
# ============================================================================
# Purpose:
#   Generate a heatmap of Z‑score normalized expression for the "viral mimicry"
#   gene set using the TNBC_CL_Epirubicin dataset, perform FRY analysis (as in
#   fig5i_5l_tnbc_cl_icd_programs_fry.R) and per‑gene limma analysis (as in
#   tnbc_cl_figure_individual_gene_panels.R).
#
#   Input files are taken from the TNBC_CL_Epirubicin dataset folder.
#   The FRY analysis mirrors the existing ICD‑program script, but applied to the
#   viral mimicry gene set as a single program.
#   The limma analysis is performed by sourcing the individual‑gene panel script.
#
#   Outputs:
#     - Results/fig5_tnbc_cl_viral_mimicry_heatmap.png (heatmap with stars)
#     - Results/fig5_tnbc_cl_viral_mimicry_fry.png (FRY plot for the program)
#     - Results/tnbc_cl_individual_gene_panels/ (per‑gene limma panels)
#
# Dependencies: data.table, dplyr, tidyr, ggplot2, limma, msigdbr, purrr,
#               tibble, stringr
# ============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)
library(msigdbr)
library(purrr)
library(tibble)
library(stringr)

# ---------------------------------------------------------------------------
# 1. Load data (TNBC_CL_Epirubicin)
# ---------------------------------------------------------------------------
input_file <- "oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epi_clustered_plot_df.csv"

clustered_plot_df <- data.table::fread(input_file, data.table = FALSE) %>%
  as.data.frame()

if ("cluster_assignments" %in% colnames(clustered_plot_df)) {
  clustered_plot_df <- clustered_plot_df %>%
    dplyr::rename(CIMIC_Cluster = cluster_assignments)
} else if (!"CIMIC_Cluster" %in% colnames(clustered_plot_df)) {
  stop("Neither 'cluster_assignments' nor 'CIMIC_Cluster' found.")
}

# Convert to factor with Dys-CIM ("2") first, then Fun-CIM ("1")
clustered_plot_df$CIMIC_Cluster <- factor(clustered_plot_df$CIMIC_Cluster, levels = c("2", "1"))
if (!all(levels(clustered_plot_df$CIMIC_Cluster) %in% c("1", "2"))) {
  stop("Unexpected CIMIC_Cluster values; expected '1' and '2'.")
}

# ---------------------------------------------------------------------------
# 2. Define the viral mimicry gene set (same as fig2 script)
# ---------------------------------------------------------------------------
viral_mimicry_set <- c(
  "CGAS","STING1","IFI16","AIM2","ZBP1",
  "RIGI","IFIH1","DHX58","TLR3","TLR7","TLR8","TLR9",
  "MAVS","TBK1","IRF3","IRF7","STAT1","STAT2",
  "IFIT1","IFIT2","IFIT3","OAS1","OAS2","OAS3","RNASEL",
  "MX1","MX2","ISG15","RSAD2","BST2","IFI44","IFI44L","IFI6","EIF2AK2",
  "B2M","HLA-A","HLA-B","HLA-C","TAP1","TAPBP","PSMB8","PSMB9","NLRC5",
  "CXCL9","CXCL10","CXCL11","CCL5"
)

# ---------------------------------------------------------------------------
# 3. Identify gene columns (metadata vs expression)
# ---------------------------------------------------------------------------
metadata_cols <- c(
  "sample_id", "CIMIC_Cluster", "base_id", "PACMAP1", "PACMAP2",
  "UMAP1", "UMAP2", "PC1", "PC2", "PC3", "V1"
)

gene_cols <- setdiff(colnames(clustered_plot_df), metadata_cols)

viral_genes_present <- intersect(viral_mimicry_set, gene_cols)
if (length(viral_genes_present) == 0) {
  stop("None of the viral mimicry genes were found in the dataset.")
}

# ---------------------------------------------------------------------------
# 4. Z‑score normalization (display only)
# ---------------------------------------------------------------------------
zscore_safe <- function(x) {
  s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(x)) else (x - m) / s
}

df_z <- clustered_plot_df %>%
  dplyr::mutate(across(all_of(viral_genes_present), zscore_safe))

# ---------------------------------------------------------------------------
# 5. Prepare data for heatmap (mean Z‑score per cluster per gene)
# ---------------------------------------------------------------------------
heatmap_df <- clustered_plot_df %>%
  dplyr::select(CIMIC_Cluster, all_of(viral_genes_present)) %>%
  tidyr::pivot_longer(cols = -CIMIC_Cluster, names_to = "Gene", values_to = "ZScore") %>%
  dplyr::group_by(CIMIC_Cluster, Gene) %>%
  dplyr::summarise(mean_z = mean(ZScore, na.rm = TRUE), .groups = "drop")

# ---------------------------------------------------------------------------
# 6. Plot heatmap (same style as fig2 script)
# ---------------------------------------------------------------------------
p_heatmap <- ggplot(heatmap_df, aes(x = CIMIC_Cluster, y = Gene, fill = mean_z)) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_gradient2(
    low = "#313695",
    mid = "white",
    high = "#A50026",
    midpoint = 0,
    name = "Mean\nΔ GE"
  ) +
  scale_x_discrete(labels = c("2" = "Dys\nCIM", "1" = "Fun\nCIM"), expand = expansion(add = c(0.5, 1.5))) +
  scale_y_discrete(expand = expansion(add = c(0, 0))) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_void(base_size = 14) +
  theme(
    axis.text.x = element_text(face = "bold", size = 16, colour = "black"),
    axis.text.y = element_text(face = "bold", size = 14, colour = "black"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(face = "bold", size = 16),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
    legend.box.spacing = unit(0, "pt")
  ) +
  coord_cartesian(clip = "off")

# ---------------------------------------------------------------------------
# 7. FRY analysis with duplicateCorrelation
# ---------------------------------------------------------------------------

# Define the viral mimicry program as a single gene set for FRY
viral_program <- list(
  ViralMimicry = viral_genes_present
)

# Expression matrix for FRY / limma (genes x samples)
expr_matrix <- clustered_plot_df %>%
  dplyr::select(
    sample_id,
    all_of(gene_cols)
  ) %>%
  tibble::column_to_rownames("sample_id") %>%
  t() %>%
  as.matrix()

storage.mode(expr_matrix) <- "numeric"

# Replace non-finite values with NA
expr_matrix[!is.finite(expr_matrix)] <- NA_real_

# ---------------------------------------------------------------------------
# Extract cell line from sample_id
#
# Example:
#   delta_BT549_EPIRUBICIN_1
#   delta_BT549_EPIRUBICIN_2
#   delta_BT549_EPIRUBICIN_3
#
# becomes:
#   BT549
# ---------------------------------------------------------------------------

clustered_plot_df$cell_line <- sub(
  "^delta_(.+?)_EPIRUBICIN_.*$",
  "\\1",
  clustered_plot_df$sample_id
)

# Check extraction
message(
  "Cell lines identified: ",
  paste(
    sort(unique(clustered_plot_df$cell_line)),
    collapse = ", "
  )
)

# Make sure extraction worked
if (
  any(
    is.na(clustered_plot_df$cell_line) |
    clustered_plot_df$cell_line == clustered_plot_df$sample_id
  )
) {
  stop(
    "Cell-line extraction failed for one or more sample IDs. ",
    "Expected format: delta_CELL_LINE_EPIRUBICIN_REPLICATE."
  )
}

# ---------------------------------------------------------------------------
# Design / contrast
#
# Cluster 1 = Dys-CIM
# Cluster 2 = Fun-CIM
#
# Contrast = Dys-CIM - Fun-CIM
# ---------------------------------------------------------------------------

cluster_factor <- factor(
  clustered_plot_df$CIMIC_Cluster,
  levels = c("2", "1")
)

design <- stats::model.matrix(
  ~ 0 + cluster_factor
)

colnames(design) <- c(
  "Dys_CIM",
  "Fun_CIM"
)

contrast <- limma::makeContrasts(
  Dys_CIM - Fun_CIM,
  levels = design
)

# ---------------------------------------------------------------------------
# duplicateCorrelation
# ---------------------------------------------------------------------------

block_factor <- factor(
  clustered_plot_df$cell_line
)

dc <- limma::duplicateCorrelation(
  expr_matrix,
  design = design,
  block = block_factor
)

message(
  sprintf(
    "duplicateCorrelation consensus correlation = %.4f",
    dc$consensus.correlation
  )
)

# ---------------------------------------------------------------------------
# FRY
# ---------------------------------------------------------------------------

viral_fry_res <- limma::fry(
  y = expr_matrix,
  index = viral_program,
  design = design,
  contrast = contrast,
  block = block_factor,
  correlation = dc$consensus.correlation
)

# Extract significance label
viral_sig <- viral_fry_res %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Program") %>%
  dplyr::mutate(
    signif_label = dplyr::case_when(
      PValue <= 0.0001 ~ "****",
      PValue <= 0.001  ~ "***",
      PValue <= 0.01   ~ "**",
      PValue <= 0.05   ~ "*",
      TRUE             ~ ""
    )
  ) %>%
  dplyr::pull(signif_label)
  
# ---------------------------------------------------------------------------
# 8. Plot FRY result – styled like the ICD program boxplots in fig5i
# ---------------------------------------------------------------------------

# Compute program scores for plotting (mean Z‑score across the viral genes)
score_program_z <- function(data, genes) {
  genes <- intersect(genes, colnames(data))
  if (length(genes) == 0) {
    stop("No genes found for program.")
  }
  data.frame(
    Score = rowMeans(data %>% dplyr::select(all_of(genes)), na.rm = TRUE)
  )
}

program_scores <- score_program_z(df_z, viral_genes_present) %>%
  dplyr::mutate(CIMIC_Cluster = clustered_plot_df$CIMIC_Cluster,
                cell_line = clustered_plot_df$cell_line)

# Define shape palette for cell‑line replicates (mirrors fig5i)
shape_values <- c(16, 17, 15, 18, 8, 3, 4, 7, 10)
cell_line_shapes <- setNames(
  rep(shape_values, length.out = length(unique(program_scores$cell_line))),
  sort(unique(program_scores$cell_line))
)

p_fry <- ggplot(program_scores,
                aes(x = CIMIC_Cluster, y = Score, fill = CIMIC_Cluster)) +
  geom_boxplot(outlier.shape = NA, colour = "black", width = 0.6, alpha = 1) +
  geom_jitter(aes(shape = cell_line),
              width = 0.15, size = 7.5, colour = "black", alpha = 1) +
  scale_shape_manual(name = "Cell line", values = cell_line_shapes) +
  scale_x_discrete(labels = c("2" = "Dys\nCIM", "1" = "Fun\nCIM"),
                   expand = expansion(add = c(0.5, 0.5))) +
  scale_fill_manual(values = c("2" = "#F8766D", "1" = "#00BFC4"),
                    guide = "none") +
  labs(title = "Viral Mimicry", x = NULL, y = "Z-Score of\nProgram Activation") +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x = element_text(face = "bold", size = 24, colour = "black"),
    axis.text.y = element_text(face = "bold", size = 24, colour = "black"),
    axis.title.y = element_text(face = "bold", size = 24, colour = "black"),
    plot.title = element_text(face = "bold", size = 30, hjust = 0.5),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  ) +
  geom_text(aes(x = 1.5,
                y = max(program_scores$Score, na.rm = TRUE) * 0.9,
                label = viral_sig),
            size = 15, fontface = "bold")

ggsave("oncoimmunology_paper/Results/fig5_tnbc_cl_viral_mimicry_fry.png",
       plot = p_fry, width = 8, height = 6, dpi = 300, bg = "white")


# ---------------------------------------------------------------------------
# 9. Per-gene limma analysis with duplicateCorrelation
# ---------------------------------------------------------------------------
#
# This replaces sourcing the individual-gene panel script.
#
# We use:
#   - duplicateCorrelation for cell-line replicates
#   - robust = TRUE empirical Bayes moderation
#   - trend = TRUE
#
# Contrast:
#   Dys-CIM - Fun-CIM
#
# Therefore:
#   positive logFC = higher in Dys-CIM
#   negative logFC = higher in Fun-CIM
# ---------------------------------------------------------------------------

# Fit the linear model with cell-line blocking
limma_fit_gene <- limma::lmFit(
  expr_matrix,
  design = design,
  block = block_factor,
  correlation = dc$consensus.correlation
)

# Apply contrast
limma_fit_gene <- limma::contrasts.fit(
  limma_fit_gene,
  contrast
)

# Robust empirical-Bayes moderation
limma_fit_gene <- limma::eBayes(
  limma_fit_gene,
  trend = TRUE,
  robust = TRUE
)

# Extract ALL genes without sorting
limma_results <- limma::topTable(
  limma_fit_gene,
  coef = 1,
  number = Inf,
  sort.by = "none"
)

limma_results$Gene <- rownames(limma_results)

# Keep only viral mimicry genes
limma_stats <- limma_results %>%
  dplyr::filter(
    Gene %in% viral_genes_present
  ) %>%
  dplyr::select(
    Gene,
    logFC,
    AveExpr,
    t,
    P.Value,
    adj.P.Val,
    B
  )

message(
  sprintf(
    "Per-gene duplicateCorrelation limma analysis completed for %d viral mimicry genes.",
    nrow(limma_stats)
  )
)

# Display the most significant genes
print(
  limma_stats %>%
    dplyr::arrange(P.Value) %>%
    dplyr::select(
      Gene,
      logFC,
      t,
      P.Value,
      adj.P.Val
    )
)


# ---------------------------------------------------------------------------
# 10. Add significance stars to heatmap
# ---------------------------------------------------------------------------

sig_labels <- limma_stats %>%
  dplyr::mutate(
    signif_label = dplyr::case_when(
      P.Value <= 0.0001 ~ "****",
      P.Value <= 0.001  ~ "***",
      P.Value <= 0.01   ~ "**",
      P.Value <= 0.05   ~ "*",
      TRUE              ~ ""
    )
  ) %>%
  dplyr::select(
    Gene,
    signif_label
  )

# Add significance information
heatmap_df <- heatmap_df %>%
  dplyr::left_join(
    sig_labels,
    by = "Gene"
  )

# Preserve gene ordering
heatmap_df$Gene <- factor(
  heatmap_df$Gene,
  levels = rev(
    sort(
      unique(
        as.character(heatmap_df$Gene)
      )
    )
  )
)

# Create a separate data frame for the stars
star_df <- heatmap_df %>%
  dplyr::distinct(
    Gene,
    signif_label
  ) %>%
  dplyr::mutate(
    y = as.numeric(Gene)
  )

# Add stars
p_heatmap <- p_heatmap +

  geom_text(
    data = star_df,
    aes(
      x = 2.55,
      y = y,
      label = signif_label
    ),
    inherit.aes = FALSE,
    size = 10,
    fontface = "bold",
    hjust = 0,
    colour = "black"
  )

# Save heatmap
heatmap_path <- "oncoimmunology_paper/Results/fig5_tnbc_cl_viral_mimicry_heatmap.png"
ggsave(filename = heatmap_path, plot = p_heatmap, width = 8, height = 10, dpi = 300, bg = "white")
message("✓ Saved heatmap: ", heatmap_path)

# ---------------------------------------------------------------------------
# End of script
# ---------------------------------------------------------------------------
