# ============================================================================
# tnbc_cl_extract_correlated_genes.R
# ============================================================================
# Purpose: Using the TNBC_CL_Epirubicin correlation analysis results, extract the genes
#          significantly correlated with the Dys‑CIM and Fun‑CIM trajectories and
#          label them in the TNBC master induction dataframe.
#
# Input files:
#   1. oncoimmunology_paper/Results/tnbc_cl_correlation_analysis.csv
#   2. oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epirubicin_master_induction_df.csv
# Output files (fixed names, no timestamps):
#   - oncoimmunology_paper/Results/tnbc_cl_fun_cim_genes.csv
#   - oncoimmunology_paper/Results/tnbc_cl_dys_cim_genes.csv
#   - oncoimmunology_paper/Results/tnbc_cl_master_labeled.csv (optional)
# ============================================================================

library(data.table)
library(dplyr)

# ---------------------------------------------------------------------------
# Define paths (relative to repository root)
# ---------------------------------------------------------------------------

correlation_path <- "oncoimmunology_paper/Results/tnbc_cl_correlation_analysis.csv"
master_path      <- "oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epirubicin_master_induction_df.csv"

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cor_df   <- fread(correlation_path)
master_df <- fread(master_path)

# Ensure the gene identifier column exists in both tables. The correlation file
# uses "gene_id"; the master file is expected to contain a matching column –
# often named "gene_id" or "gene". We'll attempt to locate it.

gene_id_col_corr <- "gene_id"

gene_id_col_master <- intersect(c("gene_id", "gene", "gene_symbol"), colnames(master_df))
if (length(gene_id_col_master) == 0) {
  stop("No gene identifier column found in TNBC master induction dataframe.")
}

gene_id_col_master <- gene_id_col_master[1]

# ---------------------------------------------------------------------------
# Filter for significant Fun‑CIM and Dys‑CIM genes
# ---------------------------------------------------------------------------
# Extract the gene identifier and trajectory label for significant genes.
signif_genes <- cor_df %>%
  filter(trajectory %in% c("Fun-CIM", "Dys-CIM")) %>%
  select(all_of(gene_id_col_corr), trajectory)

# ---------------------------------------------------------------------------
# Merge with master dataframe to obtain full gene information and label
# ---------------------------------------------------------------------------
# Merge the full correlation data (including r_rb, p_value, padj, trajectory)
merged_df <- master_df %>%
  left_join(cor_df, by = setNames(gene_id_col_corr, gene_id_col_master))

# ---------------------------------------------------------------------------
# Write separate files for each trajectory (optional combined file)
# ---------------------------------------------------------------------------
fun_cim_df <- merged_df %>% filter(trajectory == "Fun-CIM")

dys_cim_df <- merged_df %>% filter(trajectory == "Dys-CIM")

# Write outputs to fixed filenames (no timestamps)
fun_cim_path <- "oncoimmunology_paper/Results/tnbc_cl_fun_cim_genes.csv"

dys_cim_path <- "oncoimmunology_paper/Results/tnbc_cl_dys_cim_genes.csv"

master_labeled_path <- "oncoimmunology_paper/Results/tnbc_cl_master_labeled.csv"

fwrite(fun_cim_df, fun_cim_path)

fwrite(dys_cim_df, dys_cim_path)

# Also write the combined labeled master dataframe for completeness
fwrite(merged_df, master_labeled_path)

message("Extraction complete. TNBC CL files written to Results folder.")
