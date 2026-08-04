# gene_sets.R
# Run once per session via source("gene_sets.R")
# Do NOT save CSVs — MSigDB versions change and your death_programs
# are hardcoded anyway. Keep everything in one place.

library(msigdbr)
library(dplyr)
library(tidyr)
library(tibble)

# ============================================================
# 1. MSigDB gene sets
# ============================================================

msig_df <- msigdbr::msigdbr(species = "Homo sapiens") %>%
  as.data.frame()

all_gene_set_names <- c(
  "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "GOBP_B_CELL_ACTIVATION",
  "GOBP_CELLULAR_RESPONSE_TO_STRESS",
  "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  "GOBP_INFLAMMASOME_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS",
  "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
  "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "GOBP_LEUKOCYTE_CHEMOTAXIS_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  "GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "GOBP_T_CELL_ACTIVATION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS",
  "REACTOME_CELLULAR_RESPONSE_TO_MITOCHONDRIAL_STRESS",
  "WP_MOLECULAR_PATHWAY_FOR_OXIDATIVE_STRESS",
  "WP_OXIDATIVE_STRESS_RESPONSE"
)

# ============================================================
# 2. Death programs (hardcoded — version stable)
# ============================================================

death_programs <- list(
  Apoptosis = c(
    "BAX", "BAK1", "BBC3", "PMAIP1", "BCL2L11",
    "APAF1", "CASP9", "CASP3", "CASP7",
    "FAS", "TNFRSF10B", "CASP8", "FADD", "BID",
    "DIABLO", "CYCS", "TP53AIP1", "BAD", "BMF", "HRK"
  ),
  Necroptosis = c(
    "RIPK1", "RIPK3", "MLKL", "ZBP1", "TICAM1",
    "TNFRSF1A", "FADD", "CASP8"
  ),
  Pyroptosis = c(
    "GSDME", "GSDMD", "CASP1", "CASP4", "CASP5",
    "NLRP3", "AIM2", "PYCARD", "IL1B", "IL18",
    "CASP3", "NLRC4"
  ),
  PANoptosis = c(
    "ZBP1", "AIM2", "RIPK3", "RIPK1", "CASP8",
    "CASP1", "FADD", "PYCARD", "IRF1"
  ),
  Ferroptosis = c(
    "ACSL4", "LPCAT3", "TFRC", "SAT1", "PTGS2"
  )
)

# ============================================================
# 3. Build MSigDB named list with prefix filtering
# ============================================================

remove_prefixes <- c(
  "IGHV", "IGHD", "IGHJ", "IGHC", "IGHM", "IGHA", "IGHE", "IGHG",
  "IGLV", "IGLC", "IGLJ", "IGKV", "IGKC", "IGKJ",
  "TRAV", "TRBV", "TRDV", "TRGV",
  "TRAJ", "TRBJ", "TRGJ", "TRDJ",
  "TRAC", "TRBC", "TRGC", "TRDC"
)

filter_pattern <- paste0("^(", paste(remove_prefixes, collapse = "|"), ")")

msig_gene_sets <- setNames(
  lapply(all_gene_set_names, function(gs) {
    genes <- msig_df[msig_df$gs_name == gs, "gene_symbol"]
    genes[!grepl(filter_pattern, genes, perl = TRUE)]
  }),
  all_gene_set_names
)

# ============================================================
# 4. Combine all gene sets into single named list
# ============================================================

all_gene_sets           <- c(death_programs, msig_gene_sets)
all_gene_set_names_cimic <- all_gene_set_names
all_gene_set_names_death <- names(death_programs)
all_unique_gene_ids      <- unique(unlist(all_gene_sets))

# ============================================================
# 5. Long and wide formats for downstream use
# ============================================================

gene_set_long <- enframe(
  all_gene_sets,
  name  = "gene_set_ID",
  value = "genes"
) %>%
  unnest_longer(col = "genes") %>%    # quote it as a string
  dplyr::rename(gene = genes) %>%
  dplyr::group_by(gene_set_ID) %>%
  dplyr::mutate(pos = row_number()) %>%
  dplyr::ungroup()

gene_set_wide <- gene_set_long %>%
  pivot_wider(
    names_from   = pos,
    values_from  = gene,
    names_prefix = "V"
  ) %>%
  dplyr::select(gene_set_ID, starts_with("V"))

cat("Gene sets loaded:",    length(all_gene_sets),    "\n")
cat("Unique genes loaded:", length(all_unique_gene_ids), "\n")