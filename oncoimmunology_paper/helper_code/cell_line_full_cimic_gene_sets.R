# =============================================================================
# cell_line_full_cimic_gene_sets.R
# -----------------------------------------------------------------------------
# GLOBAL gene-set definition for the TNBC cell-line CIMIC analysis.
# Run ONCE (source it) before running CIMIC for any drug; it is shared by all
# five agents, so it lives in the analysis root rather than in a drug folder.
#
# It is adapted from:
#   * human_full_cimic_gene_sets.R           (structure, filtering, death programs)
#   * TNBC_CL_Data_Analysis_Anthracyclines.Rmd lines 638-659
#                                             (the reduced cell-line CIM gene-set
#                                              list actually used for clustering)
#
# What it builds (left in the calling environment AND saved to an .rds):
#   cell_line_cimic_gene_sets  named list, MSigDB CIM sets -> gene-symbol vectors
#                              (immunoglobulin / TCR gene families removed).
#                              ** This is the object passed to CIMIC as
#                                 `all_gene_sets` in STEP 4. **
#   death_programs             named list of cell-death program gene sets
#                              (hard-coded, version stable).
#   all_gene_sets              death_programs + cell_line_cimic_gene_sets.
#   cell_line_gene_set_long / _wide   tidy + wide views for downstream use.
#   all_unique_gene_ids        unique gene symbols across cell_line_cimic_gene_sets.
#   all_gene_set_names_cimic / all_gene_set_names_death   the set-name vectors.
#
# Output file:
#   cell_line_cimic_gene_sets.rds   (a single list bundling all of the above)
# =============================================================================

suppressMessages({
  library(msigdbr)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# Analysis root (so the .rds lands next to the master data).
base_dir <- "D:/OneDrive - University of Florida/Gbadamosi Lab/Mohammed Gbadamosi/9_Projects/cimic_tnbc_cl_analysis"
if (nzchar(Sys.getenv("CIMIC_BASE_DIR"))) base_dir <- Sys.getenv("CIMIC_BASE_DIR")

message("== Building cell-line CIMIC gene sets ==")

# ============================================================
# 1. MSigDB reference (Homo sapiens)
# ============================================================
msig_df <- msigdbr::msigdbr(species = "Homo sapiens") %>% as.data.frame()
stopifnot("Expected gs_name + gene_symbol in msigdbr output" =
            all(c("gs_name", "gene_symbol") %in% colnames(msig_df)))

# ============================================================
# 2. Cell-line CIM gene-set names
#    (the active, uncommented list from Rmd lines 638-659 — a focused subset of
#     the full manuscript panel, tuned for the cell-line deltas)
# ============================================================
cell_line_cimic_gene_set_names <- c(
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "GOBP_CELLULAR_RESPONSE_TO_STRESS",
  "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  "GOBP_INFLAMMASOME_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS",
  "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
  "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS",
  "REACTOME_CELLULAR_RESPONSE_TO_MITOCHONDRIAL_STRESS",
  "WP_MOLECULAR_PATHWAY_FOR_OXIDATIVE_STRESS",
  "WP_OXIDATIVE_STRESS_RESPONSE"
)

# Guard: every requested set must exist in this MSigDB build.
missing_sets <- setdiff(cell_line_cimic_gene_set_names, unique(msig_df$gs_name))
if (length(missing_sets))
  stop("Gene-set names not found in msigdbr: ", paste(missing_sets, collapse = ", "))

# Build a named list: each element = the gene-symbol vector for that set.
cell_line_cimic_gene_sets <- setNames(
  lapply(cell_line_cimic_gene_set_names, function(gs) {
    unique(msig_df[msig_df$gs_name == gs, "gene_symbol"])
  }),
  cell_line_cimic_gene_set_names
)

# ============================================================
# 3. Remove immunoglobulin / T-cell-receptor variable-gene families
#    (these are hypervariable and not informative for CIM clustering)
# ============================================================
remove_prefixes <- c(
  "IGHV", "IGHD", "IGHJ", "IGHC", "IGHM", "IGHA", "IGHE", "IGHG",
  "IGLV", "IGLC", "IGLJ", "IGKV", "IGKC", "IGKJ",
  "TRAV", "TRBV", "TRDV", "TRGV",
  "TRAJ", "TRBJ", "TRGJ", "TRDJ",
  "TRAC", "TRBC", "TRGC", "TRDC"
)
filter_pattern <- paste0("^(", paste(remove_prefixes, collapse = "|"), ")")

cell_line_cimic_gene_sets <- lapply(cell_line_cimic_gene_sets, function(gene_vec) {
  gene_vec[!grepl(filter_pattern, gene_vec, perl = TRUE)]
})

# ============================================================
# 4. Cell-death programs (hard-coded; version stable across MSigDB releases)
# ============================================================
death_programs <- list(
  Apoptosis = c(
    "BAX", "BAK1", "BBC3", "PMAIP1", "BCL2L11",
    "APAF1", "CASP9", "CASP3", "CASP7",
    "FAS", "TNFRSF10B", "CASP8", "FADD", "BID",
    "DIABLO", "CYCS", "TP53AIP1", "BAD", "BMF", "HRK"
  ),
  Necroptosis = c(
    "RIPK1", "RIPK3", "MLKL", "ZBP1", "TICAM1", "TNFRSF1A", "FADD", "CASP8"
  ),
  Pyroptosis = c(
    "GSDME", "GSDMD", "CASP1", "CASP4", "CASP5",
    "NLRP3", "AIM2", "PYCARD", "IL1B", "IL18", "CASP3", "NLRC4"
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
# 5. Combine + helper indices
# ============================================================
all_gene_sets            <- c(death_programs, cell_line_cimic_gene_sets)
all_gene_set_names_cimic <- cell_line_cimic_gene_set_names
all_gene_set_names_death <- names(death_programs)
all_unique_gene_ids      <- unique(unlist(cell_line_cimic_gene_sets))

# ============================================================
# 6. Long + wide views (one row per set/gene pair; then pivoted)
# ============================================================
cell_line_gene_set_long <- enframe(
  cell_line_cimic_gene_sets,
  name  = "gene_set_ID",
  value = "genes"
) %>%
  tidyr::unnest_longer(col = genes) %>%
  dplyr::rename(gene = genes) %>%
  dplyr::group_by(gene_set_ID) %>%
  dplyr::mutate(pos = dplyr::row_number()) %>%
  dplyr::ungroup()

cell_line_gene_set_wide <- cell_line_gene_set_long %>%
  pivot_wider(names_from = pos, values_from = gene, names_prefix = "V") %>%
  dplyr::select(gene_set_ID, starts_with("V"))

# ============================================================
# 7. Final sanity checks
# ============================================================
empty_sets <- names(cell_line_cimic_gene_sets)[
  vapply(cell_line_cimic_gene_sets, length, integer(1)) == 0
]
if (length(empty_sets))
  warning("Gene sets that are empty after filtering: ", paste(empty_sets, collapse = ", "))

cat("Cell-line CIMIC gene sets :", length(cell_line_cimic_gene_sets), "\n")
cat("Death programs            :", length(death_programs), "\n")
cat("Unique CIM gene symbols   :", length(all_unique_gene_ids), "\n")

# ============================================================
# 8. Persist everything in one bundle for fast reuse by STEP 4
# ============================================================
gene_set_bundle <- list(
  cell_line_cimic_gene_sets      = cell_line_cimic_gene_sets,
  death_programs                 = death_programs,
  all_gene_sets                  = all_gene_sets,
  all_gene_set_names_cimic       = all_gene_set_names_cimic,
  all_gene_set_names_death       = all_gene_set_names_death,
  all_unique_gene_ids            = all_unique_gene_ids,
  cell_line_gene_set_long        = cell_line_gene_set_long,
  cell_line_gene_set_wide        = cell_line_gene_set_wide
)

saveRDS(gene_set_bundle, file.path(base_dir, "cell_line_cimic_gene_sets.rds"))
message("Saved gene-set bundle -> cell_line_cimic_gene_sets.rds")
message("== Gene-set build complete ==")
