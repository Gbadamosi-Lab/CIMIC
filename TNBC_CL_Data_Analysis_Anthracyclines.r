---
title: "TNBC_Cell_Line_Data_Analysis"
output: html_document
date: "2026-02-10"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```





# Read in Human Cell Line Data and Prepare for data analysis

```{r}
# -------------------------------------------------
# 1. Load required packages (once at the top of the file)
# -------------------------------------------------
library(data.table) # fast fread()
library(dplyr) # pipe, mutate, select
library(tidyr) # row_to_names()
library(janitor) # clean column names, row_to_names()
library(stringr) # string helpers
library(tibble) # rownames_to_column()


# Set WD
setwd(
  "D:/OneDrive - University of Florida/Gbadamosi Lab/Mohammed Gbadamosi/9_Projects/CIM_TNBC_CL_Analysis"
)


# Define the filename
file_path <- "tnbc_processed_data.rds"

if (file.exists(file_path)) {
  
  data_list <- readRDS(file_path)
  tnbc_raw <- data_list$tnbc_raw
  tnbc_lognorm <- data_list$tnbc_lognorm
  tnbc_delta_final <- data_list$tnbc_delta_final

} else {
  message("RDS file not found. Processing data...")

  # -------------------------------------------------
  # 2. Read the raw TPM file and make gene IDs unique
  # -------------------------------------------------
  tnbc_raw <- fread(
    "D:/OneDrive - University of Florida/Gbadamosi Lab/Mohammed Gbadamosi/9_Projects/CIM_TNBC_CL_Analysis/TNBC_Cell_Line_Data.csv"
  ) %>%
    as.data.frame() %>% # convert to plain data.frame for downstream dplyr verbs
    mutate(gene_id = make.unique(gene_id))

  # -------------------------------------------------
  # 3. Log‑normalise all expression columns (keep gene_id untouched)
  # -------------------------------------------------
  tnbc_lognorm <- tnbc_raw %>%
    mutate(across(
      .cols = -gene_id, # all columns except gene_id
      .fns = ~ log2(.x + 1) # log2(TPM+1)
    ))


  # -------------------------------------------------
  # 4. Function: compute delta (and optional fold‑change) vs. correct control
  # -------------------------------------------------
  calc_deltas <- function(df, lognorm = FALSE) {
    # -------------------------------------------------------------------------
    #   df        : data.frame with a 'gene_id' column and sample columns
    #   lognorm   : if TRUE, apply log2(x+1) to the numeric columns first
    # -------------------------------------------------------------------------
    df <- df %>% mutate(gene_id = as.character(gene_id))

    # (optional) log‑normalise inside the function
    if (lognorm) {
      df <- df %>% mutate(across(-gene_id, ~ log2(.x + 1)))
    }

    # -------------------------------------------------
    #   Parse sample meta‑information from column names
    #   Expected pattern:  <cellLine>_<treatment>_<replicate>
    # -------------------------------------------------
    col_info <- tibble(colname = names(df)[-1]) %>% # drop gene_id
      mutate(
        cell_line = str_extract(colname, "^[^_]+"),
        treatment = str_extract(colname, "(?<=_)[^_]+(?=_)")
      )

    # -------------------------------------------------
    #   Loop over each cell line and compute deltas
    # -------------------------------------------------
    for (cell in unique(col_info$cell_line)) {
      # control columns for this cell line
      DMSO_cols <- col_info %>%
        filter(cell_line == cell, treatment == "DMSO") %>%
        pull(colname)
      dmso_cols <- col_info %>%
        filter(cell_line == cell, treatment == "DMSO") %>%
        pull(colname)

      # drug columns (exclude controls)
      drug_cols <- col_info %>%
        filter(cell_line == cell, !treatment %in% c("DMSO", "DMSO")) %>%
        pull(colname)
      if (length(drug_cols) == 0) {
        next
      }

      # pre‑compute control means (row‑wise)
      DMSO_avg <- if (length(DMSO_cols) > 0) {
        rowMeans(df[, DMSO_cols, drop = FALSE], na.rm = TRUE)
      } else {
        NULL
      }
      dmso_avg <- if (length(dmso_cols) > 0) {
        rowMeans(df[, dmso_cols, drop = FALSE], na.rm = TRUE)
      } else {
        NULL
      }

      # -----------------------------------------------------------------
      #   For each drug column, decide which control to use and store delta/fc
      # -----------------------------------------------------------------
      for (drug_col in drug_cols) {
        drug_name <- str_extract(drug_col, "(?<=_)[^_]+(?=_)")
        control_avg <- if (drug_name %in% c("EPI", "DOX")) {
          DMSO_avg
        } else {
          dmso_avg
        }
        if (is.null(control_avg)) {
          next
        } # safety – skip if no matching control

        delta_name <- paste0("delta_", drug_col)
        fc_name <- paste0("fc_", drug_col)

        df[[delta_name]] <- df[[drug_col]] - control_avg
        df[[fc_name]] <- ifelse(
          control_avg != 0,
          df[[drug_col]] / control_avg,
          NA_real_
        )
      }
    }

    df
  }

  # -------------------------------------------------
  # 5. Apply the delta calculation (no extra intermediate objects)
  # -------------------------------------------------
  tnbc_delta <- calc_deltas(tnbc_lognorm, lognorm = FALSE)

  # -------------------------------------------------
  # 6. Keep only the delta columns (plus gene_id) and reshape to long format
  # -------------------------------------------------
  delta_wide <- tnbc_delta %>%
    select(gene_id, starts_with("delta_")) # keep only deltas

  # Transpose so each sample becomes a row (matches original script’s intent)
  delta_long <- delta_wide %>%
    column_to_rownames("gene_id") %>% # gene_id → row names
    t() %>% # transpose
    as.data.frame() %>%
    rownames_to_column("sample_id")

  # -------------------------------------------------
  # 7. Clean column names (make them unique) and ensure numeric type
  # -------------------------------------------------
  colnames(delta_long) <- make.unique(colnames(delta_long))

  delta_long[-1] <- lapply(delta_long[-1], as.numeric) # convert all expression columns to numeric

  # -------------------------------------------------
  # 8. Final tidy data frame (sample_id column already present)
  # -------------------------------------------------
  tnbc_delta_final <- delta_long

  # -------------------------------------------------
  # 9. (Optional) Save the tidy object for later use
  # -------------------------------------------------
  # Option A: Save as a single list (cleanest for RDS)
  data_list <- list(
    tnbc_raw = tnbc_raw,
    tnbc_lognorm = tnbc_lognorm,
    tnbc_delta_final = tnbc_delta_final
  )
  saveRDS(data_list, "tnbc_processed_data.rds")
}

```


# Subset on drug and cells on interest


```{r}
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)


# Define the vectors you want to keep
keep_cells <- c(
  "BT549",
  "DU4475",
  "HCC1395",
  "HCC1806",
  "HCC38",
  "Hs578t",
  "MDAMB231",
  "MDAMB468",
  "MDAMB157"
)
keep_drugs <- c("EPI")


# -------------------------------------------------
# 1. Core filter function
# -------------------------------------------------
filter_delta_df <- function(
  delta_df,
  keep_cells,
  keep_drugs,
  return_sanity = FALSE
) {
  # -----------------------------------------------------------------
  # Subset the Δ‑matrix by cell line AND drug
  # -----------------------------------------------------------------
  subset_df <- delta_df %>%
    filter(
      str_detect(sample_id, str_c(keep_cells, collapse = "|")),
      str_detect(sample_id, str_c(keep_drugs, collapse = "|"))
    )

  # -----------------------------------------------------------------
  # Optional sanity‑check tables (returns a list element)
  # -----------------------------------------------------------------
  if (return_sanity) {
    cell_tbl <- table(
      str_extract(
        subset_df$sample_id,
        paste0(
          "(",
          paste(keep_cells, collapse = "|"),
          ")"
        )
      )
    )
    drug_tbl <- table(
      str_extract(
        subset_df$sample_id,
        paste0(
          "(",
          paste(keep_drugs, collapse = "|"),
          ")"
        )
      )
    )
    return(list(
      df = subset_df,
      cell_counts = cell_tbl,
      drug_counts = drug_tbl
    ))
  } else {
    return(subset_df)
  }
}

# -------------------------------------------------
# 2. Convenience wrappers (optional)
# -------------------------------------------------
# a) Return only the subset (default behaviour)
get_delta_subset <- function(delta_df, keep_cells, keep_drugs) {
  filter_delta_df(delta_df, keep_cells, keep_drugs, FALSE)
}

# b) Return subset + sanity tables
get_delta_subset_with_checks <- function(delta_df, keep_cells, keep_drugs) {
  filter_delta_df(delta_df, keep_cells, keep_drugs, TRUE)
}
library(dplyr)
library(stringr)




## 2️⃣  Filtered matrix plus quick sanity tables
delta_subset_res <- get_delta_subset_with_checks(
  tnbc_delta_final,
  keep_cells,
  keep_drugs
)

# filtered data frame
delta_subset_wide <- delta_subset_res$df %>% as.data.frame()

# sanity‑check tables (same output you had in the original script)
print(delta_subset_res$cell_counts) # counts per cell line
print(delta_subset_res$drug_counts) # counts per drug


library(dplyr)
library(stringr)

filter_with_matched_controls <- function(
  lognorm_df,
  keep_cells,
  keep_drugs,
  return_sanity = FALSE
) {
  # ------------------------------------------
  # Define drug → vehicle mapping
  # ------------------------------------------
  saline_drugs <- c("CARBO", "CIS")
  dmso_drugs <- c("DOX", "EPI", "DOC", "PAC")

  # Determine which vehicles are needed
  vehicles_needed <- c()

  if (any(keep_drugs %in% saline_drugs)) {
    vehicles_needed <- c(vehicles_needed, "SALINE")
  }
  if (any(keep_drugs %in% dmso_drugs)) {
    vehicles_needed <- c(vehicles_needed, "DMSO")
  }

  # ------------------------------------------
  # Base filter: cell lines
  # ------------------------------------------
  base_filtered <- lognorm_df %>%
    dplyr::select(
      gene_id,
      matches(str_c(keep_cells, collapse = "|"))
    )

  # ------------------------------------------
  # Post-treatment subset
  # ------------------------------------------
  post_treat_subset <- base_filtered %>%
    dplyr::select(
      gene_id,
      matches(str_c(keep_drugs, collapse = "|"))
    )

  # ------------------------------------------
  # Pre-treatment subset (matched vehicle)
  # ------------------------------------------
  pre_treat_subset <- base_filtered %>%
    dplyr::select(
      gene_id,
      matches(str_c(vehicles_needed, collapse = "|"))
    )

  # ------------------------------------------
  # Optional sanity tables
  # ------------------------------------------
  if (return_sanity) {
    post_cols <- colnames(post_treat_subset)[-1]
    pre_cols <- colnames(pre_treat_subset)[-1]

    return(list(
      post_treat_subset = post_treat_subset,
      pre_treat_subset = pre_treat_subset,
      post_n = length(post_cols),
      pre_n = length(pre_cols),
      vehicles_used = vehicles_needed
    ))
  } else {
    return(list(
      post_treat_subset = post_treat_subset,
      pre_treat_subset = pre_treat_subset
    ))
  }
}


subset_res <- filter_with_matched_controls(
  lognorm_df = tnbc_lognorm,
  keep_cells = keep_cells,
  keep_drugs = keep_drugs,
  return_sanity = TRUE
)

post_treat_subset_long <- subset_res$post_treat_subset
pre_treat_subset_long <- subset_res$pre_treat_subset


# Add in key metrics 


mhc_genes <- c("HLA-A", "HLA-B", "HLA-C")

anti_tumor <- c("TNFSF13", "CXCL10", "CXCL9", "CXCL11",
                "CCL5", "CX3CL1", "LTA")

pro_tumor <- c("CCL28", "CCL2", "CSF1", "CXCL8",
               "CCL20", "CCL22", "IL1RN", "IL6")


# ---- Keep only genes that exist in this object ----
mhc_present  <- intersect(mhc_genes, tnbc_raw$gene_id)
anti_present <- intersect(anti_tumor, tnbc_raw$gene_id)
pro_present  <- intersect(pro_tumor, tnbc_raw$gene_id)


# Gene sets (already defined)
mhc_genes   <- c("HLA-A", "HLA-B", "HLA-C")
anti_tumor <- c("TNFSF13", "CXCL10", "CXCL9", "CXCL11",
                "CCL5", "CX3CL1", "LTA")
pro_tumor  <- c("CCL28", "CCL2", "CSF1", "CXCL8",
                "CCL20", "CCL22", "IL1RN", "IL6")

# Helper to create a summary row -------------------------------------------------
make_avg_row <- function(df, genes, new_id) {
  df %>%
    filter(gene_id %in% genes) %>%          # keep only the genes of interest
    dplyr::select(-gene_id) %>%                    # drop the identifier column
    summarise(across(everything(),
                     ~ mean(.x, na.rm = TRUE))) %>%  # column‑wise mean (one value per sample)
    mutate(gene_id = new_id) %>%            # add the new identifier
    relocate(gene_id)                       # put gene_id back in first column
}

# ------------------------------------------------------------------------------

# Compute the three summary rows
avg_mhc1_row   <- make_avg_row(pre_treat_subset_long, mhc_genes,   "avg_mhc1")
avg_anti_row   <- make_avg_row(pre_treat_subset_long, anti_tumor, "avg_antitumor_cim_mediator_score")
avg_pro_row    <- make_avg_row(pre_treat_subset_long, pro_tumor,  "avg_protumor_cim_mediator_score")

# Append them to the original long table
pre_treat_subset_long <- bind_rows(pre_treat_subset_long,
                                   avg_mhc1_row,
                                   avg_anti_row,
                                   avg_pro_row)



# Compute the three summary rows
avg_mhc1_row   <- make_avg_row(post_treat_subset_long, mhc_genes,   "avg_mhc1")
avg_anti_row   <- make_avg_row(post_treat_subset_long, anti_tumor, "avg_antitumor_cim_mediator_score")
avg_pro_row    <- make_avg_row(post_treat_subset_long, pro_tumor,  "avg_protumor_cim_mediator_score")

# Append them to the original long table
post_treat_subset_long <- bind_rows(post_treat_subset_long,
                                   avg_mhc1_row,
                                   avg_anti_row,
                                   avg_pro_row)



library(dplyr)


# Safe row mean function
safe_rowmean <- function(df, genes) {
  if (length(genes) == 0) return(rep(NA, nrow(df)))
  rowMeans(dplyr::select(df, all_of(genes)), na.rm = TRUE)
}

# Add avg_MHC
delta_subset_wide <- delta_subset_wide %>%
  mutate(
    avg_mhc1 = safe_rowmean(., mhc_present),
    avg_antitumor_cim_mediator_score = safe_rowmean(., anti_present),
    avg_protumor_cim_mediator_score = safe_rowmean(., pro_present)
  )



pre_treat_subset_wide <- pre_treat_subset_long %>%
  column_to_rownames("gene_id") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")


post_treat_subset_wide <- post_treat_subset_long %>%
  column_to_rownames("gene_id") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")


delta_subset_long <- delta_subset_wide %>% as.data.frame() %>%
column_to_rownames("sample_id") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")


subset_res$vehicles_used

library(dplyr)
library(stringr)
library(tibble)

# --- 1. Mapping table ----------------------------------------------------
cluster_map <- tibble(
  cell_line = c("DU4475","HCC1806","HCC1395","Hs578t",
                "MDAMB157","MDAMB231","MDAMB468","BT549","HCC38"),
  cimic_cluster = c(1,1,1,2,2,2,1,2,2)
)

# --- 2. Helper that works for both pre/post and delta sample IDs ---------
add_cell_cluster <- function(df) {
  df %>%
    mutate(
      cell_line = sample_id,
      cell_line = str_remove(cell_line, "^delta_"),
      cell_line = str_extract(cell_line, "^[^_]+")
    ) %>%
    left_join(cluster_map, by = "cell_line")
}

# --- 3. Apply ------------------------------------------------------------
pre_treat_subset_wide  <- add_cell_cluster(pre_treat_subset_wide)
post_treat_subset_wide <- add_cell_cluster(post_treat_subset_wide)
delta_subset_wide      <- add_cell_cluster(delta_subset_wide)





library(dplyr)
library(stringr)
library(tibble)

# -------------------------------------------------
# 1. Mapping table (your cluster definition)
# -------------------------------------------------
cluster_map <- tibble(
  cell_line = c('DU4475','HCC1806','HCC1395','Hs578t',
                'MDAMB157','MDAMB231','MDAMB468','BT549','HCC38'),
  cimic_cluster   = c(1,1,1,2,2,2,1,2,2)
)

# -------------------------------------------------
# 2. Build a sample‑metadata lookup table
# -------------------------------------------------
make_sample_meta <- function(wide_df) {
  sample_names <- colnames(wide_df)[-1L]          # drop gene_id
  tibble(
    sample_id  = sample_names,
    cell_line  = str_extract(sample_names, "^[^_]+"),
    cimic_cluster    = cluster_map$cimic_cluster[match(
                  str_extract(sample_names, "^[^_]+"),
                  cluster_map$cell_line)]
  )
}

# -------------------------------------------------
# 3. Attach the metadata as an attribute
# -------------------------------------------------
attach_meta <- function(wide_df) {
  attr(wide_df, "sample_meta") <- make_sample_meta(wide_df)
  wide_df
}

# -------------------------------------------------
# 4. Apply to the three data sets
# -------------------------------------------------
pre_treat_subset_long  <- attach_meta(pre_treat_subset_long)
post_treat_subset_long <- attach_meta(post_treat_subset_long)
delta_subset_long      <- attach_meta(delta_subset_long)

# -------------------------------------------------
# 5. Helper to retrieve the metadata
# -------------------------------------------------
get_sample_meta <- function(wide_df) {
  attr(wide_df, "sample_meta")
}

# -------------------------------------------------
# Example: view the metadata for the pre‑treat set
# -------------------------------------------------
get_sample_meta(pre_treat_subset_long)


combined_TPM_lognorm <- cbind (pre_treat_subset_long, post_treat_subset_long)

```

# If need to re-run CIMIC and other analysis 

```{r}

# get clustering matrix 
clustering_matrix <- delta_subset_wide[, 1:78933] %>% 
  column_to_rownames("sample_id") %>% 
  as.matrix()



# check batch effect of deltas 

library(PCAtools)

delta_pca_mat <- as.data.frame(t(clustering_matrix))
names_vector <- colnames(delta_pca_mat)

# Identify which start with OB or NEO
cell_line_names <- str_match(names_vector, "^delta_([^_]+)")[,2]


# set up batch data 
delta_batch_data <- data.frame(
  Batch = cell_line_names,
  row.names = names_vector,
  stringsAsFactors = FALSE
)


colnames(delta_batch_data)[1] <- 'Batch'
pca_deltas <- PCAtools::pca(delta_pca_mat[apply(delta_pca_mat, 1, function(x) all(is.finite(x)) & sd(x) != 0), ], 
                  metadata = delta_batch_data, scale = TRUE,  center = TRUE)

biplot(pca_deltas,
    lab = NULL,
    colby = 'Batch',
    #shape = 'condition',
    hline = 0, vline = 0,
    legendPosition = 'right',
    title = 'PCA on deltas',
    subtitle = 'PC1 versus PC2',)

# Get Cell Line Specific Programs ----------------------


msig_df <- msigdbr::msigdbr(species = "Homo sapiens") %>% as.data.frame()

    # Default CIM Genesets from manuscript
    cell_line_cimic_gene_set_names <- c(
      #"GOBP_ADAPTIVE_IMMUNE_RESPONSE",
      "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
      #"GOBP_B_CELL_ACTIVATION",
      "GOBP_CELLULAR_RESPONSE_TO_STRESS",
      "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
      "GOBP_INFLAMMASOME_MEDIATED_SIGNALING_PATHWAY",
      "GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS",
      "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
      "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
      #"GOBP_LEUKOCYTE_CHEMOTAXIS_INVOLVED_IN_INFLAMMATORY_RESPONSE",
      #"GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
      #"GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
      #"GOBP_T_CELL_ACTIVATION",
      "HALLMARK_INFLAMMATORY_RESPONSE",
      "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
      "REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS",
      "REACTOME_CELLULAR_RESPONSE_TO_MITOCHONDRIAL_STRESS",
      "WP_MOLECULAR_PATHWAY_FOR_OXIDATIVE_STRESS",
      "WP_OXIDATIVE_STRESS_RESPONSE"
    )
    
  # Build a named list: each element is the vector of gene symbols for that set
    cell_line_cimic_gene_sets <- setNames(
      lapply(cell_line_cimic_gene_set_names, function(gs) {
        gs_df <- msig_df[msig_df$gs_name == gs, , drop = FALSE]
        gs_df$gene_symbol
      }),
      cell_line_cimic_gene_set_names

    )
  
 remove_prefixes <- c(
  # B‑cell (heavy, light, constant, joining, diversity, …)
  "IGHV", "IGHD", "IGHJ", "IGHC", "IGHM", "IGHA", "IGHE", "IGHG",
  "IGLV", "IGLC", "IGLJ",
  "IGKV", "IGKC", "IGKJ",
  # T‑cell receptor chains
  "TRAV", "TRBV", "TRDV", "TRGV",
  "TRAJ", "TRBJ", "TRGJ", "TRDJ",
  "TRAC", "TRBC", "TRGC", "TRDC"
)

# Build a single regex that matches any of the prefixes at the start
pattern <- paste0("^(", paste(remove_prefixes, collapse = "|"), ")")

# -----------------------------------------------------------------
# 2. Apply the filter to the *named list* `cell_line_cimic_gene_sets`
#    (each element is a character vector of gene symbols)
# -----------------------------------------------------------------


  cell_line_cimic_gene_sets <- lapply(cell_line_cimic_gene_sets, function(gene_vec) {

  # keep everything else
  gene_vec[!grepl(pattern, gene_vec, perl = TRUE)]
})



  # cell_line_cimic_gene_sets is a named list:
#   names(cell_line_cimic_gene_sets) → gene‑set IDs (e.g. "GOBP_ADAPTIVE_IMMUNE_RESPONSE")
#   each element          → character vector of gene symbols
# -------------------------------------------------
# 2) Convert the list → long format (one row per gene‑set / gene pair)
# -------------------------------------------------
# 2) Convert the list → long format (one row per gene‑set / gene pair)
cell_line_gene_set_long <- enframe(
  cell_line_cimic_gene_sets,
  name = "gene_set_ID",
  value = "genes"
) %>%
  tidyr::unnest_longer(col = genes) %>%
  dplyr::rename(gene = genes)

# -------------------------------------------------
# 3) Add a position index inside each gene‑set (1,2,3,…)
# -------------------------------------------------
cell_line_gene_set_long <- cell_line_gene_set_long %>% 
  group_by(gene_set_ID) %>% 
  mutate(pos = row_number()) %>%    # position of the gene inside the set
  ungroup()

# -------------------------------------------------
# 4) Wide‑reshape: one row per gene‑set, columns ID1, ID2, …
# -------------------------------------------------
cell_line_gene_set_wide <- cell_line_gene_set_long %>% 
  pivot_wider(names_from   = pos,
              values_from  = gene,
              names_prefix = "V") %>%      # → ID1, ID2, …
  # optional: keep column order tidy
  dplyr::select(gene_set_ID, starts_with("V"))


all_unique_gene_ids <- unique(unlist(cell_line_cimic_gene_sets))



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

all_gene_sets <- c(death_programs, cell_line_cimic_gene_sets)
all_gene_set_names_cimic <- cell_line_cimic_gene_set_names
all_gene_set_names_death <- names(death_programs)

    
    
# Set up for CIMIC -------------------

# Choose one clustering algorithm
clustering_alg <- "hc"  # Choose either "km" for K-means or "hc" for hierarchical clustering

#Choose clustering metrics
library(gtools)

metrics <- c("pac", "silhouette_dim_reduce_space", "silhouette_ge_space", "silhouette_combined_avg", "item_cluster_consensus")
all_perms <- list()

for (k in 1:length(metrics)) {
  perms <- permutations(n = length(metrics), r = k, v = metrics)
  all_perms[[k]] <- perms
}
# Bind into a single dataframe, pad with NA for shorter ones
max_len <- length(metrics)
df_list <- lapply(all_perms, function(mat) {
  if (ncol(mat) < max_len) {
    pad <- matrix(NA, nrow = nrow(mat), ncol = max_len - ncol(mat))
    cbind(mat, pad)
  } else {
    mat
  }
})
clust_met_options <- as.data.frame(do.call(rbind, df_list), stringsAsFactors = FALSE)
colnames(clust_met_options) <- paste0("metric", 1:max_len)
rownames(clust_met_options) <- NULL

clustering_metrics <- as.character(unlist(clust_met_options[53, ])) # silhouette combine, pac, item-cluster-consensus    
# clustering_metrics <- as.character(unlist(clust_met_options[9, ])) # "silhouette", "pac"

clustering_metrics <- clustering_metrics[!is.na(clustering_metrics)]
clustering_metrics


# Set up pacmap if needed -------------


# Only install if not already available
if (!py_module_available("pacmap")) {
  # Then come back to R to execute this code 
library(reticulate)


# This must be the exact path to your conda.exe
Sys.setenv(RETICULATE_CONDA = "C:/Users/mgbadamosi/AppData/Local/miniconda3/condabin/conda.bat")


# Set up for PacMAP uses
use_condaenv("pacmap_env", required = TRUE)
pacmap <- import("pacmap")
py_config()

}




# Run CIMIC ------------------



  
  



code_dir <- "D:/OneDrive - University of Florida/Gbadamosi Lab/Mohammed Gbadamosi/CIM_TNBC_Scripts"
source(file.path(code_dir, "CIMIC_Release_1.0.0.R"))



CIM_feature_selection_pacmap <- CIM_feature_selection_by_gene_set_pacmap(
  clustering_matrix = clustering_matrix,
  remove_immune_variable_genes = TRUE,
  all_gene_sets = cell_line_cimic_gene_sets,
  clustering_alg = clustering_alg,
  max_k = 5,
  CCP_iter = 1000,
  adj_pval_thresh = 0.05,
  max_pipeline_iter = 50,
  seed = 2024L,
  filter_approach = c("app_one", "app_two", "app_three"), 
  clustering_metrics = clustering_metrics,
  pacmap_dimensions = 2L,
  pacmap_args = NULL,               
  pacmap_guardrails = TRUE, 
  verbose = TRUE,
  working_dir = working_dir
)







# Extract out data set from clustering analyses --------------

#Plot feature importance using anova F stat
library(ggplot2)
library(EnvStats)
library(tibble)
library(grid)
library(gridtext)
library(circlize)

# Change feature selection for clustering approach
# this is the the way the input features were presented

# Must be one of: umap, pacmap, pca, individual_genes

# "umap"              # Features are reduced coordinates from UMAP
# "pacmap"            # Features are reduced coordinates from PaCMAP
# "pca"               # Features are reduced coordinates from PCA
# "individual_genes"  # Features are the original individual genes (no dimensionality reduction)
clustering_approach <- "pacmap"


# this is the way noise was removed must be one of the following below.
# app_one: Denoising within each gene set.
#   For each gene set, cluster samples and iteratively remove genes that are not significantly different between clusters.
#   Repeat for each gene set until no more genes can be removed. Combine resulting informative genes for final clustering.

# app_two: Global denoising.
#   Cluster samples using all genes together, removing genes that are not significantly different between clusters, regardless of gene set.
#   Repeat until stable gene set is reached; use this set for final clustering.

# app_three: Hierarchical denoising.
#   First apply app_one to filter genes within each gene set, then apply app_two on the combined informative genes from all sets to further denoise.
#   Final clustering is performed on the genes that pass both filters.
denoise_approach <- "app_three"

# Construct the list key based on the selected approaches
key <- paste0("final_df_", denoise_approach)


# Extract out based on approach used 
# Slick all-in-one extraction of the data frame
list_name <- paste0("CIM_feature_selection_", clustering_approach)


if (exists(list_name, envir = .GlobalEnv)) {
  results_list <- get(list_name, envir = .GlobalEnv)
  if (key %in% names(results_list)) {
    clustered_plot_df <- results_list[[key]]
  } else {
    stop(sprintf("Data frame '%s' not found in list '%s'.", key, list_name))
  }
} else {
  stop(sprintf("List '%s' does not exist.", list_name))
}

# Map denoise_approach to the corresponding gene set variable
gene_set_key <- switch(denoise_approach,
  "app_one"   = "iterated_by_gene_sets",
  "app_two"   = "iterated_over_all_genes",
  "app_three" = "iterated_use_final_gene_set_from_approach_one",
  stop("Unknown denoise_approach: ", denoise_approach)
)

# Now extract from your results list (e.g., results_list) 
if (gene_set_key %in% names(results_list)) {
  clustering_gene_set <- results_list[[gene_set_key]]
} else {
  stop(paste("Gene set", gene_set_key, "not found in results list."))
}


# read in if already ran ------------------ 



clustered_plot_df <- fread("D:/OneDrive - University of Florida/Gbadamosi Lab/Mohammed Gbadamosi/9_Projects/CIM_TNBC_CL_Analysis/EPI/CIM_states_results_hc_pacmap_app_three/clustered_samples_app_three.csv") %>% as.data.frame


# Create master induction df ------------------------

library(circlize)
library(ComplexHeatmap)
library(dunn.test)
library(colorRamp2)

clustered_plot_df_genes <- colnames(clustered_plot_df)[-c(1:10)]

# extract out optimal_k 
optimal_k <- length(unique(clustered_plot_df$cluster_assignments))

# Prepare a list to collect rows
# Predefine output matrix for efficiency
n_genes <- length(clustered_plot_df_genes)
optimal_k <- length(unique(clustered_plot_df$cluster_assignments))

# Setup column names: Mean, SD, SEM per cluster, then pairwise Dunn p-values
n_pairs <- optimal_k * (optimal_k - 1) / 2

if (optimal_k > 2){
pairs_mat <- combn(seq_len(optimal_k), 2)
pairwise_names <- apply(pairs_mat, 2, function(x) {
  sprintf("cluster%d_cluster%d_pval", x[1], x[2])
})
} else if (optimal_k == 2){
  
 pairwise_names <-  "cluster1_cluster2_pval"
}


cluster_names <- as.vector(rbind(
  paste0("cluster", seq_len(optimal_k), "_mean"),
  paste0("cluster", seq_len(optimal_k), "_sd"),
  paste0("cluster", seq_len(optimal_k), "_sem")
))

penetrance_names <- paste0("cluster", seq_len(optimal_k), "_penetrance")

out_mat <- matrix(
  NA,
  nrow = n_genes,
  ncol = 3 + length(cluster_names) + length(penetrance_names) + n_pairs
)

colnames(out_mat) <- c(
  "gene_id", "H_value", "p_value",
  cluster_names,
  penetrance_names,
  pairwise_names
)

# Loop
for (i in seq_len(n_genes)) {
  
  
  gene_name <- clustered_plot_df_genes[i]
  
  if(gene_name %in% colnames(clustered_plot_df)){
  gene_vals <- clustered_plot_df[[gene_name]]
  clusters  <- as.factor(clustered_plot_df$cluster_assignments)
  clusters  <- droplevels(clusters)  # just in case there are unused levels

  # Per-cluster stats
  cluster_means <- tapply(gene_vals, clusters, mean)
  cluster_sds   <- tapply(gene_vals, clusters, stats::sd)
  cluster_ns    <- tapply(gene_vals, clusters, function(x) sum(!is.na(x)))
  cluster_sems  <- cluster_sds / sqrt(cluster_ns)

  # Fill in full-length per-cluster vectors (in order 1:optimal_k)
  cluster_means_full <- rep(NA, optimal_k)
  cluster_sds_full   <- rep(NA, optimal_k)
  cluster_sems_full  <- rep(NA, optimal_k)
  idx <- as.integer(names(cluster_means))
  cluster_means_full[idx] <- as.numeric(cluster_means)
  cluster_sds_full[idx]   <- as.numeric(cluster_sds)
  cluster_sems_full[idx]  <- as.numeric(cluster_sems)
  
  # TRUE penetrance (% samples > 0) 
  cluster_penetrance <- tapply(gene_vals, clusters, function(x)
    mean(x > 0, na.rm = TRUE))
  
  # Fill full-length penetrance vector
  cluster_penetrance_full <- rep(NA, optimal_k)
  idx <- as.integer(names(cluster_penetrance))
  cluster_penetrance_full[idx] <- as.numeric(cluster_penetrance)
  
  # Interleave mean, sd, sem for each cluster
  combined_stats <- as.vector(rbind(cluster_means_full, cluster_sds_full, cluster_sems_full))

  # Initialize post-hoc p-value vector
  dunn_pvals <- rep(NA, n_pairs)

  # Statistical test selection
  if (nlevels(clusters) == 2) {
    # Two groups → Wilcoxon rank-sum
    test_res <- wilcox.test(gene_vals ~ clusters)
    kw_H_value <- as.numeric(test_res$statistic)
    kw_p_value <- as.numeric(test_res$p.value)
    # No Dunn’s needed for two groups
  } else if (nlevels(clusters) > 2) {
    # Kruskal–Wallis
    df_kw <- kruskal.test(gene_vals ~ clusters)
    kw_H_value <- as.numeric(df_kw$statistic)
    kw_p_value <- as.numeric(df_kw$p.value)

    # Dunn’s post-hoc
    dt <- dunn.test(gene_vals, clusters, kw = FALSE, table = FALSE, method = "bh")
    if (!is.null(dt$P.adjusted)) {
      dunn_pvals[1:length(dt$P.adjusted)] <- dt$P.adjusted
    }
  } else {
    # Less than two groups — no test possible
    kw_H_value <- NA
    kw_p_value <- NA
  }

  # Output row
  out_mat[i, ] <- c(
    gene_name,
    kw_H_value,
    kw_p_value,
    combined_stats,
    cluster_penetrance_full,
    dunn_pvals
  )
}
}

# colnames(out_mat) <- c("gene_id", "H_value", "p_value", 
#                        cluster_names, 
#                        pairwise_names)


master_induction_df <- as.data.frame(out_mat, stringsAsFactors = FALSE)

#if to remove post-hoc test if two group
 if (nlevels(clustered_plot_df$cluster_assignments) == 2) {
   
   master_induction_df <- master_induction_df %>% dplyr::select(-all_of(pairwise_names))
}

# Optionally convert numeric columns back to numeric:
num_cols <- setdiff(names(master_induction_df), "gene_id")
master_induction_df[num_cols] <- lapply(master_induction_df[num_cols], as.numeric)
master_induction_df$stat_test_adj_p <- p.adjust(master_induction_df$p_value, method = "BH")

master_induction_df <- master_induction_df[, c("gene_id", "H_value", "p_value", "stat_test_adj_p", setdiff(names(master_induction_df), c("gene_id", "H_value", "p_value", "stat_test_adj_p")))]

# Suppose optimal_k and master_induction_df are already defined
for (i in 1:(optimal_k - 1)) {
  for (j in (i + 1):optimal_k) {
    col_i_mean <- paste0("cluster", i, "_mean")
    col_j_mean <- paste0("cluster", j, "_mean")
    col_i_sd <- paste0("cluster", i, "_sd")
    col_j_sd <- paste0("cluster", j, "_sd")
    diff_col1 <- paste0("cluster", i, "_cluster", j, "_differential")
    diff_col2 <- paste0("cluster", j, "_cluster", i, "_differential")
    sd_diff_col1 <- paste0("cluster", i, "_cluster", j, "_sd")
    sd_diff_col2 <- paste0("cluster", j, "_cluster", i, "_sd")
    
    master_induction_df[[diff_col1]] <- as.numeric(master_induction_df[[col_i_mean]]) - as.numeric(master_induction_df[[col_j_mean]])
    master_induction_df[[diff_col2]] <- as.numeric(master_induction_df[[col_j_mean]]) - as.numeric(master_induction_df[[col_i_mean]])
    master_induction_df[[sd_diff_col1]] <- sqrt((as.numeric(master_induction_df[[col_i_sd]]))^2 + (as.numeric(master_induction_df[[col_j_sd]]))^2)
    master_induction_df[[sd_diff_col2]] <- master_induction_df[[sd_diff_col1]]
  }
}


# Add cimic clusters #########

clustered_plot_df <- clustered_plot_df %>% 
   dplyr::mutate(
    CIMIC_Cluster = dplyr::case_when(
      cluster_assignments == 1 ~ "Fun-CIM",
      cluster_assignments == 2 ~ "Dys-CIM",
      TRUE ~ NA_character_
    )
  )


# add avgMHC col

mhc_cols <- c("HLA-A", "HLA-B", "HLA-C")

clustered_plot_df <- clustered_plot_df %>%
  dplyr::mutate(
    Avg_MHC1 = rowMeans(dplyr::across(all_of(mhc_cols)), na.rm = TRUE)
  )
# recreate pacmap plot ---------------

# Define color palette for clusters
cols <- hue_pal()(length(unique(clustered_plot_df$cluster_assignments)))

clustered_plot_df$base_id <-  str_match(clustered_plot_df$sample_id, "^delta_([^_]+)")[,2]

# Plot using GGPlot Color by cluster
final_cluster_plot_pacmap <- ggplot(
  clustered_plot_df,
  aes(
    x = PACMAP1,
    y = PACMAP2,
    label = base_id,  # Your patient/sample label column
    color = as.factor(cluster_assignments)
  )
) +
  scale_color_manual(values = cols, name = "Cluster") +
  scale_fill_manual(values = cols, name = "Cluster") +
  ggforce::geom_mark_hull(
    aes(group = as.factor(cluster_assignments),
        fill = as.factor(cluster_assignments), label = NULL),
    concavity = 2, expand = unit(2, "mm"), alpha = 0.2
  ) +
  geom_label_repel(size = 4, fontface = "bold", label.size = 0.3, fill = "white", 
                   alpha = 1, box.padding = 0.25, max.overlaps = 100) +
  labs(x = "PacMAP 1", y = "PacMAP 2",
       fill = "Cluster",
       color = "Cluster"
       ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    axis.text = element_text(colour = "black", size = 24),
    axis.text.x = element_text(colour = "black", size = 24),
    axis.text.y = element_text(colour = "black", size = 24),
    axis.title = element_text(colour = "black", size = 30),
    axis.ticks = element_line(size = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )


final_cluster_plot_pacmap


# Volc Plot -------------------





library(ggplot2)
library(dplyr)
library(ggrepel)

# Define the two clusters of interest and other params
clusters_oi <- c("cluster1", "cluster2")
slope_cutoff <- log2(1.5)
pval_cutoff  <- 0.05


ind_cutoff <- log2(1.50)


penetrance_target_cutoff <- 0.20   # ≥20% induced in target cluster
penetrance_other_cutoff  <- 0.20   # ≤20% induced in other cluster

# Build the column names dynamically
diff_col      <- paste0(clusters_oi[1], "_", clusters_oi[2], "_differential")
diff_col_rev  <- paste0(clusters_oi[2], "_", clusters_oi[1], "_differential")
pval_col      <- "stat_test_adj_p"
mean_col	  <- paste0(clusters_oi[1], "_mean")
mean_col_rev	<- paste0(clusters_oi[2], "_mean")

pen_target <- paste0(clusters_oi[1], "_penetrance")
pen_other  <- paste0(clusters_oi[2], "_penetrance")

C1_class <- "Fun-CIM"
C2_class <- "Dys-CIM"


ind_plot_category_1 <- paste0("Induced in ",       C1_class, " only")
ind_plot_category_2 <- paste0("Induced in ",       C2_class, " only")
ind_plot_category_3 <- paste0("Highly induced in ", C1_class)
ind_plot_category_4 <- paste0("Highly induced in ", C2_class)
ind_plot_category_5 <- paste0("Highly repressed in ", C1_class)
ind_plot_category_6 <- paste0("Highly repressed in ", C2_class)
ind_plot_category_NS <- "NS"



plot_data <- master_induction_df %>%
    dplyr::select(gene_id,
                  C1 = all_of(mean_col),
                  C2 = all_of(mean_col_rev),
                  pval = all_of(pval_col),
                  diff = all_of(diff_col)) %>%
    dplyr::mutate(
        significant = pval < pval_cutoff & abs(diff) > slope_cutoff,
        category = dplyr::case_when(
            !significant                                      ~ ind_plot_category_NS,
            C1 >  0 & C2 <  0                                ~ ind_plot_category_1,
            C1 <  0 & C2 >  0                                ~ ind_plot_category_2,
            C1 >  0 & C2 >  0 & (C1 - C2) > slope_cutoff    ~ ind_plot_category_3,
            C1 >  0 & C2 >  0 & (C2 - C1) > slope_cutoff    ~ ind_plot_category_4,
            C1 <  0 & C2 <  0 & (C2 - C1) > slope_cutoff    ~ ind_plot_category_5,
            C1 <  0 & C2 <  0 & (C1 - C2) > slope_cutoff    ~ ind_plot_category_6,
            TRUE                                              ~ ind_plot_category_NS
        )
    )


category_counts <- plot_data %>%
    dplyr::count(category) %>%
    dplyr::mutate(label_with_n = paste0(category, " (N = ", n, ")"))

# Build a named vector to recode legend labels
n_labels <- setNames(category_counts$label_with_n, category_counts$category)

n_genes <- 15
available_genes <- unique(unlist(all_gene_sets[all_gene_set_names_cimic]))

top_up <- plot_data %>%
    dplyr::filter(gene_id %in% available_genes, significant,
                  category %in% c(ind_plot_category_1, ind_plot_category_3)) %>%
    dplyr::mutate(dist_from_diag = C1 - C2) %>%
    dplyr::arrange(dplyr::desc(dist_from_diag)) %>%
    dplyr::slice_head(n = n_genes) %>%
    dplyr::pull(gene_id)

top_down <- plot_data %>%
    dplyr::filter(gene_id %in% available_genes, significant,
                  category %in% c(ind_plot_category_2, ind_plot_category_4)) %>%
    dplyr::mutate(dist_from_diag = C2 - C1) %>%
    dplyr::arrange(dplyr::desc(dist_from_diag)) %>%
    dplyr::slice_head(n = n_genes) %>%
    dplyr::pull(gene_id)


c1 <- c("FOS", "FOSB", "EGR1", "ISG15", "ISG20", "MX1", "IFI35", "DHX58", "XAF1", "OAS1", "OAS2", "OAS3", "P2Y11", "ATF4")


c2 <- c("VIM", "ZEB1", "CDK6", "HMGA1", "LIFE", "XRCC4", "EIF2AK4","EIF3E", "EIF3H", "HSPA5", "UVRAG", "ATF6")

top_genes <- c(c1, c2)


# top_genes <- c(top_up, top_down)

# Add label column to plot_data
plot_data <- plot_data %>%
    dplyr::mutate(label = dplyr::if_else(gene_id %in% top_genes, gene_id, NA_character_))


size_n <- 18



ind_graph <- ggplot(plot_data, aes(x = C2, y = C1)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1.5) +
    geom_abline(slope = 1, intercept =  slope_cutoff, linetype = "dashed", color = "steelblue", linewidth = 1.5) +
    geom_abline(slope = 1, intercept = -slope_cutoff, linetype = "dashed", color = "firebrick", , linewidth = 1.5) +
    #geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "grey40", linewidth = 0.5) +
    geom_point(aes(color = category, shape = significant), size = 1.5, alpha = 0.9) +
    ggrepel::geom_label_repel(
        data = subset(plot_data, !is.na(label)),
        aes(label = label),
        size          = 4.5,
        fontface      = "bold",
        color         = "black",
        fill          = "white",
        label.size    = 0.5,
        box.padding   = 0.4,
        point.padding = 0.3,
        segment.color = "black",
        segment.size  = 0.4,
        max.overlaps  = Inf
    ) +
    scale_shape_manual(
        name   = "Significance",
        values = c("TRUE" = 16, "FALSE" = 17),
        labels = c("TRUE" = "p < 0.05", "FALSE" = "NS")
    ) +
    scale_color_manual(
        name = "Category",
        values = setNames(
            c("black", "#E41A1C", "#377EB8", "#FF7F00", "#984EA3", "#4DAF4A", "#A65628"),
            c(ind_plot_category_NS, ind_plot_category_1, ind_plot_category_2,
              ind_plot_category_3,  ind_plot_category_4, ind_plot_category_5,
              ind_plot_category_6)
        ),
        labels = n_labels
    ) +
    labs(
        x = paste0("Dys-CIM mean induction\n(Δlog2[TPM+1])"),
        y = paste0("Fun-CIM mean induction\n(Δlog2[TPM+1])")
    ) +
    ggplot2::coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3)) +
    ggplot2::scale_x_continuous(breaks = seq(-3, 3, 1)) +
    ggplot2::scale_y_continuous(breaks = seq(-3, 3, 1)) +
    theme_classic(base_size = 14) +
    theme(
        axis.title.x = element_text(face = "bold", size = size_n, color = "black"),
        axis.title.y = element_text(face = "bold", size = size_n, color = "black"),
        axis.text.x  = element_text(face = "bold", size = size_n, color = "black"),
        axis.text.y  = element_text(face = "bold", size = size_n, color = "black"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0),
        axis.line    = element_blank()
    )

print(ind_graph)

# 
# ind_graph <- ind_graph +
#   annotate("label",
#     x = 2.5,  y = -2.5, label = ind_plot_category_1,
#     fontface = "bold", size = 3.5, fill = "#E41A1C", color = "white", alpha = 0.8
#   ) +
#   annotate("label",
#     x = -2.5, y = 2.5,  label = ind_plot_category_2,
#     fontface = "bold", size = 3.5, fill = "#377EB8", color = "white", alpha = 0.8
#   ) +
#   annotate("label",
#     x = 2.5,  y = 1.5,  label = ind_plot_category_3,
#     fontface = "bold", size = 3.5, fill = "#FF7F00", color = "white", alpha = 0.8
#   ) +
#   annotate("label",
#     x = -2.5, y = -1.5, label = ind_plot_category_4,
#     fontface = "bold", size = 3.5, fill = "#984EA3", color = "white", alpha = 0.8
#   ) +
#   annotate("label",
#     x = -2.5, y = -2.5, label = ind_plot_category_5,
#     fontface = "bold", size = 3.5, fill = "#4DAF4A", color = "white", alpha = 0.8
#   ) +
#   annotate("label",
#     x = 2.5,  y = 2.5,  label = ind_plot_category_6,
#     fontface = "bold", size = 3.5, fill = "#A65628", color = "white", alpha = 0.8
#   )

induced_C1_only <- plot_data %>%
    dplyr::filter(category == ind_plot_category_1)

induced_C2_only <- plot_data %>%
    dplyr::filter(category == ind_plot_category_2)











# top heat map -----------------------
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)


ind_cutoff <- log2(1.25)
pval_cutoff <- 0.05
penetrance_target_cutoff <- 0.25
cluster_induction_cutoff <- log2(0.1)

cluster1_genes <- master_induction_df %>%
    dplyr::filter(
        !is.na(stat_test_adj_p),
        cluster1_cluster2_differential >= ind_cutoff,
        cluster1_penetrance >= penetrance_target_cutoff,
        cluster1_mean >= cluster_induction_cutoff, 
        stat_test_adj_p <= pval_cutoff
    )

cluster2_genes <- master_induction_df %>%
    dplyr::filter(
        !is.na(stat_test_adj_p),
        cluster2_cluster1_differential >= ind_cutoff,
        cluster2_penetrance >= penetrance_target_cutoff,
        cluster2_mean >= cluster_induction_cutoff,
        stat_test_adj_p <= pval_cutoff
    )

n = 100


# TOP GENES

top25_c1 <- cluster1_genes %>%
    dplyr::arrange(dplyr::desc(cluster1_cluster2_differential)) %>%
    dplyr::slice_head(n = n ) %>%
    dplyr::pull(gene_id)

top25_c2 <- cluster2_genes %>%
    dplyr::arrange(dplyr::desc(cluster2_cluster1_differential)) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::pull(gene_id)


heatmap_genes <- unique(c(top25_c1, top25_c2))

# BUILD MATRIX (samples x genes)

meta_cols <- c("sample_id","UMAP1","UMAP2","PACMAP1","PACMAP2",
               "PC1","PC2","PC3","base_id","cluster_assignments")

gene_cols <- setdiff(colnames(clustered_plot_df), meta_cols)
genes_present <- intersect(heatmap_genes, gene_cols)

mat <- clustered_plot_df %>%
    dplyr::select(sample_id, dplyr::all_of(genes_present))

rownames(mat) <- mat$sample_id
mat$sample_id <- NULL


# ORDER SAMPLES BY CIMIC

ordered_samples <- clustered_plot_df %>%
    dplyr::arrange(cluster_assignments) %>%
    dplyr::pull(sample_id)

ordered_samples <- intersect(ordered_samples, rownames(mat))

mat <- mat[ordered_samples, , drop = FALSE]


# TRANSPOSE → genes x samples

heatmap_mat <- scale(t(mat))


# GENE GROUPING (FOR SPLIT)

gene_group <- dplyr::case_when(
    rownames(heatmap_mat) %in% top25_c1 ~ "Cluster1",
    rownames(heatmap_mat) %in% top25_c2 ~ "Cluster2",
    TRUE ~ "Other"
)

gene_group <- factor(gene_group, levels = c("Cluster1", "Cluster2", "Other"))


# SAMPLE ANNOTATION

annotation_df <- data.frame(
    CIM_cluster = clustered_plot_df$cluster_assignments
)

rownames(annotation_df) <- clustered_plot_df$sample_id
annotation_df <- annotation_df[ordered_samples, , drop = FALSE]

ha <- ComplexHeatmap::HeatmapAnnotation(df = annotation_df)
col_split <- annotation_df$CIM_cluster


# 🔥 GENE HIGHLIGHTING (ADD THIS)

highlight_genes <- unique(c("FOS", "EGR1", "IFIT1", "IFIT3", "ISG20", "FOSB", "OAS3",
                            
                             
                            "FAF1", "EXT1", "SND1", "UBE2E2", "HMGA1", "LIF", "ZEB1", "CDK6", "AKT3", "NT5E", "ADAM10", "VIM", "INHBA"
                            ))




row_ha <- ComplexHeatmap::rowAnnotation(
    Highlight = ComplexHeatmap::anno_mark(
        at = which(rownames(heatmap_mat) %in% highlight_genes),
        labels = intersect(highlight_genes, rownames(heatmap_mat)),
        labels_gp = grid::gpar(fontsize = 20, fontface = "bold", col = "black")
    )
)


# HEATMAP


pdf("CIMIC_heatmap.pdf", width = 10, height = 12)


ht <- ComplexHeatmap::Heatmap(
    heatmap_mat,
    name = "ΔGE",
    
    col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), 
    
    top_annotation = ha,
    right_annotation = row_ha,   # 🔥 ADDED HERE
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    show_column_names = FALSE,
    
    row_gap = grid::unit(2, "mm"),
    row_title_gp = grid::gpar(fontsize = 12, fontface = "bold", col = "black"),
    
    column_split = col_split,
    row_split = gene_group,
    
    show_row_names = FALSE,
    row_names_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    
    rect_gp = grid::gpar(col = "grey80", lwd = 0.5)
)

ComplexHeatmap::draw(ht)
dev.off()











# ORA ANALYSIS --------------------------------------


# overrepresentation analyses 
  
library(msigdbr)
library(clusterProfiler)
library(openxlsx)
library(dplyr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)


library(paletteer)
library(enrichplot)
library(patchwork)
library(tidyverse)
  library(stringr)


# # Make plot
# library(patchwork)
# 
# # # Combine the plots using patchwork
# # combined_plot <- wrap_plots(ORA_induction_pic_list, ncol = length(1:optimal_k))  # Arrange in 2 columns
# #  print(combined_plot)
#   
#   
# 
  

ORA_induction_pic_list <- list()

cluster_cols <- paste0("cluster", 1:optimal_k, "_mean")
padj_cutoff <- 0.15
differential_threshold <- 0
induction_threshold <- log2(1.50)
penetrance_threshold <- (0.20)
universe_genes <- master_induction_df$gene_id

# Generate all non-empty, non-full subsets
cluster_subsets <- unlist(
    lapply(1:(optimal_k - 1), function(m) combn(1:optimal_k, m, simplify = FALSE)),
    recursive = FALSE
)

for (subset in cluster_subsets) {
    
    cluster_name <- paste0("cluster_", paste(subset, collapse = "_"))
    cluster_gene_set_name <- paste0(cluster_name, "_specific_induced_genes")
    
    # Work on a temp copy
    df_tmp <- master_induction_df
    
    # Combined mean for subset
    target_means <- paste0("cluster", subset, "_mean")
    penetrance_name <- paste0("cluster", subset, "_penetrance")
    
    df_tmp <- df_tmp %>%
        rowwise() %>%
        mutate(combined_mean = mean(c_across(all_of(target_means)), na.rm = TRUE)) %>%
        ungroup()
    
    # Remaining clusters
    other_clusters <- setdiff(1:optimal_k, subset)
    other_means <- paste0("cluster", other_clusters, "_mean")
    
    # Differential columns: only subset vs other
    cluster_diff_cols <- c()
    for (pos in subset) {
        cluster_diff_cols <- c(
            cluster_diff_cols,
            paste0("cluster", pos, "_cluster", other_clusters, "_differential")
        )
    }
    cluster_diff_cols <- unique(cluster_diff_cols)
    
    # Filtering
    # Now add the penetrance filter
    filtered_genes <- df_tmp %>%
      filter(
        stat_test_adj_p < 0.05,
        combined_mean > induction_threshold,
        if_all(all_of(other_means), ~ . < 0),
        if_all(all_of(cluster_diff_cols), ~ . > differential_threshold),
        # Penetrance for all clusters in subset must exceed threshold
        if_all(all_of(paste0(
          "cluster", subset, "_penetrance"
        )), ~ . >= penetrance_threshold),
      )
    
    assign(cluster_gene_set_name, filtered_genes, envir = .GlobalEnv)
    message(cluster_gene_set_name, " in progress...")
    
    
    # ORA
    ORA_induction <- enrichGO(
        gene          = filtered_genes$gene_id,
        OrgDb         = org.Hs.eg.db,
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = padj_cutoff,
        qvalueCutoff  = padj_cutoff,
        universe      = universe_genes,
        keyType       = "SYMBOL"
    )
    
    ORA_name <- paste0(cluster_gene_set_name, "_ORA_analysis")
    ORA_name_simplified <- paste0("simplified_", cluster_gene_set_name, "_ORA_analysis")
        
    ORA_simplified <- clusterProfiler::simplify(
      ORA_induction,
      cutoff = 0.7,   # similarity threshold (0.7 is common, range 0–1)
      by = "p.adjust", # keep the term with the lowest adjusted p-value
      select_fun = min,
      measure = "Wang" # semantic similarity measure; Wang is common for GO
    )

    
    if (!is.null(ORA_induction) && nrow(ORA_induction@result) > 0) {
        
        Check_ORA_induction <- rownames_to_column(ORA_induction@result)
        assign(ORA_name, Check_ORA_induction, envir = .GlobalEnv)
        
        simp_Check_ORA_induction <- rownames_to_column(ORA_simplified@result)
        assign(ORA_name_simplified, simp_Check_ORA_induction, envir = .GlobalEnv)
        
        
         n_count <- 15
        
        # Top terms
        top_per_ontology <- ORA_simplified@result %>%
            # dplyr::filter(ONTOLOGY == "BP") %>%
            slice_min(order_by = p.adjust, n = n_count, with_ties = FALSE) %>%
            mutate(Description = str_wrap(Description, width = 20))
        


        
        
        p1 <- ggplot(top_per_ontology, aes(x = Description, y = Count)) +
    geom_col(aes(fill = p.adjust), color = "black") +
    coord_flip() +
    #facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
    labs(title = paste0(ORA_name, "\nTop ", n_count, " Enriched Terms"),
         x = NULL, y = "Gene Count", fill = "Adj_p" ) +
    scale_fill_gradient(low = "red", high = "blue", trans = "reverse", 
                        labels = scales::label_scientific(digits = 2)) +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text.y = element_text(size = 18, angle = 0, face = "bold"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18)
        #plot.title = element_text(hjust = 0.5)
    )
    #       p1 <- ggplot(top_per_ontology, aes(x = Description, y = Count)) +
    # geom_col(aes(fill = p.adjust), color = "black") +
    # geom_text(
    #     aes(label = Description, y = 0),   # put label halfway inside the bar
    #     color = "black", size = 4, fontface = "bold", hjust = 0
    # ) +
    # coord_flip() +
    # facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
    # labs(title = paste0(ORA_name, "\nTop ", n_count, " Enriched Terms per ONTOLOGY"),
    #      x = NULL, y = "Gene Count", fill = "Adjusted p-value") +
    # scale_fill_gradient(low = "red", high = "lightcyan", trans = "reverse") +
    # theme_minimal(base_size = 14) +
    # theme(
    #     panel.grid.major = element_blank(),
    #     panel.grid.minor = element_blank(),
    #     panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    #     strip.background = element_rect(fill = "white", colour = "black"),
    #     strip.text.y = element_text(size = 18, angle = 0, face = "bold"),
    #     axis.text.y = element_blank(),   # remove y-axis text (since it’s inside bars now)
    #     axis.text.x = element_text(size = 16, colour = "black"),
    #     plot.title = element_text(hjust = 0.5)
    # )

          
        ORA_induction_pic_list[[ORA_name]] <- p1
        
    } else if (is.null(ORA_induction))  {
  
    assign(ORA_name, NULL, envir = .GlobalEnv)
    print(paste0(ORA_name, " is empty will put as NULL in environment"))
    
}
}


# ploting simplifed ora (using BP and N = 10; changeable)

n_count <- 15

 # Top terms
 top_res_c1 <- simplified_cluster_1_specific_induced_genes_ORA_analysis %>%
   #dplyr::filter(ONTOLOGY == "BP") %>%
   slice_min(order_by = p.adjust,
             n = n_count,
             with_ties = FALSE) %>%
   mutate(Description = str_wrap(Description, width = 50))

 
 
 p1_c1_simplifed_ora <- ggplot(top_res_c1, aes(x = Description, y = Count)) +
   geom_col(aes(fill = p.adjust), color = "black") +
   coord_flip() +
   # facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
   labs(
     title = paste0("cluster 1", "\nTop ", n_count, " Enriched Terms"),
     x = NULL,
     y = "Gene Count",
     fill = "Adj_p"
   ) +
   scale_fill_gradient(low = "red",
                       high = "blue",
                       trans = "reverse",
                        labels = scales::label_scientific(digits = 2)) +
   theme_minimal(base_size = 14) +
   theme(
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_rect(
       colour = "black",
       fill = NA,
       linewidth = 1
     ),
     strip.background = element_rect(fill = "white", colour = "black"),
     strip.text.y = element_text(size = 18, angle = 0, face = "bold"),
     axis.text.y = element_text(size = 16, colour = "black"),
     axis.text.x = element_text(size = 18, colour = "black"),
     legend.text = element_text(size = 18),
     legend.title = element_text(size = 18)
     #plot.title = element_text(hjust = 0.5)
   )

 
 
  
  # Top terms
 top_res_c2 <- simplified_cluster_2_specific_induced_genes_ORA_analysis %>%
   #dplyr::filter(ONTOLOGY == "BP") %>%
   slice_min(order_by = p.adjust,
             n = n_count,
             with_ties = FALSE) %>%
   mutate(Description = str_wrap(Description, width = 50))
 

 
  p2_c2_simplifed_ora <- ggplot(top_res_c2, aes(x = Description, y = Count)) +
   geom_col(aes(fill = p.adjust), color = "black") +
   coord_flip() +
   #facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
   labs(
     title = paste0("cluster 2", "\nTop ", n_count, " Enriched Terms per ONTOLOGY"),
     x = NULL,
     y = "Gene Count",
     fill = "Adj_p"
   ) +
   scale_fill_gradient(low = "red",
                       high = "blue",
                       trans = "reverse") +
   theme_minimal(base_size = 14) +
   theme(
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_rect(
       colour = "black",
       fill = NA,
       linewidth = 1
     ),
     strip.background = element_rect(fill = "white", colour = "black"),
     strip.text.y = element_text(size = 18, angle = 0, face = "bold"),
     axis.text.y = element_text(size = 16, colour = "black"),
     axis.text.x = element_text(size = 18, colour = "black"),
     legend.text = element_text(size = 18),
     legend.title = element_text(size = 18)
     #plot.title = element_text(hjust = 0.5)
   )
  
# Correlation Analysis ------------
  
  
  

  
  
  
  
  
# Z-score Analysis ------------------------------------------------- 

  library(rstatix)
  delta_combined_TPM_lognormnorm <- delta_subset_long
  
  
  all_programs <- all_gene_sets
  
  names(all)
df <- clustered_plot_df
gene_names <- delta_combined_TPM_lognormnorm$gene_id[1:78932]



all_cols <- colnames(df)
all_programs_present <- lapply(all_programs, function(genes) intersect(genes, all_cols))
print(all_programs_present)




gene_cols <- gene_names 




zscore_safe <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(x)) else (x - m) / s
}

df_z <- df %>%
  mutate(across(all_of(gene_cols), zscore_safe))



score_program_z <- function(data, genes, program_name) {

  # ensure genes actually exist in data
  genes <- intersect(genes, colnames(data))

  if (length(genes) == 0) {
    return(
      data.frame(
        sample_id = data$sample_id,
        cimic_cluster = data$cluster_assignment,
        Program = paste0(program_name, "_z"),
        Score = rep(NA_real_, nrow(data))
      )
    )
  }
  
  score_vec <- rowMeans(as.matrix(data[, genes, drop = FALSE]), na.rm = TRUE)
  
  data.frame(
    sample_id = data$sample_id,
    cimic_cluster = data$cluster_assignment,
    Program = paste0(program_name, "_z"),
    Score = score_vec
  )
}

program_scores_all <- purrr::imap_dfr(all_programs_present, ~score_program_z(df_z, .x, .y))



library(tidyr)
library(dplyr)

program_scores_wide <- program_scores_all %>%
  dplyr::select(sample_id, Program, Score) %>%
  tidyr::pivot_wider(
    names_from = Program,
    values_from = Score
  )

clustered_plot_df_for_merge <- clustered_plot_df %>%
  left_join(program_scores_wide, by = "sample_id")












# Boxplot Z-score -------------------



library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)




library(dplyr)

# 1. Summary stats per cluster
program_summary <- program_scores_all %>%
  mutate(cimic_cluster = factor(cimic_cluster)) %>%
  group_by(cimic_cluster, Program) %>%
  summarise(
    mean_score   = mean(Score, na.rm = TRUE),
    median_score = median(Score, na.rm = TRUE),
    sd_score     = sd(Score, na.rm = TRUE),
    n            = dplyr::n(),
    .groups = "drop"
  )

# 2. P-values per program (handles k=2 or k>=3)
program_pvals <- program_scores_all %>%
  mutate(cimic_cluster = factor(cimic_cluster)) %>%
  group_by(Program) %>%
  summarise(
    n_clusters = n_distinct(cimic_cluster),
    
    p_value = tryCatch({
      if (n_clusters == 2) {
        wilcox.test(Score ~ cimic_cluster)$p.value
      } else if (n_clusters > 2) {
        kruskal.test(Score ~ cimic_cluster)$p.value
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_),
    
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(p_value, method = "fdr")
  )

# 3. Merge
program_summary_final <- program_summary %>%
  left_join(program_pvals, by = "Program")


# program names ----------------

names_of_interest <- all_gene_set_names_death[-1]

program_scores_plot_df <- program_scores_all %>%

  dplyr::filter(Program %in% paste0(names_of_interest, "_z")) %>%
  
  mutate(
    cimic_cluster = dplyr::case_when(
      cimic_cluster == 1 ~ "Fun-CIM",
      cimic_cluster == 2 ~ "Dys-CIM"
    ),
    
    # clean labels FIRST (as character)
    Program = Program %>%
      str_remove_all("GOBP_|REACTOME_|HALLMARK_|_z") %>%
      str_replace_all("_", " "),
    
    # then convert to factor with correct ordering
    Program = factor(
      Program,
      levels = unique(
        paste0(names_of_interest, "_z") %>%
          str_remove_all("GOBP_|REACTOME_|HALLMARK_|_z") %>%
          str_replace_all("_", " ")
      )
    ), 
    cimic_cluster = factor(cimic_cluster, levels = c("Dys-CIM", "Fun-CIM"))
)

pval_df_all <- program_scores_plot_df %>%
  group_by(Program) %>%
  summarise(
    n_groups = n_distinct(cimic_cluster),
    
    y.position = max(Score, na.rm = TRUE) * 1.1,
    
    pval = tryCatch({
      if (n_groups == 2) {
        wilcox.test(Score ~ cimic_cluster)$p.value
      } else if (n_groups > 2) {
        kruskal.test(Score ~ cimic_cluster)$p.value
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_),
    
    p_label = ifelse(is.na(pval), "p = NA",
                     ifelse(pval < 0.001, "p < 0.001",
                            paste0("p = ", signif(pval, 3)))),
    
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(pval, method = "fdr"),
    
    fdr_label = ifelse(is.na(FDR), "FDR = NA",
                        ifelse(FDR < 0.001, "FDR < 0.001",
                               paste0("FDR = ", signif(FDR, 3))))
  ) %>% add_significance("FDR")


p_all_programs_p <- ggplot(program_scores_plot_df, aes(x = cimic_cluster, y = Score, fill = cimic_cluster)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 1) +
  geom_text(
    data = pval_df_all,
    aes(x = 1.5, y = y.position, label = FDR.signif),
    inherit.aes = FALSE,
    size = 12,
    fontface = "bold"
  ) +
  facet_wrap(~ Program, scales = "free_y") +
  #scale_fill_manual(values = cimic_cluster_COLORS) +
  labs(
    title = "Program scores by CIMiC cluster",
    y = "Z-score of\nProgram Activation",
    x = NULL
  ) +
  theme_bw(base_line_size = 2) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 18),
    axis.text = element_text(colour = "black", size = 24),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_text(colour = "black", size = 24),
    axis.title.y = element_text(colour = "black", size = 24),
    axis.title = element_text(colour = "black", size = 24),
    legend.text = element_text(colour = "black", size = 24),
    legend.title = element_text(colour = "black", size = 24),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 2)
  )



p_all_programs_p

# Keep only significant programs
sig_programs <- pval_df_all %>%
  filter(!is.na(FDR), FDR < 0.05) %>%
  pull(Program)

program_scores_plot_df_sig <- program_scores_plot_df %>%
  filter(Program %in% sig_programs)

pval_df_sig <- pval_df_all %>%
  filter(Program %in% sig_programs)

p_all_programs_sig <- ggplot(program_scores_plot_df_sig, aes(x = cimic_cluster, y = Score, fill = cimic_cluster)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 1) +
  geom_text(
    data = pval_df_sig,
    aes(x = 1.5, y = y.position, label = FDR.signif),
    inherit.aes = FALSE,
    size = 12,
    fontface = "bold"
  ) +
  facet_wrap(~ Program, scales = "free_y") +
  #scale_fill_manual(values = cimic_cluster_COLORS) +
  labs(
    title = "Significant program scores by CIMiC cluster",
    y = "Z-score of\nProgram Activation",
    x = NULL
  ) +
  theme_bw(base_line_size = 2) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 18),
    axis.text = element_text(colour = "black", size = 24),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_text(colour = "black", size = 24),
    axis.title.y = element_text(colour = "black", size = 24),
    axis.title = element_text(colour = "black", size = 24),
    legend.text = element_text(colour = "black", size = 24),
    legend.title = element_text(colour = "black", size = 24),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 2)
  )

p_all_programs_sig



# -------------------------------------------------
# 1. Choose the gene‑set you want to plot

# Use any name that is present in `all_gene_set_names_death`
# (e.g. "Apoptosis", "Ferroptosis", "Necroptosis","Apoptosis"   "Necroptosis" "Pyroptosis"  "PANoptosis"  "Ferroptosis")
single_gene_set <- "PANoptosis"

# indivdual z-score plot --------------
program_scores_plot_df_sig <- program_scores_plot_df %>%
  filter(Program %in% single_gene_set)

pval_df_sig <- pval_df_all %>%
  filter(Program %in% single_gene_set)



y_pos_one <- max(plot_df_one$ES, na.rm = TRUE) + 0.1



p_program_individual <- ggplot(program_scores_plot_df_sig, aes(x = cimic_cluster, y = Score, fill = cimic_cluster)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 1) +
  geom_text(
    data = pval_df_sig,
    aes(x = 1.5, y = y.position, label = FDR.signif),
    inherit.aes = FALSE,
    size = 16,
    fontface = "bold"
  ) +
  facet_wrap(~ Program, scales = "free_y") +
  #scale_fill_manual(values = cimic_cluster_COLORS) +
  labs(
    y = "Z-score of\nProgram Activation",
    x = NULL
  ) +
  theme_bw(base_line_size = 2) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 32),
    axis.text = element_text(colour = "black", size = 24),
    axis.title.x = element_blank(),
    axis.text.x = element_text(colour = "black", size = 24, face = "bold"), 
    axis.text.y = element_text(colour = "black", size = 24, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 24, face = "bold"),
    axis.title = element_text(colour = "black", size = 24),
    legend.text = element_text(colour = "black", size = 24),
    legend.title = element_text(colour = "black", size = 24),
    legend.position = "none",
    #plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 2)
  )


p_program_individual

# averaged heatmap ----------- always check order of dys vs fun CIM
str_wrap_length = 40


# Averaged Heatmap with stars above columns-------------



genesets_of_interest <- paste0(names(all_gene_sets[cell_line_cimic_gene_set_names]), "_z")


cluster_heatmap_df <- clustered_plot_df_for_merge %>%
    dplyr::select(cluster_assignments, all_of(genesets_of_interest)) %>%
    pivot_longer(
        cols = -cluster_assignments,
        names_to = "Program",
        values_to = "Score"
    ) %>%
    mutate(
        Score = as.numeric(Score),
        cluster_assignments = factor(cluster_assignments, levels = c("2", "1")),
        Program = Program %>%
            str_remove_all("GOBP_|REACTOME_|HALLMARK_|_z") %>%
            str_replace_all("_", " ") %>%
            str_wrap(width = str_wrap_length)
    ) %>%
    group_by(cluster_assignments, Program) %>%
    summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop")

pval_df <- clustered_plot_df_for_merge %>%
    dplyr::select(cluster_assignments, all_of(genesets_of_interest)) %>%
    pivot_longer(
        cols = -cluster_assignments,
        names_to = "GeneSet",
        values_to = "Score"
    ) %>%
    group_by(GeneSet) %>%
    summarise(
        p =  wilcox.test(Score ~ cluster_assignments)$p.value,
        .groups = "drop"
    ) %>%
    mutate(
        p.adj = p.adjust(p, method = "BH"),
        p.adj.signif = case_when(
            p.adj <= 0.0001 ~ "****",
            p.adj <= 0.001 ~ "***",
            p.adj <= 0.01 ~ "**",
            p.adj <= 0.05 ~ "*",
            TRUE ~ ""
        )
    )

boxplot_pvalues <- pval_df

if ("p.adj.signif" %in% colnames(boxplot_pvalues)) {
    sig_df <- boxplot_pvalues %>%
        mutate(
            Program = GeneSet %>%
                str_remove_all("GOBP_|REACTOME_|HALLMARK_|_z") %>%
                str_replace_all("_", " ") %>%
                str_wrap(width = str_wrap_length),
            signif_label = p.adj.signif
        ) %>%
        dplyr::select(Program, signif_label)
} else if ("p.signif" %in% colnames(boxplot_pvalues)) {
    sig_df <- boxplot_pvalues %>%
        mutate(
            Program = GeneSet %>%
                str_remove_all("GOBP_|REACTOME_|HALLMARK_|_z") %>%
                str_replace_all("_", " ") %>%
                str_wrap(width = str_wrap_length),
            signif_label = p.signif
        ) %>%
        dplyr::select(Program, signif_label)
} else {
    sig_df <- boxplot_pvalues %>%
        mutate(
            Program = GeneSet %>%
                str_remove_all("GOBP_|REACTOME_|HALLMARK_|_z") %>%
                str_replace_all("_", " ") %>%
                str_wrap(width = str_wrap_length),
            signif_label = case_when(
                p.adj <= 0.0001 ~ "****",
                p.adj <= 0.001 ~ "***",
                p.adj <= 0.01 ~ "**",
                p.adj <= 0.05 ~ "*",
                TRUE ~ ""
            )
        ) %>%
        dplyr::select(Program, signif_label)
}

program_order <- cluster_heatmap_df %>%
    distinct(Program) %>%
    pull(Program)

cluster_heatmap_df <- cluster_heatmap_df %>%
    mutate(Program = factor(Program, levels = rev(program_order)))

sig_df <- sig_df %>%
    mutate(Program = factor(Program, levels = rev(program_order)))

star_df <- star_df <- cluster_heatmap_df %>%
    distinct(Program) %>%
    left_join(sig_df, by = "Program") %>%
    mutate(x = 2.6)


p_cluster_heatmap <- ggplot(
    cluster_heatmap_df,
    aes(x = cluster_assignments, y = Program, fill = mean_score)
) +
    geom_tile(color = "black", linewidth = 0.8, width = 1) +
    annotate(
        "text",
        x = 2.85,
        y = seq_along(levels(cluster_heatmap_df$Program)),
        label = rev(star_df$signif_label),
        size = 20,
        fontface = "bold"
    ) +
    scale_fill_gradient2(
        low = "#313695",
        mid = "white",
        high = "#A50026",
        midpoint = 0,
        name = "Mean\nZ-score\nΔ"
    ) +
    scale_x_discrete(
        limits = c("2", "1"),
        labels = c("2" = "Dys\nCIM", "1" = "Fun\nCIM"),
        expand = expansion(add = c(0.3, 1.2))
    ) +
    # scale_y_discrete(position = "right") +
    labs(
        x = NULL,
        y = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 24) +
    theme(
        axis.text.x = element_text(face = "bold", size = 15),
        axis.text.y = element_text(
            face = "bold",
            size = 10,
            margin = margin(r = 0)
        ),
        legend.text = element_text(face = "bold", size = 24),
        axis.line = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 20, b = 10, l = 20)
    )

p_cluster_heatmap






clustered_plot_df_for_merge <- clustered_plot_df_for_merge %>% 
   dplyr::mutate(
    CIMIC_Cluster = dplyr::case_when(
      cluster_assignments == 1 ~ "Fun-CIM",
      cluster_assignments == 2 ~ "Dys-CIM",
      TRUE ~ NA_character_
    )
  )



















# ssGSEA Analysis ------------------
# Make ssGSEA Function -------------


# ssGSEA Function
ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
    row_names = rownames(X)
    num_genes = nrow(X)
    gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})

    # Ranks for genes
    R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')

    # Calculate enrichment score (es) for each sample (column)
    es = apply(R, 2, function(R_col) {
        gene_ranks = order(R_col, decreasing = TRUE)

        # Calc es for each gene set
        es_sample = sapply(gene_sets, function(gene_set_idx) {
            # pos: match (within the gene set)
            # neg: non-match (outside the gene set)
            indicator_pos = gene_ranks %in% gene_set_idx
            indicator_neg = !indicator_pos

            rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha

            step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
            step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)

            step_cdf_diff = step_cdf_pos - step_cdf_neg

            # Normalize by gene number
            if (scale) step_cdf_diff = step_cdf_diff / num_genes

            # Use ssGSEA or not
            if (single) {
                sum(step_cdf_diff)
            } else {
                step_cdf_diff[which.max(abs(step_cdf_diff))]
            }
        })
        unlist(es_sample)
    })

    if (length(gene_sets) == 1) es = matrix(es, nrow = 1)

    # Normalize by absolute diff between max and min
    if (norm) es = es / diff(range(es))

    # Prepare output
    rownames(es) = names(gene_sets)
    colnames(es) = colnames(X)
    return(es)
}








# Calc SSGSEA ------


all_gene_sets <- all_gene_sets

nki_smc_ssgsea_mat <- combined_TPM_lognorm %>%
  remove_rownames() %>%
  column_to_rownames("gene_id") %>%
  as.matrix() %>%
  {
    matrix(as.numeric(.), nrow = nrow(.), dimnames = dimnames(.))
  }


# Extract columns with "_B_" in their names (Baseline)
nki_smc_ssgsea_B <- nki_smc_ssgsea_mat[, grepl(
  "_DMSO",
  base::colnames(nki_smc_ssgsea_mat)
)]

# Extract columns with "_S_" in their names (Stimulated/Post)
nki_smc_ssgsea_S <- nki_smc_ssgsea_mat[, grepl(
  paste0("_", keep_drugs),
  base::colnames(nki_smc_ssgsea_mat)
)]

# get results for baseline and surgery
system.time(assign(
  'ssgsea_es_nki_smc_baseline',
  ssgsea(nki_smc_ssgsea_B, all_gene_sets, scale = TRUE, norm = FALSE)
))

# get results for baseline and surgery
system.time(assign(
  'ssgsea_es_nki_smc_surgery',
  ssgsea(nki_smc_ssgsea_S, all_gene_sets, scale = TRUE, norm = FALSE)
))


# Get Deltas
library(dplyr)
library(stringr)
library(purrr)


get_cell_line <- function(x) {
    str_extract(x, "^[^_]+")   # gets BT549, HCC1395, etc.
}


baseline_df <- as.data.frame(ssgsea_es_nki_smc_baseline)

# transpose to work easier (samples as rows)
baseline_t <- t(baseline_df) %>% as.data.frame()
baseline_t$sample_id <- rownames(baseline_t)

baseline_t <- baseline_t %>%
    mutate(cell_line = get_cell_line(sample_id))

# average across replicates per cell line
baseline_avg <- baseline_t %>%
    group_by(cell_line) %>%
    summarise(across(-sample_id, mean, na.rm = TRUE))

# convert back to matrix format (genesets x cell lines)
baseline_avg_mat <- baseline_avg %>%
    column_to_rownames("cell_line") %>%
    t()


surgery_df <- as.data.frame(ssgsea_es_nki_smc_surgery)

surgery_t <- t(surgery_df) %>% as.data.frame()
surgery_t$sample_id <- rownames(surgery_t)

surgery_t <- surgery_t %>%
    mutate(cell_line = get_cell_line(sample_id))


delta_list <- lapply(1:nrow(surgery_t), function(i) {
  
  sample_row <- surgery_t[i, ]
  cl <- sample_row$cell_line
  
  # get matching baseline (vector of gene sets)
  baseline_vals <- baseline_avg_mat[, cl]
  
  # extract ONLY gene set columns from sample_row
  gene_vals <- as.numeric(sample_row[, rownames(baseline_avg_mat)])
  
  # subtract
  delta_vals <- gene_vals - baseline_vals
  
  return(delta_vals)
})

# bind results
delta_mat <- do.call(rbind, delta_list)
colnames(delta_mat) <- rownames(baseline_avg_mat)
rownames(delta_mat) <- surgery_t$sample_id

# transpose back → genesets x samples
delta_mat <- t(delta_mat)


delta_ssgsea <- delta_mat

ssgsea_es_nki_smc_delta_pre_post <- delta_ssgsea

# # Z Score if needed
# # rows = gene sets, cols = patients; just for visual in this case
# ssgsea_es_nki_smc_delta_zcol <- t(apply(ssgsea_es_nki_smc_delta_pre_post, 2, function(x) {
#   (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
# }))
#
# # Convert to data frame
# ssgsea_es_nki_smc_delta_zcol <- as.data.frame(t(ssgsea_es_nki_smc_delta_zcol))

# Row-wise z-score: standardize each gene set across patients; gives a sort of GSVA in this case
ssgsea_es_nki_smc_delta_zrow <- t(apply(
  ssgsea_es_nki_smc_delta_pre_post,
  1,
  function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }
))


# Convert to data frame
ssgsea_es_nki_smc_delta_zrow <- as.data.frame(ssgsea_es_nki_smc_delta_zrow)
ssgsea_es_nki_smc_delta_pre_post <- as.data.frame(
  ssgsea_es_nki_smc_delta_pre_post
)


# Set cluster as factor to preserve order
ssgsea_nki_smc_CIM_res_df_zscore <- ssgsea_es_nki_smc_delta_zrow %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  inner_join(
    clustered_plot_df %>%
      dplyr::select(sample_id, cluster_assignments) %>%
      dplyr::mutate(sample_id = stringr::str_remove(sample_id, "^delta_")),
    
    by = "sample_id"
  ) %>%
  dplyr::mutate(
    cluster_assignments = factor(cluster_assignments, levels = sort(unique(
      cluster_assignments
    ))),
    CIMIC_Cluster = dplyr::case_when(
      cluster_assignments == 1 ~ "Fun-CIM",
      cluster_assignments == 2 ~ "Dys-CIM",
      TRUE ~ NA_character_
    )
  )

# No Z socre
ssgsea_nki_smc_CIM_res_df <- ssgsea_es_nki_smc_delta_pre_post %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  inner_join(
    clustered_plot_df %>%
      dplyr::select(sample_id, cluster_assignments) %>%
      dplyr::mutate(sample_id = stringr::str_remove(sample_id, "^delta_")),
    
    by = "sample_id"
  ) %>%
  dplyr::mutate(
    cluster_assignments = factor(cluster_assignments, levels = sort(unique(
      cluster_assignments
    ))),
    CIMIC_Cluster = dplyr::case_when(
      cluster_assignments == 1 ~ "Fun-CIM",
      cluster_assignments == 2 ~ "Dys-CIM",
      TRUE ~ NA_character_
    )
  )


ssgsea_nki_smc_es_matrix <- ssgsea_nki_smc_CIM_res_df_zscore %>%
  dplyr::select(-c(cluster_assignments, CIMIC_Cluster)) %>% # drop helper columns
  tibble::column_to_rownames("sample_id") %>% # rows = samples
  t() %>% # transpose → rows = pathways
  as.matrix()



ssgsea_nki_smc_es_long <- ssgsea_nki_smc_es_matrix %>% # ES matrix: pathways × samples
  as.data.frame() %>%
  mutate(pathway = rownames(.)) %>% # keep pathway name
  pivot_longer(
    cols = -pathway,
    names_to = "sample",
    values_to = "ES"
  ) %>%
  # Join the cluster information (CIMIC_Cluster) for each sample
  left_join(
    ssgsea_nki_smc_CIM_res_df %>% # contains sample_id, CIMIC_Cluster, etc.
      dplyr::select(sample_id, CIMIC_Cluster),
    by = c("sample" = "sample_id")
  )


ssgsea_nki_smc_es_long <- as.data.frame(ssgsea_nki_smc_es_matrix) %>%
  tibble::rownames_to_column("pathway") %>% # keep pathway name
  tidyr::pivot_longer(
    cols = -pathway, # all sample columns
    names_to = "sample",
    values_to = "ES"
  ) %>%
  dplyr::left_join(
    ssgsea_nki_smc_CIM_res_df %>%
      dplyr::select(sample_id, CIMIC_Cluster),
    by = c("sample" = "sample_id")
  )


ssgsea_nki_smc_avg_es_per_cluster <- ssgsea_nki_smc_es_long %>%
  group_by(pathway, CIMIC_Cluster) %>%
  summarise(
    avg_ES = mean(ES, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(CIMIC_Cluster = paste0(CIMIC_Cluster, "_mean")) %>% # e.g. "Dys‑CIM_mean"
  pivot_wider(
    names_from = CIMIC_Cluster,
    values_from = avg_ES
  ) %>%
  column_to_rownames("pathway")


n_clusters <- length(unique(ssgsea_nki_smc_es_long$CIMIC_Cluster))


if (n_clusters == 2) {
  stats <- ssgsea_nki_smc_es_long %>%
    group_by(pathway) %>%
    rstatix::wilcox_test(ES ~ CIMIC_Cluster) %>%
    dplyr::select(pathway, p = p)
} else {
  stats <- ssgsea_nki_smc_es_long %>%
    group_by(pathway) %>%
    rstatix::kruskal_test(ES ~ CIMIC_Cluster) %>%
    dplyr::select(pathway, p = p)
}


ssgsea_nki_smc_avg_es_per_cluster <- ssgsea_nki_smc_avg_es_per_cluster %>%
  rownames_to_column("pathway") %>%
  left_join(stats, by = "pathway") %>%
  column_to_rownames("pathway")

ssgsea_nki_smc_avg_es_per_cluster$adj_p <-
  p.adjust(ssgsea_nki_smc_avg_es_per_cluster$p, method = "fdr")


# get dataset merged
ssgsea_nki_smc_CIM_res_df_for_merge <- ssgsea_nki_smc_CIM_res_df_zscore %>% 
  dplyr::mutate(
  sample_id = paste0("delta_", sample_id)
      )


# Merge
clustered_plot_df_for_merge <- merge(
  x = clustered_plot_df,
  y = ssgsea_nki_smc_CIM_res_df_for_merge,
  by = "sample_id"
)

clustered_plot_df_for_merge$cluster_assignments <- clustered_plot_df_for_merge$cluster_assignments.x



clustered_plot_df_for_merge <- clustered_plot_df_for_merge %>%
  mutate(
    Avg_MHC1 = rowMeans(
      dplyr::select(., `HLA-A`, `HLA-B`, `HLA-C`), # keep only the three HLA columns
      na.rm = TRUE
    )
  )











# Boxplots of genesets ssGSEA  -----------------


library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
library(tidyr)
library(stringr) # <-- provides str_wrap()
library(rstatix)


# Gene sets of interest TOGGLE here 
genesets_of_interest <- c(all_gene_set_names_cimic, "vm")


clustered_plot_df_for_merge$CIMIC_Cluster <- factor(
  clustered_plot_df_for_merge$CIMIC_Cluster,
  levels = c("Dys-CIM", "Fun-CIM")
)


# Prepare data
plot_df <- clustered_plot_df_for_merge %>%
  dplyr::select(all_of(c("CIMIC_Cluster", genesets_of_interest))) %>%
  pivot_longer(-CIMIC_Cluster, names_to = "GeneSet", values_to = "ES") %>%
  mutate(
    ES = as.numeric(ES),
    GeneSet = str_remove_all(GeneSet, "GOBP_|REACTOME_|HALLMARK_") %>%
      str_replace_all("_", " ") %>%
      str_wrap(width = 25)
  )

# Compute max per GeneSet for stat_compare_means positioning
y_positions <- plot_df %>%
  group_by(GeneSet) %>%
  summarise(y_pos = max(ES, na.rm = TRUE) + 0.1)

# Generate all pairwise comparisons
comparisons <- combn(unique(plot_df$CIMIC_Cluster), 2, simplify = FALSE)

# Faceted boxplot with stat_compare_means

# Faceted boxplot with stat_compare_means
ggplot(plot_df, aes(x = CIMIC_Cluster, y = ES, fill = CIMIC_Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, lwd = 1.2) +
  geom_jitter(width = 0.2, size = 0.75, alpha = 0.7,) +
  facet_wrap(~GeneSet, scales = "free_y", ncol = 4) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    size = 15,
    fontface = "bold",
    hide.ns = TRUE,
    vjust = 1.2
  ) +
  labs(
    x = "Cluster",
    y = "Z-Score of ΔssGSEA"
    #title = "Gene Set Enrichment Across Subtypes with Significance"
  ) +
  theme_bw(base_line_size = 2) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 18),
    axis.text = element_text(colour = "black", size = 24),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_text(colour = "black", size = 24),
    axis.title.y = element_text(colour = "black", size = 24),
    axis.title = element_text(colour = "black", size = 24),
    legend.text = element_text(colour = "black", size = 24),
    legend.title = element_text(colour = "black", size = 24),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 2)
  )


library(dplyr)
library(tidyr)
library(rstatix)

# Gene sets of interest

# Prepare long dataframe (same as plotting)
plot_df <- clustered_plot_df_for_merge %>%
  dplyr::select(all_of(c("CIMIC_Cluster", genesets_of_interest))) %>%
  pivot_longer(-CIMIC_Cluster, names_to = "GeneSet", values_to = "ES") %>%
  mutate(ES = as.numeric(ES))

# Perform pairwise Wilcoxon tests per GeneSet
pval_table <- plot_df %>% as.data.frame() %>%
  group_by(GeneSet) %>%
  wilcox_test(ES ~ CIMIC_Cluster) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  ungroup()

# Clean gene set labels (optional, matching your plot)
pval_table <- pval_table %>%
  mutate(
    GeneSet = str_remove_all(GeneSet, "GOBP_|REACTOME_|HALLMARK_") %>%
      str_replace_all("_", " ")
  )

# Save object for downstream use
boxplot_pvalues <- pval_table

boxplot_pvalues



# Use any name that is present in `all_gene_set_names_death`
# (e.g. "Apoptosis", "Ferroptosis", "Necroptosis","Apoptosis"   "Necroptosis" "Pyroptosis"  "PANoptosis"  "Ferroptosis")
single_gene_set <- "GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS"

single_gene_set_label <- single_gene_set %>%
  str_remove_all("GOBP_|REACTOME_|HALLMARK_") %>%
  str_replace_all("_", " ") %>%
  str_wrap(width = 25)


# 2. Prepare a tidy data frame for that gene-set only
# -------------------------------------------------
plot_df_one <- clustered_plot_df_for_merge %>%
  dplyr::select(all_of(c("CIMIC_Cluster", single_gene_set))) %>%
  pivot_longer(
    -CIMIC_Cluster,
    names_to = "GeneSet",
    values_to = "ES"
  ) %>%
  mutate(
    ES = as.numeric(ES),
    GeneSet = single_gene_set_label
  )


y_pos_one <- max(plot_df_one$ES, na.rm = TRUE) + 0.1


ggplot(plot_df_one,
       aes(x = CIMIC_Cluster, y = ES, fill = CIMIC_Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, lwd = 2.25, width = 0.75) +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.7) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    size = 30,
    fontface = "bold",
    hide.ns = TRUE,
    vjust = 1.2
  ) +
  labs(
    x = "Cluster",
    y = "Z‑Score of ΔssGSEA",
    title = single_gene_set_label
  ) +
  theme_bw(base_line_size = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "bold", colour = "black", size = 45),
    axis.text.y = element_text(face = "bold", colour = "black", size = 45),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(face = "bold", colour = "black", size = 45),
    axis.title = element_text(colour = "black", size = 45),
    legend.text = element_text(colour = "black", size = 24),
    legend.title = element_text(colour = "black", size = 24),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 45),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 2)
  )


















# Boxplots for cell death programs ssGSEA ------------------




library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
library(tidyr)
library(stringr) # <-- provides str_wrap()
library(rstatix)


# Gene sets of interest TOGGLE here 
genesets_of_interest <- all_gene_set_names_death 


clustered_plot_df_for_merge$CIMIC_Cluster <- factor(
  clustered_plot_df_for_merge$CIMIC_Cluster,
  levels = c("Dys-CIM", "Fun-CIM")
)


# Prepare data
plot_df <- clustered_plot_df_for_merge %>%
  dplyr::select(all_of(c("CIMIC_Cluster", genesets_of_interest))) %>%
  pivot_longer(-CIMIC_Cluster, names_to = "GeneSet", values_to = "ES") %>%
  mutate(
    ES = as.numeric(ES),
    GeneSet = str_remove_all(GeneSet, "GOBP_|REACTOME_|HALLMARK_") %>%
      str_replace_all("_", " ") %>%
      str_wrap(width = 25)
  ) %>% as.data.frame()

# Compute max per GeneSet for stat_compare_means positioning
y_positions <- plot_df %>%
  group_by(GeneSet) %>%
  summarise(y_pos = max(ES, na.rm = TRUE) + 0.1)

# Generate all pairwise comparisons
comparisons <- combn(unique(plot_df$CIMIC_Cluster), 2, simplify = FALSE)

# Faceted boxplot with stat_compare_means

# Faceted boxplot with stat_compare_means
ggplot(plot_df, aes(x = CIMIC_Cluster, y = ES, fill = CIMIC_Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, lwd = 1.2) +
  geom_jitter(width = 0.2, size = 0.75, alpha = 0.7,) +
  facet_wrap(~GeneSet, scales = "free_y", ncol = 4) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    size = 15,
    fontface = "bold",
    hide.ns = TRUE,
    vjust = 1.2
  ) +
  labs(
    x = "Cluster",
    y = "Z-Score of ΔssGSEA"
    #title = "Gene Set Enrichment Across Subtypes with Significance"
  ) +
  theme_bw(base_line_size = 2) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 18),
    axis.text = element_text(colour = "black", size = 24),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_text(colour = "black", size = 24),
    axis.title.y = element_text(colour = "black", size = 24),
    axis.title = element_text(colour = "black", size = 24),
    legend.text = element_text(colour = "black", size = 24),
    legend.title = element_text(colour = "black", size = 24),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 2)
  )


library(dplyr)
library(tidyr)
library(rstatix)

# Gene sets of interest

# Prepare long dataframe (same as plotting)
plot_df <- clustered_plot_df_for_merge %>%
  dplyr::select(all_of(c("CIMIC_Cluster", genesets_of_interest))) %>%
  pivot_longer(-CIMIC_Cluster, names_to = "GeneSet", values_to = "ES") %>%
  mutate(ES = as.numeric(ES)) %>% as.data.frame()

# Perform pairwise Wilcoxon tests per GeneSet
pval_table <- plot_df %>%
  group_by(GeneSet) %>%
  wilcox_test(ES ~ CIMIC_Cluster) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  ungroup()

# Clean gene set labels (optional, matching your plot)
pval_table <- pval_table %>%
  mutate(
    GeneSet = str_remove_all(GeneSet, "GOBP_|REACTOME_|HALLMARK_") %>%
      str_replace_all("_", " ")
  )

# Save object for downstream use
boxplot_pvalues <- pval_table

boxplot_pvalues

# Individual box

library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
library(tidyr)
library(rstatix)

# -------------------------------------------------
# 1. Choose the gene‑set you want to plot
# -------------------------------------------------
# Use any name that is present in `all_gene_set_names_death`
# (e.g. "Apoptosis", "Ferroptosis", "Necroptosis","Apoptosis"   "Necroptosis" "Pyroptosis"  "PANoptosis"  "Ferroptosis")
single_gene_set <- "Pyroptosis"

# -------------------------------------------------
# 2. Prepare a tidy data frame for that gene‑set only
# -------------------------------------------------
plot_df_one <- clustered_plot_df_for_merge %>%
  dplyr::select(all_of(c("CIMIC_Cluster", single_gene_set))) %>%      # keep cluster + the chosen column
  pivot_longer(-CIMIC_Cluster,
               names_to = "GeneSet",
               values_to = "ES") %>%
  mutate(
    ES = as.numeric(ES),
    # clean the label to match the style used in the faceted plots
    GeneSet = str_remove_all(GeneSet, "GOBP_|REACTOME_|HALLMARK_") %>%
              str_replace_all("_", " ")
  ) %>% as.data.frame()

# -------------------------------------------------
# 3. Compute y‑position for the significance label
# -------------------------------------------------
y_pos_one <- max(plot_df_one$ES, na.rm = TRUE) + 0.1

# -------------------------------------------------
# 4. Plot – single panel, Wilcoxon annotation
# -------------------------------------------------
ggplot(plot_df_one,
       aes(x = CIMIC_Cluster, y = ES, fill = CIMIC_Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, lwd = 2.25, width = 0.75) +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.7) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    size = 30,
    fontface = "bold",
    hide.ns = TRUE,
    vjust = 1.2
  ) +
  labs(
    x = "Cluster",
    y = "Z‑Score of ΔssGSEA",
    title = single_gene_set
  ) +
  theme_bw(base_line_size = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "bold", colour = "black", size = 45),
    axis.text.y = element_text(face = "bold", colour = "black", size = 45),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(face = "bold", colour = "black", size = 45),
    axis.title = element_text(colour = "black", size = 45),
    legend.text = element_text(colour = "black", size = 24),
    legend.title = element_text(colour = "black", size = 24),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 45),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 2)
  )



























# Jaccard Index -------------------------------


# 2. Run on deltas

jacc_delta_mat <- delta_subset_long %>% column_to_rownames("gene_id") %>% as.matrix()


# 1. Define the Jaccard function (slightly modified for easier looping)
calc_jaccard_single <- function(delta_vec, pathway, top_p = 0.05) {
  threshold <- quantile(delta_vec, 1 - top_p, na.rm = TRUE)
  top_genes <- names(delta_vec[delta_vec >= threshold])
  
  inter <- length(intersect(top_genes, pathway))
  uni <- length(union(top_genes, pathway))
  return(inter / uni)
}

# 2. Loop through every gene set in your list
# This creates a matrix where Rows = Pathways and Cols = Patients
jaccard_delta_results <- lapply(all_gene_sets, function(gs) {
  apply(jacc_delta_mat, 2, calc_jaccard_single, pathway = gs)
})


# 3. Convert the list into a clean matrix
jaccard_res_mat <- do.call(rbind, jaccard_delta_results)

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

# 1. Transform the matrix to long format
jaccard_long <- as.data.frame(jaccard_res_mat) %>%
  tibble::rownames_to_column("Pathway") %>%
  pivot_longer(-Pathway, names_to = "sample_id", values_to = "Jaccard_Index") %>%
  mutatz

# 2. Join with your metadata (Ensure you have a 'CIMIC_Cluster' column)
# If your sample IDs in the matrix have 'delta_' prefix, make sure they match your meta
plot_data <- jaccard_long %>%
  inner_join(ssgsea_nki_smc_CIM_res_df %>% mutate(sample_id = paste0("delta_", sample_id))
               , by = "sample_id")

# 3. Create the Faceted Plot
ggplot(plot_data, aes(x = CIMIC_Cluster, y = Jaccard_Index, fill = CIMIC_Cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~Pathway, scales = "free_y", ncol = 3) + # 'free_y' is key as stress has higher baseline
  stat_compare_means(method = "t.test", label = "p.format", label.x = 1.5) + 
  scale_fill_manual(values = c("Fun-CIM" = "#74add1", "Dys-CIM" = "#f46d43")) +
  theme_minimal() +
  labs(title = "Delta Jaccard Index: Pathway Responsiveness",
       subtitle = "Top 5% of Induced Genes",
       y = "Jaccard Index (Overlap / Union)",
       x = "CIMIC Cluster") +
  theme(strip.text = element_text(face = "bold", size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1))



# Correlation CLuster Analysis --------------------


# Correlation with Cluster Assignment #################

library(dplyr)
library(purrr)
library(broom)
library(ggplot2)
library(tidyr)

# Gene set
human_CIM_genes <- unique(unlist(cell_line_cimic_gene_sets))

# Initialize lists
all_cor_results <- list()
sig_cor_results <- list()

# choose your dataset
ds_data <- clustered_plot_df   # <- change if needed

# Ensure numeric cluster
ds_data <- ds_data %>%
    mutate(cluster_numeric = as.numeric(cluster_assignments))

# ---- Gene correlation ----
expr_mat <- ds_data %>% dplyr::select(any_of(human_CIM_genes))

cor_df <- purrr::map_dfr(
    colnames(expr_mat),
    function(gene){
        test <- cor.test(expr_mat[[gene]], ds_data$cluster_numeric, method = "spearman")
        
        tibble(
            gene = gene,
            cor = test$estimate,
            cor_pval = test$p.value,
            direction = ifelse(test$estimate > 0,
                               "Fun-CIM",
                               "Dys-CIM")
        )
    }
)

# Adjust p-values
cor_df <- cor_df %>%
    mutate(cor_padj = p.adjust(cor_pval, method = "fdr")) 


# -----------------------------
# 2. MEAN EXPRESSION BY CLUSTER
# -----------------------------
ds_data$group <- factor(ds_data$cluster_numeric)
groups <- levels(ds_data$group)

cluster_means <- map_dfr(colnames(expr_mat), function(gene) {
    
    expr <- ds_data[[gene]]
    
    # ---- means per group (dynamic) ----
    mean_vals <- tapply(expr, ds_data$group, mean, na.rm = TRUE)
    
    mean_df <- as.data.frame(as.list(mean_vals))
    names(mean_df) <- paste0("mean_cluster", names(mean_vals))
    
    # ensure consistency for 2-group downstream ops
    mean_cluster1 <- if ("1" %in% names(mean_vals)) mean_vals[["1"]] else NA
    mean_cluster2 <- if ("2" %in% names(mean_vals)) mean_vals[["2"]] else NA
    
    # ---- differential ----
    cluster2_vs_1 <- mean_cluster2 - mean_cluster1
    cluster1_vs_2 <- mean_cluster1 - mean_cluster2
    
    # ---- test selection ----
    df_tmp <- data.frame(expr = expr, group = ds_data$group)
    
    if (length(groups) == 2) {
        test <- wilcox.test(expr ~ group, data = df_tmp)
        test_name <- "wilcox_pval"
    } else {
        test <- kruskal.test(expr ~ group, data = df_tmp)
        test_name <- "kruskal_pval"
    }
    
    tibble(
        gene = gene,
        !!!mean_df,
        cluster2_vs_1 = cluster2_vs_1,
        cluster1_vs_2 = cluster1_vs_2,
        expr_pval = test$p.value
    )
})


# -----------------------------
# 4. FDR CORRECTION
# -----------------------------
cluster_means <- cluster_means %>%
    mutate(
        expr_padj = p.adjust(expr_pval, method = "fdr")
    )


final_cor_df <- cor_df %>%
    left_join(cluster_means, by = "gene")


# ---- Significant subset ----

corr_thresh <- 0.6
ge_thresh <- 0

sig_cor_df <- final_cor_df %>%
    filter(cor_padj < 0.05, abs(cor) >= corr_thresh, abs(cluster2_vs_1) >= log2(ge_thresh)) %>%
    mutate(direction = ifelse(cor > 0, "CIM Functional", "CIM Dysfunctional"))






# ---- Visualization ----

library(dplyr)

top50_each_direction <- sig_cor_df %>%
    
    # define direction explicitly from biology
    mutate(
        direction = case_when(
            cluster2_vs_1 > 0 ~ "Dys-CIM",
            cluster2_vs_1 < 0 ~ "Fun-CIM"
        )
    ) %>%
    
    group_by(direction) %>%
    
    # rank within each direction
    # rank by effect size first
    dplyr::arrange(desc(abs(cluster2_vs_1)), .by_group = TRUE) %>%
    
    slice_head(n = 50) %>%
    
    ungroup()

top50_each_direction <- top50_each_direction %>%
    dplyr::group_by(direction) %>%
    dplyr::arrange(gene, .by_group = TRUE) %>%
    dplyr::ungroup()

ggplot(top50_each_direction, aes(x = gene, y = abs(cor), fill = direction)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    coord_flip() +
    # --- Add the horizontal line here ---
    geom_hline(yintercept = 0.5, linetype = "longdash", color = "black", linewidth = 1) +
    facet_wrap(~direction, scales = "free_y") +
    labs(
        # x = "Gene",
        y = "Absolute Correlation (|r|)", # Updated label for clarity
        # fill = "Dataset",
    ) +
    theme_classic(base_size = 18) +
    theme(
        axis.text.y = element_text(size = 18, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 18, colour = "black", face = "bold"),
        strip.text.x = element_text(size = 18, face = "bold") 
    )


# Extraction ###############


# Separate positive vs negative
positive_genes <- sig_cor_df %>%
    filter(cor > 0) %>%
    distinct(gene) %>%
    pull(gene)

negative_genes <- sig_cor_df %>%
    filter(cor < 0) %>%
    distinct(gene) %>%
    pull(gene)

# Check
length(positive_genes)
length(negative_genes)

# ORA ANALYSIS ###############


# overrepresentation analyses 

library(msigdbr)
library(clusterProfiler)
library(openxlsx)
library(dplyr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)


library(paletteer)
library(enrichplot)
library(patchwork)
library(tidyverse)
library(stringr)


combined_shared_genes <- master_induction_df$gene_id



universe_genes <- combined_shared_genes
padj_cutoff <- 0.05


# ORA
shared_neg_ORA <- enrichGO(
    gene          = negative_genes,
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = padj_cutoff,
    qvalueCutoff  = padj_cutoff,
    universe      = universe_genes,
    keyType       = "SYMBOL"
)

# ORA
shared_pos_ORA <- enrichGO(
    gene          = positive_genes,
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = padj_cutoff,
    qvalueCutoff  = padj_cutoff,
    universe      = universe_genes,
    keyType       = "SYMBOL"
)



# reactome 
library(ReactomePA)


# Load necessary packages
library(clusterProfiler)
library(org.Hs.eg.db)

# Convert gene symbols to Entrez IDs
universe_genes_entrez <- bitr(
    universe_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
)


# Convert gene symbols to Entrez IDs
negative_genes_entrez <- bitr(
    negative_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
)



shared_neg_reactome <- enrichPathway(
    negative_genes_entrez$ENTREZID,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    universe = universe_genes_entrez$ENTREZID,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE
)


# Convert gene symbols to Entrez IDs
positive_genes_entrez <- bitr(
    positive_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
)


shared_pos_reactome <- enrichPathway(
    positive_genes_entrez$ENTREZID,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    universe = universe_genes_entrez$ENTREZID,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE
)


# VISUALIZE ORA ##########


pos_ORA_simplified <- clusterProfiler::simplify(
    shared_pos_ORA,
    cutoff = 0.7,   # similarity threshold (0.7 is common, range 0–1)
    by = "p.adjust", # keep the term with the lowest adjusted p-value
    select_fun = min,
    measure = "Wang" # semantic similarity measure; Wang is common for GO
)



neg_ORA_simplified <- clusterProfiler::simplify(
    shared_neg_ORA,
    cutoff = 0.7,   # similarity threshold (0.7 is common, range 0–1)
    by = "p.adjust", # keep the term with the lowest adjusted p-value
    select_fun = min,
    measure = "Wang" # semantic similarity measure; Wang is common for GO
)


n_count <- 15

# Top terms
top_res_c1_pos_go<- pos_ORA_simplified@result %>%
    dplyr::filter(ONTOLOGY == "BP") %>%
    slice_min(order_by = p.adjust,
              n = n_count,
              with_ties = FALSE) %>%
    mutate(Description = str_wrap(Description, width = 50))




p1_pos_ora_go <- ggplot(top_res_c1_pos_go, aes(x = Description, y = Count)) +
    geom_col(aes(fill = p.adjust), color = "black") +
    coord_flip() +
    #facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
    labs(
        #title = paste0("positive", "\nTop ", n_count, " Enriched Terms per ONTOLOGY"),
        x = NULL,
        y = "Gene Count",
        fill = "Adj_p"
    ) +
    scale_fill_gradient(low = "red",
                        high = "blue",
                        trans = "reverse") +
    theme_minimal(base_size = 14) +
    scale_fill_gradient(low = "red", high = "blue", trans = "reverse", 
                        labels = scales::label_scientific(digits = 2)) +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text.y = element_text(size = 18, angle = 0, face = "bold"),
        axis.text.y = element_text(size = 16, colour = "black", , face = "bold"),
        axis.text.x = element_text(size = 18, colour = "black", face = "bold"),
        axis.title.x = element_text(size = 18, colour = "black", face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18)
        #plot.title = element_text(hjust = 0.5)
    )


# highlight_terms <- c(
#     "protein folding",
#     "positive regulation of execution phase of apoptosis",
#     "proteasome assembly",
#     "intrinsic apoptotic signaling pathway",
#     "protein-DNA complex organization",
#     "RNA catabolic process"
# )

top_neg_go <- neg_ORA_simplified@result %>%
    dplyr::filter(ONTOLOGY == "BP")

# top 15 by significance
top15 <- top_neg_go %>%
      # slice_min(order_by = p.adjust, n = n_count-length(highlight_terms), with_ties = FALSE)
      slice_min(order_by = p.adjust, n = n_count, with_ties = FALSE)
    
# forced terms (even if outside top 15)
# forced <- top_neg_go %>%
#     filter(Description %in% highlight_terms)

# combine + remove duplicates
top_res_c1_neg_go <- bind_rows(top15, forced) %>%
    distinct(Description, .keep_all = TRUE) %>%
    mutate(Description = str_wrap(Description, width = 40))



p1_neg_ora_go <- ggplot(top_res_c1_neg_go, aes(x = Description, y = Count)) +
    geom_col(aes(fill = p.adjust), color = "black") +
    coord_flip() +
    #facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
    labs(
        #title = paste0("negative", "\nTop ", n_count, " Enriched Terms per ONTOLOGY"),
        x = NULL,
        y = "Gene Count",
        fill = "Adj_p"
    ) +
    scale_fill_gradient(low = "red", high = "blue", trans = "reverse", 
                        labels = scales::label_scientific(digits = 2)) +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text.y = element_text(size = 18, angle = 0, face = "bold"),
        axis.text.y = element_text(size = 16, colour = "black", , face = "bold"),
        axis.text.x = element_text(size = 18, colour = "black", face = "bold"),
        axis.title.x = element_text(size = 18, colour = "black", face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18)
        #plot.title = element_text(hjust = 0.5)
    )

# Reactome VISUALIZATION ###############

# Top terms
top_res_c1_pos_react <- shared_pos_reactome@result %>%
    slice_min(order_by = p.adjust,
              n = n_count,
              with_ties = FALSE) %>%
    mutate(Description = str_wrap(Description, width = 50))




p1_pos_ora_react <- ggplot(top_res_c1_pos_react, aes(x = Description, y = Count)) +
    geom_col(aes(fill = p.adjust), color = "black") +
    coord_flip() +
    # facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
    labs(
        #title = paste0("positive", "\nTop ", n_count, " Enriched Terms per ONTOLOGY"),
        x = NULL,
        y = "Gene Count",
        fill = "Adj_p"
    ) +
    scale_fill_gradient(low = "red",
                        high = "blue",
                        trans = "reverse") +
    theme_minimal(base_size = 14) +
    scale_fill_gradient(low = "red", high = "blue", trans = "reverse", 
                        labels = scales::label_scientific(digits = 2)) +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text.y = element_text(size = 18, angle = 0, face = "bold"),
        axis.text.y = element_text(size = 16, colour = "black", , face = "bold"),
        axis.text.x = element_text(size = 18, colour = "black", face = "bold"),
        axis.title.x = element_text(size = 18, colour = "black", face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18)
        #plot.title = element_text(hjust = 0.5)
    )



# Top terms
top_res_c1_neg_react <- shared_neg_reactome@result %>%
    slice_min(order_by = p.adjust,
              n = n_count,
              with_ties = FALSE) %>%
    mutate(Description = str_wrap(Description, width = 50))




p1_neg_ora_react <- ggplot(top_res_c1_neg_react, aes(x = Description, y = Count)) +
    geom_col(aes(fill = p.adjust), color = "black") +
    coord_flip() +
    # facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
    labs(
        #title = paste0("negative", "\nTop ", n_count, " Enriched Terms per ONTOLOGY"),
        x = NULL,
        y = "Gene Count",
        fill = "Adj_p"
    ) +
    scale_fill_gradient(low = "red", high = "blue", trans = "reverse", 
                        labels = scales::label_scientific(digits = 2)) +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text.y = element_text(size = 18, angle = 0, face = "bold"),
        axis.text.y = element_text(size = 16, colour = "black", , face = "bold"),
        axis.text.x = element_text(size = 18, colour = "black", face = "bold"),
        axis.title.x = element_text(size = 18, colour = "black", face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18)
        #plot.title = element_text(hjust = 0.5)
    )

go_pos <- pos_ORA_simplified@result %>% dplyr::filter(ONTOLOGY == "BP")
go_neg <- neg_ORA_simplified@result %>% dplyr::filter(ONTOLOGY == "BP")

react_pos <- shared_pos_reactome@result
react_neg <- shared_neg_reactome@result










# Analyzing individual genes induction ----------------



library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)      # ggarrange / annotate_figure
library(rstatix)     # wilcox_test, kruskal_test, etc.
library(stringr)     # str_wrap
library(purrr)




compare_induction_boxplot <- function(df, gene) {

  gene_col <- df %>%
    dplyr::select(dplyr::matches(paste0("^", gene, "(\\.x|\\.y)?$"))) %>%
    names()
  if (length(gene_col) == 0) {
    stop("Gene column '", gene, "' not found in the data frame")
  }
  gene_col <- gene_col[[1]]   # use the first match


  df2 <- df %>%
    dplyr::mutate(
      Expression = as.numeric(as.character(.data[[gene_col]]))
    ) %>%
    dplyr::select(CIMIC_Cluster, Expression)   # keep only what we need


  pval <- df2 %>%
    rstatix::wilcox_test(Expression ~ CIMIC_Cluster) %>%
    dplyr::pull(p)

  ptxt <- paste0("P = ", format.pval(pval, digits = 4))


  ggplot(df2,
         aes(x = CIMIC_Cluster,
             y = Expression,
             fill = CIMIC_Cluster)) +
    geom_boxplot(width = 0.6,
                 outlier.shape = NA,
                 alpha = 0.8) +
    geom_jitter(width = 0.15,
                size = 1.5,
                alpha = 0.4) +
    annotate("text",
             x = 1.5,
             y = max(df2$Expression, na.rm = TRUE) * 1.1,
             label = ptxt,
             size = 8,
             fontface = "bold") +
    scale_fill_manual(values = c("Dys-CIM" = "#2980B9",
                                 "Fun-CIM" = "#C0392B")) +
    labs(title = paste0(gene),
         y = "Δ(log2[TPM+1])",
         x = NULL) +
    theme_classic(base_size = 14) +
    theme(
      plot.title      = element_text(hjust = 0.5,
                                     face = "bold",
                                     size = 16),
      legend.position = "none",
      axis.text.x     = element_text(face = "bold", size = 20),
      axis.text.y     = element_text(face = "bold", size = 20),
      axis.title.y    = element_text(face = "bold", size = 20),
      panel.border    = element_rect(colour = "black",
                                     fill = NA,
                                     linewidth = 1.5),
      axis.line       = element_line(colour = "black",
                                     linewidth = 1.5)
    )
}




genes_of_interest <- c("PDCD1", "CD274", "IRF1")
ind_plots <- purrr::map(genes_of_interest,
                    ~ compare_induction_boxplot(clustered_plot_df_for_merge, .x))






# Analyze indivdual gene penetrance ------------------






library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)      # ggarrange / annotate_figure
library(rstatix)     # wilcox_test, kruskal_test, etc.
library(stringr)     # str_wrap
library(purrr)





compare_penetrance_plot <- function(df, gene) {
  # -------------------------------------------------
  # 1.1  Find the exact column that holds the gene values
  # -------------------------------------------------
  gene_col <- df %>%
    dplyr::select(dplyr::matches(paste0("^", gene, "(\\.x|\\.y)?$"))) %>%
    names()
  if (length(gene_col) == 0) {
    stop("Gene column '", gene, "' not found in the data frame")
  }
  gene_col <- gene_col[[1]]

  # -------------------------------------------------
  # 1.2  Build a tidy frame with a binary “Induction” flag
  # -------------------------------------------------
  df2 <- df %>%
    dplyr::mutate(
      Expression = as.numeric(as.character(.data[[gene_col]])),
      Induction  = ifelse(Expression > 0,
                          "Positive Induction",
                          "Negative/No Induction")
    ) %>%
    dplyr::select(CIMIC_Cluster, Induction)

  # -------------------------------------------------
  # 1.3  Fisher or χ² test on the 2×2 table
  # -------------------------------------------------
  tab <- df2 %>%
    dplyr::count(CIMIC_Cluster, Induction) %>%
    tidyr::pivot_wider(names_from = Induction,
                       values_from = n,
                       values_fill = 0)

  mat <- as.matrix(tab[ , c("Negative/No Induction", "Positive Induction")])
  pval <- if (any(mat < 5) || sum(mat) < 50) {
    stats::fisher.test(mat)$p.value
  } else {
    stats::chisq.test(mat)$p.value
  }
  ptxt <- paste0("P = ", format.pval(pval, digits = 4))

  # -------------------------------------------------
  # 1.4  Prepare data for the stacked‑bar plot
  # -------------------------------------------------
  bar_df <- df2 %>%
    dplyr::group_by(CIMIC_Cluster, Induction) %>%
    dplyr::summarise(n = n(), .groups = "drop") %>%
    dplyr::group_by(CIMIC_Cluster) %>%
    dplyr::mutate(
      Total = sum(n),
      Prop  = n / Total
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(
      CIMIC_Cluster,
      factor(Induction,
             levels = c("Negative/No Induction", "Positive Induction"))
    )

  # -------------------------------------------------
  # 1.5  Plot -------------------------------------------------
  # -------------------------------------------------
  ggplot(bar_df,
         aes(x = CIMIC_Cluster,
             y = Prop,
             fill = Induction)) +
    geom_bar(stat = "identity",
             position = "fill",
             width = 0.7,
             aes(fill = factor(Induction,
                               levels = c("Negative/No Induction",
                                          "Positive Induction")))) +
    geom_text(aes(y = Prop,
                  label = ifelse(Prop > 0,
                                 paste0(round(Prop * 100), "%"),
                                 "")),
              position = position_stack(vjust = 0.5),
              colour = "black",
              size = 8,
              fontface = "bold") +
    geom_text(aes(x = 1.5, y = 1.06, label = ptxt),
              inherit.aes = FALSE,
              size = 8,
              fontface = "bold") +
    geom_text(aes(x = CIMIC_Cluster,
                  y = -0.05,
                  label = paste0("N = ", Total)),
              inherit.aes = FALSE,
              size = 6,
              fontface = "bold") +
    labs(title = paste0(gene),
         y = "Percentage",
         fill = "Induction") +
    scale_fill_manual(values = c("Positive Induction" = "#3fb949ff",
                                 "Negative/No Induction" = "#f76661ff")) +
    scale_y_continuous(labels = scales::percent,
                       limits = c(-0.1, 1.12),
                       breaks = seq(0, 1, 0.2)) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    theme_classic(base_size = 18) +
    theme(
      plot.title      = element_text(hjust = 0.5,
                                     face = "bold",
                                     size = 18),
      legend.position = "top",
      panel.border    = element_rect(colour = "black",
                                     fill = NA,
                                     linewidth = 1.5),
      axis.line.x     = element_line(colour = "black",
                                     linewidth = 1.5),
      axis.line.y     = element_line(colour = "black",
                                     linewidth = 1.5),
      axis.text.x     = element_text(face = "bold",
                                     size = 20),
      axis.text.y     = element_text(face = "bold",
                                     size = 20),
      axis.title.y    = element_text(face = "bold",
                                     size = 20),
      legend.text     = element_text(size = 18),
      legend.title    = element_text(size = 18,
                                     face = "bold")
    )
}

genes_of_interest <- c("PDCD1", "CD274", "IRF1")
pen_plots <- purrr::map(genes_of_interest,
                    ~ compare_penetrance_plot(clustered_plot_df_for_merge, .x))





# Analyze PBSC for individual Genes -------------




compare_pointbiserial_plot <- function(df, gene) {

  gene_col <- df %>%
    dplyr::select(dplyr::matches(paste0("^", gene, "(\\.x|\\.y)?$"))) %>%
    names()
  if (length(gene_col) == 0) {
    stop("Gene column '", gene, "' not found in the data frame")
  }
  gene_col <- gene_col[[1]]


  df2 <- df %>%
    dplyr::mutate(
      Expression = as.numeric(as.character(.data[[gene_col]])),
      Cluster_bin = ifelse(CIMIC_Cluster == "Fun-CIM", 1, 0)
    ) %>%
    dplyr::select(CIMIC_Cluster, Cluster_bin, Expression)


  cor_res <- cor.test(df2$Cluster_bin, df2$Expression,
                      method = "pearson")
  rtxt <- sprintf("r = %.3f\nP = %s",
                  cor_res$estimate,
                  format.pval(cor_res$p.value, digits = 4))


  ggplot(df2,
         aes(x = Cluster_bin,
             y = Expression,
             colour = CIMIC_Cluster)) +
    geom_jitter(width = 0.1,
                size = 2,
                alpha = 0.6) +
    geom_smooth(method = "lm",
                se = FALSE,
                colour = "black",
                linewidth = 1.2) +
    scale_x_continuous(breaks = c(0, 1),
                       labels = c("Dys-CIM", "Fun-CIM")) +
    scale_colour_manual(values = c("Dys-CIM" = "#C0392B",
                                   "Fun-CIM" = "#2980B9")) +
    annotate("text",
             x = 0.5,
             y = max(df2$Expression, na.rm = TRUE) * 1.1,
             label = rtxt,
             size = 5,
             fontface = "bold") +
    labs(title = paste0(gene),
         x = NULL,
         y = "Δ(log2[TPM+1])") +
    theme_classic(base_size = 14) +
    theme(
      plot.title      = element_text(hjust = 0.5,
                                     face = "bold",
                                     size = 18),
      legend.position = "none",
      axis.text.x     = element_text(face = "bold",
                                     size = 18),
      axis.text.y     = element_text(face = "bold",
                                     size = 18),
      axis.title.y    = element_text(face = "bold",
                                     size = 18),
      panel.border    = element_rect(colour = "black",
                                     fill = NA,
                                     linewidth = 1.5),
      axis.line       = element_line(colour = "black",
                                     linewidth = 1.5)
    )
}

genes_of_interest <- c("PDCD1", "CD274",  "PDCD1LG2",   # PD-L2
  "CTLA4",
  "LAG3",
  "HAVCR2",     # TIM-3
  "TIGIT",       # VISTA
  "BTLA")
pbs_plots <- purrr::map(genes_of_interest,
                    ~ compare_pointbiserial_plot(clustered_plot_df_for_merge, .x))











# Individual Gene three plot -------------------------



genes_of_interest <- c("PDCD1", "CD274", "IRF1")



# optional: combine the three panels for a single gene
make_three_panel <- function(gene) {
  ind  <- compare_induction_boxplot(clustered_plot_df_for_merge, gene)
  pen  <- compare_penetrance_plot(clustered_plot_df_for_merge, gene)
  pb   <- compare_pointbiserial_plot(clustered_plot_df_for_merge, gene)

  ggpubr::ggarrange(ind, pen, pb,
                    ncol = 3, nrow = 1,
                    common.legend = FALSE) %>%
    ggpubr::annotate_figure(
      top = ggpubr::text_grob(paste0(gene, " – Summary"),
                              face = "bold", size = 18)
    )
}

# Example: three‑panel view for PDCD1
make_three_panel("PDCD1")








# Individual Gene two plot -------------------






# optional: combine the three panels for a single gene
make_two_panel <- function(gene) {
  ## 1️⃣  Create the two plots (they still carry their titles)
  ind <- compare_induction_boxplot(clustered_plot_df_for_merge, gene)
  pen <- compare_penetrance_plot(clustered_plot_df_for_merge, gene)

  ## 2️⃣  Remove the titles from each plot
  ind <- ind + theme(plot.title = element_blank())
  pen <- pen + theme(plot.title = element_blank(), 
    axis.title.x = element_blank(), legend.position = "none")

  ## 3️⃣  Combine them side‑by‑side
  ggpubr::ggarrange(
    ind, pen,
    ncol = 2, nrow = 1,
    common.legend = FALSE          # keep the legend if you like it
  ) %>%
    ggpubr::annotate_figure(
      top = ggpubr::text_grob(
        paste0(gene),
        face = "bold", size = 18
      )
    )
}

# Example: 2‑panel view for PDCD1
# make_two_panel("PDCD1")

genes_of_interest <- c("PDCD1", "CD274", "Avg_MHC1", "CALR", "NLRC5", "STAT1")
two_panel_plots <- purrr::map(
  genes_of_interest,
  ~ make_two_panel(.x)      # .x is the gene name
)


classical_tcell_inhibitory_checkpoints <- c(
  "PDCD1",      # PD-1
  "CD274",      # PD-L1
  "PDCD1LG2",   # PD-L2
  "CTLA4",      
  "LAG3",
  "HAVCR2",     # TIM-3
  "TIGIT",
  "BTLA",
  "CD160",
  "CD244"
)


two_panel_plots <- purrr::map(
  classical_tcell_inhibitory_checkpoints,
  ~ make_two_panel(.x)      # .x is the gene name
)


classical_icd_genes <- c(
  "CALR",     # Calreticulin (ecto-CALR exposure signal)
  "HMGB1",    # High mobility group box 1 (DAMP released during ICD)
  "HSP90AA1", # Heat shock protein 90
  "HSP90AB1",
  "HSPA1A",   # HSP70
  "HSPA1B",
  "HSPA5",    # ER stress / unfolded protein response
  "ANXA1",    # Annexin A1 (DC recruitment signal)
  "P2RX7",    # ATP sensing receptor
  "PANX1",    # ATP release channel
  "TLR4",     # HMGB1 receptor on dendritic cells
  "EIF2AK3",  # PERK (ER stress signaling)
  "ATF4",     # Integrated stress response
  "DDIT3",    # CHOP
  "CXCL10",   # Type I IFN chemokine
  "IFNB1",     # Type I interferon
  "IFNAR1",
  "IFNAR2",
  "IFNA1",
  "IFNG",
  "IFNGR1",
  "IFNGR2"
)


two_panel_plots_icd <- purrr::map(
  classical_icd_genes,
  ~ make_two_panel(.x)      # .x is the gene name
)





















# Extra stuff for CIMIC Paper IC30, etc ---------


df <- df %>%
    mutate(
        Subtype_Group = case_when(
            str_detect(TNBC_Subtype, "BL") ~ "BL-like",
            TNBC_Subtype %in% c("M", "MSL") ~ "M-like",
            TRUE ~ "Other"
        )
    )

bar_df <- df %>%
    dplyr::count(CIMIC_Cluster, Subtype_Group) %>%
    group_by(CIMIC_Cluster) %>%
    mutate(
        prop = n / sum(n),
        label = paste0(round(prop * 100, 1), "%")
    )


tbl <- table(df$CIMIC_Cluster, df$Subtype_Group)

if(any(tbl < 5)) {
    test <- fisher.test(tbl)
    test_name <- "Fisher"
} else {
    test <- chisq.test(tbl)
    test_name <- "Chi-square"
}

pval <- signif(test$p.value, 3)


p_bar <- ggplot(bar_df, aes(x = CIMIC_Cluster, y = prop, fill = Subtype_Group)) +
    geom_bar(stat = "identity", color = "black") +
    geom_text(
        aes(label = label),
        position = position_stack(vjust = 0.5),
        size = 5,
        fontface = "bold"
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
        y = "Percentage",
        x = NULL,
        fill = "Subtype",
        title = paste0("Subtype Distribution (", test_name, " p = ", pval, ")")
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title      = element_text(hjust = 0.5,
                                       face = "bold",
                                       size = 16),
        legend.position = "none",
        axis.text.x     = element_text(face = "bold", size = 20),
        axis.text.y     = element_text(face = "bold", size = 20),
        axis.title.y    = element_text(face = "bold", size = 20),
        panel.border    = element_rect(colour = "black",
                                       fill = NA,
                                       linewidth = 1.5),
        axis.line       = element_line(colour = "black",
                                       linewidth = 1.5)
    )



p_box <- ggplot(df, aes(x = CIMIC_Cluster, y = IC30, fill = CIMIC_Cluster)) +
    geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 3) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    labs(
        y = "IC30 (µM)",
        x = NULL,
        title = "IC30 Comparison"
    ) +
     theme_classic(base_size = 14) +
    theme(
        plot.title      = element_text(hjust = 0.5,
                                       face = "bold",
                                       size = 16),
        legend.position = "none",
        axis.text.x     = element_text(face = "bold", size = 20),
        axis.text.y     = element_text(face = "bold", size = 20),
        axis.title.y    = element_text(face = "bold", size = 20),
        panel.border    = element_rect(colour = "black",
                                       fill = NA,
                                       linewidth = 1.5),
        axis.line       = element_line(colour = "black",
                                       linewidth = 1.5)
    )


p_bar
p_box



# PAM 50 Analysis -------------


library(data.table)
library(dplyr)
library(genefu)
library(org.Hs.eg.db)
library(AnnotationDbi)


# Pam 50 analysis

# Load PAM50 data
data(pam50)
data(pam50.robust)

  expr_mat <- combined_TPM_lognorm %>% column_to_rownames("gene_id") %>% as.matrix()

  # Transpose: samples as rows, genes as columns
  texp_mat <- t(expr_mat)
  texp_mat_numeric <- as.matrix(texp_mat)
  class(texp_mat_numeric) <- "numeric"
  if(any(is.na(texp_mat_numeric))) {
    print("Warning: NAs detected. Replacing with row means or 0.")
    texp_mat_numeric[is.na(texp_mat_numeric)] <- 0 
}
  
  
  gene_info <- data.frame(Genes = colnames(texp_mat))
  rownames(gene_info) <- gene_info$Genes
  


  # PAM50 classification
  pam50_predictions <- molecular.subtyping(
    sbt.model = "pam50",
    data = texp_mat_numeric,
    annot = gene_info,
    do.mapping = FALSE
  )
  
  # Build output
  out <- data.frame(
    Sample = rownames(texp_mat),
    PAM50_Subtype = pam50_predictions$subtype,
    stringsAsFactors = FALSE
  )
  
  
library(tidyverse)


pam50_df <- out %>%
  dplyr::filter(Sample != "gene_id") %>% 
  dplyr::select(Sample, PAM50_Subtype)


cluster_df <- clustered_plot_df %>%
  dplyr::select(sample_id, cluster_assignments) %>%
  dplyr::mutate(
    Sample = stringr::str_remove(sample_id, "^delta_"),
    
    CIMIC_Cluster = dplyr::case_when(
      cluster_assignments == 1 ~ "Fun-CIM",
      cluster_assignments == 2 ~ "Dys-CIM"
    )
  ) %>%
  dplyr::select(Sample, CIMIC_Cluster)


pam50_df <- pam50_df %>%
  mutate(
    Treatment = case_when(
      str_detect(Sample, "_DMSO") ~ "Baseline",
      str_detect(Sample, "_EPI")  ~ "Post-Treatment"
    )
  )

plot_df <- pam50_df %>%
  left_join(cluster_df, by = "Sample")

baseline_df_pam_50 <- plot_df %>%
  dplyr::filter(Treatment == "Baseline")

post_df_pam50 <- plot_df %>%
  dplyr::filter(Treatment == "Post-Treatment")

create_pam50_plot <- function(plot_df, condition_name) {
  
  # -----------------------------------
  # FORMAT DATA
  # -----------------------------------
  plot_data <- plot_df %>%
    dplyr::group_by(CIMIC_Cluster, PAM50_Subtype) %>%
    dplyr::summarise(total_n = n(), .groups = "drop") %>%
    dplyr::group_by(CIMIC_Cluster) %>%
    dplyr::mutate(
      total_samples = sum(total_n),
      percentage = (total_n / total_samples) * 100
    ) %>%
    dplyr::ungroup()
  
  # -----------------------------------
  # SAMPLE COUNTS
  # -----------------------------------
  sample_sizes <- plot_data %>%
    dplyr::group_by(CIMIC_Cluster) %>%
    dplyr::summarise(N = sum(total_n), .groups = "drop")
  
  labels_df <- data.frame(
    x = sample_sizes$CIMIC_Cluster,
    y = -8,
    label = paste0("N = ", sample_sizes$N)
  )
  
  # -----------------------------------
  # FISHER TEST
  # -----------------------------------
  contingency_table <- plot_data %>%
    dplyr::select(CIMIC_Cluster, PAM50_Subtype, total_n) %>%
    tidyr::pivot_wider(
      names_from = PAM50_Subtype,
      values_from = total_n,
      values_fill = 0
    )
  
  mat <- as.matrix(contingency_table[, -1])
  rownames(mat) <- contingency_table$CIMIC_Cluster
  
  fisher_p <- fisher.test(mat, simulate.p.value = TRUE)$p.value
  
  fish_stats_P <- paste0("P = ", signif(fisher_p, 3))
  
  # -----------------------------------
  # COLORS
  # -----------------------------------
  subtype_colors <- c(
    "Basal"  = "#C44E52",
    "Her2"   = "#DD8DA6",
    "LumA"   = "#0072B2",
    "LumB"   = "#5C8C45",
    "Normal" = "#9467BD"
  )
  
  # -----------------------------------
  # PLOT
  # -----------------------------------
  p <- ggplot(
    plot_data,
    aes(
      x = factor(CIMIC_Cluster, levels = c("Dys-CIM", "Fun-CIM")),
      y = percentage,
      fill = PAM50_Subtype
    )
  ) +
    
    geom_bar(
      stat = "identity",
      width = 0.85,
      color = "black",
      linewidth = 1
    ) +
    
    geom_text(
      aes(
        label = ifelse(
          percentage > 0,
          paste0(round(percentage, 1), "%"),
          ""
        )
      ),
      position = position_stack(vjust = 0.5),
      size = 6,
      fontface = "bold",
      color = "black"
    ) +
    
    geom_text(
      data = labels_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 6,
      fontface = "bold"
    ) +
    
    annotate(
      "text",
      x = 1.5,
      y = 110,
      label = fish_stats_P,
      size = 6,
      fontface = "bold"
    ) +
    
    scale_fill_manual(values = subtype_colors) +
    
    scale_y_continuous(
      limits = c(-15, 117),
      breaks = seq(0, 100, 20)
    ) +
    
    scale_x_discrete(
      expand = expansion(mult = c(0.55, 0.55))
    ) +
    
    labs(
      title = paste0("PAM50 Subtypes — ", condition_name),
      x = "",
      y = "Percentage",
      fill = ""
    ) +
    
    theme_classic(base_size = 18) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 22,
        face = "bold"
      ),
      
      axis.text = element_text(
        colour = "black",
        size = 18,
        face = "bold"
      ),
      
      axis.title = element_text(
        colour = "black",
        size = 20,
        face = "bold"
      ),
      
      axis.ticks = element_line(linewidth = 1.5),
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1.2
      ),
      
      legend.position = "right",
      
      legend.text = element_text(
        size = 16,
        face = "bold"
      ),
      
      legend.title = element_text(
        size = 16,
        face = "bold"
      )
    ) +
    
    guides(
      fill = guide_legend(
        title = NULL,
        keywidth = unit(1, "cm"),
        keyheight = unit(0.8, "cm")
      )
    )
  
  return(p)
}

create_pam50_plot(baseline_df_pam_50, "Pre-Treatment")
create_pam50_plot(post_df_pam50, "Post-Treatment")


# Viral Mimicry Z-score  ----------------


# Boxplot of Z-score



library(rstatix)
delta_combined_TPM_lognormnorm <- delta_subset_long


all_programs <- all_gene_sets

viral_mimicry_set <- c(

# dsDNA sensing
"CGAS","STING1","IFI16","AIM2","ZBP1",

# dsRNA sensing
"RIGI","IFIH1","DHX58","TLR3","TLR7","TLR8","TLR9",

# adaptor/signaling
"MAVS","TBK1","IRF3","IRF7",
"STAT1","STAT2",

# antiviral effectors
"IFIT1","IFIT2","IFIT3",
"OAS1","OAS2","OAS3","RNASEL",
"MX1","MX2",
"ISG15","RSAD2","BST2",
"IFI44","IFI44L","IFI6",
"EIF2AK2",

# antigen presentation
"B2M","HLA-A","HLA-B","HLA-C",
"TAP1","TAPBP","PSMB8","PSMB9","NLRC5",

# inflammatory outputs
"CXCL9","CXCL10","CXCL11","CCL5"

)



df <- clustered_plot_df
gene_names <- delta_combined_TPM_lognormnorm$gene_id[1:78932]

all_cols <- colnames(df)
viral_mimicry_set <-  intersect(viral_mimicry_set, all_cols)
print(viral_mimicry_set)

gene_cols <- gene_names 




zscore_safe <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(x)) else (x - m) / s
}

df_z <- df %>%
  mutate(across(all_of(gene_cols), zscore_safe))



score_program_z <- function(data, genes, program_name) {

  # ensure genes actually exist in data
  genes <- intersect(genes, colnames(data))

  if (length(genes) == 0) {
    return(
      data.frame(
        sample_id = data$sample_id,
        cimic_cluster = data$cluster_assignment,
        Program = paste0(program_name, "_z"),
        Score = rep(NA_real_, nrow(data))
      )
    )
  }
  
  score_vec <- rowMeans(as.matrix(data[, genes, drop = FALSE]), na.rm = TRUE)
  
  data.frame(
    sample_id = data$sample_id,
    cimic_cluster = data$cluster_assignment,
    Program = paste0(program_name, "_z"),
    Score = score_vec
  )
}

viral_mimicry_set_zscore <- score_program_z(
  data = df_z,
  genes = viral_mimicry_set,
  program_name = "Viral_Mimicry"
)




library(tidyr)
library(dplyr)

viral_mimicry_set_zscore_wide <- viral_mimicry_set_zscore %>%
  dplyr::select(sample_id, Program, Score) %>%
  tidyr::pivot_wider(
    names_from = Program,
    values_from = Score
  )

clustered_plot_df_for_merge <- clustered_plot_df_for_merge %>%
  left_join(viral_mimicry_set_zscore_wide, by = "sample_id")




viral_mimicry_plot_df <- viral_mimicry_set_zscore %>%
  
  mutate(
    cimic_cluster = dplyr::case_when(
      cimic_cluster == 1 ~ "Fun-CIM",
      cimic_cluster == 2 ~ "Dys-CIM"
    ),
    
    # clean labels FIRST (as character)
    Program = Program %>%
      str_remove_all("GOBP_|REACTOME_|HALLMARK_|_z") %>%
      str_replace_all("_", " "),
    
    # then convert to factor with correct ordering
    Program = factor(
      Program,
      levels = unique(
        paste0("Viral_Mimicry_z") %>%
          str_remove_all("GOBP_|REACTOME_|HALLMARK_|_z") %>%
          str_replace_all("_", " ")
      )
    ), 
    cimic_cluster = factor(cimic_cluster, levels = c("Dys-CIM", "Fun-CIM"))
)

pval_df_all <- viral_mimicry_plot_df %>%
  group_by(Program) %>%
  summarise(
    n_groups = n_distinct(cimic_cluster),
    
    y.position = max(Score, na.rm = TRUE) * 1.1,
    
    pval = tryCatch({
      if (n_groups == 2) {
        wilcox.test(Score ~ cimic_cluster)$p.value
      } else if (n_groups > 2) {
        kruskal.test(Score ~ cimic_cluster)$p.value
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_),
    
    p_label = ifelse(is.na(pval), "p = NA",
                     ifelse(pval < 0.001, "p < 0.001",
                            paste0("p = ", signif(pval, 3)))),
    
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(pval, method = "fdr"),
    
    fdr_label = ifelse(is.na(FDR), "FDR = NA",
                        ifelse(FDR < 0.001, "FDR < 0.001",
                               paste0("FDR = ", signif(FDR, 3))))
  ) %>% add_significance("FDR")


p_vm <- ggplot(viral_mimicry_plot_df, aes(x = cimic_cluster, y = Score, fill = cimic_cluster)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 1) +
  geom_text(
    data = pval_df_all,
    aes(x = 1.5, y = y.position, label = FDR.signif),
    inherit.aes = FALSE,
    size = 12,
    fontface = "bold"
  ) +
  facet_wrap(~ Program, scales = "free_y") +
  #scale_fill_manual(values = cimic_cluster_COLORS) +
  labs(
    y = "Z-score of\nProgram Activation",
    x = NULL
  ) +
    theme_bw(base_line_size = 2) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 32),
    axis.text = element_text(colour = "black", size = 24),
    axis.title.x = element_blank(),
    axis.text.x = element_text(colour = "black", size = 24, face = "bold"), 
    axis.text.y = element_text(colour = "black", size = 24, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 24, face = "bold"),
    axis.title = element_text(colour = "black", size = 24),
    legend.text = element_text(colour = "black", size = 24),
    legend.title = element_text(colour = "black", size = 24),
    legend.position = "none",
    #plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 2)
  )



p_vm


# Heatmap of individual Viral Mimi gemnes ------------


library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(matrixStats)
library(tibble)
library(rstatix)





viral_mimicry_gene_categories <- list(

# dsDNA sensing
"CGAS","STING1","IFI16","AIM2","ZBP1",

# dsRNA sensing
"RIGI","IFIH1","DHX58","TLR3","TLR7","TLR8","TLR9",

# adaptor/signaling
"MAVS","TBK1","IRF3","IRF7",
"STAT1","STAT2",

# antiviral effectors
"IFIT1","IFIT2","IFIT3",
"OAS1","OAS2","OAS3","RNASEL",
"MX1","MX2",
"ISG15","RSAD2","BST2",
"IFI44","IFI44L","IFI6",
"EIF2AK2",

# antigen presentation
"B2M","HLA-A","HLA-B","HLA-C",
"TAP1","TAPBP","PSMB8","PSMB9","NLRC5",

# inflammatory outputs
"CXCL9","CXCL10","CXCL11","CCL5"

)


genesets_of_interest <- unlist(viral_mimicry_gene_categories)

genesets_of_interest <- genesets_of_interest[genesets_of_interest %in% names(clustered_plot_df_for_merge)]


plot_df <- clustered_plot_df_for_merge %>%
  dplyr::select(all_of(c("CIMIC_Cluster", genesets_of_interest))) %>%
  pivot_longer(
    cols = -CIMIC_Cluster,
    names_to = "Gene",
    values_to = "Score"
  ) %>%
  mutate(
    Score = as.numeric(Score),
    CIMIC_Cluster = factor(CIMIC_Cluster, levels = c("Dys-CIM", "Fun-CIM")),
    Gene = Gene %>%
      str_replace_all("_", " ")
  ) %>%
  group_by(CIMIC_Cluster, Gene) %>%
  summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop")


pval_df <- clustered_plot_df_for_merge %>%
  dplyr::select(all_of(c("CIMIC_Cluster", genesets_of_interest))) %>%
  pivot_longer(-CIMIC_Cluster, names_to = "Gene", values_to = "Score") %>%
  group_by(Gene) %>%
  summarise(
    p = wilcox.test(Score ~ CIMIC_Cluster)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    signif = case_when(
      p.adj <= 0.0001 ~ "****",
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01 ~ "**",
      p.adj <= 0.05 ~ "*",
      TRUE ~ ""
    ),
    Gene = str_replace_all(Gene, "_", " ")
  )


program_order <- plot_df %>%
  distinct(Gene) %>%
  pull(Gene)

plot_df <- plot_df %>%
  mutate(Gene = factor(Gene, levels = rev(program_order)))

pval_df <- pval_df %>%
  mutate(Gene = factor(Gene, levels = rev(program_order)))

star_df <- plot_df %>%
  distinct(Gene) %>%
  left_join(pval_df, by = "Gene") %>%
  mutate(x = 2.6)

p <- ggplot(plot_df, aes(x = CIMIC_Cluster, y = Gene, fill = mean_score)) +
  
  geom_tile(color = "black", linewidth = 0.8) +
  
  annotate(
    "text",
    x = 2.65,
    y = seq_along(levels(plot_df$Gene)),
    label = rev(star_df$signif),
    size = 10,
    fontface = "bold"
  ) +
  
  scale_fill_gradient2(
    low = "#313695",
    mid = "white",
    high = "#A50026",
    midpoint = 0,
    name = "Mean\nZ-score\nΔ"
  ) +
  
  scale_x_discrete(
    limits = c("Dys-CIM", "Fun-CIM"),
    labels = c("Dys-CIM" = "Dys\nCIM", "Fun-CIM" = "Fun\nCIM"),
    expand = expansion(add = c(0.3, 1.2))
  ) +
  
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  
  theme_void(base_size = 18) +
  theme(
    axis.text.x = element_text(face = "bold", size = 24),
    axis.text.y = element_text(face = "bold", size = 16),
    legend.position = "right",
    legend.text = element_text(face = "bold", size = 20),
    plot.margin = margin(10, 20, 10, 20)
  )

p



```







# Baseline mhc1 induction by drug by cell line

```{r}


library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(ggpubr)

# ---- 1. Extract CellLine and Drug from sample_id ----
pre_treat_subset_long <- pre_treat_subset_wide %>%
  mutate(
    parts = strsplit(as.character(sample_id), "_"),
    CellLine = sapply(parts, `[`, 1),
    Drug     = sapply(parts, `[`, 2)
  ) %>%
  select(-parts)



# Combine with original delta_long
plot_df <- pre_treat_subset_long

# ---- 4. Plot box + jitter ----
ggplot(plot_df, aes(x = CellLine, y = avg_mhc1, fill = Drug, color = Drug)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 2) +
  theme_classic() +
  labs(
    y = "Average MHC-I Induction",
    x = "Cell Line",
    title = "Average MHC-I Across Cell Lines by Drug"
  ) +
      # Perform Global ANOVA to see if ANY cell line differs significantly
  stat_compare_means(method = "anova", label.y = max(plot_df$avg_mhc1) * 1.1) + 
  
  # Optional: Compare all cell lines against the global mean
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
  scale_fill_manual(values = c("DMSO" = "firebrick", "SALINE" = "steelblue", "Both" = "darkgreen")) +
  scale_color_manual(values = c("DMSO" = "firebrick", "SALINE" = "steelblue", "Both" = "darkgreen"))


```


# MHC - I Induction by drug by cell ----------

```{r}

library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(ggpubr)

# ---- 1. Extract CellLine and Drug from sample_id ----
delta_long <- delta_subset_wide %>%
  mutate(
    parts = strsplit(as.character(sample_id), "_"),
    CellLine = sapply(parts, `[`, 2),
    Drug     = sapply(parts, `[`, 3)
  ) %>%
  select(-parts)


# ---- 3. Create "Both" drug rows ----
both_df <- delta_long %>%
  group_by(CellLine) %>%
  summarise(
    avg_mhc1 = avg_mhc1, # keep all replicates
    sample_id = sample_id,
    Drug = "Both",
    .groups = "drop"
  )

# Combine with original delta_long
plot_df <- bind_rows(delta_long, both_df)

# ---- 4. Plot box + jitter ----
ggplot(plot_df, aes(x = CellLine, y = avg_mhc1, fill = Drug, color = Drug)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 2) +
  theme_classic() +
  labs(
    y = "Average MHC-I Induction",
    x = "Cell Line",
    title = "Average MHC-I Across Cell Lines by Drug"
  ) +
  scale_fill_manual(values = c("EPI" = "firebrick", "DOX" = "steelblue", "Both" = "darkgreen")) +
  scale_color_manual(values = c("EPI" = "firebrick", "DOX" = "steelblue", "Both" = "darkgreen"))

```


# Does MHC-I Induction Differ by DOX vs EPI in cell lines? 

```{r}
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(ggpubr)

# ---- 1. Extract CellLine and Drug from sample_id ----
delta_long <- delta_subset_wide %>%
  mutate(
    parts = strsplit(as.character(sample_id), "_"),
    CellLine = sapply(parts, `[`, 2),
    Drug     = sapply(parts, `[`, 3)
  ) %>%
  select(-parts)

# ---- 3. Create "Both" drug rows ----
both_df <- delta_long %>%
  group_by(CellLine) %>%
  summarise(
    avg_mhc1 = avg_mhc1, # keep all replicates
    sample_id = sample_id,
    Drug = "Both",
    .groups = "drop"
  )

# Combine with original delta_long
plot_df <- bind_rows(delta_long, both_df)

library(ggpubr)

# Define the specific drug pairs you want to compare
my_comparisons <- list( c("EPI", "DOX"))

ggplot(plot_df, aes(x = Drug, y = avg_mhc1, fill = Drug)) +
  # 1. Boxplot and Jitter
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, aes(color = Drug)) +
  
  # 2. Add Stats (Now it works because Drug is on the x-axis)
  stat_compare_means(
    comparisons = my_comparisons,
    label = "p.signif",
    method = "t.test",
    step.increase = 0.12 # Adds space between the brackets vertically
  ) +
  
  # 3. Separate by Cell Line
  facet_wrap(~CellLine) + 
  
  # 4. Styling
  theme_bw() + # Clean grid helps compare across facets
  theme(strip.background = element_rect(fill = "gray90"), # Labels for cell lines
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    y = "Average MHC-I Induction",
    title = "Drug Effect Significance per Cell Line"
  ) +
  scale_fill_manual(values = c("EPI" = "firebrick", "DOX" = "steelblue", "Both" = "darkgreen")) +
  scale_color_manual(values = c("EPI" = "firebrick", "DOX" = "steelblue", "Both" = "darkgreen"))


```

# Does MHC-I Induction Differ by cell line when treated with DOX? 

```{r}
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)

# ---- 1. Data Preparation ----
# Filter for DOX only, since the question is about differences in DOX response
DOX_data <- delta_subset_wide %>%
  mutate(
    parts = strsplit(as.character(sample_id), "_"),
    CellLine = sapply(parts, `[`, 2),
    Drug     = sapply(parts, `[`, 3)
  ) %>%
  filter(Drug == "DOX") # Focus strictly on DOXplatin

# ---- 3. Visualization & Statistics ----
# We move CellLine to the X-axis to compare them directly
ggplot(DOX_data, aes(x = CellLine, y = avg_mhc1, fill = CellLine)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  
  # Perform Global ANOVA to see if ANY cell line differs significantly
  stat_compare_means(method = "anova", label.y = max(DOX_data$avg_mhc1) * 1.1) + 
  
  # Optional: Compare all cell lines against the global mean
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
  
  theme_minimal() +
  labs(
    title = "MHC-I Induction by Cell Line (DOXplatin Treatment)",
    subtitle = "Comparing average HLA-A/B/C expression across cell lines",
    y = "Average MHC-I Expression (Delta)",
    x = "Cell Line"
  ) +
  theme(legend.position = "none")

```



# Does MHC-I Induction Differ by cell line when treated with EPI? 

```{r}
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)

# ---- 1. Data Preparation ----
# Filter for EPI only, since the question is about differences in EPI response
EPI_data <- delta_subset_wide %>%
  mutate(
    parts = strsplit(as.character(sample_id), "_"),
    CellLine = sapply(parts, `[`, 2),
    Drug     = sapply(parts, `[`, 3)
  ) %>%
  filter(Drug == "EPI") # Focus strictly on EPIplatin


# ---- 3. Visualization & Statistics ----
# We move CellLine to the X-axis to compare them directly
ggplot(EPI_data, aes(x = CellLine, y = avg_mhc1, fill = CellLine)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  
  # Perform Global ANOVA to see if ANY cell line differs significantly
  stat_compare_means(method = "anova", label.y = max(EPI_data$avg_mhc1) * 1.1) + 
  
  # Optional: Compare all cell lines against the global mean
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
  
  theme_minimal() +
  labs(
    title = "MHC-I Induction by Cell Line (EPIplatin Treatment)",
    subtitle = "Comparing average HLA-A/B/C expression across cell lines",
    y = "Average MHC-I Expression (Delta)",
    x = "Cell Line"
  ) +
  theme(legend.position = "none")

```





# PD-L1 Induction by drug by cell ----------


```{r}

library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)

# ---- 1. Extract CellLine and Drug from sample_id ----
delta_long <- delta_subset_wide %>%
  mutate(
    parts = strsplit(as.character(sample_id), "_"),
    CellLine = sapply(parts, `[`, 2),
    Drug     = sapply(parts, `[`, 3)
  ) %>%
  select(-parts)

# ---- 2. Compute avg_pdl1 ----
pdl1_genes <- c("CD274")
pdl1_genes <- pdl1_genes[pdl1_genes %in% colnames(delta_long)]

delta_long <- delta_long %>%
  rowwise() %>%
  mutate(avg_pdl1 = mean(c_across(all_of(pdl1_genes)), na.rm = TRUE)) %>%
  ungroup()

# ---- 3. Create "Both" drug rows ----
both_df <- delta_long %>%
  group_by(CellLine) %>%
  summarise(
    avg_pdl1 = avg_pdl1, # keep all replicates
    sample_id = sample_id,
    Drug = "Both",
    .groups = "drop"
  )

# Combine with original delta_long
plot_df <- bind_rows(delta_long, both_df)

# ---- 4. Plot box + jitter ----
ggplot(plot_df, aes(x = CellLine, y = avg_pdl1, fill = Drug, color = Drug)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 2) +
  theme_classic() +
  labs(
    y = "Average PD-L1 Induction",
    x = "Cell Line",
    title = "Average PD-L1 Induction Cell Lines by Drug"
  ) +
  scale_fill_manual(values = c("EPI" = "firebrick", "DOX" = "steelblue", "Both" = "darkgreen")) +
  scale_color_manual(values = c("EPI" = "firebrick", "DOX" = "steelblue", "Both" = "darkgreen"))

```




# Does PD-L1 Induction Differ by cell line when treated with DOX? 

```{r}
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)

# ---- 1. Data Preparation ----
# Filter for DOX only, since the question is about differences in DOX response
DOX_data <- delta_subset_wide %>%
  mutate(
    parts = strsplit(as.character(sample_id), "_"),
    CellLine = sapply(parts, `[`, 2),
    Drug     = sapply(parts, `[`, 3)
  ) %>%
  filter(Drug == "DOX") # Focus strictly on DOXplatin

# ---- 2. Compute avg_pdl1 ----
pdl1_genes <- c("CD274")
pdl1_genes <- pdl1_genes[pdl1_genes %in% colnames(DOX_data)]

DOX_data <- DOX_data %>%
  rowwise() %>%
  mutate(avg_pdl1 = mean(c_across(all_of(pdl1_genes)), na.rm = TRUE)) %>%
  ungroup()

# ---- 3. Visualization & Statistics ----
# We move CellLine to the X-axis to compare them directly
ggplot(DOX_data, aes(x = CellLine, y = avg_pdl1, fill = CellLine)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  
  # Perform Global ANOVA to see if ANY cell line differs significantly
  stat_compare_means(method = "anova", label.y = max(DOX_data$avg_pdl1) * 1.1) + 
  
  # Optional: Compare all cell lines against the global mean
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
  
  theme_minimal() +
  labs(
    title = "PD-L1 Induction by Cell Line (DOXplatin Treatment)",
    y = "Average PD-L1 Expression (Delta)",
    x = "Cell Line"
  ) +
  theme(legend.position = "none")

```



# Does PD-L1 Induction Differ by cell line when treated with EPI? 

```{r}
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)

# ---- 1. Data Preparation ----
# Filter for EPI only, since the question is about differences in EPI response
EPI_data <- delta_subset_wide %>%
  mutate(
    parts = strsplit(as.character(sample_id), "_"),
    CellLine = sapply(parts, `[`, 2),
    Drug     = sapply(parts, `[`, 3)
  ) %>%
  filter(Drug == "EPI") # Focus strictly on EPIplatin

# ---- 2. Compute avg_pdl1 ----
pdl1_genes <- c("CD274")
pdl1_genes <- pdl1_genes[pdl1_genes %in% colnames(EPI_data)]

EPI_data <- EPI_data %>%
  rowwise() %>%
  mutate(avg_pdl1 = mean(c_across(all_of(pdl1_genes)), na.rm = TRUE)) %>%
  ungroup()

# ---- 3. Visualization & Statistics ----
# We move CellLine to the X-axis to compare them directly
ggplot(EPI_data, aes(x = CellLine, y = avg_pdl1, fill = CellLine)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  
  # Perform Global ANOVA to see if ANY cell line differs significantly
  stat_compare_means(method = "anova", label.y = max(EPI_data$avg_pdl1) * 1.1) + 
  
  # Optional: Compare all cell lines against the global mean
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
  
  theme_minimal() +
  labs(
    title = "PD-L1 Induction by Cell Line (EPIplatin Treatment)",
    y = "Average PD-L1 Expression (Delta)",
    x = "Cell Line"
  ) +
  theme(legend.position = "none")

```



# Heatmap of key CIM programs by cell Line for DOX 

```{r}

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(tibble)

# -------------------------------------------------
# 0. Gene‑set definition (unchanged)
# -------------------------------------------------
gene_categories <- list(
  "Summary"          = c("avg_mhc1"),
  "MHC-I"            = c("HLA-A", "HLA-B", "HLA-C"),
  "PD-L1"            = c("CD274"),
  "IFN-1/3"          = c("IFNA1", "IFNA2", "IFNB1", "IFNL1", "IFNL2", "IFNL3"),
  "STING Pathway"    = c("CGAS", "STING1", "IRF3"),
  "Chemokines (Pro)" = c("CXCL9", "CXCL10", "CXCL13", "CXCL14", "CXCL6", "CCL5"),
  "Chemokines (Supp)"= c("IL6", "CXCL8", "CCL2", "CCL22", "CCL28"),
  "Pro-inflammatory" = c("TNF", "IL1B", "IL6", "IL12A", "IFNG"),
  "ICD/Cytokines"    = c("CALR", "HMGB1", "TNF", "IFNG"),
  "Ferroptosis (Pro)"= c("ACSL4","LPCAT3","PTGS2","ALOX5","ALOX15","CHAC1"),
  "Ferroptosis (Anti)"=c("GPX4","SLC7A11","AIFM2","GSS","GCLC"),
  "ATP Release/Sense"=c("ATP5F1","PANX1","VNUT","P2RX7","P2RY2")
)

target_genes <- unlist(gene_categories)

# -------------------------------------------------
# 1. Prepare DOX‑only data 
# -------------------------------------------------
DOX_samples <- delta_subset_wide %>%
  mutate(
    parts    = strsplit(as.character(sample_id), "_"),
    CellLine = sapply(parts, `[`, 2),
    Drug     = sapply(parts, `[`, 3)
  ) %>%
  filter(Drug == "DOX") %>%
  rowwise() %>%
  ungroup()

# -------------------------------------------------
# 2. Build a numeric matrix (genes × samples)
# -------------------------------------------------
mat <- DOX_samples %>%
  select(sample_id, any_of(intersect(target_genes, colnames(.)))) %>%   # keep only genes we care about
  column_to_rownames("sample_id") %>%       # rows = samples
  as.matrix() %>%
  t()                                        # rows = genes, columns = samples

# -------------------------------------------------
# 3. Coerce to numeric and drop rows that are completely NA
# -------------------------------------------------
mat_rownames <- rownames(mat)
mat <- apply(mat, 2, as.numeric)            # force numeric
rownames(mat) <- mat_rownames


keep_rows <- rowSums(!is.na(mat)) > 0        # keep genes that have at least one non‑NA value
if (!any(keep_rows)) {
  stop("No genes with non‑missing values remain after conversion. ",
       "Check the gene list and the DOX subset.")
}
mat <- mat[keep_rows, , drop = FALSE]

# -------------------------------------------------
# 4. Build row_info AFTER the matrix is final
# -------------------------------------------------
row_info <- stack(gene_categories) %>%
  rename(gene = values, category = ind) %>%
  filter(gene %in% rownames(mat))

# Re‑order the matrix to follow the order in row_info
mat <- mat[row_info$gene, , drop = FALSE]

# -------------------------------------------------
# 5. Colour function (robust 99 % limit)
# -------------------------------------------------
limit   <- quantile(abs(mat), probs = 0.99, na.rm = TRUE)
col_fun <- colorRamp2(c(-limit, 0, limit), c("blue", "white", "red"))

# -------------------------------------------------
# 6. Column annotation (CellLine)
# -------------------------------------------------

# 6. Column annotation (CellLine) – corrected color assignment
celllines   <- unique(DOX_samples$CellLine)
n_celllines <- length(celllines)

# use a palette that can produce at least n_celllines colors
col_vec <- if (n_celllines <= 12) {
  RColorBrewer::brewer.pal(n = n_celllines, name = "Set3")
} else {
  # for >12 groups fall back to a continuous palette
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_celllines)
}

col_ha <- HeatmapAnnotation(
  CellLine = DOX_samples$CellLine,
  col = list(CellLine = setNames(
    RColorBrewer::brewer.pal(min(8, length(unique(DOX_samples$CellLine))), "Set3"),
    unique(DOX_samples$CellLine)
  )),
  show_annotation_name = TRUE
)

# -------------------------------------------------
# 7. Draw the heat‑map – corrected `row_names_gp`
# -------------------------------------------------
Heatmap(
  mat,
  name            = "Delta Exp",
  col             = col_fun,
  row_split       = row_info$category,          # now length == nrow(mat)
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  top_annotation  = col_ha,
  rect_gp         = gpar(col = "white", lwd = 0.5),
  row_names_gp    = gpar(
    fontsize = 14,
    fontface = if ("avg_mhc1" %in% rownames(mat)) "bold" else "plain"
  ),
  row_title_rot   = 0,
  column_names_gp = gpar(fontsize = 14),
  column_title    = "DOXplatin‑Induced Gene Expression Changes per Replicate"
)


```


# Heatmap of key CIM programs by cell Line for EPI 

```{r}
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(tibble)

# -------------------------------------------------
# 0. Updated Gene Categories
# -------------------------------------------------
gene_categories <- list(
  "Summary"          = c("avg_mhc1"), # Our calculated row
  "MHC-I"            = c("HLA-A", "HLA-B", "HLA-C"),
  "PD-L1"            = c("CD274"),
  "IFN-1/3"          = c("IFNA1", "IFNA2", "IFNB1", "IFNL1", "IFNL2", "IFNL3"),
  "STING Pathway"    = c("CGAS", "STING1", "IRF3"),
  "Chemokines (Pro)" = c("CXCL9", "CXCL10", "CXCL13", "CXCL14", "CXCL6", "CCL5"),
  "Chemokines (Supp)"= c("IL6", "CXCL8", "CCL2", "CCL22", "CCL28"),
  "Pro-inflammatory"  = c("TNF", "IL1B", "IL6", "IL12A", "IFNG"),
  "ICD/Cytokines"    = c("CALR", "HMGB1", "TNF", "IFNG"),
  "Ferroptosis (Pro)"= c("ACSL4", "LPCAT3", "PTGS2", "ALOX5", "ALOX15", "CHAC1"),
  "Ferroptosis (Anti)"= c("GPX4", "SLC7A11", "AIFM2", "GSS", "GCLC"),
  "ATP Release/Sense"= c("ATP5F1", "PANX1", "VNUT", "P2RX7", "P2RY2")
)

target_genes <- unlist(gene_categories)

# -------------------------------------------------
# 1. Prepare EPI-only data
# -------------------------------------------------
EPI_samples <- delta_subset_wide %>%
  mutate(
    parts    = strsplit(as.character(sample_id), "_"),
    CellLine = sapply(parts, `[`, 2),
    Drug     = sapply(parts, `[`, 3)
  ) %>% 
  filter(Drug == "EPI") %>%
  # --- NEW: Calculate Individual Sample avg_mhc1 ---
  rowwise() %>%
  ungroup() %>%
  # ------------------------------------------------
  select(sample_id, CellLine, any_of(intersect(target_genes, colnames(.))))

# -------------------------------------------------
# 2. Build numeric matrix
# -------------------------------------------------
mat <- EPI_samples %>% 
  column_to_rownames("sample_id") %>% 
  select(-CellLine) %>% 
  as.matrix() %>% 
  t() 

# -------------------------------------------------
# 3. Align gene names with categories
# -------------------------------------------------
row_info <- stack(gene_categories) %>% 
  rename(gene = values, category = ind) %>% 
  filter(gene %in% rownames(mat))

# Re-order matrix and handle "avg_mhc1" bolding later
mat <- mat[row_info$gene, , drop = FALSE]

# -------------------------------------------------
# 4. Color function & Annotation
# -------------------------------------------------
# Use a fixed symmetric limit for easier comparison between plots
limit <- quantile(abs(mat), 0.99, na.rm = TRUE) # Robust limit
col_fun <- colorRamp2(c(-limit, 0, limit), c("blue", "white", "red"))

col_ha <- HeatmapAnnotation(
  CellLine = EPI_samples$CellLine,
  col = list(CellLine = setNames(
    RColorBrewer::brewer.pal(min(8, length(unique(EPI_samples$CellLine))), "Set3"),
    unique(EPI_samples$CellLine)
  )),
  show_annotation_name = TRUE
)

# -------------------------------------------------
# 5. Draw the Heatmap
# -------------------------------------------------
Heatmap(
  mat,
  name            = "Delta Exp",
  col             = col_fun,
  row_split       = row_info$category,
  cluster_rows    = FALSE,
  cluster_columns = FALSE, # Clusters similar responding replicates
  top_annotation  = col_ha,
  
  # Styling individual boxes
  rect_gp         = gpar(col = "white", lwd = 0.5),
  
  # Row Labeling (Make avg_mhc1 Bold)
  row_names_gp = gpar(
    fontsize = 14, 
    fontface = ifelse(rownames(mat) == "avg_mhc1", "bold", "plain")
  ),
  
  row_title_rot   = 0,
  column_names_gp = gpar(fontsize = 14),
  column_title    = "EPIplatin-Induced Gene Expression Changes per Replicate"
)


```


# Heatmap for CIM Mediatiors by cell line for DOX 
```{r}
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(RColorBrewer)

# 0. Gene‑set definition (unchanged) ……………………………………… (omitted)

# 1. Prepare DOX‑only data 
DOX_samples <- delta_subset_wide %>%
  mutate(
    parts    = strsplit(as.character(sample_id), "_"),
    CellLine = sapply(parts, `[`, 2),
    Drug     = sapply(parts, `[`, 3)
  ) %>%
  filter(Drug == "DOX") %>%
  rowwise() %>%
  ungroup() %>%
  select(sample_id, CellLine, any_of(intersect(target_genes, colnames(.))))

# 2. Build numeric matrix (genes × samples)
mat <- DOX_samples %>%
  column_to_rownames("sample_id") %>%   # rows = samples
  select(-CellLine) %>%                 # drop the annotation column
  as.matrix() %>%
  t()                                   # now rows = genes

# 3. Align rows with gene‑set categories
row_info <- stack(gene_categories) %>%
  rename(gene = values, category = ind) %>%
  filter(gene %in% rownames(mat))

mat <- mat[row_info$gene, , drop = FALSE]

# 4. Colour function (robust 99 % limit)
limit   <- quantile(abs(mat), probs = 0.99, na.rm = TRUE)
col_fun <- colorRamp2(c(-limit, 0, limit), c("blue", "white", "red"))

# 5. Column annotation – **fixed colour mapping**
cell_lines      <- unique(DOX_samples$CellLine)
cell_line_cols  <- setNames(
  brewer.pal(n = length(cell_lines), name = "Set3"),
  cell_lines
)

col_ha <- HeatmapAnnotation(
  CellLine = DOX_samples$CellLine,
  col = list(CellLine = cell_line_cols),
  show_annotation_name = TRUE
)

# 6. Draw the heat‑map
Heatmap(
  mat,
  name            = "Delta Exp",
  col             = col_fun,
  row_split       = row_info$category,
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  top_annotation  = col_ha,
  rect_gp         = gpar(col = "white", lwd = 0.5),
  row_names_gp    = gpar(
    fontsize = 14,
    fontface = if ("avg_mhc1" %in% rownames(mat)) "bold" else "plain"
  ),
  row_title_rot   = 0,
  column_names_gp = gpar(fontsize = 14),
  column_title    = "DOXplatin‑Induced Gene Expression Changes per Replicate"
)



```



# Get all Gene Sets of Interest
```{r}

library(msigdbr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)


msig_df <- msigdbr::msigdbr(species = "Homo sapiens") %>% as.data.frame()

    # Default CIM Genesets from manuscript
    all_gene_set_names <- c(
"GOBP_TRANSLATIONAL_INITIATION", 
    "GOBP_RIBOSOME_BIOGENESIS", 
    "GOBP_PROTEIN_FOLDING",
    "REACTOME_RRNA_PROCESSING",
    "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
    "GOBP_PRODUCTION_OF_MOLECULAR_MEDIATOR_OF_IMMUNE_RESPONSE",
    "GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
    "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
    "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
    "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
    "GOBP_CELLULAR_RESPONSE_TO_STRESS",
    "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
    "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
    "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
    "GOBP_T_CELL_ACTIVATION",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS",
    "WP_OXIDATIVE_STRESS_RESPONSE",
    "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
    "GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS",
    "GOBP_LEUKOCYTE_CHEMOTAXIS_INVOLVED_IN_INFLAMMATORY_RESPONSE",
    "GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "GOBP_B_CELL_ACTIVATION"
    )
    

all_gene_set_names  <- c(all_gene_set_names, unique(msig_df$gs_name[startsWith(msig_df$gs_name, "HALLMARK")]))

  # Build a named list: each element is the vector of gene symbols for that set
    all_gene_sets <- setNames(
      lapply(all_gene_set_names, function(gs) {
        gs_df <- msig_df[msig_df$gs_name == gs, , drop = FALSE]
        gs_df$gene_symbol
      }),
      all_gene_set_names

    )
  

  # all_gene_sets is a named list:
#   names(all_gene_sets) → gene‑set IDs (e.g. "GOBP_ADAPTIVE_IMMUNE_RESPONSE")
#   each element          → character vector of gene symbols
# -------------------------------------------------
# 2) Convert the list → long format (one row per gene‑set / gene pair)
# -------------------------------------------------
# 2) Convert the list → long format (one row per gene‑set / gene pair)
gene_set_long <- enframe(
  all_gene_sets,               # list element → vector of genes
  name = "gene_set_ID",
  value = "genes"
) %>% 
  # 2) Expand the list‑column; give the column name directly
  unnest_longer(genes) %>%     # <- correct call, no `cols=` argument
  # 3) Rename the expanded column to the name you use later
  rename(gene = genes)

# -------------------------------------------------
# 3) Add a position index inside each gene‑set (1,2,3,…)
# -------------------------------------------------
gene_set_long <- gene_set_long %>% 
  group_by(gene_set_ID) %>% 
  mutate(pos = row_number()) %>%    # position of the gene inside the set
  ungroup()

# -------------------------------------------------
# 4) Wide‑reshape: one row per gene‑set, columns ID1, ID2, …
# -------------------------------------------------
gene_set_wide <- gene_set_long %>% 
  pivot_wider(names_from   = pos,
              values_from  = gene,
              names_prefix = "V") %>%      # → ID1, ID2, …
  # optional: keep column order tidy
  select(gene_set_ID, starts_with("V"))




cimic_selected_gene_sets <- all_gene_sets[c(
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
)]


  # all_gene_sets is a named list:
#   names(all_gene_sets) → gene‑set IDs (e.g. "GOBP_ADAPTIVE_IMMUNE_RESPONSE")
#   each element          → character vector of gene symbols
# -------------------------------------------------
# 2) Convert the list → long format (one row per gene‑set / gene pair)
# -------------------------------------------------
# 2) Convert the list → long format (one row per gene‑set / gene pair)
cimic_gene_set_long <- enframe(
  cimic_selected_gene_sets,               # list element → vector of genes
  name = "gene_set_ID",
  value = "genes"
) %>% 
  # 2) Expand the list‑column; give the column name directly
  unnest_longer(genes) %>%     # <- correct call, no `cols=` argument
  # 3) Rename the expanded column to the name you use later
  rename(gene = genes)

# -------------------------------------------------
# 3) Add a position index inside each gene‑set (1,2,3,…)
# -------------------------------------------------
cimic_gene_set_long <- cimic_gene_set_long %>% 
  group_by(gene_set_ID) %>% 
  mutate(pos = row_number()) %>%    # position of the gene inside the set
  ungroup()

# -------------------------------------------------
# 4) Wide‑reshape: one row per gene‑set, columns ID1, ID2, …
# -------------------------------------------------
cimic_gene_set_wide <- cimic_gene_set_long %>% 
  pivot_wider(names_from   = pos,
              values_from  = gene,
              names_prefix = "V") %>%      # → ID1, ID2, …
  # optional: keep column order tidy
  select(gene_set_ID, starts_with("V"))



```


# Make ssGSEA Function

```{r}

# ssGSEA Function
ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
    row_names = rownames(X)
    num_genes = nrow(X)
    gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})

    # Ranks for genes
    R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')

    # Calculate enrichment score (es) for each sample (column)
    es = apply(R, 2, function(R_col) {
        gene_ranks = order(R_col, decreasing = TRUE)

        # Calc es for each gene set
        es_sample = sapply(gene_sets, function(gene_set_idx) {
            # pos: match (within the gene set)
            # neg: non-match (outside the gene set)
            indicator_pos = gene_ranks %in% gene_set_idx
            indicator_neg = !indicator_pos

            rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha

            step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
            step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)

            step_cdf_diff = step_cdf_pos - step_cdf_neg

            # Normalize by gene number
            if (scale) step_cdf_diff = step_cdf_diff / num_genes

            # Use ssGSEA or not
            if (single) {
                sum(step_cdf_diff)
            } else {
                step_cdf_diff[which.max(abs(step_cdf_diff))]
            }
        })
        unlist(es_sample)
    })

    if (length(gene_sets) == 1) es = matrix(es, nrow = 1)

    # Normalize by absolute diff between max and min
    if (norm) es = es / diff(range(es))

    # Prepare output
    rownames(es) = names(gene_sets)
    colnames(es) = colnames(X)
    return(es)
}




```



# Make Deltas for Gene programs using ssGSEA for DOX and EPI

```{r}

library(dplyr)
library(tibble)
library(GSVA)  # for ssgsea
library(tidyr)

# ---- 1. Convert to numeric matrix ----
tnbc_mat <- tnbc_raw %>%
  remove_rownames() %>%
  column_to_rownames("gene_id") %>%
  as.matrix() %>%
  { matrix(as.numeric(.), nrow = nrow(.), dimnames = dimnames(.)) }

# ---- 2. Function to extract columns for a given drug / control ----
get_drug_matrix <- function(mat, drug_name, control_name = "DMSO") {
  
  drug_cols    <- grep(paste0("_", drug_name, "_"), colnames(mat), value = TRUE)
  control_cols <- grep(paste0("_", control_name, "_"), colnames(mat), value = TRUE)
  
  list(
    drug    = mat[, drug_cols, drop = FALSE],
    control = mat[, control_cols, drop = FALSE]
  )
}


# ---- 3. Function to calculate delta ssGSEA ----
compute_delta_ssgsea <- function(drug_mat, control_mat, gene_sets, scale = TRUE, norm = FALSE) {
  
  # ssGSEA for drug
  es_drug <- ssgsea(drug_mat, gene_sets, scale = scale, norm = norm)
  
  # ssGSEA for control
  es_control <- ssgsea(control_mat, gene_sets, scale = scale, norm = norm)
  
  # delta
  es_delta <- es_drug - es_control
  
  # optional: row-wise z-score
  es_delta_zrow <- t(apply(es_delta, 1, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
  
  list(
    delta_raw = as.data.frame(es_delta),
    delta_zrow = as.data.frame(es_delta_zrow),
    drug_es = as.data.frame(es_drug),
    control_es = as.data.frame(es_control)
  )
}

# ---- 4. Example: EPI vs DMSO ----
EPI_data <- get_drug_matrix(tnbc_mat, drug_name = "EPI", control_name = "DMSO")
ssgsea_EPI_delta <- compute_delta_ssgsea(EPI_data$drug, EPI_data$control, all_gene_sets)

# Access results
delta_EPI_raw   <- ssgsea_EPI_delta$delta_raw
delta_EPI_zrow  <- ssgsea_EPI_delta$delta_zrow
EPI_es_drug     <- ssgsea_EPI_delta$drug_es
EPI_es_control  <- ssgsea_EPI_delta$control_es



# ---- 5. Example: DOX vs DMSO ----
DOX_data <- get_drug_matrix(tnbc_mat, drug_name = "DOX", control_name = "DMSO")
ssgsea_DOX_delta <- compute_delta_ssgsea(DOX_data$drug, DOX_data$control, all_gene_sets)


# Access results
delta_DOX_raw   <- ssgsea_DOX_delta$delta_raw
delta_DOX_zrow  <- ssgsea_DOX_delta$delta_zrow
DOX_es_drug     <- ssgsea_DOX_delta$drug_es
DOX_es_control  <- ssgsea_DOX_delta$control_es


```


# Box plots of gene_sets 

```{r}
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(tibble)

# -------------------------------------------------
# 1. Prepare Tidy Data from delta_DOX_zrow
# -------------------------------------------------
# We extract metadata directly from the column names of your provided dataframe
plot_df_DOX <- delta_DOX_zrow %>% 
  rownames_to_column("GeneSet") %>% 
  pivot_longer(
    cols      = -GeneSet,
    names_to  = "sample_id",
    values_to = "Delta"
  ) %>%
  mutate(
    # Extracts "BT549" from "BT549_DOX_R1"
    CellLine = sapply(strsplit(sample_id, "_"), `[`, 1),
    # Extracts "DOX" from "BT549_DOX_R1"
    Drug     = sapply(strsplit(sample_id, "_"), `[`, 2)
  ) %>%
  filter(Drug == "DOX")

# Define the list of gene sets to loop over if not already defined
all_gene_set_names <- unique(plot_df_DOX$GeneSet)

# -------------------------------------------------
# 2. Loop over GeneSets for Plots and Stats
# -------------------------------------------------
all_plots  <- list()
stats_list <- list()

for (gs in all_gene_set_names) {
  current_data <- plot_df_DOX %>% filter(GeneSet == gs)

  if (nrow(current_data) < 2) next

  # --- Plot with Wilcox and Kruskal ---
  p <- ggplot(current_data, aes(x = CellLine, y = Delta, fill = CellLine)) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2) +
    # Global non-parametric test
    stat_compare_means(method = "kruskal.test", label.y.npc = "top") +
    # Comparison of each cell line against the grand median
    stat_compare_means(method = "wilcox.test", 
                       ref.group = ".all.", 
                       label = "p.signif") +
    theme_bw() +
    labs(title = gs, y = "Delta ssGSEA", x = "Cell Line") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

  all_plots[[gs]] <- p

  # --- Statistics Collection ---
  gs_stats <- tryCatch({
      compare_means(Delta ~ CellLine, 
                    data = current_data, 
                    method = "wilcox.test", 
                    ref.group = ".all.") %>%
        mutate(GeneSet = gs)
    }, error = function(e) {
      tibble(GeneSet = gs, p = NA_real_, method = "wilcox.test", group1 = ".all.", group2 = NA_character_)
    })
  
  stats_list[[gs]] <- gs_stats
}

# -------------------------------------------------
# 3. Combine and Adjust P-values
# -------------------------------------------------
all_stats_df <- bind_rows(stats_list) %>%
  # Ensure p is numeric for adjustment
  mutate(p = as.numeric(p)) %>% 
  # Apply Benjamini-Hochberg adjustment across all tests performed
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  # Re-assign significance stars based on adjusted p-values
  mutate(p.adj.signif = case_when(
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01  ~ "**",
    p.adj <= 0.05  ~ "*",
    TRUE           ~ "ns"
  )) %>%
  select(GeneSet, group2, p, p.adj, p.adj.signif, method) %>%
  rename(CellLine = group2) %>%
  arrange(p.adj)

# -------------------------------------------------
# 4. Final Review
# -------------------------------------------------
print(head(all_stats_df))

# Example: Display the plot for Immune Response mediators
# print(all_plots[["GOBP_PRODUCTION_OF_MOLECULAR_MEDIATOR_OF_IMMUNE_RESPONSE"]])
```

# Heatmap of Gene sets 

```{r}

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

# 1. Define Categorized List
gene_categories <- list(
  "MHC-I / Antigen Presentation" = c(
    "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
    "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
    "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"
  ),
  "Immune Activation" = c(
    "GOBP_ADAPTIVE_IMMUNE_RESPONSE", "GOBP_T_CELL_ACTIVATION", "GOBP_B_CELL_ACTIVATION",
    "GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL", "GOBP_PRODUCTION_OF_MOLECULAR_MEDIATOR_OF_IMMUNE_RESPONSE",
    "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_COMPLEMENT"
  ),
  "Inflammatory Signaling" = c(
    "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
    "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
    "GOBP_LEUKOCYTE_CHEMOTAXIS_INVOLVED_IN_INFLAMMATORY_RESPONSE"
  ),
  "Translation & Proteostasis" = c(
    "GOBP_TRANSLATIONAL_INITIATION", "GOBP_RIBOSOME_BIOGENESIS", "REACTOME_RRNA_PROCESSING",
    "GOBP_PROTEIN_FOLDING", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_PROTEIN_SECRETION"
  ),
  "Stress & Apoptosis" = c(
    "GOBP_CELLULAR_RESPONSE_TO_STRESS", "WP_OXIDATIVE_STRESS_RESPONSE",
    "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY", "HALLMARK_APOPTOSIS", 
    "HALLMARK_P53_PATHWAY", "HALLMARK_DNA_REPAIR"
  ),
  "Proliferation" = c(
    "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MITOTIC_SPINDLE", 
    "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MTORC1_SIGNALING"
  )
)

# 2. Filter Matrix to these selected groups
target_gs <- unlist(gene_categories)
heatmap_mat <- delta_DOX_raw[rownames(delta_DOX_raw) %in% target_gs, ]

# 3. Create Row Split Mapping
row_info <- stack(gene_categories) %>%
  rename(GeneSet = values, Category = ind) %>%
  filter(GeneSet %in% rownames(heatmap_mat))

heatmap_mat <- heatmap_mat[row_info$GeneSet, ]

# 4. Column Annotation (Cell Lines)
col_meta <- data.frame(sample_id = colnames(heatmap_mat)) %>%
  mutate(CellLine = sapply(strsplit(sample_id, "_"), `[`, 1))

column_ha <- HeatmapAnnotation(
  CellLine = col_meta$CellLine,
  show_annotation_name = TRUE
)

# 5. Plot
col_fun <- colorRamp2(c(-.1, 0, .1), c("blue", "white", "red"))

Heatmap(as.matrix(heatmap_mat), 
        name = "Delta ssGSEA",
        col = col_fun,
        row_split = row_info$Category,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        top_annotation = column_ha,
        rect_gp = gpar(col = "white", lwd = 0.5),
        row_names_gp = gpar(fontsize = 8),
        column_title = "DOXplatin-Induced Pathway Shifts")

```

# Investigation : What Genes Are most Negatively Correlated with MHC-I for DOX
```{r}



library(dplyr)
library(purrr)
library(tibble)
library(dplyr)

# -------------------------------------------------
# 1. Prepare DOX‑only data 
# -------------------------------------------------
DOX_samples <- delta_subset_wide %>%                     # delta_subset_wide is already filtered
  mutate(
    parts    = strsplit(as.character(sample_id), "_"),
    CellLine = sapply(parts, `[`, 2),
    Drug     = sapply(parts, `[`, 3)
  ) %>% 
  filter(Drug == "DOX") %>%                         # keep only DOX samples
  rowwise() %>% 
  ungroup()

# -------------------------------------------------
# 2. Build a **numeric** matrix of gene values
# -------------------------------------------------
gene_matrix <- DOX_samples %>% 
  # drop identifier columns
  select(-sample_id, -CellLine) %>%             
  # keep only columns that are (or can be) numeric
  mutate(across(everything(),
                ~ as.numeric(as.character(.x)))) %>% 
  # finally keep only the truly numeric columns (removes any remaining NA‑only cols)
  select(where(~ !all(is.na(.x))))              

# -------------------------------------------------
# 3. Correlate each gene with avg_mhc1
# -------------------------------------------------
DOX_mhc1_gene_cor_results <- lapply(colnames(gene_matrix), function(gene) {
  cor.test(gene_matrix[[gene]],
           DOX_samples$avg_mhc1,          # numeric vector
           method = "spearman")
})

# -------------------------------------------------
# 4. Extract estimates and p‑values
# -------------------------------------------------
DOX_mhc1_gene_cor_df <- data.frame(
  gene        = colnames(gene_matrix),
  correlation = sapply(DOX_mhc1_gene_cor_results, function(x) x$estimate),
  p_value     = sapply(DOX_mhc1_gene_cor_results, function(x) x$p.value)
) %>% 
  # optional: drop self‑correlation if it ever appears
  filter(!gene %in% c("avg_mhc1", mhc1_genes)) %>% 
  arrange(desc(correlation)) %>% 
  mutate(p_adj = p.adjust(p_value, method = "fdr"))




# Top 25
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ---- 1. Select top 25 positive and negative correlations ----
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

# -------------------------------------------------
# 1. Use the tidy results data frame (already built)
# -------------------------------------------------
# DOX_mhc1_gene_cor_df  # columns: gene, correlation, p_value, p_adj

# -------------------------------------------------
# 2. Top 25 positive correlations
# -------------------------------------------------
top_pos <- DOX_mhc1_gene_cor_df %>%          # <-- use the data frame
  arrange(desc(correlation)) %>%             # highest rho first
  slice_head(n = 25)                         # keep first 25 rows

# -------------------------------------------------
# 3. Top 25 negative correlations
# -------------------------------------------------
top_neg <- DOX_mhc1_gene_cor_df %>%
  arrange(correlation) %>%                   # most negative rho first
  slice_head(n = 25)

# -------------------------------------------------
# 4. Combine and keep only unique genes (a gene could appear in both lists)
# -------------------------------------------------
top50_cor <- bind_rows(top_pos, top_neg) %>%
  distinct(gene, .keep_all = TRUE)

gene_list <- top50_cor$gene   # 50 gene symbols for the heat‑map

# -------------------------------------------------
# 5. Build the heat‑map matrix (genes × samples)
# -------------------------------------------------
heatmap_mat <- delta_subset_wide %>%                     # full delta data
  select(sample_id, all_of(gene_list)) %>%          # keep only the 50 genes
  column_to_rownames("sample_id") %>%               # rows = samples
  t() %>%                                            # transpose → genes as rows
  as.matrix()

# -------------------------------------------------
# 6. Row annotation: positive vs. negative direction
# -------------------------------------------------
gene_direction <- c(
  rep("Pos_MHC1", nrow(top_pos)),
  rep("Neg_MHC1", nrow(top_neg))
)

row_ha <- rowAnnotation(
  Direction = gene_direction,
  col = list(
    Direction = c(
      Pos_MHC1 = "firebrick",
      Neg_MHC1 = "steelblue"
    )
  )
)

# -------------------------------------------------
# 7. Plot the heat‑map (no Z‑scoring, with cell borders)
# -------------------------------------------------
Heatmap(
  heatmap_mat,
  name            = "Expression",
  col             = colorRamp2(
    c(min(heatmap_mat), median(heatmap_mat), max(heatmap_mat)),
    c("blue", "white", "red")
  ),
  cluster_rows    = TRUE,
  cluster_columns = TRUE,
  show_row_names  = TRUE,
  show_column_names = FALSE,
  left_annotation = row_ha,
  rect_gp         = gpar(col = "black", lwd = 0.5),   # cell borders
  row_names_gp    = gpar(fontsize = 14)
)


# library(ggplot2)
# library(ggrepel)

# # Prepare data for plotting
# plot_data <- DOX_mhc1_gene_cor_df %>%
#   mutate(neg_log_p = -log10(p_value),
#          label = ifelse(p_adj < 0.05 & abs(correlation) > 0.5, gene, ""))

# ggplot(plot_data, aes(x = correlation, y = neg_log_p)) +
#   geom_point(aes(color = correlation > 0), alpha = 0.6) +
#   geom_text_repel(aes(label = label), max.overlaps = 15, box.padding = 0.5) +
#   scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue")) +
#   theme_minimal() +
#   labs(title = "Gene Correlation with avg_mhc1 (DOXplatin)",
#        x = "Spearman Correlation (rho)",
#        y = "-log10(p-value)",
#        color = "Positive Correlation") +
#   geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey") +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey")


library(ggplot2)
library(ggpubr)

# Define a function to plot correlation for any gene
plot_gene_correlation <- function(gene_name, data_matrix, meta_data) {
  
  # Create a temporary data frame for plotting
  df_plot <- data.frame(
    GeneExpression = data_matrix[[gene_name]],
    avg_mhc1 = meta_data$avg_mhc1,
    CellLine = meta_data$CellLine
  )
  
  ggplot(df_plot, aes(x = GeneExpression, y = avg_mhc1)) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
    geom_point(aes(color = CellLine), size = 3, alpha = 0.8) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    theme_minimal() +
    labs(
      title = paste("Correlation:", gene_name, "vs avg_mhc1"),
      subtitle = "Treatment: DOXplatin",
      x = paste(gene_name, "(Delta Expression)"),
      y = "Avg MHC-I (Delta Expression)"
    ) +
    scale_color_brewer(palette = "Set1")
}

# Example Usage:
# plot_gene_correlation("CD274", gene_matrix, DOX_samples)
# plot_gene_correlation("STING1", gene_matrix, DOX_samples)
```




# Investigation : What Genes Are most Negatively Correlated with PD-L1 for DOX
```{r}



library(dplyr)
library(purrr)
library(tibble)
library(dplyr)

# -------------------------------------------------
# 1. Prepare DOX‑only data and compute Avg_PDL1 (numeric)
# -------------------------------------------------
DOX_samples <- delta_subset_wide %>%                     # delta_subset_wide is already filtered
  mutate(
    parts    = strsplit(as.character(sample_id), "_"),
    CellLine = sapply(parts, `[`, 2),
    Drug     = sapply(parts, `[`, 3)
  ) %>% 
  filter(Drug == "DOX") %>%                         # keep only DOX samples
  rowwise() %>% 
  mutate(Avg_PDL1 = mean(c_across(any_of(c("CD274"))),
                         na.rm = TRUE)) %>%      # numeric result
  ungroup()

# -------------------------------------------------
# 2. Build a **numeric** matrix of gene values
# -------------------------------------------------
gene_matrix <- DOX_samples %>% 
  # drop identifier columns
  select(-sample_id, -CellLine) %>%             
  # keep only columns that are (or can be) numeric
  mutate(across(everything(),
                ~ as.numeric(as.character(.x)))) %>% 
  # finally keep only the truly numeric columns (removes any remaining NA‑only cols)
  select(where(~ !all(is.na(.x))))              

# -------------------------------------------------
# 3. Correlate each gene with Avg_PDL1
# -------------------------------------------------
DOX_PDL1_gene_cor_results <- lapply(colnames(gene_matrix), function(gene) {
  cor.test(gene_matrix[[gene]],
           DOX_samples$Avg_PDL1,          # numeric vector
           method = "spearman")
})

# -------------------------------------------------
# 4. Extract estimates and p‑values
# -------------------------------------------------
DOX_PDL1_gene_cor_df <- data.frame(
  gene        = colnames(gene_matrix),
  correlation = sapply(DOX_PDL1_gene_cor_results, function(x) x$estimate),
  p_value     = sapply(DOX_PDL1_gene_cor_results, function(x) x$p.value)
) %>%
  # optional: drop the target gene itself
  filter(!gene %in% c("Avg_PDL1", PDL1_genes)) %>%
  arrange(desc(correlation)) %>%          # now arrange works
  mutate(p_adj = p.adjust(p_value, method = "fdr"))




# Top 25
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ---- 1. Select top 25 positive and negative correlations ----

top_pos <- DOX_PDL1_gene_cor_results %>%
  arrange(desc(correlation)) %>%
  slice(1:25)

top_neg <- DOX_PDL1_gene_cor_results %>%
  arrange(correlation) %>%
  slice(1:25)

top50_cor <- bind_rows(top_pos, top_neg) %>%
  distinct(gene, .keep_all = TRUE)

gene_list <- top50_cor$gene

# ---- 2. Build heatmap matrix from original data ----

heatmap_mat <- delta_subset_wide %>%
  select(sample_id, all_of(gene_list)) %>%
  column_to_rownames("sample_id") %>%
  t() %>%                 # genes as rows
  as.matrix()

# Optional: log transform if needed
# heatmap_mat <- log2(heatmap_mat + 1)

# ---- 3. Row annotation (positive vs negative) ----

gene_direction <- c(rep("Pos_PDL1", 25), rep("Neg_PDL1", 25))

row_ha <- rowAnnotation(
  Direction = gene_direction,
  col = list(Direction = c("Pos_PDL1" = "firebrick",
                           "Neg_PDL1" = "steelblue"))
)

# ---- 4. Plot heatmap (NON-Z-scored, with borders) ----

Heatmap(
  heatmap_mat,
  name = "Expression",
  col = colorRamp2(
    c(min(heatmap_mat), median(heatmap_mat), max(heatmap_mat)),
    c("blue", "white", "red")
  ),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  left_annotation = row_ha,
  rect_gp = gpar(col = "black", lwd = 0.5),  # borders
  row_names_gp = gpar(fontsize = 14)
)

```



# Investigation : What Gene Programs are most negatively correlated with MHC-I for DOX
```{r}

library(tidyr)
library(tibble)

mhc_genes <- c("HLA-A", "HLA-B", "HLA-C")
 
pathway_df <- delta_DOX_raw %>%
  rownames_to_column("geneset") %>%
  pivot_longer(-geneset, names_to = "sample_id", values_to = "pathway_score") %>%
  pivot_wider(names_from = geneset, values_from = pathway_score)


formerge_tnbc_delta_final <- tnbc_delta_final %>%
  mutate(sample_id_clean = gsub("^delta_", "", sample_id)) %>%
  select(sample_id_clean, avg_mhc1)


merged_df <- pathway_df %>%
  inner_join(formerge_tnbc_delta_final, by = c("sample_id" = "sample_id_clean"))


library(purrr)

pathway_matrix <- merged_df %>%
  select(-sample_id, where(is.numeric)) %>%                     # drop the non‑numeric identifier
  as.data.frame() 



# -----------------------------------------------------------------
# 1. Keep only numeric columns (already done, but repeat for safety)
# -----------------------------------------------------------------
pathway_matrix <- merged_df %>%
  select(-sample_id, where(is.numeric)) %>%   # numeric scores + avg_mhc1
  as.data.frame()

# -----------------------------------------------------------------
# 2. Compute column variances, ignoring NA values
# -----------------------------------------------------------------
col_var <- apply(pathway_matrix, 2, var, na.rm = TRUE)

# -----------------------------------------------------------------
# 3. Build a logical mask: variance is finite and > 0
# -----------------------------------------------------------------
keep_cols <- !is.na(col_var) & col_var != 0

# -----------------------------------------------------------------
# 4. Sub‑set the matrix safely
# -----------------------------------------------------------------
pathway_matrix <- pathway_matrix[, keep_cols, drop = FALSE]

# -----------------------------------------------------------------
# 5. (Optional) sanity check
# -----------------------------------------------------------------
if (ncol(pathway_matrix) == 0) {
  warning("All columns were removed – check that the data contain variance.")
}


mhc1_gene_pathways_cor_results <- lapply(colnames(pathway_matrix)[
  colnames(pathway_matrix) != "avg_mhc1"
], function(gs) {
  
  cor.test(pathway_matrix[[gs]],
           pathway_matrix$avg_mhc1,
           method = "spearman")
})

mhc1_gene_pathways_cor_df <- data.frame(
  geneset = colnames(pathway_matrix)[
    colnames(pathway_matrix) != "avg_mhc1"
  ],
  correlation = sapply(mhc1_gene_pathways_cor_results, function(x) x$estimate),
  p_value = sapply(mhc1_gene_pathways_cor_results, function(x) x$p.value)
)

mhc1_gene_pathways_cor_df <- mhc1_gene_pathways_cor_df %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  arrange(desc(correlation))

```


# Investigation : What Gene Programs are most negatively correlated with PD-L1 for DOX
```{r}

library(tidyr)
library(tibble)
library(dplyr)
library(purrr)

# 1. Transform Pathway Data
pathway_df <- delta_DOX_raw %>%
  rownames_to_column("geneset") %>%
  pivot_longer(-geneset, names_to = "sample_id", values_to = "pathway_score") %>%
  pivot_wider(names_from = geneset, values_from = pathway_score)

# 2. Extract PD-L1 (CD274) induction instead of MHC-I
# We assume 'CD274' is the column name for PD-L1 in your tnbc_delta_final df
formerge_tnbc_delta_final <- tnbc_delta_final %>%
  mutate(sample_id_clean = gsub("^delta_", "", sample_id)) %>%
  # If you have multiple PD-L1 probes/isoforms, you can use rowMeans, 
  # but usually, just select the CD274 gene:
  mutate(target_PDL1 = CD274) %>% 
  select(sample_id_clean, target_PDL1)

# 3. Merge
merged_df <- pathway_df %>%
  inner_join(formerge_tnbc_delta_final, by = c("sample_id" = "sample_id_clean"))

# 4. Prepare Matrix and Remove Non-Variable Genes
pathway_matrix <- merged_df %>%
  select(-sample_id, where(is.numeric)) %>% 
  as.data.frame()

col_var <- apply(pathway_matrix, 2, var, na.rm = TRUE)
keep_cols <- !is.na(col_var) & col_var != 0
pathway_matrix <- pathway_matrix[, keep_cols, drop = FALSE]

if (ncol(pathway_matrix) == 0) {
  warning("All columns were removed – check variance.")
}

# 5. Run Spearman Correlation against PD-L1
# We iterate through all columns EXCEPT our target_PDL1
pdl1_gene_pathways_cor_results <- lapply(colnames(pathway_matrix)[
  colnames(pathway_matrix) != "target_PDL1"
], function(gs) {
  
  cor.test(pathway_matrix[[gs]],
           pathway_matrix$target_PDL1,
           method = "spearman",
           exact = FALSE) # exact=FALSE avoids warnings with ties
})

# 6. Final Results Dataframe
pdl1_gene_pathways_cor_df <- data.frame(
  geneset = colnames(pathway_matrix)[
    colnames(pathway_matrix) != "target_PDL1"
  ],
  correlation = sapply(pdl1_gene_pathways_cor_results, function(x) x$estimate),
  p_value = sapply(pdl1_gene_pathways_cor_results, function(x) x$p.value)
)

pdl1_gene_pathways_cor_df <- pdl1_gene_pathways_cor_df %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  arrange(desc(correlation))
```



# Investigation: which cytokines and chemokines are coexpressed


```{r}


library(dplyr)
library(purrr)
library(tibble)
library(dplyr)

# -------------------------------------------------
# 1. Prepare DOX‑only data and compute avg_mhc1 (numeric)
# -------------------------------------------------
DOX_samples <- delta_subset_wide %>%                     # delta_subset_wide is already filtered
  mutate(
    parts    = strsplit(as.character(sample_id), "_"),
    CellLine = sapply(parts, `[`, 2),
    Drug     = sapply(parts, `[`, 3)
  ) %>% 
  filter(Drug == "DOX") %>%                         # keep only DOX samples
  rowwise() %>%
  ungroup()

# -------------------------------------------------
# 2. Build a **numeric** matrix of gene values
# -------------------------------------------------
gene_matrix <- DOX_samples %>% 
  # drop identifier columns
  select(-sample_id, -CellLine) %>%             
  # keep only columns that are (or can be) numeric
  mutate(across(everything(),
                ~ as.numeric(as.character(.x)))) %>% 
  # finally keep only the truly numeric columns (removes any remaining NA‑only cols)
  select(where(~ !all(is.na(.x))))              

# Define your panels based on Gene Symbols
panel_genes <- c(
  "CXCL10", "CXCL9", "CXCL11", "HMGB1", "CCL28", "CXCL8", 
  "IL6", "TNF", "CCL20", "CX3CL1", "CXCL1", "CCL5",
  "CSF1", "IL1RN", "TNFSF13", "LTA"
)



# Subset your existing gene_matrix to include these + avg_mhc1
panel_mat <- gene_matrix[, colnames(gene_matrix) %in% c(panel_genes)]

# Recompute BOTH from the same numeric matrix, then align
mat <- as.matrix(panel_mat)

# (Optional but recommended) drop columns that are all NA or have zero variance
keep <- apply(mat, 2, function(x) !all(is.na(x)) && stats::sd(x, na.rm = TRUE) > 0)
mat  <- mat[, keep, drop = FALSE]


library(corrplot)


# ----------------------------
# 1) Correlation + p-value matrices
#    (Hmisc::rcorr gives both r and P)
# ----------------------------
rc <- Hmisc::rcorr(mat, type = "spearman")
cor_mat  <- rc$r
pval_mat <- rc$P


# Plot using p-value matrix to blank out non-sig cells
corrplot(cor_mat,
         method = "color",
         type = "upper",
         order = "hclust",
         tl.col = "black",
         tl.srt = 45,
         p.mat = pval_mat,
         sig.level = c(0.001, 0.01, 0.05),
         insig = "label_sig",   # prints significance labels
         pch.cex = 1.2,
         pch.col = "black")

         
# Optional: save/inspect p-values directly
pval_mat
```



# Investigation: linear models for association with DOX Considering that Baseline MHC is different

```{r}


library(dplyr)
library(stringr)
library(purrr)
library(broom)

## 1. Parse the delta identifiers (wide)
delta_df <- delta_subset_wide %>%
  mutate(
    base_sample_id = str_remove(sample_id, "^delta_"),
    parts = str_split_fixed(base_sample_id, "_", 3),
    cell_line = parts[, 1],
    drug      = parts[, 2],
    replicate = parts[, 3]
  ) %>% 
  select(-parts) %>%
  rename_with(
    .fn = ~ paste0(.x, "_delta"),
    .cols = -c(cell_line, drug, replicate)   # **do not rename replicate**
  )

## -------------------------------------------------
## 1. Baseline (pre‑treatment) – keep replicate info,
##    then average across replicates
## -------------------------------------------------
pre_df <- pre_treat_subset_wide %>%
  mutate(
    base_sample_id = sample_id,
    parts = str_split_fixed(sample_id, "_", 3),
    cell_line = parts[, 1],
    vehicle   = parts[, 2],
    replicate = parts[, 3]
  ) %>%
  select(-parts) %>%
  rename_with(
    .fn = ~ paste0(.x, "_baseline"),
    .cols = -c(cell_line, vehicle, replicate)   # do not rename the key columns
  )

# average across replicates for each (cell_line, vehicle) pair
pre_avg <- pre_df %>%
  group_by(cell_line, vehicle) %>%                 # *ignore* replicate here
  summarise(
    across(ends_with("_baseline"), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

## -------------------------------------------------
## 2. Post‑treatment – keep replicate info,
##    then average across replicates
## -------------------------------------------------
post_df <- post_treat_subset_wide %>%
  mutate(
    base_sample_id = sample_id,
    parts = str_split_fixed(sample_id, "_", 3),
    cell_line = parts[, 1],
    drug      = parts[, 2],
    replicate = parts[, 3]
  ) %>%
  select(-parts) %>%
  rename_with(
    .fn = ~ paste0(.x, "_post"),
    .cols = -c(cell_line, drug, replicate)        # keep key columns
  )

## 4. Join the three tables
full_df <- delta_df %>%
  left_join(pre_avg,  by = c("cell_line")) %>%
  left_join(post_df, by = c("cell_line", "drug", "replicate"))

# Ensure numeric columns are really numeric
full_df <- full_df %>%
  mutate(
    avg_mhc1_delta   = as.numeric(avg_mhc1_delta),
    avg_mhc1_baseline = as.numeric(avg_mhc1_baseline)
  )

full_df <- full_df %>%
  mutate(
    avg_mhc1_delta   = as.numeric(avg_mhc1_delta),
    avg_mhc1_baseline = as.numeric(avg_mhc1_baseline)
  )

## 1. Identify which columns are numeric in `full_df`
numeric_cols <- names(which(sapply(full_df, is.numeric)))

## 2. Keep only the numeric columns that end with “_delta”
##    (and drop the summary column `avg_mhc1_delta`)
delta_gene_cols <- numeric_cols[grepl("_delta$", numeric_cols) &
                               numeric_cols != "avg_mhc1_delta"]
# ------------------------------------------------------------------
#  Run a linear model for each gene
# ------------------------------------------------------------------
library(dplyr)
library(purrr)
library(lme4)
library(lmerTest)  # gives p-values


# ---- 1. Build the list of per‑gene results (same loop as before) ----
results_list <- list()
counter <- 1

for (delta_col in delta_gene_cols) {
  # Skip if column not present
  if (!delta_col %in% colnames(full_df)) next

  tmp_df <- full_df %>%
    dplyr::select(
      avg_mhc1_delta,
      avg_mhc1_baseline,
      cell_line,
      all_of(delta_col)
    ) %>%
    na.omit()

  # Require variation
  if (nrow(tmp_df) == 0 ||
      length(unique(tmp_df[[delta_col]])) < 2) next

  # Rename for modeling clarity
  colnames(tmp_df)[4] <- "gene_delta"

  fit <- lm(
    avg_mhc1_delta ~ gene_delta + avg_mhc1_baseline + cell_line,
    data = tmp_df
  )

  co <- summary(fit)$coefficients

  results_list[[delta_col]] <- tibble(
    gene      = sub("_delta$", "", delta_col),
    estimate  = co["gene_delta", "Estimate"],
    std.error = co["gene_delta", "Std. Error"],
    statistic = co["gene_delta", "t value"],
    p.value   = co["gene_delta", "Pr(>|t|)"]
  )

  counter <- counter + 1
}

results_df <- bind_rows(results_list)

# ---- 2. Remove empty slots and bind into a tibble -----------------
delta_mhc1_linear_model_results <- results_list %>% 
  compact() %>%          # drop NULL entries
  bind_rows() %>%        # creates a data.frame/tibble
  mutate(FDR = p.adjust(p.value, method = "BH")) %>% 
  arrange(p.value)       # equivalent to order(p.value)

# ---- 3. Keep only significant results --------------------------------
sig_delta_mhc1_linear_model_results <-
  delta_mhc1_linear_model_results %>% 
  filter(FDR <= 0.05)
```













# Defining CIM axes




```{r}
# Define ER-Stress / ISR module
er_stress_genes <- c(
  "ATF4", # Master ISR transcription factor
  "DDIT3", # CHOP, pro-apoptotic ER stress effector
  "EIF2AK3", # PERK, ER stress sensor kinase
  "XBP1", # UPR transcription factor
  "ATF6", # ER stress sensor
  "ERN1", # IRE1, splices XBP1
  "HSPA5" # GRP78/BiP chaperone
)

icd_damps <- c(
  "CALR", # Calreticulin – surface exposure promotes phagocytosis
  "HMGB1", # nuclear protein released during ICD, activates TLR4
  "ANXA1", # Annexin A1 – contributes to apoptotic signaling
  "PANX1", # pannexin channel involved in ATP release
  "VNUT",
  "PDIA3"
)

icd_sensors <- c(
  "TLR4", # HMGB1 receptor; Recognizes HMGB1 released by dying cells; triggers dendritic cell activation
  "P2RX7", # ATP receptor; ATP receptor; senses extracellular ATP released via PANX1 channels; activates inflammasome
  "CGAS", # cGAS; Cytosolic DNA sensor; part of DNA sensing → STING → IFN pathway
  "STING1", # STING
  "DDX58", # RIG-I; RNA sensor, triggers type I IFN signaling
  "IFIH1", # MDA5; RNA sensor, triggers type I IFN signaling
  "IRF3", # transcription factor; Downstream transcription factors that respond to cGAS/STING or RIG-I/IFIH1 signaling
  "IRF7" # transcription factor; Downstream transcription factors that respond to cGAS/STING or RIG-I/IFIH1 signaling
)

# ifn_module <- c("IFNB1", "IFNA1") # usually co-expressed so this should capture

ifn_sensor_and_signaling <- c(
  "IFNAR1",
  "IFNAR2",

  # Downstream ISGs
  "ISG15",
  "MX1",
  "OAS1",
  "RSAD2",
  "IRF9"
)


# Antigen presentation
mhc1_antigen_presentation = c(
  "HLA-A",
  "HLA-B",
  "HLA-C",
  "B2M",
  "TAP1",
  "TAP2"
)

# Immune checkpoints
immune_checkpoints = c(
  "CD274", # PD-L1
  "PDCD1LG2" # PD-L2
  # "CD47",      # “Don’t eat me” signal
  # "CD276",     # B7-H3
  # "VTCN1"      # B7-H4
)

# Anti-tumoral (Th1 recruitment + effector trafficking)
chemokines_anti_tumoral <- c(
  "CXCL9",
  "CXCL10",
  "CXCL11",
  "CCL5",
  "CX3CL1",
  "CXCL16",
  "IL1B"
)

# Pro-tumoral / suppressive
chemokines_pro_tumoral <- c(
  "CCL2",
  "CSF1",
  "CXCL8",
  "CCL20",
  "CCL22",
  "IL6",
  "IL1RN",
  "TGFB1",
  "IL10"
)


gene_sets_user_defined <- list(
  ER_Stress = er_stress_genes,
  ICD_DAMPs = icd_damps,
  ICD_Sensors = icd_sensors,
  #IFN = ifn_module,
  IFN_Signaling = ifn_sensor_and_signaling,
  mhc1_antigen_presentation = mhc1_antigen_presentation,
  Immune_Checkpoints = immune_checkpoints,
  Chemokines_AntiTumoral = chemokines_anti_tumoral,
  Chemokines_ProTumoral = chemokines_pro_tumoral
)
```




# scaled Averages for user_cim genesets

```{r}

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales) # for nice colors

# Assume your list of modules is `cimic_modules` and your data is `delta_subset_wide`
cimic_modules <- gene_sets_user_defined


# 1. Remove non-gene columns (like sample_id) for scaling
gene_cols <- setdiff(colnames(delta_subset_wide), "sample_id")
delta_scaled <- delta_subset_wide

# 2. Scale genes across all samples to give equal weight
delta_scaled[gene_cols] <- scale(delta_scaled[gene_cols], center = TRUE, scale = TRUE)

# 3. Function to compute module score (average of scaled genes per module)
compute_module_score <- function(df, module_genes, module_name){
  # Keep only genes in the module that exist in df
  genes_in_df <- intersect(module_genes, colnames(df))
  if(length(genes_in_df) == 0){
    warning(paste("No genes from module", module_name, "found in dataset"))
    return(NULL)
  }
  
  # Average across scaled genes
  score <- rowMeans(df[, genes_in_df, drop = FALSE], na.rm = TRUE)
  return(score)
}

# 4. Compute module scores for all modules
module_scores <- data.frame(sample_id = delta_scaled$sample_id)

for(mod_name in names(cimic_modules)){
  module_scores[[mod_name]] <- compute_module_score(delta_scaled, cimic_modules[[mod_name]], mod_name)
}

# 5. Add metadata for plotting (cell line and treatment)
module_scores <- module_scores %>%
  mutate(
    cell_line = str_extract(sample_id, "(?<=delta_)[^_]+"),
    treatment = str_extract(sample_id, "(DOX|EPI)")
  )


# 6. Convert to long format for ggplot
module_scores_long <- module_scores %>%
  pivot_longer(cols = names(cimic_modules), names_to = "module", values_to = "score")


# 7. Boxplot by module and cell line
library(dplyr)
library(ggplot2)
library(ggpubr)

# Filter to DOX samples only
module_scores_DOX <- module_scores_long %>%
  filter(treatment == "DOX")

# Plot
ggplot(module_scores_DOX, aes(x = cell_line, y = score, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.8) +
  facet_wrap(~module, scales = "free_y") +
  theme_bw() +
  labs(x = "Cell line", y = "Scaled module score") +
  
  # Stats for DOX only
  stat_compare_means(method = "kruskal.test", label.y.npc = "top") +
  stat_compare_means(method = "wilcox.test", ref.group = ".all.", label = "p.signif") +
  
  scale_fill_manual(values = c("DOX" = "#E41A1C")) +  # only DOX color
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


library(dplyr)
library(ggplot2)
library(ggpubr)

# Filter to EPI samples only
module_scores_EPI <- module_scores_long %>%
  filter(treatment == "EPI")

# Plot
ggplot(module_scores_EPI, aes(x = cell_line, y = score, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.8) +
  facet_wrap(~module, scales = "free_y") +
  theme_bw() +
  labs(x = "Cell line", y = "Scaled module score") +
  
  # Stats for EPI only
  stat_compare_means(method = "kruskal.test", label.y.npc = "top") +
  stat_compare_means(method = "wilcox.test", ref.group = ".all.", label = "p.signif") +
  
  scale_fill_manual(values = c("EPI" = "#E41A1C")) +  # only EPI color
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# CIM SCORE 

# Example: scale modules
module_scores_by_module <- module_scores_long %>%
  group_by(module) %>%
  mutate(module_scores = score) %>%
  ungroup()

# Define pro-CIM modules
pro_cim_modules <- c("ER_Stress", "ICD_DAMPs", "ICD_Sensors", "IFN", "IFN_Signaling",
                     "Antigen_Presentation", "Chemokines_AntiTumoral", "Immune_Checkpoints", "Chemokines_ProTumoral")


# Approach 1 -- General CIM scores 
# Compute CIM score
cim_scores <- module_scores_by_module %>%
  filter(module %in% pro_cim_modules) %>%
  group_by(sample_id, cell_line, treatment) %>%
  summarise(CIM_score = mean(module_scores, na.rm = TRUE))


ggplot(cim_scores, aes(x = cell_line, y = CIM_score, fill = treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.8) +
  theme_bw() +
  labs(x = "Cell line", y = "CIM score") +
  scale_fill_manual(values = c("DOX" = "#E41A1C", "EPI" = "#377EB8")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# --- DOX only ---
cim_scores_DOX <- cim_scores %>% filter(treatment == "DOX")

ggplot(cim_scores_DOX, aes(x = cell_line, y = CIM_score, fill = treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.8) +
  theme_bw() +
  labs(x = "Cell line", y = "CIM score") +
  scale_fill_manual(values = c("DOX" = "#E41A1C")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- EPI only ---
cim_scores_EPI <- cim_scores %>% filter(treatment == "EPI")

ggplot(cim_scores_EPI, aes(x = cell_line, y = CIM_score, fill = treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.8) +
  theme_bw() +
  labs(x = "Cell line", y = "CIM score") +
  scale_fill_manual(values = c("EPI" = "#377EB8")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Approach 2 by +/- system

library(dplyr)

cim_summary <- module_scores_long %>%
  group_by(cell_line, treatment, module) %>%
  summarize(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

cim_summary <- cim_summary %>%
  group_by(module, treatment) %>%
  mutate(
    percentile = percent_rank(mean_score),
    induction_rating = case_when(
      # --- Top Tier (Top 25% + Positive) ---
      mean_score > 0 & percentile >= 0.75 ~ "++",
      
      # --- High-Mid Tier (50-75% + Positive) ---
      mean_score > 0 & percentile >= 0.50 ~ "+",
      
      # --- The Neutral Zone (The "Zero" catch-all) ---
      # Catches: 1) Negative scores in the top 50% OR 2) Low positive scores (< 50%)
      mean_score < 0 & percentile >= 0.50 ~ "0",
      mean_score >= 0 & percentile < 0.50 ~ "0",
      
      # --- Low-Mid Tier (25-50% + Negative) ---
      mean_score < 0 & percentile >= 0.25 ~ "-",
      
      # --- Bottom Tier (Bottom 25% + Negative) ---
      TRUE ~ "--"
    )
  ) %>%
  ungroup()

library(tidyr)

cim_table <- cim_summary %>%
  select(cell_line, treatment, module, induction_rating) %>%
  pivot_wider(names_from = module, values_from = induction_rating)


cim_table_DOX <- cim_table %>% filter(treatment == "DOX")
cim_table_EPI <- cim_table %>% filter(treatment == "EPI")


library(dplyr)
library(tidyr)
library(pheatmap)

# 1️⃣ Map induction ratings to numeric values
rating_to_numeric <- c("++" = 2, "+" = 1, "0" = 0, "-" = -1, "--" = -2)  # or adjust scale if you like

cim_table_numeric <- cim_table %>%
  mutate(across(-c(cell_line, treatment), ~ rating_to_numeric[.])) %>%
  filter(treatment == "DOX") %>%  # example for DOX
  select(-treatment)

# 2️⃣ Convert to matrix for heatmap
cim_matrix <- cim_table_numeric %>%
  column_to_rownames("cell_line") %>%  # rownames = cell lines
  as.matrix()

# 3️⃣ Generate heatmap with clustering
pheatmap(
  cim_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "CIM Module Induction (DOX)"
)

library(dplyr)
library(ggplot2)
library(factoextra)  # for nice PCA plots


# 1️⃣ Prepare wide matrix for PCA
pca_data <- module_scores %>%
  select(sample_id, names(cimic_modules)) %>%
  column_to_rownames("sample_id")

# 2️⃣ Add treatment info
pca_meta <- module_scores %>%
  select(sample_id, cell_line, treatment) %>%
  distinct()

# 3️⃣ Function to run PCA for a subset
run_pca_plot <- function(subset_treatment) {
  # Filter samples
  subset_ids <- pca_meta %>% filter(treatment == subset_treatment) %>% pull(sample_id)
  subset_matrix <- pca_data[subset_ids, ]
  
  # PCA
  pca_res <- prcomp(subset_matrix, scale. = TRUE)
  
  # Scores + metadata
  pca_scores <- as.data.frame(pca_res$x) %>%
    mutate(sample_id = rownames(.)) %>%
    left_join(pca_meta, by = "sample_id")
  
  # Plot
  ggplot(pca_scores, aes(x = PC1, y = PC2, color = cell_line)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_bw() +
    labs(
      x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2,1],1), "%)"),
      y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2,2],1), "%)"),
      title = paste("PCA of CIM Module Scores:", subset_treatment)
    )
}

# 4️⃣ Run separately for DOX and EPI
pca_DOX_plot <- run_pca_plot("DOX")
pca_EPI_plot <- run_pca_plot("EPI")

# Display plots
pca_DOX_plot
pca_EPI_plot

# # by CIMIC clusteres

# cluster_df <- cimic_clusters %>%
#   distinct(cell_lines, cimic_cluster_DOX, cimic_cluster_EPI)

# cluster_df <- cluster_df %>%
#   rename(cell_line = cell_lines)

# module_scores_long <- module_scores_long %>%
#   left_join(cluster_df, by = "cell_line")

# ggplot(module_scores_long %>% filter(treatment=="DOX"),
#        aes(x = cimic_cluster_DOX, y = score, fill = factor(cimic_cluster_DOX))) +
#   geom_boxplot() +
#   facet_wrap(~module)


# ggplot(module_scores_long %>% filter(treatment=="EPI"),
#        aes(x = cimic_cluster_EPI, y = score, fill = factor(cimic_cluster_EPI))) +
#   geom_boxplot() +
#   facet_wrap(~module)

```



# ssGSEA for user_cim genesets 

```{r}



# ssGSEA Function
ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
    row_names = rownames(X)
    num_genes = nrow(X)
    gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})

    # Ranks for genes
    R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')

    # Calculate enrichment score (es) for each sample (column)
    es = apply(R, 2, function(R_col) {
        gene_ranks = order(R_col, decreasing = TRUE)

        # Calc es for each gene set
        es_sample = sapply(gene_sets, function(gene_set_idx) {
            # pos: match (within the gene set)
            # neg: non-match (outside the gene set)
            indicator_pos = gene_ranks %in% gene_set_idx
            indicator_neg = !indicator_pos

            rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha

            step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
            step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)

            step_cdf_diff = step_cdf_pos - step_cdf_neg

            # Normalize by gene number
            if (scale) step_cdf_diff = step_cdf_diff / num_genes

            # Use ssGSEA or not
            if (single) {
                sum(step_cdf_diff)
            } else {
                step_cdf_diff[which.max(abs(step_cdf_diff))]
            }
        })
        unlist(es_sample)
    })

    if (length(gene_sets) == 1) es = matrix(es, nrow = 1)

    # Normalize by absolute diff between max and min
    if (norm) es = es / diff(range(es))

    # Prepare output
    rownames(es) = names(gene_sets)
    colnames(es) = colnames(X)
    return(es)
}



library(dplyr)
library(tibble)
library(GSVA)  # for ssgsea
library(tidyr)

# ---- 1. Convert to numeric matrix ----
tnbc_mat <- tnbc_raw %>%
  remove_rownames() %>%
  column_to_rownames("gene_id") %>%
  as.matrix() %>%
  { matrix(as.numeric(.), nrow = nrow(.), dimnames = dimnames(.)) }

# ---- 2. Function to extract columns for a given drug / control ----
get_drug_matrix <- function(mat, drug_name, control_name = "DMSO") {
  
  drug_cols    <- grep(paste0("_", drug_name, "_"), colnames(mat), value = TRUE)
  control_cols <- grep(paste0("_", control_name, "_"), colnames(mat), value = TRUE)
  
  list(
    drug    = mat[, drug_cols, drop = FALSE],
    control = mat[, control_cols, drop = FALSE]
  )
}


# ---- 3. Function to calculate delta ssGSEA ----
compute_delta_ssgsea <- function(drug_mat, control_mat, gene_sets, scale = TRUE, norm = FALSE) {
  
  # ssGSEA for drug
  es_drug <- ssgsea(drug_mat, gene_sets, scale = scale, norm = norm)
  
  # ssGSEA for control
  es_control <- ssgsea(control_mat, gene_sets, scale = scale, norm = norm)
  
  # delta
  es_delta <- es_drug - es_control
  
  # optional: row-wise z-score
  es_delta_zrow <- t(apply(es_delta, 1, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
  
  list(
    delta_raw = as.data.frame(es_delta),
    delta_zrow = as.data.frame(es_delta_zrow),
    drug_es = as.data.frame(es_drug),
    control_es = as.data.frame(es_control)
  )
}

# ---- 4. Example: EPI vs DMSO ----
EPI_data <- get_drug_matrix(tnbc_mat, drug_name = "EPI", control_name = "SALINE")
ssgsea_EPI_delta <- compute_delta_ssgsea(EPI_data$drug, EPI_data$control, gene_sets_user_defined)

# Access results
delta_EPI_raw   <- ssgsea_EPI_delta$delta_raw
delta_EPI_zrow  <- ssgsea_EPI_delta$delta_zrow
EPI_es_drug     <- ssgsea_EPI_delta$drug_es
EPI_es_control  <- ssgsea_EPI_delta$control_es



# ---- 5. Example: DOX vs DMSO ----
DOX_data <- get_drug_matrix(tnbc_mat, drug_name = "DOX", control_name = "SALINE")
ssgsea_DOX_delta <- compute_delta_ssgsea(DOX_data$drug, DOX_data$control, gene_sets_user_defined)

# Access results
delta_DOX_raw   <- ssgsea_DOX_delta$delta_raw
delta_DOX_zrow  <- ssgsea_DOX_delta$delta_zrow
DOX_es_drug     <- ssgsea_DOX_delta$drug_es
DOX_es_control  <- ssgsea_DOX_delta$control_es

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(tibble)

# -------------------------------------------------
# 1. Prepare Tidy Data from delta_DOX_zrow
# -------------------------------------------------
# We extract metadata directly from the column names of your provided dataframe
plot_df_DOX <- delta_DOX_raw %>% 
  rownames_to_column("GeneSet") %>% 
  pivot_longer(
    cols      = -GeneSet,
    names_to  = "sample_id",
    values_to = "Delta"
  ) %>%
  mutate(
    # Extracts "BT549" from "BT549_DOX_R1"
    CellLine = sapply(strsplit(sample_id, "_"), `[`, 1),
    # Extracts "DOX" from "BT549_DOX_R1"
    Drug     = sapply(strsplit(sample_id, "_"), `[`, 2)
  ) %>%
  filter(Drug == "DOX")

# Define the list of gene sets to loop over if not already defined
all_gene_set_names <- unique(plot_df_DOX$GeneSet)

# -------------------------------------------------
# 2. Loop over GeneSets for Plots and Stats
# -------------------------------------------------
all_plots  <- list()
stats_list <- list()

for (gs in all_gene_set_names) {
  current_data <- plot_df_DOX %>% filter(GeneSet == gs)

  if (nrow(current_data) < 2) next

  # --- Plot with Wilcox and Kruskal ---
  p <- ggplot(current_data, aes(x = CellLine, y = Delta, fill = CellLine)) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2) +
    # Global non-parametric test
    stat_compare_means(method = "kruskal.test", label.y.npc = "top") +
    # Comparison of each cell line against the grand median
    stat_compare_means(method = "wilcox.test", 
                       ref.group = ".all.", 
                       label = "p.signif") +
    theme_bw() +
    labs(title = gs, y = "Delta ssGSEA", x = "Cell Line") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

  all_plots[[gs]] <- p

  # --- Statistics Collection ---
  gs_stats <- tryCatch({
      compare_means(Delta ~ CellLine, 
                    data = current_data, 
                    method = "wilcox.test", 
                    ref.group = ".all.") %>%
        mutate(GeneSet = gs)
    }, error = function(e) {
      tibble(GeneSet = gs, p = NA_real_, method = "wilcox.test", group1 = ".all.", group2 = NA_character_)
    })
  
  stats_list[[gs]] <- gs_stats
}

# -------------------------------------------------
# 3. Combine and Adjust P-values
# -------------------------------------------------
all_stats_df <- bind_rows(stats_list) %>%
  # Ensure p is numeric for adjustment
  mutate(p = as.numeric(p)) %>% 
  # Apply Benjamini-Hochberg adjustment across all tests performed
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  # Re-assign significance stars based on adjusted p-values
  mutate(p.adj.signif = case_when(
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01  ~ "**",
    p.adj <= 0.05  ~ "*",
    TRUE           ~ "ns"
  )) %>%
  select(GeneSet, group2, p, p.adj, p.adj.signif, method) %>%
  rename(CellLine = group2) %>%
  arrange(p.adj)

# -------------------------------------------------
# 4. Final Review
# -------------------------------------------------
print(head(all_stats_df))

# Example: Display the plot for Immune Response mediators
# print(all_plots[["GOBP_PRODUCTION_OF_MOLECULAR_MEDIATOR_OF_IMMUNE_RESPONSE"]])
```


# Figuring out death programs at work


```{r}

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)

# -----------------------------
# 1. Define gene sets
# -----------------------------
death_programs <- list(
  Apoptosis = c(
    "BAX", "BAK1", "BBC3", "PMAIP1", "BCL2L11",
    "APAF1", "CASP9", "CASP3", "CASP7",
    "FAS", "TNFRSF10B", "CASP8", "FADD", "BID",
    "DIABLO", "CYCS", "TP53AIP1"
  ),
  Necroptosis = c(
    "RIPK1", "RIPK3", "MLKL", "ZBP1", "TICAM1"
  ),
  Pyroptosis = c(
    "GSDME", "GSDMD", "CASP1", "CASP4", "CASP5",
    "NLRP3", "AIM2", "PYCARD", "IL1B", "IL18"
  ),
  PANoptosis = c(
    "ZBP1", "AIM2", "RIPK3", "RIPK1", "CASP8",
    "CASP1", "FADD", "IRF1", "PYCARD"
  ),
  Ferroptosis = c(
    "ACSL4", "LPCAT3", "TFRC", "SAT1", "PTGS2"
  )
)

other_programs <- list(
  Senescence = c(
    "CDKN1A", "CDKN2A", "CDKN2B", "GLB1", "SERPINE1",
    "IL6", "CXCL8", "CXCL1", "CCL2", "MMP1", "MMP3"
  ),
  ICD_Markers = c(
    "CALR", "HMGB1", "ANXA1", "PANX1",
    "CXCL10", "CXCL11", "IFNB1", "CCL5"
  ),
  Cytosolic_DNA_Sensing = c(
    "CGAS", "MB21D1", "TMEM173", "STING1", "TBK1", "IKBKE",
    "IRF3", "IRF7", "IFI16", "AIM2", "ZBP1",
    "IFNB1", "CXCL10", "CXCL11", "CCL5", "ISG15", "IFIT1", "IFIT3", "MX1", "OAS1"
  ),
  ER_Stress_ISR = c(
    "DDIT3", "ATF4", "ATF3", "XBP1", "ERN1", "EIF2AK3",
    "HSPA5", "DNAJB9", "PPP1R15A", "EDEM1", "HERPUD1", "TRIB3"
  ),
  Proteostasis = c(
    "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA8", "HSPB1",
    "DNAJA1", "DNAJB1", "BAG3", "STIP1",
    "HSPA5", "CANX", "CALR", "PDIA4", "PDIA6", "ERP44",
    "HERPUD1", "EDEM1", "SEL1L", "SYVN1",
    "PSMB5", "PSMD14", "UBC", "SQSTM1",
    "MAP1LC3B", "ATG5", "ATG7", "BECN1", "LAMP1", "CTSD",
    "HMOX1", "NQO1", "TXN", "TXNRD1", "GCLC"
  ),
  Immune_Checkpoints = c(
    "CD274", "PDCD1LG2", "CD47", "LGALS9", "PVR", "NECTIN2",
    "IDO1", "VTCN1", "CD276", "HHLA2"
  ),
  AntiTumoral_Cytokines = c(
    "CXCL9", "CXCL10", "CXCL11", "CCL5",
    "IFNB1", "IFNA1", "IFNG", "TNF",
    "IL12A", "IL12B", "IL15", "CSF2"
  ),
  ProTumoral_Cytokines = c(
    "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL8",
    "CCL20", "CCL22", "CCL28",
    "IL6", "IL1B", "TGFB1", "CSF1", "VEGFA"
  )
)

all_programs <- c(death_programs, other_programs)

# -----------------------------
# 2. Metadata
# -----------------------------
df <- delta_subset_wide %>%
  mutate(
    sample_id = as.character(sample_id),
    cell_line = str_match(sample_id, "^delta_([^_]+)")[,2],
    treatment = str_match(sample_id, "^delta_[^_]+_([^_]+)")[,2]
  )

# -----------------------------
# 3. Keep genes present
# -----------------------------
all_gene_cols <- colnames(df)
all_programs_present <- lapply(all_programs, function(genes) intersect(genes, all_gene_cols))
print(all_programs_present)

# -----------------------------
# 4. Gene-wise z-score across samples
# -----------------------------
meta_cols <- c("sample_id", "cimic_cluster", "cell_line", "treatment")
gene_cols <- setdiff(colnames(df), meta_cols)

zscore_safe <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(x)) else (x - m) / s
}

df_z <- df %>%
  mutate(across(all_of(gene_cols), zscore_safe))

# -----------------------------
# 5. Score each program
# -----------------------------
score_program_z <- function(data, genes, program_name) {
  if (length(genes) == 0) {
    return(
      data.frame(
        sample_id = data$sample_id,
        cimic_cluster = data$cimic_cluster,
        cell_line = data$cell_line,
        treatment = data$treatment,
        Program = program_name,
        Score = NA_real_
      )
    )
  }
  
  score_vec <- rowMeans(as.matrix(data[, genes, drop = FALSE]), na.rm = TRUE)
  
  data.frame(
    sample_id = data$sample_id,
    cimic_cluster = data$cimic_cluster,
    cell_line = data$cell_line,
    treatment = data$treatment,
    Program = program_name,
    Score = score_vec
  )
}


program_scores_all <- purrr::imap_dfr(all_programs_present, ~score_program_z(df_z, .x, .y))

# -----------------------------
# 6. Summary
# -----------------------------
program_summary <- program_scores_all %>%
  group_by(cell_line, cimic_cluster, Program) %>%
  summarise(
    mean_score = mean(Score, na.rm = TRUE),
    median_score = median(Score, na.rm = TRUE),
    sd_score = sd(Score, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )

print(program_summary)

cimic_cluster_COLORS <- c("1" = "#1b9e77", "2" = "#d95f02")

# -----------------------------
# 7. Faceted boxplot by cluster
# -----------------------------
# -----------------------------
# 7. Faceted boxplot by cluster
# -----------------------------
program_scores_plot_df <- program_scores_all %>%
  mutate(
    cimic_cluster = factor(cimic_cluster),

    # Use unique names for the factor levels to avoid duplicates
    Program = factor(
      Program,
      levels = unique(names(all_programs))   # ← deduplicated level list
    )
  )

pval_df_all <- program_scores_plot_df %>%
  group_by(Program) %>%
  summarise(
    y.position = max(Score, na.rm = TRUE) * 1.12,
    pval = tryCatch(kruskal.test(Score ~ cimic_cluster)$p.value, error = function(e) NA_real_),
    p_label = ifelse(is.na(pval), "p = NA",
                     ifelse(pval < 0.001, "p < 0.001", paste0("p = ", signif(pval, 3)))),
    .groups = "drop"
  )

p_all_programs_p <- ggplot(program_scores_plot_df, aes(x = cimic_cluster, y = Score, fill = cimic_cluster)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.4) +
  geom_text(
    data = pval_df_all,
    aes(x = 1.5, y = y.position, label = p_label),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold"
  ) +
  facet_wrap(~ Program, scales = "free_y") +
  scale_fill_manual(values = cimic_cluster_COLORS) +
  labs(
    title = "Program scores by CIMiC cluster",
    y = "Program z-score",
    x = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    axis.line = element_line(colour = "black", linewidth = 1.5)
  )

print(p_all_programs_p)

# Keep only significant programs
sig_programs <- pval_df_all %>%
  filter(!is.na(pval), pval < 0.05) %>%
  pull(Program)

program_scores_plot_df_sig <- program_scores_plot_df %>%
  filter(Program %in% sig_programs)

pval_df_sig <- pval_df_all %>%
  filter(Program %in% sig_programs)

p_all_programs_sig <- ggplot(program_scores_plot_df_sig, aes(x = cimic_cluster, y = Score, fill = cimic_cluster)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.4) +
  geom_text(
    data = pval_df_sig,
    aes(x = 1.5, y = y.position, label = p_label),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold"
  ) +
  facet_wrap(~ Program, scales = "free_y") +
  scale_fill_manual(values = cimic_cluster_COLORS) +
  labs(
    title = "Significant program scores by CIMiC cluster",
    y = "Program z-score",
    x = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    axis.line = element_line(colour = "black", linewidth = 1.5)
  )

print(p_all_programs_sig)

# -----------------------------
# 8. Heatmap by cell line + cluster
# -----------------------------
program_heatmap_df <- program_scores_all %>%
  group_by(cell_line, cimic_cluster, Program) %>%
  summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    cell_cluster = paste0(cell_line, "_C", cimic_cluster),
    # Use only unique names for the factor levels
    Program = factor(Program, levels = unique(names(all_programs)))
  )



p_heatmap <- ggplot(program_heatmap_df, aes(x = Program, y = cell_cluster, fill = mean_score)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_gradient2(
    low = "#313695",
    mid = "white",
    high = "#A50026",
    midpoint = 0,
    name = "Mean/nz-score"
  ) +
  labs(
    title = "Program landscape by cell line and CIMiC cluster",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 11),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)
  )

print(p_heatmap)
library(dplyr)
library(ggplot2)

# -------------------------------------------------
# 9. Cluster heatmap (same issue)
# -------------------------------------------------
cluster_heatmap_df <- program_scores_all %>%
  group_by(cimic_cluster, Program) %>%
  summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    cimic_cluster = factor(cimic_cluster),
    Program = factor(Program, levels = unique(names(all_programs)))   # ← deduped levels
  )


# 2. Per-program Wilcoxon p-values and FDR correction
program_pvals <- program_scores_all %>%
  mutate(cimic_cluster = factor(cimic_cluster)) %>%
  group_by(Program) %>%
  summarise(
    p_value = tryCatch(wilcox.test(Score ~ cimic_cluster)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(p_value, method = "fdr"),
    fdr_label = case_when(
      is.na(FDR) ~ "FDR=NA",
      FDR < 0.001 ~ "FDR<0.001",
      TRUE ~ paste0("FDR=", signif(FDR, 2))
    )
  )

# 3. Append FDR to x-axis labels
program_labels <- setNames(
  paste0(program_pvals$Program, "/n", program_pvals$fdr_label),
  program_pvals$Program
)

# 4. Plot
p_cluster_heatmap <- ggplot(cluster_heatmap_df, aes(x = Program, y = cimic_cluster, fill = mean_score)) +
  geom_tile(color = "white", linewidth = 0.8) +
  #geom_text(aes(label = round(mean_score, 2)), size = 4, fontface = "bold") +
  scale_fill_gradient2(
    low = "#313695",
    mid = "white",
    high = "#A50026",
    midpoint = 0,
    name = "Mean/nz-score"
  ) +
  # scale_x_discrete(labels = program_labels) +
  labs(
    title = "Mean program z-score by CIMiC cluster",
    x = NULL,
    y = "CIMiC cluster"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)
  )

print(p_cluster_heatmap)
```


# checking individual genes for apoptosis 

```{r}
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# Apoptosis genes
apoptosis_genes <- c(
  "BAX", "BAK1", "BBC3", "PMAIP1", "APAF1", "DIABLO", "CYCS",
  "CASP9", "CASP3", "FAS", "FASLG", "TNFRSF10B", "TNFRSF10A",
  "TNFSF10", "CASP8", "BID", "FADD", "CFLAR", "BCL2L11", "TP53AIP1"
)

# Keep only genes present
apoptosis_genes_present <- intersect(apoptosis_genes, colnames(delta_subset_wide))

# Long format
gene_df <- delta_subset_wide %>%
  mutate(
    sample_id = as.character(sample_id),
    cell_line = stringr::str_remove(sample_id, "^delta_"),
    cell_line = stringr::str_extract(cell_line, "^[^_]+"),
    cimic_cluster = factor(cimic_cluster)
  ) %>%
  select(sample_id, cell_line, cimic_cluster, all_of(apoptosis_genes_present)) %>%
  pivot_longer(
    cols = all_of(apoptosis_genes_present),
    names_to = "Gene",
    values_to = "Expression"
  )

# Box colors for clusters
cimic_cluster_COLORS <- c(
  "1" = "#1b9e77",
  "2" = "#d95f02"
)

# Point colors for cell lines
cell_line_colors <- c(
  "DU4475"   = "#1f78b4",
  "HCC1806"  = "#33a02c",
  "HCC1395"  = "#e31a1c",
  "Hs578t"   = "#ff7f00",
  "MDAMB157" = "#6a3d9a",
  "MDAMB231" = "#b15928",
  "MDAMB468" = "#a6cee3",
  "BT549"    = "#fb9a99",
  "HCC38"    = "#b2df8a"
)

# p-value positions per facet
pval_df <- gene_df %>%
  group_by(Gene) %>%
  summarise(
    y.position = max(Expression, na.rm = TRUE) * 1.12,
    .groups = "drop"
  )

# Plot
p_apoptosis <- ggplot(gene_df, aes(x = cimic_cluster, y = Expression)) +
  geom_boxplot(
    aes(fill = cimic_cluster),
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    aes(color = cell_line),
    width = 0.15,
    size = 2,
    alpha = 0.7
  ) +
  stat_compare_means(
    comparisons = list(c("1", "2")),
    method = "wilcox.test",
    label = "p.format",
    label.y = pval_df$y.position
  ) +
  scale_fill_manual(values = cimic_cluster_COLORS) +
  scale_color_manual(values = cell_line_colors) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(
    title = "Apoptosis genes by CIMiC cluster",
    y = "delta(log2(TPM+1))",
    x = NULL,
    color = "Cell line"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    axis.line = element_line(colour = "black", linewidth = 1.5)
  )

print(p_apoptosis)

```


# Checking out individual genes for panoptosis

```{r}

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# Panoptosis genes
panoptosis_genes <- c("ZBP1", "AIM2", "RIPK3", "RIPK1", "CASP8", "CASP1", "FADD", "IRF1", "ADAR",
"MX1", "CGAS", "STING1")

# Keep only genes present
panoptosis_genes_present <- intersect(panoptosis_genes, colnames(delta_subset_wide))

# Long format
gene_df <- delta_subset_wide %>%
  mutate(
    sample_id = as.character(sample_id),
    cell_line = stringr::str_remove(sample_id, "^delta_"),
    cell_line = stringr::str_extract(cell_line, "^[^_]+"),
    cimic_cluster = factor(cimic_cluster)
  ) %>%
  select(sample_id, cell_line, cimic_cluster, all_of(panoptosis_genes_present)) %>%
  pivot_longer(
    cols = all_of(panoptosis_genes_present),
    names_to = "Gene",
    values_to = "Expression"
  )

# Box colors for clusters
cimic_cluster_COLORS <- c(
  "1" = "#1b9e77",
  "2" = "#d95f02"
)

# Point colors for cell lines
cell_line_colors <- c(
  "DU4475"   = "#1f78b4",
  "HCC1806"  = "#33a02c",
  "HCC1395"  = "#e31a1c",
  "Hs578t"   = "#ff7f00",
  "MDAMB157" = "#6a3d9a",
  "MDAMB231" = "#b15928",
  "MDAMB468" = "#a6cee3",
  "BT549"    = "#fb9a99",
  "HCC38"    = "#b2df8a"
)

# p-value positions per facet
pval_df <- gene_df %>%
  group_by(Gene) %>%
  summarise(
    y.position = max(Expression, na.rm = TRUE) * 1.12,
    .groups = "drop"
  )

# Plot
p_panoptosis <- ggplot(gene_df, aes(x = cimic_cluster, y = Expression)) +
  geom_boxplot(
    aes(fill = cimic_cluster),
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    aes(color = cell_line),
    width = 0.15,
    size = 2,
    alpha = 0.7
  ) +
  stat_compare_means(
    comparisons = list(c("1", "2")),
    method = "wilcox.test",
    label = "p.format",
    label.y = pval_df$y.position
  ) +
  scale_fill_manual(values = cimic_cluster_COLORS) +
  scale_color_manual(values = cell_line_colors) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(
    title = "Panoptosis genes by CIMiC cluster",
    y = "delta(log2(TPM+1))",
    x = NULL,
    color = "Cell line"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    axis.line = element_line(colour = "black", linewidth = 1.5)
  )

print(p_panoptosis)

```


# checking out individual genes for ferroptosis

```{r}




library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# Ferroptosis genes
ferroptosis_genes <- c("GPX4", "SLC7A11", "SLC3A2", "SAT1", "TFRC", "ACSL4", "LPCAT3", "PTGS2")

# Keep only genes present
ferroptosis_genes_present <- intersect(ferroptosis_genes, colnames(delta_subset_wide))

# Long format
gene_df <- delta_subset_wide %>%
  mutate(
    sample_id = as.character(sample_id),
    cell_line = stringr::str_remove(sample_id, "^delta_"),
    cell_line = stringr::str_extract(cell_line, "^[^_]+"),
    cimic_cluster = factor(cimic_cluster)
  ) %>%
  select(sample_id, cell_line, cimic_cluster, all_of(ferroptosis_genes_present)) %>%
  pivot_longer(
    cols = all_of(ferroptosis_genes_present),
    names_to = "Gene",
    values_to = "Expression"
  )

# Box colors for clusters
cimic_cluster_COLORS <- c(
  "1" = "#1b9e77",
  "2" = "#d95f02"
)

# Point colors for cell lines
cell_line_colors <- c(
  "DU4475"   = "#1f78b4",
  "HCC1806"  = "#33a02c",
  "HCC1395"  = "#e31a1c",
  "Hs578t"   = "#ff7f00",
  "MDAMB157" = "#6a3d9a",
  "MDAMB231" = "#b15928",
  "MDAMB468" = "#a6cee3",
  "BT549"    = "#fb9a99",
  "HCC38"    = "#b2df8a"
)

# p-value positions per facet
pval_df <- gene_df %>%
  group_by(Gene) %>%
  summarise(
    y.position = max(Expression, na.rm = TRUE) * 1.12,
    .groups = "drop"
  )

# Plot
p_ferroptosis <- ggplot(gene_df, aes(x = cimic_cluster, y = Expression)) +
  geom_boxplot(
    aes(fill = cimic_cluster),
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    aes(color = cell_line),
    width = 0.15,
    size = 2,
    alpha = 0.7
  ) +
  stat_compare_means(
    comparisons = list(c("1", "2")),
    method = "wilcox.test",
    label = "p.format",
    label.y = pval_df$y.position
  ) +
  scale_fill_manual(values = cimic_cluster_COLORS) +
  scale_color_manual(values = cell_line_colors) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(
    title = "Ferroptosis genes by CIMiC cluster",
    y = "delta(log2(TPM+1))",
    x = NULL,
    color = "Cell line"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    axis.line = element_line(colour = "black", linewidth = 1.5)
  )

print(p_ferroptosis)
```


# checking out individual genes for ICD 

```{r}


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# ICD marker genes
icd_genes <- icd_genes <- c("CALR", "HMGB1", "ANXA1", "PANX1", "ENTPD1", "NT5E", "CXCL10", "CXCL11")

# Keep only genes present
icd_genes_present <- intersect(icd_genes, colnames(delta_subset_wide))

# Long format
gene_df <- delta_subset_wide %>%
  mutate(
    sample_id = as.character(sample_id),
    cell_line = stringr::str_remove(sample_id, "^delta_"),
    cell_line = stringr::str_extract(cell_line, "^[^_]+"),
    cimic_cluster = factor(cimic_cluster)
  ) %>%
  select(sample_id, cell_line, cimic_cluster, all_of(icd_genes_present)) %>%
  pivot_longer(
    cols = all_of(icd_genes_present),
    names_to = "Gene",
    values_to = "Expression"
  )

# Box colors for clusters
cimic_cluster_COLORS <- c(
  "1" = "#1b9e77",
  "2" = "#d95f02"
)

# Point colors for cell lines
cell_line_colors <- c(
  "DU4475"   = "#1f78b4",
  "HCC1806"  = "#33a02c",
  "HCC1395"  = "#e31a1c",
  "Hs578t"   = "#ff7f00",
  "MDAMB157" = "#6a3d9a",
  "MDAMB231" = "#b15928",
  "MDAMB468" = "#a6cee3",
  "BT549"    = "#fb9a99",
  "HCC38"    = "#b2df8a"
)

# p-value positions per facet
pval_df <- gene_df %>%
  group_by(Gene) %>%
  summarise(
    y.position = max(Expression, na.rm = TRUE) * 1.12,
    .groups = "drop"
  )

# Plot
p_icd <- ggplot(gene_df, aes(x = cimic_cluster, y = Expression)) +
  geom_boxplot(
    aes(fill = cimic_cluster),
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    aes(color = cell_line),
    width = 0.15,
    size = 2,
    alpha = 0.7
  ) +
  stat_compare_means(
    comparisons = list(c("1", "2")),
    method = "wilcox.test",
    label = "p.format",
    label.y = pval_df$y.position
  ) +
  scale_fill_manual(values = cimic_cluster_COLORS) +
  scale_color_manual(values = cell_line_colors) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(
    title = "ICD markers by CIMiC cluster",
    y = "delta(log2(TPM+1))",
    x = NULL,
    color = "Cell line"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    axis.line = element_line(colour = "black", linewidth = 1.5)
  )

print(p_icd)

```



# checking out individual genes for pyroptosis 
```{r}
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# Pyroptosis genes
pyroptosis_genes <- c("GSDME", "GSDMD", "CASP1", "CASP4", "CASP5", "NLRP3", "AIM2", "PYCARD", "IL1B", "IL18")

# Keep only genes present
pyroptosis_genes_present <- intersect(pyroptosis_genes, colnames(delta_subset_wide))

# Long format
gene_df <- delta_subset_wide %>%
  mutate(
    sample_id = as.character(sample_id),
    cell_line = stringr::str_remove(sample_id, "^delta_"),
    cell_line = stringr::str_extract(cell_line, "^[^_]+"),
    cimic_cluster = factor(cimic_cluster)
  ) %>%
  select(sample_id, cell_line, cimic_cluster, all_of(pyroptosis_genes_present)) %>%
  pivot_longer(
    cols = all_of(pyroptosis_genes_present),
    names_to = "Gene",
    values_to = "Expression"
  )

# Box colors for clusters
cimic_cluster_COLORS <- c(
  "1" = "#1b9e77",
  "2" = "#d95f02"
)

# Point colors for cell lines
cell_line_colors <- c(
  "DU4475"   = "#1f78b4",
  "HCC1806"  = "#33a02c",
  "HCC1395"  = "#e31a1c",
  "Hs578t"   = "#ff7f00",
  "MDAMB157" = "#6a3d9a",
  "MDAMB231" = "#b15928",
  "MDAMB468" = "#a6cee3",
  "BT549"    = "#fb9a99",
  "HCC38"    = "#b2df8a"
)

# p-value positions per facet
pval_df <- gene_df %>%
  group_by(Gene) %>%
  summarise(
    y.position = max(Expression, na.rm = TRUE) * 1.12,
    .groups = "drop"
  )

# Plot
p_pyroptosis <- ggplot(gene_df, aes(x = cimic_cluster, y = Expression)) +
  geom_boxplot(
    aes(fill = cimic_cluster),
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    aes(color = cell_line),
    width = 0.15,
    size = 2,
    alpha = 0.7
  ) +
  stat_compare_means(
    comparisons = list(c("1", "2")),
    method = "wilcox.test",
    label = "p.format",
    label.y = pval_df$y.position
  ) +
  scale_fill_manual(values = cimic_cluster_COLORS) +
  scale_color_manual(values = cell_line_colors) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(
    title = "Pyroptosis genes by CIMiC cluster",
    y = "delta(log2(TPM+1))",
    x = NULL,
    color = "Cell line"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    axis.line = element_line(colour = "black", linewidth = 1.5)
  )

print(p_pyroptosis)


```






# checking out individual genes for ER stress


```{r}


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# ER stress / UPR genes
er_stress_genes <- c(
  "DDIT3", "ATF4", "ATF3", "XBP1", "ERN1", "EIF2AK3", "HSPA5",
  "DNAJB9", "PPP1R15A", "EDEM1", "HERPUD1", "GADD45A", "TRIB3", 
  "EIF2AK1"
)

# Keep only genes present
er_stress_genes_present <- intersect(er_stress_genes, colnames(delta_subset_wide))

# Long format
gene_df <- delta_subset_wide %>%
  mutate(
    sample_id = as.character(sample_id),
    cell_line = stringr::str_remove(sample_id, "^delta_"),
    cell_line = stringr::str_extract(cell_line, "^[^_]+"),
    cimic_cluster = factor(cimic_cluster)
  ) %>%
  select(sample_id, cell_line, cimic_cluster, all_of(er_stress_genes_present)) %>%
  pivot_longer(
    cols = all_of(er_stress_genes_present),
    names_to = "Gene",
    values_to = "Expression"
  )

# Box colors for clusters
cimic_cluster_COLORS <- c(
  "1" = "#1b9e77",
  "2" = "#d95f02"
)

# Point colors for cell lines
cell_line_colors <- c(
  "DU4475"   = "#1f78b4",
  "HCC1806"  = "#33a02c",
  "HCC1395"  = "#e31a1c",
  "Hs578t"   = "#ff7f00",
  "MDAMB157" = "#6a3d9a",
  "MDAMB231" = "#b15928",
  "MDAMB468" = "#a6cee3",
  "BT549"    = "#fb9a99",
  "HCC38"    = "#b2df8a"
)

# p-value positions per facet
pval_df <- gene_df %>%
  group_by(Gene) %>%
  summarise(
    y.position = max(Expression, na.rm = TRUE) * 1.12,
    .groups = "drop"
  )

# Plot
p_er_stress <- ggplot(gene_df, aes(x = cimic_cluster, y = Expression)) +
  geom_boxplot(
    aes(fill = cimic_cluster),
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    aes(color = cell_line),
    width = 0.15,
    size = 2,
    alpha = 0.7
  ) +
  stat_compare_means(
    comparisons = list(c("1", "2")),
    method = "wilcox.test",
    label = "p.format",
    label.y = pval_df$y.position
  ) +
  scale_fill_manual(values = cimic_cluster_COLORS) +
  scale_color_manual(values = cell_line_colors) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(
    title = "ER stress / UPR genes by CIMiC cluster",
    y = "delta(log2(TPM+1))",
    x = NULL,
    color = "Cell line"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    axis.line = element_line(colour = "black", linewidth = 1.5)
  )

print(p_er_stress)
```







# checking individual genes for proteostasis 

```{r}



library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# Proteostasis / protein quality control genes
proteostasis_genes <- c(
  # Cytosolic chaperones
  "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA8", "HSPB1",
  # Co-chaperones
  "DNAJA1", "DNAJB1", "BAG3", "STIP1",
  # ER proteostasis / folding
  "HSPA5", "CANX", "CALR", "PDIA4", "PDIA6", "ERP44",
  # ERAD / misfolded protein handling
  "HERPUD1", "EDEM1", "SEL1L", "SYVN1",
  # Proteasome / ubiquitin stress
  "PSMB5", "PSMD14", "UBC", "SQSTM1",
  # Autophagy / lysosome
  "MAP1LC3B", "ATG5", "ATG7", "BECN1", "LAMP1", "CTSD",
  # Redox / proteotoxic stress buffering
  "HMOX1", "NQO1", "TXN", "TXNRD1", "GCLC"
)

# Keep only genes present
proteostasis_genes_present <- intersect(proteostasis_genes, colnames(delta_subset_wide))

# Long format
gene_df <- delta_subset_wide %>%
  mutate(
    sample_id = as.character(sample_id),
    cell_line = stringr::str_remove(sample_id, "^delta_"),
    cell_line = stringr::str_extract(cell_line, "^[^_]+"),
    cimic_cluster = factor(cimic_cluster)
  ) %>%
  select(sample_id, cell_line, cimic_cluster, all_of(proteostasis_genes_present)) %>%
  pivot_longer(
    cols = all_of(proteostasis_genes_present),
    names_to = "Gene",
    values_to = "Expression"
  )

# Box colors for clusters
cimic_cluster_COLORS <- c(
  "1" = "#1b9e77",
  "2" = "#d95f02"
)

# Point colors for cell lines
cell_line_colors <- c(
  "DU4475"   = "#1f78b4",
  "HCC1806"  = "#33a02c",
  "HCC1395"  = "#e31a1c",
  "Hs578t"   = "#ff7f00",
  "MDAMB157" = "#6a3d9a",
  "MDAMB231" = "#b15928",
  "MDAMB468" = "#a6cee3",
  "BT549"    = "#fb9a99",
  "HCC38"    = "#b2df8a"
)

# p-value positions per facet
pval_df <- gene_df %>%
  group_by(Gene) %>%
  summarise(
    y.position = max(Expression, na.rm = TRUE) * 1.12,
    .groups = "drop"
  )

# Plot
p_proteostasis <- ggplot(gene_df, aes(x = cimic_cluster, y = Expression)) +
  geom_boxplot(
    aes(fill = cimic_cluster),
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    aes(color = cell_line),
    width = 0.15,
    size = 2,
    alpha = 0.7
  ) +
  stat_compare_means(
    comparisons = list(c("1", "2")),
    method = "wilcox.test",
    label = "p.format",
    label.y = pval_df$y.position
  ) +
  scale_fill_manual(values = cimic_cluster_COLORS) +
  scale_color_manual(values = cell_line_colors) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(
    title = "Proteostasis / protein quality control genes by CIMiC cluster",
    y = "delta(log2(TPM+1))",
    x = NULL,
    color = "Cell line"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    axis.line = element_line(colour = "black", linewidth = 1.5)
  )

print(p_proteostasis)

```



# checking individual genes for DNA sensing




```{r}



library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# Cytosolic DNA sensing + downstream signaling genes
cytosolic_dna_genes <- c(
  "CGAS", "MB21D1", "TMEM173", "STING1", "TBK1", "IKBKE", "IRF3", "IRF7",
  "IFI16", "AIM2", "ZBP1", "TREX1",
  "IFNB1", "CXCL10", "CXCL11", "CCL5", "ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "TLR9"
)

# Keep only genes present
cytosolic_dna_genes_present <- intersect(cytosolic_dna_genes, colnames(delta_subset_wide))

# Long format
gene_df <- delta_subset_wide %>%
  mutate(
    sample_id = as.character(sample_id),
    cell_line = stringr::str_remove(sample_id, "^delta_"),
    cell_line = stringr::str_extract(cell_line, "^[^_]+"),
    cimic_cluster = factor(cimic_cluster)
  ) %>%
  select(sample_id, cell_line, cimic_cluster, all_of(cytosolic_dna_genes_present)) %>%
  pivot_longer(
    cols = all_of(cytosolic_dna_genes_present),
    names_to = "Gene",
    values_to = "Expression"
  )

# Box colors for clusters
cimic_cluster_COLORS <- c(
  "1" = "#1b9e77",
  "2" = "#d95f02"
)

# Point colors for cell lines
cell_line_colors <- c(
  "DU4475"   = "#1f78b4",
  "HCC1806"  = "#33a02c",
  "HCC1395"  = "#e31a1c",
  "Hs578t"   = "#ff7f00",
  "MDAMB157" = "#6a3d9a",
  "MDAMB231" = "#b15928",
  "MDAMB468" = "#a6cee3",
  "BT549"    = "#fb9a99",
  "HCC38"    = "#b2df8a"
)

# p-value positions per facet
pval_df <- gene_df %>%
  group_by(Gene) %>%
  summarise(
    y.position = max(Expression, na.rm = TRUE) * 1.12,
    .groups = "drop"
  )

# Plot
p_cytosolic_dna <- ggplot(gene_df, aes(x = cimic_cluster, y = Expression)) +
  geom_boxplot(
    aes(fill = cimic_cluster),
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    aes(color = cell_line),
    width = 0.15,
    size = 2,
    alpha = 0.7
  ) +
  stat_compare_means(
    comparisons = list(c("1", "2")),
    method = "wilcox.test",
    label = "p.format",
    label.y = pval_df$y.position
  ) +
  scale_fill_manual(values = cimic_cluster_COLORS) +
  scale_color_manual(values = cell_line_colors) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(
    title = "Cytosolic DNA sensing and downstream signaling by CIMiC cluster",
    y = "delta(log2(TPM+1))",
    x = NULL,
    color = "Cell line"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    axis.line = element_line(colour = "black", linewidth = 1.5)
  )

print(p_cytosolic_dna)
```





# Create Master_induction df

```{r}

library(circlize)
library(ComplexHeatmap)
library(dunn.test)
library(colorRamp2)

delta_subset_wide_genes <- setdiff(
  colnames(delta_subset_wide),
  c("sample_id", "cimic_cluster", "cell_line")
)

# extract out optimal_k 
optimal_k <- length(unique(delta_subset_wide$cimic_cluster))

# Prepare a list to collect rows
# Predefine output matrix for efficiency
n_genes <- length(delta_subset_wide_genes)
optimal_k <- length(unique(delta_subset_wide$cimic_cluster))

# Setup column names: Mean, SD, SEM per cluster, then pairwise Dunn p-values
n_pairs <- optimal_k * (optimal_k - 1) / 2

if (optimal_k > 2){
pairs_mat <- combn(seq_len(optimal_k), 2)
pairwise_names <- apply(pairs_mat, 2, function(x) {
  sprintf("cluster%d_cluster%d_pval", x[1], x[2])
})
} else if (optimal_k == 2){
  
 pairwise_names <-  "cluster1_cluster2_pval"
}


cluster_names <- as.vector(rbind(
  paste0("cluster", seq_len(optimal_k), "_mean"),
  paste0("cluster", seq_len(optimal_k), "_sd"),
  paste0("cluster", seq_len(optimal_k), "_sem")
))

penetrance_names <- paste0("cluster", seq_len(optimal_k), "_penetrance")

out_mat <- matrix(
  NA,
  nrow = n_genes,
  ncol = 3 + length(cluster_names) + length(penetrance_names) + n_pairs
)

colnames(out_mat) <- c(
  "gene_id", "H_value", "p_value",
  cluster_names,
  penetrance_names,
  pairwise_names
)

# Loop
for (i in seq_len(n_genes)) {
  
  
  gene_name <- delta_subset_wide_genes[i]
  
  if(gene_name %in% colnames(delta_subset_wide)){
  gene_vals <- delta_subset_wide[[gene_name]]
  clusters  <- as.factor(delta_subset_wide$cimic_cluster)
  clusters  <- droplevels(clusters)  # just in case there are unused levels

  # Per-cluster stats
  cluster_means <- tapply(gene_vals, clusters, mean)
  cluster_sds   <- tapply(gene_vals, clusters, stats::sd)
  cluster_ns    <- tapply(gene_vals, clusters, function(x) sum(!is.na(x)))
  cluster_sems  <- cluster_sds / sqrt(cluster_ns)

  # Fill in full-length per-cluster vectors (in order 1:optimal_k)
  cluster_means_full <- rep(NA, optimal_k)
  cluster_sds_full   <- rep(NA, optimal_k)
  cluster_sems_full  <- rep(NA, optimal_k)
  idx <- as.integer(names(cluster_means))
  cluster_means_full[idx] <- as.numeric(cluster_means)
  cluster_sds_full[idx]   <- as.numeric(cluster_sds)
  cluster_sems_full[idx]  <- as.numeric(cluster_sems)
  
  # TRUE penetrance (% samples > 0) 
  cluster_penetrance <- tapply(gene_vals, clusters, function(x)
    mean(x > 0, na.rm = TRUE))
  
  # Fill full-length penetrance vector
  cluster_penetrance_full <- rep(NA, optimal_k)
  idx <- as.integer(names(cluster_penetrance))
  cluster_penetrance_full[idx] <- as.numeric(cluster_penetrance)
  
  # Interleave mean, sd, sem for each cluster
  combined_stats <- as.vector(rbind(cluster_means_full, cluster_sds_full, cluster_sems_full))

  # Initialize post-hoc p-value vector
  dunn_pvals <- rep(NA, n_pairs)

  # Statistical test selection
  if (nlevels(clusters) == 2) {
    # Two groups → Wilcoxon rank-sum
    test_res <- wilcox.test(gene_vals ~ clusters)
    kw_H_value <- as.numeric(test_res$statistic)
    kw_p_value <- as.numeric(test_res$p.value)
    # No Dunn’s needed for two groups
  } else if (nlevels(clusters) > 2) {
    # Kruskal–Wallis
    df_kw <- kruskal.test(gene_vals ~ clusters)
    kw_H_value <- as.numeric(df_kw$statistic)
    kw_p_value <- as.numeric(df_kw$p.value)

    # Dunn’s post-hoc
    dt <- dunn.test(gene_vals, clusters, kw = FALSE, table = FALSE, method = "bh")
    if (!is.null(dt$P.adjusted)) {
      dunn_pvals[1:length(dt$P.adjusted)] <- dt$P.adjusted
    }
  } else {
    # Less than two groups — no test possible
    kw_H_value <- NA
    kw_p_value <- NA
  }

  # Output row
  out_mat[i, ] <- c(
    gene_name,
    kw_H_value,
    kw_p_value,
    combined_stats,
    cluster_penetrance_full,
    dunn_pvals
  )
}
}

# colnames(out_mat) <- c("gene_id", "H_value", "p_value", 
#                        cluster_names, 
#                        pairwise_names)


master_induction_df <- as.data.frame(out_mat, stringsAsFactors = FALSE)

#if to remove post-hoc test if two group
 if (nlevels(delta_subset_wide$cimic_cluster) == 2) {
   
   master_induction_df <- master_induction_df %>% dplyr::select(-all_of(pairwise_names))
}

# Optionally convert numeric columns back to numeric:
num_cols <- setdiff(names(master_induction_df), "gene_id")
master_induction_df[num_cols] <- lapply(master_induction_df[num_cols], as.numeric)
master_induction_df$stat_test_adj_p <- p.adjust(master_induction_df$p_value, method = "BH")

master_induction_df <- master_induction_df[, c("gene_id", "H_value", "p_value", "stat_test_adj_p", setdiff(names(master_induction_df), c("gene_id", "H_value", "p_value", "stat_test_adj_p")))]

# Suppose optimal_k and master_induction_df are already defined
for (i in 1:(optimal_k - 1)) {
  for (j in (i + 1):optimal_k) {
    col_i_mean <- paste0("cluster", i, "_mean")
    col_j_mean <- paste0("cluster", j, "_mean")
    col_i_sd <- paste0("cluster", i, "_sd")
    col_j_sd <- paste0("cluster", j, "_sd")
    diff_col1 <- paste0("cluster", i, "_cluster", j, "_differential")
    diff_col2 <- paste0("cluster", j, "_cluster", i, "_differential")
    sd_diff_col1 <- paste0("cluster", i, "_cluster", j, "_sd")
    sd_diff_col2 <- paste0("cluster", j, "_cluster", i, "_sd")
    
    master_induction_df[[diff_col1]] <- as.numeric(master_induction_df[[col_i_mean]]) - as.numeric(master_induction_df[[col_j_mean]])
    master_induction_df[[diff_col2]] <- as.numeric(master_induction_df[[col_j_mean]]) - as.numeric(master_induction_df[[col_i_mean]])
    master_induction_df[[sd_diff_col1]] <- sqrt((as.numeric(master_induction_df[[col_i_sd]]))^2 + (as.numeric(master_induction_df[[col_j_sd]]))^2)
    master_induction_df[[sd_diff_col2]] <- master_induction_df[[sd_diff_col1]]
  }
}



```

# Volcano plot by cimic cluster

```{r}


# Volcano plot (modify as needed ) ###################################



library(dplyr)
library(EnhancedVolcano)


# Define the two clusters of interest
clusters_oi <- c("cluster1", "cluster2")
pval_cutoff <- 0.05
ind_cutoff <- 0.25

penetrance_target_cutoff <- 0.20   # ≥20% induced in target cluster
penetrance_other_cutoff  <- 0.20   # ≤20% induced in other cluster

# Build the column names dynamically
diff_col      <- paste0(clusters_oi[1], "_", clusters_oi[2], "_differential")
diff_col_rev  <- paste0(clusters_oi[2], "_", clusters_oi[1], "_differential")
pval_col      <- "stat_test_adj_p"

pen_target <- paste0(clusters_oi[1], "_penetrance")
pen_other  <- paste0(clusters_oi[2], "_penetrance")


n_genes <- 20

# 1. Define the pool of available genes from your filtered set
available_genes <- unique(gene_set_long$gene)

# 2. Get the Top 10 genes based on the MAGNITUDE of difference
# We take the absolute value so we catch the biggest movers in either direction
top_genes <- master_induction_df %>%
  # Only consider genes that are actually in your filtered plot data
  dplyr::filter(gene_id %in% available_genes) %>%
  # Only consider significant genes
  dplyr::filter(.data[[pval_col]] < pval_cutoff) %>%
  # Calculate absolute difference to find the "biggest" changes
  dplyr::mutate(abs_diff = abs(.data[[diff_col]])) %>%
  # Sort by magnitude
  dplyr::arrange(dplyr::desc(abs_diff)) %>%
  # Take the top 10
  dplyr::slice_head(n = n_genes) %>%
  dplyr::pull(gene_id)

# Alternative: 5 up, 5 down
top_up <- master_induction_df %>% 
  filter(gene_id %in% available_genes, .data[[pval_col]] < pval_cutoff) %>%
  arrange(desc(.data[[diff_col]])) %>% slice_head(n = n_genes) %>% pull(gene_id)

top_down <- master_induction_df %>% 
  filter(gene_id %in% available_genes, .data[[pval_col]] < pval_cutoff) %>%
  arrange(.data[[diff_col]]) %>% slice_head(n = n_genes) %>% pull(gene_id)

top_genes <- c(top_up, top_down)

# # Top 10 where cluster1 > cluster2
# top_up <- master_induction_df %>%
#   dplyr::filter(.data[[pval_col]] < pval_cutoff) %>%
#   dplyr::arrange(dplyr::desc(.data[[diff_col]])) %>%
#   dplyr::slice_head(n = n_genes) %>%
#   dplyr::pull(gene_id)
# 
# # Top 10 where cluster2 > cluster1
# top_down <- master_induction_df %>%
#   dplyr::filter(.data[[pval_col]] < pval_cutoff) %>%
#   dplyr::arrange(dplyr::desc(.data[[diff_col_rev]])) %>%
#   dplyr::slice_head(n = n_genes) %>%
#   dplyr::pull(gene_id)
# 
# # Combine
# top_genes <- union(top_up, top_down) 

# Volcano plot
ind_graph <- EnhancedVolcano(
  master_induction_df,
  lab = master_induction_df$gene_id,
  x = diff_col_rev,
  y = pval_col,
  xlab = bquote("Induction differential (" * Delta * " log"[2]*"(TPM+1))"),
  selectLab = top_genes,
  boxedLabels = TRUE,
  labSize = 6.0,
  labCol = 'black',
  labFace = 'bold',
  pCutoff = pval_cutoff,
  FCcutoff = ind_cutoff,
  pointSize = 2.0,
  colAlpha = 4/5,
  col = c('black', 'black', 'black', 'red'),
  legendPosition = 'none',
  legendLabSize = 16,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  border = 'full',
  borderColour = 'black',
  borderWidth = 1.0,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  title = paste0(clusters_oi[1], " vs ", clusters_oi[2], " induction"),
  max.overlaps = 30,
  maxoverlapsConnectors = Inf
) +
  ggplot2::coord_cartesian(xlim=c(-5, 5)) +
  ggplot2::scale_x_continuous(breaks = seq(-5, 5, 1)) +
  ggplot2::coord_cartesian(ylim = c(0, 8))

ind_graph



```




# ORA for cluster specific induction 

```{r}
# overrepresentation analyses 
  
library(msigdbr)
library(clusterProfiler)
library(openxlsx)
library(dplyr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)


library(paletteer)
library(enrichplot)
library(patchwork)
library(tidyverse)
  library(stringr)


# # Make plot
# library(patchwork)
# 
# # # Combine the plots using patchwork
# # combined_plot <- wrap_plots(ORA_induction_pic_list, ncol = length(1:optimal_k))  # Arrange in 2 columns
# #  print(combined_plot)
#   
#   
# 
  

ORA_induction_pic_list <- list()

cluster_cols <- paste0("cluster", 1:optimal_k, "_mean")
padj_cutoff <- 0.05
differential_threshold <- 0
induction_threshold <- 0.25
penetrance_threshold <- (0.20)
universe_genes <- master_induction_df$gene_id

# Generate all non-empty, non-full subsets
cluster_subsets <- unlist(
    lapply(1:(optimal_k - 1), function(m) combn(1:optimal_k, m, simplify = FALSE)),
    recursive = FALSE
)

for (subset in cluster_subsets) {
    
    cluster_name <- paste0("cluster_", paste(subset, collapse = "_"))
    cluster_gene_set_name <- paste0(cluster_name, "_specific_induced_genes")
    
    # Work on a temp copy
    df_tmp <- master_induction_df
    
    # Combined mean for subset
    target_means <- paste0("cluster", subset, "_mean")
    penetrance_name <- paste0("cluster", subset, "_penetrance")
    
    df_tmp <- df_tmp %>%
        rowwise() %>%
        mutate(combined_mean = mean(c_across(all_of(target_means)), na.rm = TRUE)) %>%
        ungroup()
    
    # Remaining clusters
    other_clusters <- setdiff(1:optimal_k, subset)
    other_means <- paste0("cluster", other_clusters, "_mean")
    
    # Differential columns: only subset vs other
    cluster_diff_cols <- c()
    for (pos in subset) {
        cluster_diff_cols <- c(
            cluster_diff_cols,
            paste0("cluster", pos, "_cluster", other_clusters, "_differential")
        )
    }
    cluster_diff_cols <- unique(cluster_diff_cols)
    
    # Filtering
    # Now add the penetrance filter
    filtered_genes <- df_tmp %>%
      filter(
        stat_test_adj_p < 0.05,
        combined_mean > induction_threshold,
        if_all(all_of(other_means), ~ . < 0),
        if_all(all_of(cluster_diff_cols), ~ . > differential_threshold),
        # Penetrance for all clusters in subset must exceed threshold
        if_all(all_of(paste0(
          "cluster", subset, "_penetrance"
        )), ~ . >= penetrance_threshold),
      )
    
    assign(cluster_gene_set_name, filtered_genes, envir = .GlobalEnv)
    message(cluster_gene_set_name, " in progress...")
    
    
    # ORA
    ORA_induction <- enrichGO(
        gene          = filtered_genes$gene_id,
        OrgDb         = org.Hs.eg.db,
        ont           = "ALL",
        pAdjustMethod = "BH",
        pvalueCutoff  = padj_cutoff,
        qvalueCutoff  = padj_cutoff,
        universe      = universe_genes,
        keyType       = "SYMBOL"
    )
    
    ORA_name <- paste0(cluster_gene_set_name, "_ORA_analysis")
    ORA_name_simplified <- paste0("simplified_", cluster_gene_set_name, "_ORA_analysis")
        
    ORA_simplified <- clusterProfiler::simplify(
      ORA_induction,
      cutoff = 0.7,   # similarity threshold (0.7 is common, range 0–1)
      by = "p.adjust", # keep the term with the lowest adjusted p-value
      select_fun = min,
      measure = "Wang" # semantic similarity measure; Wang is common for GO
    )

    
    if (!is.null(ORA_induction) && nrow(ORA_induction@result) > 0) {
        
        Check_ORA_induction <- rownames_to_column(ORA_induction@result)
        assign(ORA_name, Check_ORA_induction, envir = .GlobalEnv)
        
        simp_Check_ORA_induction <- rownames_to_column(ORA_simplified@result)
        assign(ORA_name_simplified, simp_Check_ORA_induction, envir = .GlobalEnv)
        
        
         n_count <- 5
        
        # Top terms
        top_per_ontology <- ORA_simplified@result %>%
            group_by(ONTOLOGY) %>%
            slice_min(order_by = p.adjust, n = n_count, with_ties = FALSE) %>%
            mutate(Description = str_wrap(Description, width = 20))
        


        
        
        p1 <- ggplot(top_per_ontology, aes(x = Description, y = Count)) +
    geom_col(aes(fill = p.adjust), color = "black") +
    coord_flip() +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
    labs(title = paste0(ORA_name, "\nTop ", n_count, " Enriched Terms per ONTOLOGY"),
         x = NULL, y = "Gene Count", fill = "Adj_p" ) +
    scale_fill_gradient(low = "red", high = "blue", trans = "reverse") +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text.y = element_text(size = 18, angle = 0, face = "bold"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18)
        #plot.title = element_text(hjust = 0.5)
    )
    #       p1 <- ggplot(top_per_ontology, aes(x = Description, y = Count)) +
    # geom_col(aes(fill = p.adjust), color = "black") +
    # geom_text(
    #     aes(label = Description, y = 0),   # put label halfway inside the bar
    #     color = "black", size = 4, fontface = "bold", hjust = 0
    # ) +
    # coord_flip() +
    # facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
    # labs(title = paste0(ORA_name, "\nTop ", n_count, " Enriched Terms per ONTOLOGY"),
    #      x = NULL, y = "Gene Count", fill = "Adjusted p-value") +
    # scale_fill_gradient(low = "red", high = "lightcyan", trans = "reverse") +
    # theme_minimal(base_size = 14) +
    # theme(
    #     panel.grid.major = element_blank(),
    #     panel.grid.minor = element_blank(),
    #     panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    #     strip.background = element_rect(fill = "white", colour = "black"),
    #     strip.text.y = element_text(size = 18, angle = 0, face = "bold"),
    #     axis.text.y = element_blank(),   # remove y-axis text (since it’s inside bars now)
    #     axis.text.x = element_text(size = 16, colour = "black"),
    #     plot.title = element_text(hjust = 0.5)
    # )

          
        ORA_induction_pic_list[[ORA_name]] <- p1
        
    } else if (is.null(ORA_induction))  {
  
    assign(ORA_name, NULL, envir = .GlobalEnv)
    print(paste0(ORA_name, " is empty will put as NULL in environment"))
    
}
}


# ploting simplifed ora (using BP and N = 10; changeable)

n_count <- 15

 # Top terms
 top_res_c1 <- simplified_cluster_1_specific_induced_genes_ORA_analysis %>%
   dplyr::filter(ONTOLOGY == "BP") %>%
   slice_min(order_by = p.adjust,
             n = n_count,
             with_ties = FALSE) %>%
   mutate(Description = str_wrap(Description, width = 50))

 
 
 p1_c1_simplifed_ora <- ggplot(top_res_c1, aes(x = Description, y = Count)) +
   geom_col(aes(fill = p.adjust), color = "black") +
   coord_flip() +
   facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
   labs(
     title = paste0("cluster 1", "\nTop ", n_count, " Enriched Terms per ONTOLOGY"),
     x = NULL,
     y = "Gene Count",
     fill = "Adj_p"
   ) +
   scale_fill_gradient(low = "red",
                       high = "blue",
                       trans = "reverse") +
   theme_minimal(base_size = 14) +
   theme(
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_rect(
       colour = "black",
       fill = NA,
       linewidth = 1
     ),
     strip.background = element_rect(fill = "white", colour = "black"),
     strip.text.y = element_text(size = 18, angle = 0, face = "bold"),
     axis.text.y = element_text(size = 16, colour = "black"),
     axis.text.x = element_text(size = 18, colour = "black"),
     legend.text = element_text(size = 18),
     legend.title = element_text(size = 18)
     #plot.title = element_text(hjust = 0.5)
   )

 
 
  
  # Top terms
 top_res_c2 <- simplified_cluster_2_specific_induced_genes_ORA_analysis %>%
   dplyr::filter(ONTOLOGY == "BP") %>%
   slice_min(order_by = p.adjust,
             n = n_count,
             with_ties = FALSE) %>%
   mutate(Description = str_wrap(Description, width = 50))
 

 
  p2_c2_simplifed_ora <- ggplot(top_res_c2, aes(x = Description, y = Count)) +
   geom_col(aes(fill = p.adjust), color = "black") +
   coord_flip() +
   facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
   labs(
     title = paste0("cluster 2", "\nTop ", n_count, " Enriched Terms per ONTOLOGY"),
     x = NULL,
     y = "Gene Count",
     fill = "Adj_p"
   ) +
   scale_fill_gradient(low = "red",
                       high = "blue",
                       trans = "reverse") +
   theme_minimal(base_size = 14) +
   theme(
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_rect(
       colour = "black",
       fill = NA,
       linewidth = 1
     ),
     strip.background = element_rect(fill = "white", colour = "black"),
     strip.text.y = element_text(size = 18, angle = 0, face = "bold"),
     axis.text.y = element_text(size = 16, colour = "black"),
     axis.text.x = element_text(size = 18, colour = "black"),
     legend.text = element_text(size = 18),
     legend.title = element_text(size = 18)
     #plot.title = element_text(hjust = 0.5)
   )

 
# Keyword pull out --------------------------

# Define keywords of interest
keywords <- c("lymphocyte", "chemotaxis", "immunity", "leukocyte", "killing", "immune", "adaptive", "toll-like")

top_per_oncology_filtered <- simplified_cluster_1_specific_induced_genes_ORA_analysis

# Build regex (case-insensitive)
pattern <- paste(keywords, collapse = "|")

# Filter descriptions that match keywords
top_per_ontology_filtered <- top_per_oncology_filtered %>%
  dplyr::filter(str_detect(tolower(Description), pattern))

# Plot
p1 <- ggplot(top_per_ontology_filtered, aes(x = Description, y = Count)) +
  geom_col(aes(fill = p.adjust), color = "black") +
  coord_flip() +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
  labs(
    x = NULL, y = "Gene Count", fill = "Adj_p"
  ) +
  scale_fill_gradient(low = "red", high = "blue", trans = "reverse") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text.y = element_text(size = 18, angle = 0, face = "bold"),
    axis.text.y = element_text(size = 16, colour = "black"),
    axis.text.x = element_text(size = 18, colour = "black"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  )

```






# ssGSEA for CIMIC

```{r}
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(matrixStats)

# ---- 1. Convert to numeric matrix ----
tnbc_mat <- tnbc_raw %>%
  remove_rownames() %>%
  column_to_rownames("gene_id") %>%
  as.matrix() %>%
  { matrix(as.numeric(.), nrow = nrow(.), dimnames = dimnames(.)) }

# ---- 2. Function to extract DMSO + one drug together ----
get_drug_subset_matrix <- function(mat, drug_name, control_name = "DMSO") {
  
  keep_cols <- grep(
    paste0("_((", control_name, ")|(", drug_name, "))_"),
    colnames(mat),
    value = TRUE
  )
  
  mat[, keep_cols, drop = FALSE]
}

# ---- 3. Function to run ssGSEA and compute matched delta ----
compute_delta_ssgsea <- function(mat_subset, gene_sets, drug_name, control_name = "DMSO",
                                 alpha = 0.25, scale = TRUE, norm = FALSE, single = TRUE) {
  
  # Run ssGSEA on all samples in this drug-specific subset
  es <- ssgsea(
    X = mat_subset,
    gene_sets = gene_sets,
    alpha = alpha,
    scale = scale,
    norm = norm,
    single = single
  )
  
  es_df <- as.data.frame(es)
  
  # Long format
  es_long <- es_df %>%
    rownames_to_column("GeneSet") %>%
    pivot_longer(
      cols = -GeneSet,
      names_to = "sample_id",
      values_to = "ssGSEA_score"
    ) %>%
    mutate(
      sample_id = as.character(sample_id),
      cell_line = str_extract(sample_id, "^[^_]+"),
      treatment = str_match(sample_id, "^[^_]+_([^_]+)")[,2],
      replicate = str_match(sample_id, "_(R[0-9]+)$")[,2]
    )
  
  # Split into control and drug
  control_df <- es_long %>%
    filter(treatment == control_name) %>%
    dplyr::select(GeneSet, cell_line, replicate, control_score = ssGSEA_score)
  
  drug_df <- es_long %>%
    filter(treatment == drug_name) %>%
    dplyr::select(GeneSet, sample_id, cell_line, treatment, replicate, drug_score = ssGSEA_score)
  
  # Matched delta
  delta_df <- drug_df %>%
    left_join(control_df, by = c("GeneSet", "cell_line", "replicate")) %>%
    mutate(delta_ssgsea = drug_score - control_score)
  
  # Row-wise z-score across samples within each GeneSet
  delta_df <- delta_df %>%
    group_by(GeneSet) %>%
    mutate(delta_ssgsea_z = as.numeric(scale(delta_ssgsea))) %>%
    ungroup()
  
  # Wide outputs
  delta_raw_wide <- delta_df %>%
    dplyr::select(GeneSet, sample_id, delta_ssgsea) %>%
    pivot_wider(names_from = sample_id, values_from = delta_ssgsea)
  
  delta_z_wide <- delta_df %>%
    dplyr::select(GeneSet, sample_id, delta_ssgsea_z) %>%
    pivot_wider(names_from = sample_id, values_from = delta_ssgsea_z)
  
  list(
    delta_long = delta_df,
    delta_raw = delta_raw_wide,
    delta_zrow = delta_z_wide,
    all_es = es_df
  )
}

# ---- 4. EPI vs DMSO ----
EPI_mat <- get_drug_subset_matrix(tnbc_mat, drug_name = "EPI", control_name = "DMSO")

ssgsea_EPI_delta <- compute_delta_ssgsea(
  mat_subset = EPI_mat,
  gene_sets = all_gene_sets,
  drug_name = "EPI",
  control_name = "DMSO",
  alpha = 0.25,
  scale = TRUE,
  norm = FALSE,
  single = TRUE
)

delta_EPI_long  <- ssgsea_EPI_delta$delta_long
delta_EPI_raw   <- ssgsea_EPI_delta$delta_raw
delta_EPI_zrow  <- ssgsea_EPI_delta$delta_zrow
EPI_all_es      <- ssgsea_EPI_delta$all_es

# ---- 5. DOX vs DMSO ----
DOX_mat <- get_drug_subset_matrix(tnbc_mat, drug_name = "DOX", control_name = "DMSO")

ssgsea_DOX_delta <- compute_delta_ssgsea(
  mat_subset = DOX_mat,
  gene_sets = all_gene_sets,
  drug_name = "DOX",
  control_name = "DMSO",
  alpha = 0.25,
  scale = TRUE,
  norm = FALSE,
  single = TRUE
)

delta_DOX_long  <- ssgsea_DOX_delta$delta_long
delta_DOX_raw   <- ssgsea_DOX_delta$delta_raw
delta_DOX_zrow  <- ssgsea_DOX_delta$delta_zrow
DOX_all_es      <- ssgsea_DOX_delta$all_es



delta_DOX_long <- delta_DOX_long %>%
  mutate(cell_line = stringr::str_extract(sample_id, "^[^_]+")) %>%
  left_join(cluster_map, by = "cell_line")

delta_EPI_long <- delta_EPI_long %>%
  mutate(cell_line = stringr::str_extract(sample_id, "^[^_]+")) %>%
  left_join(cluster_map, by = "cell_line")

# cluster_map should contain:
# cell_line and cimic_cluster

cluster_map <- tibble::tibble(
  cell_line = c("DU4475","HCC1806","HCC1395","Hs578t",
                "MDAMB157","MDAMB231","MDAMB468","BT549","HCC38"),
  cimic_cluster = c(1,1,1,2,2,2,1,2,2)
)

# Join CIMIC clusters to delta_DOX_long
delta_DOX_long <- delta_DOX_long %>%
  left_join(cluster_map, by = "cell_line")

# Join CIMIC clusters to delta_EPI_long
delta_EPI_long <- delta_EPI_long %>%
  left_join(cluster_map, by = "cell_line")

```


# heatmap for dox and epi by CIMIC

```{r}

library(dplyr)
library(tidyr)
library(ggplot2)

plot_df <- delta_DOX_zrow %>%
  left_join(cluster_map, by = "cell_line") %>%
  mutate(
    CIMIC_Cluster = factor(cimic_cluster)
  )

# 1. Mean z-score by cluster and gene set
heatmap_df <- plot_df %>%
  group_by(GeneSet, CIMIC_Cluster) %>%
  summarise(
    mean_ES = mean(delta_ssgsea , na.rm = TRUE),
    .groups = "drop"
  )

# 2. Per-gene-set stats across clusters
stats_df <- plot_df %>%
  group_by(GeneSet) %>%
  summarise(
    p_value = tryCatch(wilcox.test(delta_ssgsea_z ~ CIMIC_Cluster)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# 3. Order gene sets by Cluster 2 - Cluster 1 difference
gene_order <- heatmap_df %>%
  tidyr::pivot_wider(names_from = CIMIC_Cluster, values_from = mean_ES) %>%
  mutate(diff = `2` - `1`) %>%
  arrange(diff) %>%
  pull(GeneSet)

heatmap_df <- heatmap_df %>%
  mutate(GeneSet = factor(GeneSet, levels = gene_order))

stats_df <- stats_df %>%
  mutate(GeneSet = factor(GeneSet, levels = gene_order))

# 4. Create significance annotation column
sig_df <- stats_df %>%
  transmute(
    GeneSet,
    CIMIC_Cluster = factor("sig", levels = c("1", "2", "sig")),
    mean_ES = NA_real_,
    label = as.character(p_label),
    p_label = as.character(p_label)
  )

# 5. Main heatmap values
heatmap_with_sig <- heatmap_df %>%
  transmute(
    GeneSet,
    CIMIC_Cluster = factor(as.character(CIMIC_Cluster), levels = c("1", "2", "sig")),
    mean_ES,
    label = as.character(round(mean_ES, 2)),
    p_label = NA_character_
  )

# 6. Combine
plot_df_final <- bind_rows(
  heatmap_with_sig,
  sig_df
) %>%
  mutate(
    GeneSet = factor(GeneSet, levels = gene_order),
    CIMIC_Cluster = factor(CIMIC_Cluster, levels = c("1", "2", "sig"))
  )

# 7. Plot
p_heatmap_stats <- ggplot(plot_df_final, aes(x = CIMIC_Cluster, y = GeneSet, fill = mean_ES)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    data = subset(plot_df_final, CIMIC_Cluster %in% c("1", "2")),
    aes(label = label),
    size = 5,
    fontface = "bold"
  ) +
  geom_text(
    data = subset(plot_df_final, CIMIC_Cluster == "sig"),
    aes(label = p_label),
    size = 7,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#313695",
    mid = "white",
    high = "#A50026",
    midpoint = 0,
    na.value = "white",
    name = "Mean/nZ-score"
  ) +
  scale_x_discrete(labels = c("1", "2", "")) +
  labs(
    x = "CIMIC Cluster",
    y = NULL,
    title = "ΔssGSEA z-score heatmap by CIMIC cluster"
  ) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
    axis.text.x = element_text(face = "bold", size = 16),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 18),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 1.5)
  )

print(p_heatmap_stats)
```



































# Looking at most correlated CIM genes per cluster


```{r}


library(dplyr)
library(purrr)
library(broom)
library(ggplot2)
library(tidyr)

# -----------------------------
# 1. Set up data
# -----------------------------
cor_df_input <- delta_subset_wide %>%
  mutate(
    cluster_numeric = as.numeric(cimic_cluster)
  )

# CIMiC gene list
cimic_genes <- unique(cimic_gene_set_long$gene)

# Metadata columns to exclude
meta_cols <- c("sample_id", "cimic_cluster", "cluster_numeric", "cell_line", "treatment")

# Keep only CIMiC genes that are actually present
gene_cols <- intersect(cimic_genes, setdiff(colnames(cor_df_input), meta_cols))

expr_mat <- cor_df_input %>%
  select(all_of(gene_cols))

# -----------------------------
# 2. Correlate each CIMiC gene with cluster
# -----------------------------
gene_cluster_cor_df <- map_dfr(
  colnames(expr_mat),
  function(gene) {
    test <- cor.test(
      cor_df_input[[gene]],
      cor_df_input$cluster_numeric,
      method = "pearson"
    )
    
    tibble(
      gene = gene,
      cor = unname(test$estimate),
      pval = test$p.value,
      direction = ifelse(test$estimate > 0, "Higher in Cluster 2", "Higher in Cluster 1")
    )
  }
)

gene_cluster_cor_df <- gene_cluster_cor_df %>%
  mutate(
    padj = p.adjust(pval, method = "fdr"),
    abs_cor = abs(cor)
  ) %>%
  arrange(padj)

# -----------------------------
# 3. Significant correlated CIMiC genes
# -----------------------------
sig_gene_cluster_cor_df <- gene_cluster_cor_df %>%
  filter(padj < 0.05)

sig_gene_cluster_cor_df_strict <- gene_cluster_cor_df %>%
  filter(padj < 0.05, abs_cor > 0.5)

# -----------------------------
# 4. Top genes for plotting
#    top 30 in each direction
# -----------------------------
top_n_genes <- 30

plot_df <- gene_cluster_cor_df %>%
  filter(padj < 0.05) %>%
  group_by(direction) %>%
  slice_max(order_by = abs_cor, n = top_n_genes, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    gene = factor(gene, levels = rev(unique(gene)))
  )

# -----------------------------
# 5. Plot
# -----------------------------
p_gene_cluster_cor <- ggplot(plot_df, aes(x = gene, y = abs_cor, fill = direction)) +
  geom_col(width = 0.7) +
  coord_flip() +
  geom_hline(
    yintercept = 0.3,
    linetype = "longdash",
    color = "black",
    linewidth = 1
  ) +
  scale_fill_manual(
    values = c(
      "Higher in Cluster 2" = "firebrick",
      "Higher in Cluster 1" = "steelblue"
    )
  ) +
  facet_wrap(~direction, scales = "free_y") +
  labs(
    x = "Gene",
    y = "Absolute correlation with CIMiC cluster",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    legend.position = "none"
  )

print(p_gene_cluster_cor)

# -----------------------------
# 6. Extract positive/negative gene lists
# -----------------------------
positive_genes <- gene_cluster_cor_df %>%
  filter(padj < 0.05, cor > 0.5) %>%
  arrange(desc(cor)) %>%
  pull(gene)

negative_genes <- gene_cluster_cor_df %>%
  filter(padj < 0.05, cor < -0.5) %>%
  arrange(cor) %>%
  pull(gene)

length(positive_genes)
length(negative_genes)

head(positive_genes, 20)
head(negative_genes, 20)

```


# ORA analysis on correlated genes


```{r}

# ORA ANALYSIS ############################

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(dplyr)
library(ggplot2)
library(stringr)

# -----------------------------------------
# 1. Gene universe from your current dataset
# -----------------------------------------
meta_cols <- c("sample_id", "cimic_cluster", "cluster_numeric", "cell_line", "treatment")
universe_genes <- intersect(
  setdiff(colnames(cor_df_input), meta_cols),
  keys(org.Hs.eg.db, keytype = "SYMBOL")
)

padj_cutoff <- 0.05

# -----------------------------------------
# 2. ORA: Cluster 1-associated genes
#    negative_genes = higher in Cluster 1
# -----------------------------------------
C1_associated_ORA <- enrichGO(
  gene          = negative_genes,
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = padj_cutoff,
  qvalueCutoff  = padj_cutoff,
  universe      = universe_genes,
  keyType       = "SYMBOL"
)

# -----------------------------------------
# 3. ORA: Cluster 2-associated genes
#    positive_genes = higher in Cluster 2
# -----------------------------------------
C2_associated_ORA <- enrichGO(
  gene          = positive_genes,
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = padj_cutoff,
  qvalueCutoff  = padj_cutoff,
  universe      = universe_genes,
  keyType       = "SYMBOL"
)

# -----------------------------------------
# 4. Convert to Entrez for Reactome
# -----------------------------------------
universe_genes_entrez <- bitr(
  universe_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

negative_genes_entrez <- bitr(
  negative_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

positive_genes_entrez <- bitr(
  positive_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# -----------------------------------------
# 5. Reactome enrichment
# -----------------------------------------
C1_associated_reactome <- enrichPathway(
  gene          = negative_genes_entrez$ENTREZID,
  organism      = "human",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.2,
  universe      = universe_genes_entrez$ENTREZID,
  minGSSize     = 10,
  maxGSSize     = 500,
  readable      = FALSE
)

C2_associated_reactome <- enrichPathway(
  gene          = positive_genes_entrez$ENTREZID,
  organism      = "human",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.2,
  universe      = universe_genes_entrez$ENTREZID,
  minGSSize     = 10,
  maxGSSize     = 500,
  readable      = FALSE
)

# -----------------------------------------
# 6. Simplify GO ORA results
# -----------------------------------------
C1_associated_ORA_simplified <- clusterProfiler::simplify(
  C1_associated_ORA,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)

C2_associated_ORA_simplified <- clusterProfiler::simplify(
  C2_associated_ORA,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)

# -----------------------------------------
# 7. Plot top GO BP terms
# -----------------------------------------
n_count <- 20

top_res_C1 <- C1_associated_ORA_simplified@result %>%
  dplyr::filter(ONTOLOGY == "BP") %>%
  slice_min(order_by = p.adjust, n = n_count, with_ties = FALSE) %>%
  mutate(Description = str_wrap(Description, width = 50))

p_C1_ora <- ggplot(top_res_C1, aes(x = Description, y = Count)) +
  geom_col(aes(fill = p.adjust), color = "black") +
  coord_flip() +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
  labs(
    title = paste0("C1_associated", "/nTop ", n_count, " Enriched Terms per ONTOLOGY"),
    x = NULL,
    y = "Gene Count",
    fill = "Adj_p"
  ) +
  scale_fill_gradient(
    low = "red",
    high = "blue",
    trans = "reverse"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text.y = element_text(size = 18, angle = 0, face = "bold"),
    axis.text.y = element_text(size = 16, colour = "black"),
    axis.text.x = element_text(size = 18, colour = "black"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  )

print(p_C1_ora)

top_res_C2 <- C2_associated_ORA_simplified@result %>%
  dplyr::filter(ONTOLOGY == "BP") %>%
  slice_min(order_by = p.adjust, n = n_count, with_ties = FALSE) %>%
  mutate(Description = str_wrap(Description, width = 50))

p_C2_ora <- ggplot(top_res_C2, aes(x = Description, y = Count)) +
  geom_col(aes(fill = p.adjust), color = "black") +
  coord_flip() +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
  labs(
    title = paste0("C2_associated", "/nTop ", n_count, " Enriched Terms per ONTOLOGY"),
    x = NULL,
    y = "Gene Count",
    fill = "Adj_p"
  ) +
  scale_fill_gradient(
    low = "red",
    high = "blue",
    trans = "reverse"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text.y = element_text(size = 18, angle = 0, face = "bold"),
    axis.text.y = element_text(size = 16, colour = "black"),
    axis.text.x = element_text(size = 18, colour = "black"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  )

print(p_C2_ora)



```





