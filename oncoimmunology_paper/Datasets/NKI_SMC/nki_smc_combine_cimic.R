#' Run ONE CIMIC-LIMMA (moderated linear model feature selection).
#' Standalone - run this by itself. 
#' -------------------------------------------------------------------------
script_dir <- "D:/OneDrive - University of Florida/Gbadamosi Lab/Mohammed Gbadamosi/Gbadamosi_Lab_GitHub/CIMIC/cimic_releases/v1"
input_csv  <- "D:/OneDrive - University of Florida/Gbadamosi Lab/Mohammed Gbadamosi/Gbadamosi_Lab_GitHub/CIMIC/oncoimmunology_paper/Datasets/NKI_SMC/nki_smc_combine_initial_clustering_mat.csv"
output_dir <- "D:/OneDrive - University of Florida/Gbadamosi Lab/Mohammed Gbadamosi/Gbadamosi_Lab_GitHub/CIMIC/oncoimmunology_paper/Results/nki_smc_cimic_results"
variant_script <- "CIMIC_1.0.0.R"




# ---- cimic settings ----
clustering_alg    <- "hc"                                   # "hc" or "km"
max_k             <- 5
CCP_iter          <- 5000                                   # lower (e.g. 500) for a quick test
adj_pval_thresh   <- 0.05
max_pipeline_iter <- 50
seedval           <- 2026L
filter_approach   <- c("app_one", "app_two", "app_three")

# ---- load input (samples in rows, genes in columns; 1st column = sample ids) ----
clustering_matrix <- as.matrix(read.csv(input_csv, row.names = 1, check.names = FALSE))
storage.mode(clustering_matrix) <- "double"
cat(sprintf("input: %d samples x %d genes\n", nrow(clustering_matrix), ncol(clustering_matrix)))


source(file.path(script_dir, variant_script))

res <- CIMIC(
  clustering_matrix = clustering_matrix,
  # One sample per patient (36 unique ids, NKI numeric + SMC OB_), so there is no
  # replicate structure to model: every sample is independent. With NULL the fit
  # is numerically identical to plain unblocked limma (CIMIC_limma_1.0.0.R).
  block_metadata = NULL,
  all_gene_sets = NULL,
  clustering_alg = clustering_alg, max_k = max_k, CCP_iter = CCP_iter,
  adj_pval_thresh = adj_pval_thresh, max_pipeline_iter = max_pipeline_iter,
  seed = seedval, filter_approach = filter_approach, verbose = TRUE, working_dir = output_dir
)
saveRDS(res, file.path(output_dir, paste0("CIMIC_result", ".rds")))
cat("Done results in ->", output_dir, "\n")
