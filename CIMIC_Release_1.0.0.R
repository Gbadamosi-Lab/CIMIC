#' CIMIC pipeline – iteratively identify important genes
#' using PacMAP dimensionality reduction, k‑means (or hierarchical)
#' clustering, consensus‑cluster validation and non‑parametric feature selection.
#'
#' 
#' 
#' 
#' 
#'#### Dependency Management ####


# ---- Dependency management: Install & load needed packages ----
needed_pkgs <- c(
  "umap",
  "cluster",
  "factoextra",
  "dplyr",
  "stringr",
  "magrittr",
  "tibble",
  "rstatix",
  "coin",
  "RColorBrewer",
  "ggforce",
  "grid",
  "gridExtra",
  "corrplot",
  "concaveman",
  "paran",
  "reticulate",
  "utils",
  "plotly",
  "scales",
  "withr",
  "ggrepel",
  "ggplot2",
  "tidyr",
  "msigdbr",
  "jsonlite"
  
)

bioc_pkgs <- c(
  "ConsensusClusterPlus",
  "ComplexHeatmap"
)

# Install missing CRAN packages
new_pkgs <- needed_pkgs[!(needed_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

# Install missing Bioconductor packages
new_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[, "Package"])]
if (length(new_bioc)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(new_bioc, ask = FALSE)
}

# Now load all
all_pkgs <- c(needed_pkgs, bioc_pkgs)
suppressPackageStartupMessages({
  lapply(all_pkgs, require, character.only = TRUE)
})




# ---- Function Purpose and Definitions  ----
#' Chemoimmunomodulation induction classifier
#' Iterative, per-gene-set marker selection and clustering pipeline
#' @param all_gene_sets List of gene sets (character vectors) **or** a character string specifying an MSigDB category, or a data frame returned by `msigdbr::msigdbr()`.
#'   If a character string is supplied, the function will query MSigDB via `msigdbr` for the given category (e.g., "H" for Hallmark, "C2" for curated gene sets) and build a named list of gene vectors.
#'   If a data frame with columns `gs_name` and `gene_symbol` is supplied, it will be converted to the required list format.
#'   Otherwise, a named list of character vectors is expected.
#' @param clustering_matrix Numeric matrix, genes in columns, samples in rows.
#' @param adj_pval_thresh FDR threshold for gene selection.
#' @param max_pipeline_iter Maximum number of refinement iterations per gene set.
#' @param max_k Maximum number of clusters to try.
#' @param clustering_alg Algorithm for clustering ("hc", etc.).
#' @param verbose Logical, print progress messages?
#' @param working_dir working_directory where CIM_states_results_folder will be stored
#' @param clustering_metrics Character vector specifying one or more cluster quality metrics to use for ranking cluster solutions.
#'        Valid options are any combination of "pac", "silhouette", and "item_cluster_consensus".
#'        The overall ranking is computed as the sum of the ranks for the selected metrics (lower is better).
#'        In the case of ties for the overall rank, the function defaults to using the last metric specified in
#'        \code{clustering_metrics} as a tiebreaker (selecting maximum value for "silhouette" or "item_cluster_consensus", minimum value for "pac").
#'        
#'        
#' @return List with:
#'   - $iterated_by_gene_sets: 
#'       (character vector) Final union of important genes selected across all gene sets.
#'   - $iterated_over_all_genes: 
#'       (character vector) Final gene set when using all genes for selection (if performed).
#'   - $iterated_use_final_gene_set_from_approach_one: 
#'       (character vector) Final gene set from the first (or primary) approach/run.
#'   - $details_per_set_analysis_by_gene_sets: 
#'       (named list) For each gene set, a list with:
#'          * final_gene_set: genes retained for that set
#'          * stats_per_iter: marker stats (e.g., p-values) for each iteration
#'          * cluster_metrics_per_iter: (if implemented) clustering metrics
#'          * optimal_k_per_iter: (if implemented) the chosen number of clusters per iteration
#'          * cluster_stability_per_iter: (if implemented) cluster stability scores per iteration
#'          * wide_cluster_stability_per_iter: (if implemented) stability scores in wide format
#'   - $details_per_set_analysis_by_all_genes: 
#'       (list) Per-iteration stats and cluster metrics when using all genes at once.
#'   - $details_per_set_analysis_by_gene_set_then_all_genes: 
#'       (list) Optionally, per-iteration stats for a two-stage process (first by gene set, then refined over all genes).
#'   - $iter_log: 
#'       (data.frame) Log of number of genes before/after each refinement iteration, by gene set.
#'       
#'       
#'
#' 


# ---- Helper Functions  ----
# function to classify stability
# ConsensusClusterPlus does not have a predefined threshold
# using here the qualitative scale for analytic hierarchy process (AHP) consensus indicator (DOI: 10.13140/RG.2.2.28745.83045)
get_stability_category <- function(score) {
  if (is.na(score)) return(NA)
  else if (score < 0.50) return("Very Low (<50%)")
  else if (score < 0.625) return("Low (50-62.5%)")
  else if (score < 0.75) return("Moderate (62.5-75%)")
  else if (score < 0.875) return("High (75-87.5%)")
  else return("Very High (>=87.5%)")
}


# --------------------------------------------------------------
# 0️⃣ Load pacmap  ------------------------------------------------
# --------------------------------------------------------------
if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate")
}
library(reticulate)

# ---- PaCMAP Loading ----

setup_conda <- function(env_name = "pacmap_env") {


  # If pacmap already loaded, skip conda setup
  if (exists("pacmap")) {
    message("pacmap already present; skipping conda setup.")
      pacmap <- import("pacmap")
  return(pacmap)
  }

  
  message("Trying to set up PaCMAP...")


# Clear out old jsons
cache_file <- file.path(Sys.getenv("HOME"),
                        ".local", "share", "r-reticulate", "conda_list.json")
if (file.exists(cache_file)) file.remove(cache_file)

  ## 1️⃣  Honor an explicit RETICULATE_CONDA value
  cp <- Sys.getenv("RETICULATE_CONDA", unset = NA)

  if (!is.na(cp) && file.exists(cp)) {
    # The conda executable is already known – we can continue without prompting.
    # Nothing else is needed here; the later code will use `cp`.
  } else {
    ## 2️⃣  Try reticulate's built‑in discovery, but catch any error
    cp <- tryCatch(
      reticulate::conda_binary(),
      error = function(e) NULL
    )
  }

 ## -----------------------------------------------------------------
  ## 3️⃣  If we still don't have a valid conda binary, fall back to the
  ##     interactive prompt (max three attempts)
  ## -----------------------------------------------------------------
  
  if (is.null(cp) || !file.exists(cp)) {
  max_tries <- 3
  install_fallback <- FALSE # <-- flag that tells us whether to install Miniconda
  for (i in seq_len(max_tries)) {
    cat("\n=== Conda not detected ===\n")
    cat(
      "Please type the **full path** to your conda executable (WITH NO QUOTES), e.g.\n"
    )
    cat("  C:/Users/you/AppData/Local/anaconda3/Scripts/conda.exe\n")
    cat("  if copy and paste for path isn't working try ctrl+shift+v (Windows users)n")
    cat("Or type **5** to let the script install Miniconda automatically.\n")
    message("\n[Info] To avoid this interactive prompt in the future, set the environment variable RETICULATE_CONDA once, e.g.:")
  message('Sys.setenv(RETICULATE_CONDA = "C:/Users/user/AppData/Local/anaconda3/Scripts/conda.exe")')
  message("You can place that line in your .Rprofile or run it before calling the pipeline.")
    user_input <- readline(prompt = "conda path (or 5): ")

    # ---- 5 = install Miniconda ------------------------------------
    if (trimws(user_input) == "5" || nzchar(user_input) == FALSE) {
      # break out of the loop – we will install Miniconda below
      install_fallback <- TRUE
      break
    }

    # ---- User supplied a path – validate it -----------------------
    if (file.exists(user_input)) {
      # **STORE THE PATH IN `cp`** so the later code can use it
      cp <- normalizePath(user_input, winslash = "/", mustWork = TRUE)
      Sys.setenv(RETICULATE_CONDA = cp) # make it visible to reticulate
      message("User provided file exists and contains valid conda…")
    break
    } else {
      message("The file you entered does not exist: ", user_input)
      if (i < max_tries) {
        message("Please try again (attempt ", i + 1, " of ", max_tries, ").")
      }
    }
  }

  ## -------------------------------------------------
  ## 4️⃣  **Install Miniconda (fallback)**   <-- *step 4 goes here*
  ## -------------------------------------------------
  if (install_fallback) {
    message("⚠️ No valid Conda path provided – installing Miniconda locally …")
    install_miniconda() # reticulate helper
    cp <- reticulate::conda_binary()
    if (is.null(cp) || !file.exists(cp)) {
      stop(
        "❌ Failed to install or locate Miniconda. Install it manually."
      )
    }
    message("✅ Miniconda installed at: ", cp)
  }
  }

  ## -------------------------------------------------
  ## 5️⃣  Create / reuse the requested conda environment
  ## -------------------------------------------------
  envs <- conda_list(conda = cp)


  if (!(env_name %in% envs$name)) {
    message(sprintf(
      "Creating conda environment '%s' and installing pacmap.",
      env_name
    ))
    conda_create(envname = env_name, conda = cp)
    conda_install(
      envname = env_name,
      packages = "pacmap",
      pip = TRUE,
      conda = cp
    )
  } else {
    # Ensure pacmap is present in an existing env
    out <- system2(cp, c("list", "-n", env_name, "--json"), stdout = TRUE)
    conda_package_list <- jsonlite::fromJSON(paste(out, collapse = "\n"))
    pkgs <- conda_package_list$name

    if (!any(pkgs == "pacmap")) {
      message("Installing pacmap into existing environment …")
      conda_install(
        envname = env_name,
        packages = "pacmap",
        pip = TRUE,
        conda = cp
      )
      message("Pacmap set up.")
    }
if (any(pkgs == "pacmap")) {
message("Environment with pacmap present in user environments. Pacmap set up.")
  
  
  }
  }
  
  ## -------------------------------------------------
  ## 6️⃣  Activate the environment and import PacMAP
  ## -------------------------------------------------
  use_condaenv(env_name, conda = cp, required = TRUE)
  pacmap <- import("pacmap")
  return(pacmap)
}



# ========================================================================
# Main function: Cluster/Iterative Marker feature selection
# ========================================================================




CIM_feature_selection_by_gene_set_pacmap <- function(clustering_matrix, # rows: samples x cols: genes
                                                     all_gene_sets = NULL,
                                                     remove_immune_variable_genes = TRUE,
                                                     clustering_alg,
                                                     max_k = 5,
                                                     CCP_iter = 5000,
                                                     adj_pval_thresh = 0.05,
                                                     max_pipeline_iter = 50,
                                                     seed,
                                                     clustering_metrics = c("pac", "silhouette_combined_avg", "item_cluster_consensus"),
                                                     filter_approach = c("app_one", "app_two", "app_three"),
                                                     pacmap_dimensions = 2L,
                                                     pacmap_args = NULL,               # user overrides (list)
                                                     pacmap_guardrails = TRUE,         # enforce allowed arg names
                                                     verbose = TRUE,
                                                     working_dir) {
  
  
  
  
  
  # Set seed 
  set.seed(2024L)


# Run on script load
pacmap <- setup_conda(env_name = "pacmap_env")

# Get all_gene_sets if not provided (DEFAULT IS CIM GENES USED IN PAPER)
# Default is to use most updated list from msig

 ## 1. If the user already supplied a list, just return it
  if (!is.null(all_gene_sets)){
  message("Using user provided genesets...")
  }


  if (is.null(all_gene_sets)) {
    message("Using default CIM manuscript genesets...")
    msig_df <- msigdbr::msigdbr(species = "Homo sapiens") %>% as.data.frame()

    # Default CIM Genesets from manuscript
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
    
  # Build a named list: each element is the vector of gene symbols for that set
    all_gene_sets <- setNames(
      lapply(all_gene_set_names, function(gs) {
        gs_df <- msig_df[msig_df$gs_name == gs, , drop = FALSE]
        gs_df$gene_symbol
      }),
      all_gene_set_names
    )
  }

 

# Clean up all gene sets by removing immune genes if user desires

if (remove_immune_variable_genes){

removed_immune_genes <- c()
  # -----------------------------------------------------------------
# 1. Define the prefixes that identify immunoglobulin (B‑cell) and
#    T‑cell receptor (TCR) genes
# -----------------------------------------------------------------
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
# 2. Apply the filter to the *named list* `all_gene_sets`
#    (each element is a character vector of gene symbols)
# -----------------------------------------------------------------


  all_gene_sets <- lapply(all_gene_sets, function(gene_vec) {
  # genes that match the immune pattern → remove
  removed_immune_hits <- gene_vec[grepl(pattern, gene_vec, perl = TRUE)]
  # accumulate them globally
  removed_immune_genes <<- c(removed_immune_genes, removed_immune_hits)

  # keep everything else
  gene_vec[!grepl(pattern, gene_vec, perl = TRUE)]
})

  
  immune_genes_removed_all <- unique(removed_immune_genes)

  
}
  
  
  # Set random state for UMAP
  custom.config <- umap.defaults
  custom.config$random_state <- 2024L
  
  # Set options for approaches 
  no_app_one <- FALSE
  no_app_two <- FALSE
  no_app_three <- FALSE
  
  # keep genes that are in clustering matrix 
  genes <- colnames(clustering_matrix)
  all_gene_sets <- lapply(all_gene_sets, function(gs) intersect(gs, genes))
  
  # Genes removed for each set because they're not in clustering_matrix
  removed_by_membership_list <- lapply(all_gene_sets, function(gs) {
    setdiff(gs, genes)
  })
  
  removed_by_membership <- unique(unlist(removed_by_membership_list))
  
  # Remove zero var genes
  
  # Get genes with zero variance in the clustering_matrix
  zero_var_genes <- colnames(clustering_matrix)[apply(clustering_matrix, 2, sd) == 0]
  
  
  # Filter each gene set in all_gene_sets
  all_gene_sets <- lapply(all_gene_sets, function(gene_vec) {
    setdiff(gene_vec, zero_var_genes)
  })
  
  # Get genes from all_gene_sets that were filtered out
  removed_genes <- lapply(all_gene_sets, function(gene_vec) {
    intersect(gene_vec, zero_var_genes)
  })
  
 # Combine the three collections of removed genes into one vector
removed_genes <- c(
  unlist(removed_genes),          # genes removed by zero‑variance filtering
  unlist(removed_by_membership),  # genes absent from the clustering matrix
  immune_genes_removed_all        # immune genes filtered out earlier
)

# (Optional) make the list unique and sorted
removed_genes <- sort(unique(removed_genes))
  
  # Validate clustering_alg
  if (!(clustering_alg %in% c("km", "hc"))) {
    stop("clustering_alg must be either 'km' (for k-means) or 'hc' (for hierarchical clustering)")
  }
  
  # validate clustering metrics
  valid_metrics <- c("pac", "silhouette_dim_reduce_space", "silhouette_ge_space", "silhouette_combined_avg", "item_cluster_consensus")
  if (any(!clustering_metrics %in% valid_metrics)) {
    stop(
      paste0(
        "clustering_metrics must singular metric or a  vector combination (i.e. c('pac', 'silhouette_dim_reduce_space')) from: ",
        paste(valid_metrics, collapse = ", ")
      )
    )
  }
  
  
# -------------------------------------------------
# 3️⃣  **Check retained gene counts**
#      – stop if a set is empty
#      – warn if the retained size is < 10 % of the original size
# -------------------------------------------------
for (gs_name in names(all_gene_sets)) {
  ## original size before any filtering
  orig_len   <- length( ## pull the raw list that was created earlier ##
    if (exists("all_gene_set_names")) {
      # when the default manuscript sets were used
      # (they are stored in `all_gene_set_names`)
      # we can safely reuse that vector
      setdiff(
        msig_df$gene_symbol[msig_df$gs_name == gs_name],
        character()
      )
    } else {
      # when the user supplied a list we already have it in `all_gene_sets`
      # (but we need the *pre‑filter* version – store it before any filtering)
      # In this script the pre‑filter version is the same object that we are
      # iterating over, so we simply use `all_gene_sets[[gs_name]]` before
      # any removal was applied – that value is kept in `orig_len` here.
      all_gene_sets[[gs_name]]
    }
  )
  orig_len   <- length(orig_len)

  ## number of genes that survived filtering
  n_retained <- length(all_gene_sets[[gs_name]])

  ## 3a.  Zero‑gene set → fatal error
  if (n_retained == 0) {
    stop(sprintf(
      "Gene set '%s' has **no** genes after filtering (zero‑variance or immune removal).",
      gs_name
    ), call. = FALSE)
  }

  ## 3b.  Low‑gene warning – 10 % of the original size (minimum of 1 gene)
  min_warn   <- max(ceiling(0.10 * orig_len), 1L)
  if (n_retained < min_warn) {
    warning(sprintf(
      "Gene set '%s' retained only %d genes (10 %% of original = %d). Consider relaxing the filter or checking the input data.",
      gs_name, n_retained, min_warn
    ), call. = FALSE)
  }
}
  
  
  # because app_three dependent on app_one
  if ("app_three" %in% filter_approach && !"app_one" %in% filter_approach) {
    filter_approach <- c(filter_approach, "app_one")
  }
  
  # Set PacMAP for Script 

  # Seed on python side (incase)
  np <- reticulate::import("numpy", convert = FALSE)
  np$random$seed <- 2024L  # reproducibility on Python side

  
  # n_neigh_pacmap code
  sample_number <- nrow(clustering_matrix) 
  n_neigh_pacmap <- if (sample_number <= 20L){
    as.integer(max(5L, floor(sample_number * 0.25)))
  } else NULL 
  
  
  # Defaults
  pm_defaults <- list(
    n_components = pacmap_dimensions,
    n_neighbors  = n_neigh_pacmap,     
    random_state = seed
  )
  
  # in case of user overrides
  if (!is.null(pacmap_args)) {
    if (pacmap_guardrails) {
      allowed <- c("n_components","n_neighbors","MN_ratio","FP_ratio","num_iters",
                   "apply_pca","random_state","PCA_dim")
      bad <- setdiff(names(pacmap_args), allowed)
      if (length(bad)) {
        warning("Ignoring unrecognized pacmap_args: ", paste(bad, collapse = ", "))
        pacmap_args[bad] <- NULL
      }
    }
    pm_defaults <- utils::modifyList(pm_defaults, pacmap_args)
  }
  
  
  # for writing csv
  pm_settings <- as.data.frame(unlist(pm_defaults))
  
  # Build embedder for analysis 
  script_pacmap_reducer <- do.call(pacmap$PaCMAP, pm_defaults)

  # Use default reducer for visualization 
  
  # create options
  pm_visualization <- pm_defaults
  pm_visualization$n_components <- 2L
  
  # create reducer 
  pm_visualization_reducer <- do.call(pacmap$PaCMAP, pm_visualization)
  
  # n_neigh code umap
  sample_number <- nrow(clustering_matrix) 
  n_neigh_umap <- if (sample_number <= 20L){
    as.integer(max(5L, floor(sample_number * 0.25)))
  } else 15L
  
  # Defaults for umap visual
  umap_option_visual <- list(
    n_components = 2L,
    n_neighbors  = n_neigh_umap,     
    random_state = 2024L
  )
  
  custom.config.visual <- utils::modifyList(umap.defaults, umap_option_visual) # Set options for UMAP
  
  # ---- Initialize outputs ----
  details_per_set_analysis_by_gene_sets <- list()
  details_per_set_analysis_by_all_genes <- list()
  details_per_set_gene_sets_then_all_genes <- list()
  
  iter_log <- data.frame(gene_set=character(), iteration=integer(),
                         before=integer(), after=integer(), stringsAsFactors=FALSE)
  new_final <- character(0)
  

  
  # ---- Approach One Iterate over all gene sets ----
  
  

  
  if ("app_one" %in% filter_approach){

    # Set seed appropriately
    withr::with_seed(seed = 2024L, {
    
  for (gene_set_i in seq_along(all_gene_sets)) {

    
    # set max_k for run
    working_max_k <- max_k
    
    gene_set_name <- names(all_gene_sets)[gene_set_i]
    if (is.null(gene_set_name)) gene_set_name <- paste0("set", gene_set_i)
    if (verbose) message(sprintf("\nGene set: %s", gene_set_name))
    
    iterated_final_gene_set <- all_gene_sets[[gene_set_i]]
    gene_set_check_1 <- 0
    gene_set_check_2 <- length(iterated_final_gene_set)
    n <- 0
    details_per_set_analysis_by_gene_sets[[gene_set_name]] <- list()
    stats_per_iter <- list()
    cluster_metrics_per_iter <- list()
    optimal_k_per_iter <- list()
    cluster_stability_per_iter <- list()
    wide_cluster_stability_per_iter <- list()
    
    # ---- Iteratively refine gene set based on feature selection after clustering ----
    while (gene_set_check_1 != gene_set_check_2 && n < max_pipeline_iter) {
      if (verbose) {
        message(sprintf(
          "%s Iteration %d Start: %d genes (%d previously)",
          gene_set_name,
          n,
          gene_set_check_2,
          gene_set_check_1
        ))
      }

      iterated_feature_df <- data.frame(
        gene_id = character(0),
        F_value = numeric(0),
        p_value = numeric(0)
      )

      # --- Run PaCMAP; attempt safely in case of too-small gene sets ---

      # Build embedder for analysis (also resets seed)
      script_pacmap_reducer <- do.call(pacmap$PaCMAP, pm_defaults)

      safe_pacmap <- tryCatch(
        {
          # Run PaCMAP
          embedding <- script_pacmap_reducer$fit_transform(
            clustering_matrix[, iterated_final_gene_set, drop = FALSE]
          )
          embedding
        },
        error = function(e) {
          message(sprintf(
            "PaCMAP failed for gene set %s iteration %d. Skipping this set. Most iterated gene set will be kept (if iter > 0)",
            gene_set_name,
            n
          ))
          NULL
        }
      )

      if (is.null(safe_pacmap)) {
        if (n >= 1) {
          # Keep original gene set if iteration > 1; i.e. some genes have been iterated out
          # Convert back to old geneset
          details_per_set_analysis_by_gene_sets[[
            gene_set_name
          ]]$final <- iterated_final_gene_set
        }
        details_per_set_analysis_by_gene_sets[[
          gene_set_name
        ]]$feature_importance_per_iter <- stats_per_iter
        details_per_set_analysis_by_gene_sets[[
          gene_set_name
        ]]$optimal_k_per_iter <- optimal_k_per_iter
        details_per_set_analysis_by_gene_sets[[
          gene_set_name
        ]]$cluster_metrics_per_iter <- cluster_metrics_per_iter
        details_per_set_analysis_by_gene_sets[[
          gene_set_name
        ]]$cluster_stability_per_iter <- cluster_stability_per_iter
        details_per_set_analysis_by_gene_sets[[
          gene_set_name
        ]]$wide_cluster_stability_per_iter <- wide_cluster_stability_per_iter

        new_final <- unique(c(new_final, iterated_final_gene_set))

        break
      }

      rownames(safe_pacmap) <- rownames(clustering_matrix)

      # ---- Prepare PaCMAP result for clustering and feature selection ----
      pacmap_plot_df <- data.frame(safe_pacmap) %>%
        tibble::rownames_to_column("sample_id") %>%
        dplyr::inner_join(
          as.data.frame(clustering_matrix) %>% rownames_to_column("sample_id"),
          by = "sample_id"
        )

      # ---- Determine optimal number of clusters via Consensus Clustering ----

      # Your input matrix: samples x features (PaCMAP embedding)
      input_matrix <- t(safe_pacmap) # transpose to get dims as rows (like you did with UMAP)
      colnames(input_matrix) <- rownames(clustering_matrix)

      # Initialize starting value for pItem (percentage of items to sample per iteration)
      pItem_val <- 0.8

      # Set maximum allowable pItem value (cannot exceed 1.0)
      max_pItem <- 1.0

      # Step size for incrementing pItem when NA values are found
      step <- 0.05

      # Flag to track presence of NA values in itemConsensus
      has_na <- TRUE

      # if in the case of working max k = 2
      if (working_max_k == 2) {
        # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
        while (has_na && pItem_val <= max_pItem) {
          # Run ConsensusClusterPlus with the current pItem_val
          consensus_results <- ConsensusClusterPlus(
            input_matrix,
            maxK = 3,
            reps = CCP_iter,
            pItem = pItem_val, # Percentage of items sampled per iteration
            pFeature = 1, # Use all features
            clusterAlg = clustering_alg,
            innerLinkage = "ward.D2",
            finalLinkage = "ward.D2",
            distance = "euclidean",
            seed = 2024L,
            plot = "none",
            verbose = FALSE,
            writeTable = FALSE
          )

          # Only extract out result for K = 2
          consensus_results <- consensus_results[1:2]

          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL) # Start invisible plotting device
          icl_results <- calcICL(
            consensus_results,
            plot = "png",
            writeTable = TRUE
          )
          dev.off() # Close device

          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]

          # Rename columns for clarity
          colnames(item_consensus_df) <- c(
            "k",
            "cluster",
            "sample_id",
            "item_consensus"
          )

          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL) # Start invisible plotting device
          icl_results <- calcICL(
            consensus_results,
            plot = "png",
            writeTable = TRUE
          )
          dev.off() # Close device

          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]

          # Rename columns for clarity
          colnames(item_consensus_df) <- c(
            "k",
            "cluster",
            "sample_id",
            "item_consensus"
          )

          # Check if any values in the item_consensus column are NA
          if (any(is.na(item_consensus_df$item_consensus))) {
            # Print message indicating NA detected and that pItem will be increased
            message(sprintf(
              "NA found in item_consensus at pItem=%.2f, increasing pItem...",
              pItem_val
            ))

            # Increase pItem by the predefined step to sample more items per iteration
            pItem_val <- pItem_val + step
          } else {
            # If no NA found, set flag to FALSE to exit the loop
            has_na <- FALSE
          }
        }
      } else {
        # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
        while (has_na && pItem_val <= max_pItem) {
          # Run ConsensusClusterPlus with the current pItem_val
          consensus_results <- ConsensusClusterPlus(
            input_matrix,
            maxK = working_max_k,
            reps = CCP_iter,
            pItem = pItem_val, # Percentage of items sampled per iteration
            pFeature = 1, # Use all features
            clusterAlg = clustering_alg,
            innerLinkage = "ward.D2",
            finalLinkage = "ward.D2",
            distance = "euclidean",
            seed = 2024L,
            plot = "none",
            verbose = FALSE,
            writeTable = FALSE
          )

          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL) # Start invisible plotting device
          icl_results <- calcICL(
            consensus_results,
            plot = "png",
            writeTable = TRUE
          )
          dev.off() # Close device

          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]

          # Rename columns for clarity
          colnames(item_consensus_df) <- c(
            "k",
            "cluster",
            "sample_id",
            "item_consensus"
          )

          # Check if any values in the item_consensus column are NA
          if (any(is.na(item_consensus_df$item_consensus))) {
            # Print message indicating NA detected and that pItem will be increased
            message(sprintf(
              "NA found in item_consensus at pItem=%.2f, increasing pItem...",
              pItem_val
            ))

            # Increase pItem by the predefined step to sample more items per iteration
            pItem_val <- pItem_val + step
          } else {
            # If no NA found, set flag to FALSE to exit the loop
            has_na <- FALSE
          }
        }
      }

      # After loop finishes, check if NA values still exist (meaning max pItem was reached)
      if (has_na) {
        # Warn the user that maximum pItem was reached but NAs persist
        message(
          "Reached max pItem but NA values still exist in item_consensus."
        )
        message("Reducing working_max_k to the maximum most stable K")

        # Find all K values with no NA item consensus
        ks_without_na <- item_consensus_df %>%
          group_by(k) %>%
          summarize(no_na = all(!is.na(item_consensus))) %>%
          filter(no_na) %>%
          pull(k)

        # if there is no solution skip it
        if (length(ks_without_na) == 0) {
          message(
            "No stable K values found. Skipping this geneset/iteration. \n Will keep most iterated if iteration > 0"
          )

          if (n >= 1) {
            # Document in iter_log
            iter_log <- rbind(
              iter_log,
              data.frame(
                gene_set = gene_set_name,
                iteration = n,
                before = length(old_iterated_final_gene_set),
                after = length(iterated_final_gene_set)
              )
            )
            # update iter log
            iter_log[nrow(iter_log), "after"] <- paste0(
              length(iterated_final_gene_set),
              "; No stable clustering produced at higher iterations. Took most iterated geneset."
            )

            break # exits the while loop
          }

          if (n == 0) {
            # Document in iter_log
            iter_log <- rbind(
              iter_log,
              data.frame(
                gene_set = gene_set_name,
                iteration = n,
                before = length(iterated_final_gene_set),
                after = 0
              )
            )
            # update iter log
            iter_log[nrow(iter_log), "after"] <- paste0(
              "0 ; No stable clustering produced at first iteration. Skipped geneset."
            )

            break # exits the while loop
          }
        } else {
          # Get maximal most stable K
          max_stable_k <- max(ks_without_na)
          message("maximum most stable K = ", max_stable_k)
        }

        # special condition for if maximum most stable K = 2
        if (max_stable_k == 2) {
          # let user know
          message("since maximum most stable K = ", max_stable_k)
          message("Will have to extract k = 2 results using k = 3 analysis")

          # set that to working_k
          working_max_k <- 2

          # Initialize starting value for pItem (percentage of items to sample per iteration)
          pItem_val <- 0.8

          # Set maximum allowable pItem value (cannot exceed 1.0)
          max_pItem <- 1.0

          # Step size for incrementing pItem when NA values are found
          step <- 0.05

          # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
          while (has_na && pItem_val <= max_pItem) {
            # Run ConsensusClusterPlus with the current pItem_val
            consensus_results <- ConsensusClusterPlus(
              input_matrix,
              maxK = 3,
              reps = CCP_iter,
              pItem = pItem_val, # Percentage of items sampled per iteration
              pFeature = 1, # Use all features
              clusterAlg = clustering_alg,
              innerLinkage = "ward.D2",
              finalLinkage = "ward.D2",
              distance = "euclidean",
              seed = 2024L,
              plot = "none",
              verbose = FALSE,
              writeTable = FALSE
            )

            # Only extract out result for K = 2
            consensus_results <- consensus_results[1:2]

            # Calculate item consensus and cluster metrics (ICL)
            pdf(file = NULL) # Start invisible plotting device
            icl_results <- calcICL(
              consensus_results,
              plot = "png",
              writeTable = TRUE
            )
            dev.off() # Close device

            # Extract the itemConsensus dataframe from ICL results
            item_consensus_df <- icl_results[["itemConsensus"]]

            # Rename columns for clarity
            colnames(item_consensus_df) <- c(
              "k",
              "cluster",
              "sample_id",
              "item_consensus"
            )

            # Calculate item consensus and cluster metrics (ICL)
            pdf(file = NULL) # Start invisible plotting device
            icl_results <- calcICL(
              consensus_results,
              plot = "png",
              writeTable = TRUE
            )
            dev.off() # Close device

            # Extract the itemConsensus dataframe from ICL results
            item_consensus_df <- icl_results[["itemConsensus"]]

            # Rename columns for clarity
            colnames(item_consensus_df) <- c(
              "k",
              "cluster",
              "sample_id",
              "item_consensus"
            )

            # Check if any values in the item_consensus column are NA
            if (any(is.na(item_consensus_df$item_consensus))) {
              # Print message indicating NA detected and that pItem will be increased
              message(sprintf(
                "NA found in item_consensus at pItem=%.2f, increasing pItem...",
                pItem_val
              ))

              # Increase pItem by the predefined step to sample more items per iteration
              pItem_val <- pItem_val + step
            } else {
              # If no NA found, set flag to FALSE to exit the loop
              has_na <- FALSE
            }
          }
        } else {
          # set that to working_k
          working_max_k <- max_stable_k

          # Initialize starting value for pItem (percentage of items to sample per iteration)
          pItem_val <- 0.8

          # Set maximum allowable pItem value (cannot exceed 1.0)
          max_pItem <- 1.0

          # Step size for incrementing pItem when NA values are found
          step <- 0.05

          # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
          while (has_na && pItem_val <= max_pItem) {
            # Run ConsensusClusterPlus with the current pItem_val
            consensus_results <- ConsensusClusterPlus(
              input_matrix,
              maxK = working_max_k,
              reps = CCP_iter,
              pItem = pItem_val, # Percentage of items sampled per iteration
              pFeature = 1, # Use all features
              clusterAlg = clustering_alg,
              innerLinkage = "ward.D2",
              finalLinkage = "ward.D2",
              distance = "euclidean",
              seed = 2024L,
              plot = "none",
              verbose = FALSE,
              writeTable = FALSE
            )

            # Calculate item consensus and cluster metrics (ICL)
            pdf(file = NULL) # Start invisible plotting device
            icl_results <- calcICL(
              consensus_results,
              plot = "png",
              writeTable = TRUE
            )
            dev.off() # Close device

            # Extract the itemConsensus dataframe from ICL results
            item_consensus_df <- icl_results[["itemConsensus"]]

            # Rename columns for clarity
            colnames(item_consensus_df) <- c(
              "k",
              "cluster",
              "sample_id",
              "item_consensus"
            )

            # Calculate item consensus and cluster metrics (ICL)
            pdf(file = NULL) # Start invisible plotting device
            icl_results <- calcICL(
              consensus_results,
              plot = "png",
              writeTable = TRUE
            )
            dev.off() # Close device

            # Extract the itemConsensus dataframe from ICL results
            item_consensus_df <- icl_results[["itemConsensus"]]

            # Rename columns for clarity
            colnames(item_consensus_df) <- c(
              "k",
              "cluster",
              "sample_id",
              "item_consensus"
            )

            # Check if any values in the item_consensus column are NA
            if (any(is.na(item_consensus_df$item_consensus))) {
              # Print message indicating NA detected and that pItem will be increased
              message(sprintf(
                "NA found in item_consensus at pItem=%.2f, increasing pItem...",
                pItem_val
              ))

              # Increase pItem by the predefined step to sample more items per iteration
              pItem_val <- pItem_val + step
            } else {
              # If no NA found, set flag to FALSE to exit the loop
              has_na <- FALSE
            }
          }
        }
      }

      # If no NA found, proceed with downstream analysis steps
      if (!has_na) {
        # Classify stability based on item_consensus values using custom function
        item_consensus_df$stability <- sapply(
          item_consensus_df$item_consensus,
          get_stability_category
        )

        # Split the dataframe into a list of dataframes, one per cluster number k
        item_consensus_by_k <- split(item_consensus_df, item_consensus_df$k)
      }

      # Calculate PAC metric (Proportion of Ambiguous Clustering) using ecdf

      Kvec = 2:working_max_k # Vector of cluster numbers to evaluate
      x1 = 0.1
      x2 = 0.9 # Thresholds defining "ambiguous" region
      PAC = rep(NA, length(Kvec)) # Placeholder for PAC values
      names(PAC) = paste("K=", Kvec, sep = "") # Label PAC entries with cluster numbers

      for (i in Kvec) {
        # Loop over cluster numbers
        M = consensus_results[[i]]$consensusMatrix # Extract consensus matrix for current K
        Fn = ecdf(M[lower.tri(M)]) # Empirical CDF of the lower triangle of consensus matrix
        PAC[i - 1] = Fn(x2) - Fn(x1) # The fraction of values between x1 and x2 (ambiguous assignments)
      }
      optK = Kvec[which.min(PAC)] # Select K with the lowest PAC (most stable clustering)

      # Convert PAC to data frame for plotting
      pac_values <- data.frame(k = Kvec, pac = unname(PAC))

      # Silhouette calculations using consensus cluster

      input_matrix <- t(input_matrix) # transpose to get dims as samples x features
      rownames(input_matrix) <- rownames(clustering_matrix)

      # Distances
      dist_space_for_clustering <- dist(input_matrix)
      dist_ge_space <- dist(clustering_matrix[,
        iterated_final_gene_set,
        drop = FALSE
      ])

      # Collect silhouettes for all k
      sil_results <- lapply(2:length(consensus_results), function(k) {
        clusters <- consensus_results[[k]]$consensusClass

        # Silhouette in gene space
        sil_space_for_clustering <- silhouette(
          clusters,
          dist_space_for_clustering
        )
        avg_space_for_clustering <- mean(sil_space_for_clustering[, 3])

        # Silhouette in PaCMAP space
        sil_ge_space <- silhouette(clusters, dist_ge_space)
        avg_ge_space <- mean(sil_ge_space[, 3])

        data.frame(
          k = k,
          avg_sil_space_for_clustering = avg_space_for_clustering,
          avg_sil_ge_space = avg_ge_space
        )
      })

      # Combine
      sil_df <- do.call(rbind, sil_results)
      sil_df$avg_sil_combined <- ((sil_df$avg_sil_space_for_clustering +
        sil_df$avg_sil_ge_space) /
        2)

      #Get cluster_consensus_df
      cluster_consensus_df <- as.data.frame(icl_results[["clusterConsensus"]])

      # Create item consensus table
      avg_cluster_consensus <- cluster_consensus_df %>%
        group_by(k) %>%
        summarize(avg_cluster_consensus = mean(clusterConsensus))

      # Create metrics summary dataframe
      metrics_summary <- data.frame(
        k = 2:working_max_k,
        silhouette_dim_reduce_space = NA, # renamed to indicate dim‑reduce (original) space
        silhouette_ge_space = NA, # silhouette computed in PacMAP (ge) space
        silhouette_combined_avg = NA, # silhouette computed in PacMAP (ge) space
        pac = NA,
        avg_cluster_consensus = NA
      )

      # Fill in the metrics
      for (i in 1:nrow(metrics_summary)) {
        k_val <- metrics_summary$k[i]

        # Silhouette in the original (dim‑reduce) space
        idx <- which(sil_df$k == k_val)
        if (length(idx) > 0) {
          metrics_summary$silhouette_dim_reduce_space[
            i
          ] <- sil_df$avg_sil_space_for_clustering[idx]

          metrics_summary$silhouette_ge_space[i] <- sil_df$avg_sil_ge_space[idx]
          metrics_summary$silhouette_combined_avg[i] <- sil_df$avg_sil_combined[
            idx
          ]
        }

        # PAC
        k_idx <- which(Kvec == k_val)
        if (length(k_idx) > 0) {
          metrics_summary$pac[i] <- PAC[k_idx]
        }

        # avg_cluster_consensus
        cc_idx <- which(avg_cluster_consensus$k == k_val)
        if (length(cc_idx) > 0) {
          metrics_summary$avg_cluster_consensus[
            i
          ] <- avg_cluster_consensus$avg_cluster_consensus[cc_idx]
        }
      }

      # Alternatively,  Optimal K can be selected using rank based compromise
      # The rank-based compromise approach assigns a rank to each k for both silhouette (descending) and PAC (ascending), then sums the ranks to identify the k with the best overall trade-off between high silhouette and low PAC. The k with the lowest combined rank is selected as the optimal choice.
      # https://dl.acm.org/doi/10.1145/371920.372165

      # Assign ranks (lowest rank = best)
      metrics_summary$rank_silhouette_dim_reduce_space <- rank(
        -metrics_summary$silhouette_dim_reduce_space,
        ties.method = "min"
      ) # higher silhouette = rank 1 (original space)
      metrics_summary$rank_silhouette_ge_space <- rank(
        -metrics_summary$silhouette_ge_space,
        ties.method = "min"
      ) # higher silhouette = rank 1 (ge space)

      metrics_summary$rank_silhouette_combined_avg <- rank(
        -metrics_summary$silhouette_combined_avg,
        ties.method = "min"
      ) # higher silhouette = rank 1 (averaged space)
      metrics_summary$rank_pac <- rank(metrics_summary$pac, ties.method = "min") # lower PAC = rank 1
      metrics_summary$rank_item_cluster_consensus <- rank(
        -metrics_summary$avg_cluster_consensus,
        ties.method = "min"
      ) # higher consensus = rank 1

      # Build the column names dynamically.
      rank_columns <- c(
        paste0("rank_", clustering_metrics) # e.g. "rank_pac", "rank_item_cluster_consensus"
      )

      # Compute the overall rank as the row sum across selected rank columns
      metrics_summary$overall_rank <- rowSums(
        metrics_summary[, rank_columns, drop = FALSE],
        na.rm = TRUE
      )

      # Find all rows with the minimum overall rank
      best_rows <- metrics_summary[
        metrics_summary$overall_rank == min(metrics_summary$overall_rank),
      ]

      # Use the last metric in clustering_metrics as the tiebreaker
      tie_metric <- clustering_metrics[length(clustering_metrics)]

      if (verbose) {
        cat(
          "Metrics being used for ranking:\n  ",
          paste(clustering_metrics, collapse = ", "),
          "\n"
        )
        cat(
          "Tie‑breaker metric (last in clustering_metrics vector):",
          tie_metric,
          "\n"
        )
      }

      # Determine if the tiebreaking metric should be maximized or minimized
      # (Assuming: max for 'silhouette' and 'item_cluster_consensus', min for 'pac')
      if (tie_metric == "silhouette_dim_reduce_space") {
        idx <- which.max(best_rows$silhouette_dim_reduce_space)
      } else if (tie_metric == "silhouette_ge") {
        idx <- which.max(best_rows$silhouette_ge_space)
      } else if (tie_metric == "item_cluster_consensus") {
        idx <- which.max(best_rows$avg_cluster_consensus)
      } else if (tie_metric == "pac") {
        idx <- which.min(best_rows$pac)
      } else {
        stop("Unknown tiebreaker metric.")
      }

      best_k_row <- best_rows[idx, ]
      best_k <- best_k_row$k

      # add cluster metrics used for evaluation and tie
      metrics_summary$metrics_used <- paste(clustering_metrics, collapse = ", ")
      metrics_summary$metrics_used_for_tie <- tie_metric

      # Based on examining the plots and metrics
      optimal_k <- best_k

      # Extract cluster assignments for optimal k
      consensus_clusters <- consensus_results[[optimal_k]]$consensusClass

      # Get stability/reliability metrics
      optimal_stability_cluster_metrics <- item_consensus_by_k[[as.character(
        optimal_k
      )]]

      # Start from your original data frame
      wide_optimal_stability_cluster_metrics <- optimal_stability_cluster_metrics %>%
        # Only keep relevant columns
        dplyr::select(sample_id, cluster, item_consensus) %>%
        # Pivot wider: each cluster makes a new column
        pivot_wider(
          names_from = cluster,
          values_from = item_consensus,
          names_prefix = "item_consensus_cluster"
        )

      # ---- assign labels from consensus plus clustering  ----

      pacmap_plot_df$kmeans_cluster <- consensus_clusters

      # ---- Feature selection for each gene ----
      # Preallocate
      Fvals <- numeric(length(iterated_final_gene_set)) # This will be effect size now
      Pvals <- numeric(length(iterated_final_gene_set))
      TestType <- character(length(iterated_final_gene_set))

      for (i in seq_along(iterated_final_gene_set)) {
        gene_name <- iterated_final_gene_set[i]
        vals <- pacmap_plot_df[[gene_name]]
        cluster <- as.factor(pacmap_plot_df[["kmeans_cluster"]])
        n_groups <- nlevels(cluster)
        df <- data.frame(vals = vals, cluster = cluster)

        if (n_groups == 2) {
          # Wilcoxon test and effect size with rstatix
          wt <- wilcox.test(vals ~ cluster)
          eff <- tryCatch(
            {
              ef <- rstatix::wilcox_effsize(df, vals ~ cluster, ci = FALSE)
              ef$effsize
            },
            error = function(e) NA
          )
          Pvals[i] <- wt$p.value
          Fvals[i] <- eff
          TestType[i] <- "wilcox"
        } else if (n_groups >= 3) {
          # Kruskal test and effect size with rstatix
          kw <- kruskal.test(vals ~ cluster)
          eff <- tryCatch(
            {
              ef <- rstatix::kruskal_effsize(df, vals ~ cluster, ci = FALSE)
              ef$effsize
            },
            error = function(e) NA
          )
          Pvals[i] <- kw$p.value
          Fvals[i] <- eff
          TestType[i] <- "kruskal"
        } else {
          Pvals[i] <- NA
          Fvals[i] <- NA
          TestType[i] <- "NA"
        }
      }

      # Adjusted P-values (FDR)
      adj_Pvals <- p.adjust(Pvals, method = "fdr")

      iterated_feature_df <- data.frame(
        gene_id = iterated_final_gene_set,
        effect_size = as.numeric(Fvals),
        p_value = as.numeric(Pvals),
        adj_p_value = as.numeric(adj_Pvals),
        test = TestType
      )

      old_iterated_feature_df <- iterated_feature_df
      iterated_feature_df <- type.convert(iterated_feature_df, as.is = TRUE)
      iterated_feature_df <- iterated_feature_df[
        iterated_feature_df$p_value <= 0.05 &
          iterated_feature_df$adj_p_value <= adj_pval_thresh,
      ]
      old_iterated_final_gene_set <- iterated_final_gene_set
      iterated_final_gene_set <- unique(sort(iterated_feature_df$gene_id))
      stats_per_iter[[n + 1]] <- iterated_feature_df
      cluster_metrics_per_iter[[n + 1]] <- metrics_summary
      optimal_k_per_iter[[n + 1]] <- optimal_k
      cluster_stability_per_iter[[n + 1]] <- optimal_stability_cluster_metrics
      wide_cluster_stability_per_iter[[
        n + 1
      ]] <- wide_optimal_stability_cluster_metrics

      iter_log <- rbind(
        iter_log,
        data.frame(
          gene_set = gene_set_name,
          iteration = n,
          before = length(old_iterated_final_gene_set),
          after = length(iterated_final_gene_set)
        )
      )

      # --- check convergence; increment loop ---
      gene_set_check_1 <- length(old_iterated_final_gene_set)
      gene_set_check_2 <- length(iterated_final_gene_set)
      if (verbose) {
        message(sprintf(
          "%s Iteration %d End: %d genes (%d previously)",
          gene_set_name,
          n,
          gene_set_check_2,
          gene_set_check_1
        ))
      }
      n <- n + 1
      if (n >= max_pipeline_iter) {
        warning(sprintf(
          "Reached max_pipeline_iter (%d) for gene set %s",
          max_pipeline_iter,
          gene_set_name
        ))
      }

      # Check if there is non-immediate convergence to zero
      if (gene_set_check_2 == 0 && gene_set_check_1 < length(gene_set_i)) {
        cat(
          "Non-immediate convergence to zero detected: reverting to most recent non-zero iteration"
        )

        # Convert back to old geneset
        iterated_final_gene_set <- old_iterated_final_gene_set

        # update last row to reflect immediate convergence iter log
        iter_log[nrow(iter_log), "after"] <- paste0(
          length(iterated_final_gene_set),
          "; non-immediate convergence reverted to most recent non-zero iteration"
        )

        # break the loop
        break
      }

      # Check if there is immediate convergence to zero
      if (gene_set_check_2 == 0 && n == 1) {
        cat(
          "Immediate convergence to zero detected: will retain the top 10% of genes by pval then eff_size. Consider relaxing adj_pval_thresh\n"
        )

        # Rank old gene set by adj_p_value (ascending) and effect_size (descending)
        ranked_genes <- old_iterated_feature_df[
          order(
            old_iterated_feature_df$adj_p_value,
            -old_iterated_feature_df$effect_size
          ),
        ]

        # Keep top 10%
        top_n <- ceiling(0.10 * nrow(ranked_genes))
        iterated_feature_df <- ranked_genes[seq_len(top_n), ]
        iterated_final_gene_set <- iterated_feature_df$gene_id

        # update last row to reflect immediate convergence iter log
        iter_log[nrow(iter_log), "after"] <- paste0(
          length(iterated_final_gene_set),
          "; immediate convergence to zero retained top 10% by pval then eff_size"
        )

        # Break the loop
        break
      }
    }
    # ---- Save per-set details ----
    
    
   
    details_per_set_analysis_by_gene_sets[[gene_set_name]]$final <- iterated_final_gene_set
    details_per_set_analysis_by_gene_sets[[gene_set_name]]$feature_importance_per_iter <- stats_per_iter
    details_per_set_analysis_by_gene_sets[[gene_set_name]]$optimal_k_per_iter <- optimal_k_per_iter
    details_per_set_analysis_by_gene_sets[[gene_set_name]]$cluster_metrics_per_iter <- cluster_metrics_per_iter
    details_per_set_analysis_by_gene_sets[[gene_set_name]]$cluster_stability_per_iter <- cluster_stability_per_iter
    details_per_set_analysis_by_gene_sets[[gene_set_name]]$wide_cluster_stability_per_iter <- wide_cluster_stability_per_iter
    
    new_final <- sort(unique(c(new_final, iterated_final_gene_set)))
  }

  # Check if Gene set is empty
  if (length(new_final) == 0) {
    cat("\033[1;31mWARNING: Final Gene Set for Approach one empty! Consider Relaxing P-val\nApproach output will not be produced\033[0m\n")
    
    # toggle Switch
    no_app_one <- TRUE
    
  }
    
      })
    
    }
  
  # ---- Approach Two iterate over pooled union of all sets ----
  
  if ("app_two" %in% filter_approach){
  
    
    # Set seed appropriately
    withr::with_seed(seed = 2024L, {
    
  if (verbose) message("\n[Second Approach: iteration over all unique genes]")
  iterated_final_gene_set <- unique(unlist(all_gene_sets))
  gene_set_check_1 <- 0
  gene_set_check_2 <- length(iterated_final_gene_set)
  n <- 0
  overall_stats_per_iter <- list()
  details_per_set_analysis_by_all_genes <- list()
  stats_per_iter <- list()
  cluster_metrics_per_iter <- list()
  optimal_k_per_iter <- list()
  cluster_stability_per_iter <- list()
  wide_cluster_stability_per_iter <- list()
  
  
  while (gene_set_check_1 != gene_set_check_2 && n < max_pipeline_iter) {
    if (verbose) message(sprintf("Analysis by all genes Iter %d Start: %d genes (%d prev)", n, gene_set_check_2, gene_set_check_1))
    
    # Set working_max k

    working_max_k <- max_k
    
    
    # --- Run PaCMAP; attempt safely in case of too-small gene sets ---
    
    # Build embedder for analysis (also resets seed)
    script_pacmap_reducer <- do.call(pacmap$PaCMAP, pm_defaults)

    
    safe_pacmap <- tryCatch({

      # Run PaCMAP
      embedding <- script_pacmap_reducer$fit_transform(
        clustering_matrix[, iterated_final_gene_set, drop = FALSE]
      )
      embedding
    }, error = function(e) {
      message(sprintf("PaCMAP failed for gene set %s iteration %d. Skipping this set. Most iterated gene set will be kept (if iter > 0)", gene_set_name, n))
      NULL
    })
    
    if (is.null(safe_pacmap)) {
      details_per_set_analysis_by_all_genes$final <- iterated_final_gene_set
      details_per_set_analysis_by_all_genes$feature_importance_per_iter <- stats_per_iter
      details_per_set_analysis_by_all_genes$optimal_k_per_iter <- optimal_k_per_iter
      details_per_set_analysis_by_all_genes$cluster_metrics_per_iter <- cluster_metrics_per_iter
      details_per_set_analysis_by_all_genes$cluster_stability_per_iter <- cluster_stability_per_iter
      details_per_set_analysis_by_all_genes$wide_cluster_stability_per_iter <- wide_cluster_stability_per_iter
      
      iterated_over_all_genes <- iterated_final_gene_set
      
      break
    }
    
    rownames(safe_pacmap) <- rownames(clustering_matrix)
    
    # ---- Prepare PaCMAP result for clustering and feature selection ----
    
    pacmap_plot_df <- data.frame(safe_pacmap) %>%
      tibble::rownames_to_column("sample_id") %>%
      dplyr::inner_join(as.data.frame(clustering_matrix) %>% rownames_to_column("sample_id"), by = "sample_id")

    
    
    
    # Your input matrix: samples x features (PaCMAP embedding)
    input_matrix <- t(safe_pacmap)  # transpose to get dims as rows (like you did with UMAP)
    colnames(input_matrix) <- rownames(clustering_matrix)
    
 

    
    # ---- Determine optimal number of clusters via Consensus Clustering ----
    
    
    # Initialize starting value for pItem (percentage of items to sample per iteration)
    pItem_val <- 0.8
    
    # Set maximum allowable pItem value (cannot exceed 1.0)
    max_pItem <- 1.0
    
    # Step size for incrementing pItem when NA values are found
    step <- 0.05
    
    # Flag to track presence of NA values in itemConsensus
    has_na <- TRUE
    
    # if in the case of working max k = 2 
    if (working_max_k == 2){
      
      # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
      while (has_na && pItem_val <= max_pItem) {
        # Run ConsensusClusterPlus with the current pItem_val
        consensus_results <- ConsensusClusterPlus(
          input_matrix, 
          maxK = 3, 
          reps = CCP_iter,            
          pItem = pItem_val,            # Percentage of items sampled per iteration
          pFeature = 1,                # Use all features
          clusterAlg = clustering_alg,      
          innerLinkage = "ward.D2",
          finalLinkage = "ward.D2",
          distance = "euclidean", 
          seed = 2024L,
          plot ="none",           
          verbose = FALSE,
          writeTable = FALSE
        )
        
        # Only extract out result for K = 2
        consensus_results <- consensus_results[1:2]
        
        # Calculate item consensus and cluster metrics (ICL)
        pdf(file = NULL)  # Start invisible plotting device
        icl_results <- calcICL(consensus_results, 
                               plot = "png", 
                               writeTable = TRUE)
        dev.off()          # Close device
        
        # Extract the itemConsensus dataframe from ICL results
        item_consensus_df <- icl_results[["itemConsensus"]]
        
        # Rename columns for clarity
        colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
        
        
        # Calculate item consensus and cluster metrics (ICL)
        pdf(file = NULL)  # Start invisible plotting device
        icl_results <- calcICL(consensus_results, 
                               plot = "png", 
                               writeTable = TRUE)
        dev.off()          # Close device
        
        # Extract the itemConsensus dataframe from ICL results
        item_consensus_df <- icl_results[["itemConsensus"]]
        
        # Rename columns for clarity
        colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
        
        # Check if any values in the item_consensus column are NA
        if (any(is.na(item_consensus_df$item_consensus))) {
          # Print message indicating NA detected and that pItem will be increased
          message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
          
          # Increase pItem by the predefined step to sample more items per iteration
          pItem_val <- pItem_val + step
        } else {
          # If no NA found, set flag to FALSE to exit the loop
          has_na <- FALSE
        }
        
      } 
      
    }
    else {
      
      # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
      while (has_na && pItem_val <= max_pItem) {
        
        # Run ConsensusClusterPlus with the current pItem_val
        consensus_results <- ConsensusClusterPlus(
          input_matrix, 
          maxK = working_max_k, 
          reps = CCP_iter,            
          pItem = pItem_val,            # Percentage of items sampled per iteration
          pFeature = 1,                # Use all features
          clusterAlg = clustering_alg,      
          innerLinkage = "ward.D2",
          finalLinkage = "ward.D2",
          distance = "euclidean", 
          seed = 2024L,
          plot ="none",           
          verbose = FALSE,
          writeTable = FALSE
        )
        
        # Calculate item consensus and cluster metrics (ICL)
        pdf(file = NULL)  # Start invisible plotting device
        icl_results <- calcICL(consensus_results, 
                               plot = "png", 
                               writeTable = TRUE)
        dev.off()          # Close device
        
        # Extract the itemConsensus dataframe from ICL results
        item_consensus_df <- icl_results[["itemConsensus"]]
        
        # Rename columns for clarity
        colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
        
        # Check if any values in the item_consensus column are NA
        if (any(is.na(item_consensus_df$item_consensus))) {
          # Print message indicating NA detected and that pItem will be increased
          message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
          
          # Increase pItem by the predefined step to sample more items per iteration
          pItem_val <- pItem_val + step
        } else {
          # If no NA found, set flag to FALSE to exit the loop
          has_na <- FALSE
        }
      }
    }
    
    # After loop finishes, check if NA values still exist (meaning max pItem was reached)
    if (has_na) {
      # Warn the user that maximum pItem was reached but NAs persist
      message("Reached max pItem but NA values still exist in item_consensus.")
      message("Reducing working_max_k to the maximum most stable K")
      
      # Find all K values with no NA item consensus
      ks_without_na <- item_consensus_df %>%
        group_by(k) %>%
        summarize(no_na = all(!is.na(item_consensus))) %>%
        filter(no_na) %>%
        pull(k)
      
      
      # if there is no solution skip it 
      if (length(ks_without_na) == 0) {
        
        message("No stable K values found. Skipping this geneset/iteration. \n Will keep most iterated if iteration > 0")
        
        if(n >= 1){
          
          # Document in iter_log
          iter_log <- rbind(iter_log,
                            data.frame(gene_set = gene_set_name, iteration = n,
                                       before = length(old_iterated_final_gene_set),
                                       after = length(iterated_final_gene_set)))
          # update iter log
          iter_log[nrow(iter_log), "after"] <- paste0(length(iterated_final_gene_set), "; No stable clustering produced at higher iterations. Took most iterated geneset.") 
          
          break  # exits the while loop
          
        }
        
        
        if(n == 0){
          
          # Document in iter_log
          iter_log <- rbind(iter_log,
                            data.frame(gene_set = gene_set_name, iteration = n,
                                       before = length(iterated_final_gene_set),
                                       after = 0))
          # update iter log
          iter_log[nrow(iter_log), "after"] <- paste0("0 ; No stable clustering produced at first iteration. Skipped geneset.") 
          
          break  # exits the while loop
          
        }
        
      } else {
        # Get maximal most stable K  
        max_stable_k <- max(ks_without_na)
        message("maximum most stable K = ", max_stable_k)
        
      }
      
      # special condition for if maximum most stable K = 2 
      if (max_stable_k == 2){
        
        
        # let user know 
        message("since maximum most stable K = ", max_stable_k)
        message("Will have to extract k = 2 results using k = 3 analysis")
        
        # set that to working_k 
        working_max_k <- 2
        
        
        # Initialize starting value for pItem (percentage of items to sample per iteration)
        pItem_val <- 0.8
        
        # Set maximum allowable pItem value (cannot exceed 1.0)
        max_pItem <- 1.0
        
        # Step size for incrementing pItem when NA values are found
        step <- 0.05
        
        
        # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
        while (has_na && pItem_val <= max_pItem) {
          # Run ConsensusClusterPlus with the current pItem_val
          consensus_results <- ConsensusClusterPlus(
            input_matrix, 
            maxK = 3, 
            reps = CCP_iter,            
            pItem = pItem_val,            # Percentage of items sampled per iteration
            pFeature = 1,                # Use all features
            clusterAlg = clustering_alg,      
            innerLinkage = "ward.D2",
            finalLinkage = "ward.D2",
            distance = "euclidean", 
            seed = 2024L,
            plot ="none",           
            verbose = FALSE,
            writeTable = FALSE
          )
          
          # Only extract out result for K = 2
          consensus_results <- consensus_results[1:2]
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          # Check if any values in the item_consensus column are NA
          if (any(is.na(item_consensus_df$item_consensus))) {
            # Print message indicating NA detected and that pItem will be increased
            message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
            
            # Increase pItem by the predefined step to sample more items per iteration
            pItem_val <- pItem_val + step
          } else {
            # If no NA found, set flag to FALSE to exit the loop
            has_na <- FALSE
          }
        }
        
        
      } else {
        
        
        # set that to working_k 
        working_max_k <- max_stable_k
        
        # Initialize starting value for pItem (percentage of items to sample per iteration)
        pItem_val <- 0.8
        
        # Set maximum allowable pItem value (cannot exceed 1.0)
        max_pItem <- 1.0
        
        # Step size for incrementing pItem when NA values are found
        step <- 0.05
        
        
        # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
        while (has_na && pItem_val <= max_pItem) {
          # Run ConsensusClusterPlus with the current pItem_val
          consensus_results <- ConsensusClusterPlus(
            input_matrix, 
            maxK = working_max_k, 
            reps = CCP_iter,            
            pItem = pItem_val,            # Percentage of items sampled per iteration
            pFeature = 1,                # Use all features
            clusterAlg = clustering_alg,      
            innerLinkage = "ward.D2",
            finalLinkage = "ward.D2",
            distance = "euclidean", 
            seed = 2024L,
            plot ="none",           
            verbose = FALSE,
            writeTable = FALSE
          )
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          # Check if any values in the item_consensus column are NA
          if (any(is.na(item_consensus_df$item_consensus))) {
            # Print message indicating NA detected and that pItem will be increased
            message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
            
            # Increase pItem by the predefined step to sample more items per iteration
            pItem_val <- pItem_val + step
          } else {
            # If no NA found, set flag to FALSE to exit the loop
            has_na <- FALSE
          }
          
        }
        
      }
      
      
      
    }
    
    # If no NA found, proceed with downstream analysis steps
    if (!has_na) {
      # Classify stability based on item_consensus values using custom function
      item_consensus_df$stability <- sapply(item_consensus_df$item_consensus, get_stability_category)
      
      # Split the dataframe into a list of dataframes, one per cluster number k
      item_consensus_by_k <- split(item_consensus_df, item_consensus_df$k)
    }
    
    
    
    
    # Calculate PAC metric (Proportion of Ambiguous Clustering) using ecdf
    
    Kvec = 2:working_max_k                                  # Vector of cluster numbers to evaluate
    x1 = 0.1
    x2 = 0.9                                        # Thresholds defining "ambiguous" region
    PAC = rep(NA, length(Kvec))                     # Placeholder for PAC values
    names(PAC) = paste("K=", Kvec, sep = "")        # Label PAC entries with cluster numbers
    
    for (i in Kvec) {
      # Loop over cluster numbers
      M = consensus_results[[i]]$consensusMatrix    # Extract consensus matrix for current K
      Fn = ecdf(M[lower.tri(M)])                    # Empirical CDF of the lower triangle of consensus matrix
      PAC[i - 1] = Fn(x2) - Fn(x1)                  # The fraction of values between x1 and x2 (ambiguous assignments)
    }
    optK = Kvec[which.min(PAC)]                     # Select K with the lowest PAC (most stable clustering)
    
    
    # Convert PAC to data frame for plotting
    pac_values <- data.frame(k = Kvec, pac = unname(PAC))
    
    # Silhouette calculations using consensus cluster 
    
    input_matrix <- t(input_matrix)  # transpose to get dims as samples x features
    rownames(input_matrix) <- rownames(clustering_matrix)
    
    # Distances
    dist_space_for_clustering <- dist(input_matrix)
    dist_ge_space <- dist(clustering_matrix[,iterated_final_gene_set, drop=FALSE])
    
    # Collect silhouettes for all k
    sil_results <- lapply(2:length(consensus_results), function(k) {
      
      clusters <- consensus_results[[k]]$consensusClass
      
      # Silhouette in gene space
      sil_space_for_clustering   <- silhouette(clusters, dist_space_for_clustering)
      avg_space_for_clustering   <- mean(sil_space_for_clustering[, 3])
      
      # Silhouette in PaCMAP space
      sil_ge_space <- silhouette(clusters, dist_ge_space)
      avg_ge_space <- mean(sil_ge_space[, 3])
      
      data.frame(
        k = k,
        avg_sil_space_for_clustering = avg_space_for_clustering,
        avg_sil_ge_space = avg_ge_space
      )
    })
    
     # Combine
      sil_df <- do.call(rbind, sil_results)
      sil_df$avg_sil_combined <- ((sil_df$avg_sil_space_for_clustering +
        sil_df$avg_sil_ge_space) /
        2)

      #Get cluster_consensus_df
      cluster_consensus_df <- as.data.frame(icl_results[["clusterConsensus"]])

      # Create item consensus table
      avg_cluster_consensus <- cluster_consensus_df %>%
        group_by(k) %>%
        summarize(avg_cluster_consensus = mean(clusterConsensus))

      # Create metrics summary dataframe
      metrics_summary <- data.frame(
        k = 2:working_max_k,
        silhouette_dim_reduce_space = NA, # renamed to indicate dim‑reduce (original) space
        silhouette_ge_space = NA, # silhouette computed in PacMAP (ge) space
        silhouette_combined_avg = NA, # silhouette computed in PacMAP (ge) space
        pac = NA,
        avg_cluster_consensus = NA
      )

      # Fill in the metrics
      for (i in 1:nrow(metrics_summary)) {
        k_val <- metrics_summary$k[i]

        # Silhouette in the original (dim‑reduce) space
        idx <- which(sil_df$k == k_val)
        if (length(idx) > 0) {
          metrics_summary$silhouette_dim_reduce_space[
            i
          ] <- sil_df$avg_sil_space_for_clustering[idx]

          metrics_summary$silhouette_ge_space[i] <- sil_df$avg_sil_ge_space[idx]
          metrics_summary$silhouette_combined_avg[i] <- sil_df$avg_sil_combined[
            idx
          ]
        }

        # PAC
        k_idx <- which(Kvec == k_val)
        if (length(k_idx) > 0) {
          metrics_summary$pac[i] <- PAC[k_idx]
        }

        # avg_cluster_consensus
        cc_idx <- which(avg_cluster_consensus$k == k_val)
        if (length(cc_idx) > 0) {
          metrics_summary$avg_cluster_consensus[
            i
          ] <- avg_cluster_consensus$avg_cluster_consensus[cc_idx]
        }
      }

      # Alternatively,  Optimal K can be selected using rank based compromise
      # The rank-based compromise approach assigns a rank to each k for both silhouette (descending) and PAC (ascending), then sums the ranks to identify the k with the best overall trade-off between high silhouette and low PAC. The k with the lowest combined rank is selected as the optimal choice.
      # https://dl.acm.org/doi/10.1145/371920.372165

      # Assign ranks (lowest rank = best)
      metrics_summary$rank_silhouette_dim_reduce_space <- rank(
        -metrics_summary$silhouette_dim_reduce_space,
        ties.method = "min"
      ) # higher silhouette = rank 1 (original space)
      metrics_summary$rank_silhouette_ge_space <- rank(
        -metrics_summary$silhouette_ge_space,
        ties.method = "min"
      ) # higher silhouette = rank 1 (ge space)

      metrics_summary$rank_silhouette_combined_avg <- rank(
        -metrics_summary$silhouette_combined_avg,
        ties.method = "min"
      ) # higher silhouette = rank 1 (averaged space)
      metrics_summary$rank_pac <- rank(metrics_summary$pac, ties.method = "min") # lower PAC = rank 1
      metrics_summary$rank_item_cluster_consensus <- rank(
        -metrics_summary$avg_cluster_consensus,
        ties.method = "min"
      ) # higher consensus = rank 1

      # Build the column names dynamically.
      rank_columns <- c(
        paste0("rank_", clustering_metrics) # e.g. "rank_pac", "rank_item_cluster_consensus"
      )

      # Compute the overall rank as the row sum across selected rank columns
      metrics_summary$overall_rank <- rowSums(
        metrics_summary[, rank_columns, drop = FALSE],
        na.rm = TRUE
      )

      # Find all rows with the minimum overall rank
      best_rows <- metrics_summary[
        metrics_summary$overall_rank == min(metrics_summary$overall_rank),
      ]

      # Use the last metric in clustering_metrics as the tiebreaker
      tie_metric <- clustering_metrics[length(clustering_metrics)]

      if (verbose) {
        cat(
          "Metrics being used for ranking:\n  ",
          paste(clustering_metrics, collapse = ", "),
          "\n"
        )
        cat(
          "Tie‑breaker metric (last in clustering_metrics vector):",
          tie_metric,
          "\n"
        )
      }

      # Determine if the tiebreaking metric should be maximized or minimized
      # (Assuming: max for 'silhouette' and 'item_cluster_consensus', min for 'pac')
      if (tie_metric == "silhouette_dim_reduce_space") {
        idx <- which.max(best_rows$silhouette_dim_reduce_space)
      } else if (tie_metric == "silhouette_ge") {
        idx <- which.max(best_rows$silhouette_ge_space)
      } else if (tie_metric == "item_cluster_consensus") {
        idx <- which.max(best_rows$avg_cluster_consensus)
      } else if (tie_metric == "pac") {
        idx <- which.min(best_rows$pac)
      } else {
        stop("Unknown tiebreaker metric.")
      }

      best_k_row <- best_rows[idx, ]
      best_k <- best_k_row$k

      # add cluster metrics used for evaluation and tie
      metrics_summary$metrics_used <- paste(clustering_metrics, collapse = ", ")
      metrics_summary$metrics_used_for_tie <- tie_metric

      # Based on examining the plots and metrics
      optimal_k <- best_k
    
    # Extract cluster assignments for optimal k
    consensus_clusters <- consensus_results[[optimal_k]]$consensusClass
    
    
    # Get stability/reliability metrics 
    optimal_stability_cluster_metrics <- item_consensus_by_k[[as.character(optimal_k)]]
    
    # Start from your original data frame
    wide_optimal_stability_cluster_metrics <- optimal_stability_cluster_metrics %>%
      # Only keep relevant columns
      dplyr::select(sample_id, cluster, item_consensus) %>%
      # Pivot wider: each cluster makes a new column
      pivot_wider(
        names_from = cluster,
        values_from = item_consensus,
        names_prefix = "item_consensus_cluster"
      )
    
    
    # ---- assign labels from consensus plus clustering  ----
    
    pacmap_plot_df$kmeans_cluster <- consensus_clusters
    
    
    # ---- Nonparametric test and effect size for all genes in unified set ----
    Fvals <- numeric(length(iterated_final_gene_set))       # Store effect size (not F)
    Pvals <- numeric(length(iterated_final_gene_set))
    TestType <- character(length(iterated_final_gene_set))
    
    for (i in seq_along(iterated_final_gene_set)) {
      gene_name <- iterated_final_gene_set[i]
      vals <- pacmap_plot_df[[gene_name]]
      cluster <- as.factor(pacmap_plot_df[["kmeans_cluster"]])
      n_groups <- nlevels(cluster)
      df <- data.frame(vals = vals, cluster = cluster)
      
      if (n_groups == 2) {
        # Wilcoxon test and effect size
        wt <- wilcox.test(vals ~ cluster)
        eff <- tryCatch({
          ef <- rstatix::wilcox_effsize(df, vals ~ cluster, ci = FALSE)
          ef$effsize
        }, error = function(e) NA)
        Pvals[i] <- wt$p.value
        Fvals[i] <- eff
        TestType[i] <- "wilcox"
      } else if (n_groups >= 3) {
        # Kruskal test and effect size
        kw <- kruskal.test(vals ~ cluster)
        eff <- tryCatch({
          ef <- rstatix::kruskal_effsize(df, vals ~ cluster, ci = FALSE)
          ef$effsize
        }, error = function(e) NA)
        Pvals[i] <- kw$p.value
        Fvals[i] <- eff
        TestType[i] <- "kruskal"
      } else {
        Pvals[i] <- NA
        Fvals[i] <- NA
        TestType[i] <- "NA"
      }
    }
    
    # Adjusted P-values (FDR)
    adj_Pvals <- p.adjust(Pvals, method = "fdr")
    
    iterated_feature_df <- data.frame(
      gene_id = iterated_final_gene_set,
      effect_size = as.numeric(Fvals),
      p_value = as.numeric(Pvals),
      adj_p_value = as.numeric(adj_Pvals),
      test = TestType
    )
    
    old_iterated_feature_df <- iterated_feature_df
    iterated_feature_df <- type.convert(iterated_feature_df, as.is = TRUE)
    iterated_feature_df <- iterated_feature_df[iterated_feature_df$p_value <= 0.05 & iterated_feature_df$adj_p_value <= adj_pval_thresh, ]
    old_iterated_final_gene_set <- iterated_final_gene_set
    iterated_final_gene_set <- unique(sort(iterated_feature_df$gene_id))
    stats_per_iter[[n + 1]] <- iterated_feature_df
    cluster_metrics_per_iter[[n+1]] <- metrics_summary
    optimal_k_per_iter[[n + 1]] <- optimal_k
    cluster_stability_per_iter[[n + 1]] <- optimal_stability_cluster_metrics
    wide_cluster_stability_per_iter[[n + 1]] <- wide_optimal_stability_cluster_metrics
    
    iter_log <- rbind(
      iter_log,
      data.frame(
        gene_set = "All Genes", iteration = n,
        before = length(old_iterated_final_gene_set),
        after = length(iterated_final_gene_set))
    )
    
    # --- check convergence / next iteration ---
    gene_set_check_1 <- length(old_iterated_final_gene_set)
    gene_set_check_2 <- length(iterated_final_gene_set)
    if (verbose) message(sprintf("Analysis by all genes Iter %d End: %d genes (%d prev)", n, gene_set_check_2, gene_set_check_1))
    
    n <- n + 1
    if (n >= max_pipeline_iter) warning(sprintf("[analysis_by_all_genes] Reached max_pipeline_iter (%d) in final iteration", max_pipeline_iter))
  
    # Check if there is non-immediate convergence to zero 
    if(gene_set_check_2 == 0 && gene_set_check_1 < length(unique(unlist(all_gene_sets)))){
      
      cat("Non-immediate convergence to zero detected: reverting to most recent non-zero iteration")
      
      # Convert back to old geneset 
      iterated_final_gene_set <- old_iterated_final_gene_set
      
      # update last row to reflect immediate convergence iter log 
      iter_log[nrow(iter_log), "after"] <- paste0(length(iterated_final_gene_set), "; non-immediate convergence reverted to most recent non-zero iteration") 
      
      # break the loop
      break
    }
    
    
    # Check if there is immediate convergence to zero 
    if (gene_set_check_2 == 0 && n == 1) {
      
      cat("Immediate convergence to zero detected: will retain the top 10% of genes by pval then eff_size. Consider relaxing adj_pval_thresh\n")
      
      # Rank old gene set by adj_p_value (ascending) and effect_size (descending)
      ranked_genes <- old_iterated_feature_df[order(old_iterated_feature_df$adj_p_value,
                                                    -old_iterated_feature_df$effect_size), ]
      
      # Keep top 10%
      top_n <- ceiling(0.10 * nrow(ranked_genes))
      iterated_feature_df <- ranked_genes[seq_len(top_n), ]
      iterated_final_gene_set <- iterated_feature_df$gene_id
      
      # update last row to reflect immediate convergence iter log 
      iter_log[nrow(iter_log), "after"] <- paste0(length(iterated_final_gene_set), "; immediate convergence to zero retained top 10% by pval then eff_size") 
      
      # Break the loop
      break
    }
    
    }
  
  
  # ---- Save per-set details ----
  details_per_set_analysis_by_all_genes$final <- iterated_final_gene_set
  details_per_set_analysis_by_all_genes$feature_importance_per_iter <- stats_per_iter
  details_per_set_analysis_by_all_genes$optimal_k_per_iter <- optimal_k_per_iter
  details_per_set_analysis_by_all_genes$cluster_metrics_per_iter <- cluster_metrics_per_iter
  details_per_set_analysis_by_all_genes$cluster_stability_per_iter <- cluster_stability_per_iter
  details_per_set_analysis_by_all_genes$wide_cluster_stability_per_iter <- wide_cluster_stability_per_iter
  
  iterated_over_all_genes <- iterated_final_gene_set
  
  # Check if Gene set is empty
  if (length(iterated_over_all_genes) == 0) {
    cat("\033[1;31mWARNING: Final Gene Set for approach two empty! Consider Relaxing P-val\nApproach output will not be produced\033[0m\n")
    
    # toggle Switch
    no_app_two <- TRUE
    
  }
  

  
  })
  
    }
  
  

  # ---- Approach three: iterate over pooled union of all sets ----
  
  
  if ("app_three" %in% filter_approach){
  
    # Set seed appropriately
    withr::with_seed(seed = 2024L, {
    
    if (verbose) message("\n[Third Approach: use final gene set from approach one ]")
  iterated_final_gene_set <- new_final
  gene_set_check_1 <- 0
  gene_set_check_2 <- length(iterated_final_gene_set)
  n <- 0
  overall_stats_per_iter <- list()
  details_per_set_analysis_by_gene_set_then_all_genes <- list()
  stats_per_iter <- list()
  cluster_metrics_per_iter <- list()
  optimal_k_per_iter <- list()
  cluster_stability_per_iter <- list()
  wide_cluster_stability_per_iter <- list()
  
  
  while (gene_set_check_1 != gene_set_check_2 && n < max_pipeline_iter) {
  
    if (verbose)
      message(
        sprintf(
          "Use final gene set from approach one analysis
                                 Iter %d Start: %d genes (%d prev)",
          n,
          gene_set_check_2,
          gene_set_check_1
        )
      )
    
    # set max_k for run
    working_max_k <- max_k
    
    # --- Run PaCMAP; attempt safely in case of too-small gene sets ---
   
    # Build embedder for analysis (also resets seed)
    script_pacmap_reducer <- do.call(pacmap$PaCMAP, pm_defaults)
    
     safe_pacmap <- tryCatch({
      
      # Run PaCMAP
      embedding <- script_pacmap_reducer$fit_transform(
        clustering_matrix[, iterated_final_gene_set, drop = FALSE]
      )
      embedding
    }, error = function(e) {
      message(sprintf("PaCMAP failed for gene set %s iteration %d. Skipping this set. Most iterated gene set will be kept (if iter > 0)", gene_set_name, n))
      NULL
    })
    
    if (is.null(safe_pacmap)) {
      
      details_per_set_analysis_by_gene_set_then_all_genes$final <- iterated_final_gene_set
      details_per_set_analysis_by_gene_set_then_all_genes$feature_importance_per_iter <- stats_per_iter
      details_per_set_analysis_by_gene_set_then_all_genes$optimal_k_per_iter <- optimal_k_per_iter
      details_per_set_analysis_by_gene_set_then_all_genes$cluster_metrics_per_iter <- cluster_metrics_per_iter
      details_per_set_analysis_by_gene_set_then_all_genes$cluster_stability_per_iter <- cluster_stability_per_iter
      details_per_set_analysis_by_gene_set_then_all_genes$wide_cluster_stability_per_iter <- wide_cluster_stability_per_iter
      
      iterated_use_final_gene_set_from_approach_one <- iterated_final_gene_set
      
      break
    }
    
    rownames(safe_pacmap) <- rownames(clustering_matrix)
    
    # ---- Prepare PaCMAP result for clustering and feature selection ----
    pacmap_plot_df <- data.frame(safe_pacmap) %>%
      tibble::rownames_to_column("sample_id") %>%
      dplyr::inner_join(as.data.frame(clustering_matrix) %>% rownames_to_column("sample_id"), by = "sample_id")
    
    
    # ---- Determine optimal number of clusters via Consensus Clustering ----
    
    # Your input matrix: samples x features (PaCMAP embedding)
    input_matrix <- t(safe_pacmap)  # transpose to get dims as rows (like you did with UMAP)
    colnames(input_matrix) <- rownames(clustering_matrix)
    
    
    
    # Initialize starting value for pItem (percentage of items to sample per iteration)
    pItem_val <- 0.8
    
    # Set maximum allowable pItem value (cannot exceed 1.0)
    max_pItem <- 1.0
    
    # Step size for incrementing pItem when NA values are found
    step <- 0.05
    
    # Flag to track presence of NA values in itemConsensus
    has_na <- TRUE
    
    # if in the case of working max k = 2 
    if (working_max_k == 2){
      
      # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
      while (has_na && pItem_val <= max_pItem) {
        # Run ConsensusClusterPlus with the current pItem_val
        consensus_results <- ConsensusClusterPlus(
          input_matrix, 
          maxK = 3, 
          reps = CCP_iter,            
          pItem = pItem_val,            # Percentage of items sampled per iteration
          pFeature = 1,                # Use all features
          clusterAlg = clustering_alg,      
          innerLinkage = "ward.D2",
          finalLinkage = "ward.D2",
          distance = "euclidean", 
          seed = 2024L,
          plot ="none",           
          verbose = FALSE,
          writeTable = FALSE
        )
        
        # Only extract out result for K = 2
        consensus_results <- consensus_results[1:2]
        
        # Calculate item consensus and cluster metrics (ICL)
        pdf(file = NULL)  # Start invisible plotting device
        icl_results <- calcICL(consensus_results, 
                               plot = "png", 
                               writeTable = TRUE)
        dev.off()          # Close device
        
        # Extract the itemConsensus dataframe from ICL results
        item_consensus_df <- icl_results[["itemConsensus"]]
        
        # Rename columns for clarity
        colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
        
        
        # Calculate item consensus and cluster metrics (ICL)
        pdf(file = NULL)  # Start invisible plotting device
        icl_results <- calcICL(consensus_results, 
                               plot = "png", 
                               writeTable = TRUE)
        dev.off()          # Close device
        
        # Extract the itemConsensus dataframe from ICL results
        item_consensus_df <- icl_results[["itemConsensus"]]
        
        # Rename columns for clarity
        colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
        
        # Check if any values in the item_consensus column are NA
        if (any(is.na(item_consensus_df$item_consensus))) {
          # Print message indicating NA detected and that pItem will be increased
          message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
          
          # Increase pItem by the predefined step to sample more items per iteration
          pItem_val <- pItem_val + step
        } else {
          # If no NA found, set flag to FALSE to exit the loop
          has_na <- FALSE
        }
        
      } 
      
    }
    else {
      
      # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
      while (has_na && pItem_val <= max_pItem) {
        
        # Run ConsensusClusterPlus with the current pItem_val
        consensus_results <- ConsensusClusterPlus(
          input_matrix, 
          maxK = working_max_k, 
          reps = CCP_iter,            
          pItem = pItem_val,            # Percentage of items sampled per iteration
          pFeature = 1,                # Use all features
          clusterAlg = clustering_alg,      
          innerLinkage = "ward.D2",
          finalLinkage = "ward.D2",
          distance = "euclidean", 
          seed = 2024L,
          plot ="none",           
          verbose = FALSE,
          writeTable = FALSE
        )
        
        # Calculate item consensus and cluster metrics (ICL)
        pdf(file = NULL)  # Start invisible plotting device
        icl_results <- calcICL(consensus_results, 
                               plot = "png", 
                               writeTable = TRUE)
        dev.off()          # Close device
        
        # Extract the itemConsensus dataframe from ICL results
        item_consensus_df <- icl_results[["itemConsensus"]]
        
        # Rename columns for clarity
        colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
        
        # Check if any values in the item_consensus column are NA
        if (any(is.na(item_consensus_df$item_consensus))) {
          # Print message indicating NA detected and that pItem will be increased
          message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
          
          # Increase pItem by the predefined step to sample more items per iteration
          pItem_val <- pItem_val + step
        } else {
          # If no NA found, set flag to FALSE to exit the loop
          has_na <- FALSE
        }
      }
    }
    
    # After loop finishes, check if NA values still exist (meaning max pItem was reached)
    if (has_na) {
      # Warn the user that maximum pItem was reached but NAs persist
      message("Reached max pItem but NA values still exist in item_consensus.")
      message("Reducing working_max_k to the maximum most stable K")
      
      # Find all K values with no NA item consensus
      ks_without_na <- item_consensus_df %>%
        group_by(k) %>%
        summarize(no_na = all(!is.na(item_consensus))) %>%
        filter(no_na) %>%
        pull(k)
      
      
      # if there is no solution skip it 
      if (length(ks_without_na) == 0) {
        
        message("No stable K values found. Skipping this geneset/iteration. \n Will keep most iterated if iteration > 0")
        
        if(n >= 1){
          
          # Document in iter_log
          iter_log <- rbind(iter_log,
                            data.frame(gene_set = gene_set_name, iteration = n,
                                       before = length(old_iterated_final_gene_set),
                                       after = length(iterated_final_gene_set)))
          # update iter log
          iter_log[nrow(iter_log), "after"] <- paste0(length(iterated_final_gene_set), "; No stable clustering produced at higher iterations. Took most iterated geneset.") 
          
          break  # exits the while loop
          
        }
        
        
        if(n == 0){
          
          # Document in iter_log
          iter_log <- rbind(iter_log,
                            data.frame(gene_set = gene_set_name, iteration = n,
                                       before = length(iterated_final_gene_set),
                                       after = 0))
          # update iter log
          iter_log[nrow(iter_log), "after"] <- paste0("0 ; No stable clustering produced at first iteration. Skipped geneset.") 
          
          break  # exits the while loop
          
        }
        
      } else {
        # Get maximal most stable K  
        max_stable_k <- max(ks_without_na)
        message("maximum most stable K = ", max_stable_k)
        
      }
      
      # special condition for if maximum most stable K = 2 
      if (max_stable_k == 2){
        
        
        # let user know 
        message("since maximum most stable K = ", max_stable_k)
        message("Will have to extract k = 2 results using k = 3 analysis")
        
        # set that to working_k 
        working_max_k <- 2
        
        
        # Initialize starting value for pItem (percentage of items to sample per iteration)
        pItem_val <- 0.8
        
        # Set maximum allowable pItem value (cannot exceed 1.0)
        max_pItem <- 1.0
        
        # Step size for incrementing pItem when NA values are found
        step <- 0.05
        
        
        # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
        while (has_na && pItem_val <= max_pItem) {
          # Run ConsensusClusterPlus with the current pItem_val
          consensus_results <- ConsensusClusterPlus(
            input_matrix, 
            maxK = 3, 
            reps = CCP_iter,            
            pItem = pItem_val,            # Percentage of items sampled per iteration
            pFeature = 1,                # Use all features
            clusterAlg = clustering_alg,      
            innerLinkage = "ward.D2",
            finalLinkage = "ward.D2",
            distance = "euclidean", 
            seed = 2024L,
            plot ="none",           
            verbose = FALSE,
            writeTable = FALSE
          )
          
          # Only extract out result for K = 2
          consensus_results <- consensus_results[1:2]
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          # Check if any values in the item_consensus column are NA
          if (any(is.na(item_consensus_df$item_consensus))) {
            # Print message indicating NA detected and that pItem will be increased
            message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
            
            # Increase pItem by the predefined step to sample more items per iteration
            pItem_val <- pItem_val + step
          } else {
            # If no NA found, set flag to FALSE to exit the loop
            has_na <- FALSE
          }
        }
        
        
      } else {
        
        
        # set that to working_k 
        working_max_k <- max_stable_k
        
        # Initialize starting value for pItem (percentage of items to sample per iteration)
        pItem_val <- 0.8
        
        # Set maximum allowable pItem value (cannot exceed 1.0)
        max_pItem <- 1.0
        
        # Step size for incrementing pItem when NA values are found
        step <- 0.05
        
        
        # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
        while (has_na && pItem_val <= max_pItem) {
          # Run ConsensusClusterPlus with the current pItem_val
          consensus_results <- ConsensusClusterPlus(
            input_matrix, 
            maxK = working_max_k, 
            reps = CCP_iter,            
            pItem = pItem_val,            # Percentage of items sampled per iteration
            pFeature = 1,                # Use all features
            clusterAlg = clustering_alg,      
            innerLinkage = "ward.D2",
            finalLinkage = "ward.D2",
            distance = "euclidean", 
            seed = 2024L,
            plot ="none",           
            verbose = FALSE,
            writeTable = FALSE
          )
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          # Check if any values in the item_consensus column are NA
          if (any(is.na(item_consensus_df$item_consensus))) {
            # Print message indicating NA detected and that pItem will be increased
            message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
            
            # Increase pItem by the predefined step to sample more items per iteration
            pItem_val <- pItem_val + step
          } else {
            # If no NA found, set flag to FALSE to exit the loop
            has_na <- FALSE
          }
          
        }
        
      }
      
      
      
    }
    
    # If no NA found, proceed with downstream analysis steps
    if (!has_na) {
      # Classify stability based on item_consensus values using custom function
      item_consensus_df$stability <- sapply(item_consensus_df$item_consensus, get_stability_category)
      
      # Split the dataframe into a list of dataframes, one per cluster number k
      item_consensus_by_k <- split(item_consensus_df, item_consensus_df$k)
    }
    
    
    
    # Calculate PAC metric (Proportion of Ambiguous Clustering) using ecdf
    
    Kvec = 2:working_max_k                                  # Vector of cluster numbers to evaluate
    x1 = 0.1
    x2 = 0.9                                        # Thresholds defining "ambiguous" region
    PAC = rep(NA, length(Kvec))                     # Placeholder for PAC values
    names(PAC) = paste("K=", Kvec, sep = "")        # Label PAC entries with cluster numbers
    
    for (i in Kvec) {
      # Loop over cluster numbers
      M = consensus_results[[i]]$consensusMatrix    # Extract consensus matrix for current K
      Fn = ecdf(M[lower.tri(M)])                    # Empirical CDF of the lower triangle of consensus matrix
      PAC[i - 1] = Fn(x2) - Fn(x1)                  # The fraction of values between x1 and x2 (ambiguous assignments)
    }
    optK = Kvec[which.min(PAC)]                     # Select K with the lowest PAC (most stable clustering)
    
    
    # Convert PAC to data frame for plotting
    pac_values <- data.frame(k = Kvec, pac = unname(PAC))
    
  
    # Silhouette calculations using consensus cluster 
    
    input_matrix <- t(input_matrix)  # transpose to get dims as samples x features
    rownames(input_matrix) <- rownames(clustering_matrix)
    
    # Distances
    dist_space_for_clustering <- dist(input_matrix)
    dist_ge_space <- dist(clustering_matrix[,iterated_final_gene_set, drop=FALSE])
    
    # Collect silhouettes for all k
    sil_results <- lapply(2:length(consensus_results), function(k) {
      
      clusters <- consensus_results[[k]]$consensusClass
      
      # Silhouette in gene space
      sil_space_for_clustering   <- silhouette(clusters, dist_space_for_clustering)
      avg_space_for_clustering   <- mean(sil_space_for_clustering[, 3])
      
      # Silhouette in PaCMAP space
      sil_ge_space <- silhouette(clusters, dist_ge_space)
      avg_ge_space <- mean(sil_ge_space[, 3])
      
      data.frame(
        k = k,
        avg_sil_space_for_clustering = avg_space_for_clustering,
        avg_sil_ge_space = avg_ge_space
      )
    })
    
     # Combine
      sil_df <- do.call(rbind, sil_results)
      sil_df$avg_sil_combined <- ((sil_df$avg_sil_space_for_clustering +
        sil_df$avg_sil_ge_space) /
        2)

      #Get cluster_consensus_df
      cluster_consensus_df <- as.data.frame(icl_results[["clusterConsensus"]])

      # Create item consensus table
      avg_cluster_consensus <- cluster_consensus_df %>%
        group_by(k) %>%
        summarize(avg_cluster_consensus = mean(clusterConsensus))

      # Create metrics summary dataframe
      metrics_summary <- data.frame(
        k = 2:working_max_k,
        silhouette_dim_reduce_space = NA, # renamed to indicate dim‑reduce (original) space
        silhouette_ge_space = NA, # silhouette computed in PacMAP (ge) space
        silhouette_combined_avg = NA, # silhouette computed in PacMAP (ge) space
        pac = NA,
        avg_cluster_consensus = NA
      )

      # Fill in the metrics
      for (i in 1:nrow(metrics_summary)) {
        k_val <- metrics_summary$k[i]

        # Silhouette in the original (dim‑reduce) space
        idx <- which(sil_df$k == k_val)
        if (length(idx) > 0) {
          metrics_summary$silhouette_dim_reduce_space[
            i
          ] <- sil_df$avg_sil_space_for_clustering[idx]

          metrics_summary$silhouette_ge_space[i] <- sil_df$avg_sil_ge_space[idx]
          metrics_summary$silhouette_combined_avg[i] <- sil_df$avg_sil_combined[
            idx
          ]
        }

        # PAC
        k_idx <- which(Kvec == k_val)
        if (length(k_idx) > 0) {
          metrics_summary$pac[i] <- PAC[k_idx]
        }

        # avg_cluster_consensus
        cc_idx <- which(avg_cluster_consensus$k == k_val)
        if (length(cc_idx) > 0) {
          metrics_summary$avg_cluster_consensus[
            i
          ] <- avg_cluster_consensus$avg_cluster_consensus[cc_idx]
        }
      }

      # Alternatively,  Optimal K can be selected using rank based compromise
      # The rank-based compromise approach assigns a rank to each k for both silhouette (descending) and PAC (ascending), then sums the ranks to identify the k with the best overall trade-off between high silhouette and low PAC. The k with the lowest combined rank is selected as the optimal choice.
      # https://dl.acm.org/doi/10.1145/371920.372165

      # Assign ranks (lowest rank = best)
      metrics_summary$rank_silhouette_dim_reduce_space <- rank(
        -metrics_summary$silhouette_dim_reduce_space,
        ties.method = "min"
      ) # higher silhouette = rank 1 (original space)
      metrics_summary$rank_silhouette_ge_space <- rank(
        -metrics_summary$silhouette_ge_space,
        ties.method = "min"
      ) # higher silhouette = rank 1 (ge space)

      metrics_summary$rank_silhouette_combined_avg <- rank(
        -metrics_summary$silhouette_combined_avg,
        ties.method = "min"
      ) # higher silhouette = rank 1 (averaged space)
      metrics_summary$rank_pac <- rank(metrics_summary$pac, ties.method = "min") # lower PAC = rank 1
      metrics_summary$rank_item_cluster_consensus <- rank(
        -metrics_summary$avg_cluster_consensus,
        ties.method = "min"
      ) # higher consensus = rank 1

      # Build the column names dynamically.
      rank_columns <- c(
        paste0("rank_", clustering_metrics) # e.g. "rank_pac", "rank_item_cluster_consensus"
      )

      # Compute the overall rank as the row sum across selected rank columns
      metrics_summary$overall_rank <- rowSums(
        metrics_summary[, rank_columns, drop = FALSE],
        na.rm = TRUE
      )

      # Find all rows with the minimum overall rank
      best_rows <- metrics_summary[
        metrics_summary$overall_rank == min(metrics_summary$overall_rank),
      ]

      # Use the last metric in clustering_metrics as the tiebreaker
      tie_metric <- clustering_metrics[length(clustering_metrics)]

      if (verbose) {
        cat(
          "Metrics being used for ranking:\n  ",
          paste(clustering_metrics, collapse = ", "),
          "\n"
        )
        cat(
          "Tie‑breaker metric (last in clustering_metrics vector):",
          tie_metric,
          "\n"
        )
      }

      # Determine if the tiebreaking metric should be maximized or minimized
      # (Assuming: max for 'silhouette' and 'item_cluster_consensus', min for 'pac')
      if (tie_metric == "silhouette_dim_reduce_space") {
        idx <- which.max(best_rows$silhouette_dim_reduce_space)
      } else if (tie_metric == "silhouette_ge") {
        idx <- which.max(best_rows$silhouette_ge_space)
      } else if (tie_metric == "item_cluster_consensus") {
        idx <- which.max(best_rows$avg_cluster_consensus)
      } else if (tie_metric == "pac") {
        idx <- which.min(best_rows$pac)
      } else {
        stop("Unknown tiebreaker metric.")
      }

      best_k_row <- best_rows[idx, ]
      best_k <- best_k_row$k

      # add cluster metrics used for evaluation and tie
      metrics_summary$metrics_used <- paste(clustering_metrics, collapse = ", ")
      metrics_summary$metrics_used_for_tie <- tie_metric

      # Based on examining the plots and metrics
      optimal_k <- best_k
    
    # Extract cluster assignments for optimal k
    consensus_clusters <- consensus_results[[optimal_k]]$consensusClass
    
    
    # Get stability/reliability metrics 
    optimal_stability_cluster_metrics <- item_consensus_by_k[[as.character(optimal_k)]]
    
    # Start from your original data frame
    wide_optimal_stability_cluster_metrics <- optimal_stability_cluster_metrics %>%
      # Only keep relevant columns
      dplyr::select(sample_id, cluster, item_consensus) %>%
      # Pivot wider: each cluster makes a new column
      pivot_wider(
        names_from = cluster,
        values_from = item_consensus,
        names_prefix = "item_consensus_cluster"
      )
    
    
    # ---- assign labels from consensus plus clustering  ----
    
    pacmap_plot_df$kmeans_cluster <- consensus_clusters
    
    
    # ---- Nonparametric test and effect size for all genes in unified set ----
    Fvals <- numeric(length(iterated_final_gene_set))       # Store effect size (not F)
    Pvals <- numeric(length(iterated_final_gene_set))
    TestType <- character(length(iterated_final_gene_set))
    
    for (i in seq_along(iterated_final_gene_set)) {
      gene_name <- iterated_final_gene_set[i]
      vals <- pacmap_plot_df[[gene_name]]
      cluster <- as.factor(pacmap_plot_df[["kmeans_cluster"]])
      n_groups <- nlevels(cluster)
      df <- data.frame(vals = vals, cluster = cluster)
      
      if (n_groups == 2) {
        # Wilcoxon test and effect size
        wt <- wilcox.test(vals ~ cluster)
        eff <- tryCatch({
          ef <- rstatix::wilcox_effsize(df, vals ~ cluster, ci = FALSE)
          ef$effsize
        }, error = function(e) NA)
        Pvals[i] <- wt$p.value
        Fvals[i] <- eff
        TestType[i] <- "wilcox"
      } else if (n_groups >= 3) {
        # Kruskal test and effect size
        kw <- kruskal.test(vals ~ cluster)
        eff <- tryCatch({
          ef <- rstatix::kruskal_effsize(df, vals ~ cluster, ci = FALSE)
          ef$effsize
        }, error = function(e) NA)
        Pvals[i] <- kw$p.value
        Fvals[i] <- eff
        TestType[i] <- "kruskal"
      } else {
        Pvals[i] <- NA
        Fvals[i] <- NA
        TestType[i] <- "NA"
      }
    }
    
    # Adjusted P-values (FDR)
    adj_Pvals <- p.adjust(Pvals, method = "fdr")
    
    iterated_feature_df <- data.frame(
      gene_id = iterated_final_gene_set,
      effect_size = as.numeric(Fvals),
      p_value = as.numeric(Pvals),
      adj_p_value = as.numeric(adj_Pvals),
      test = TestType
    )
    
    old_iterated_feature_df <- iterated_feature_df
    iterated_feature_df <- type.convert(iterated_feature_df, as.is = TRUE)
    iterated_feature_df <- iterated_feature_df[iterated_feature_df$p_value <= 0.05 & iterated_feature_df$adj_p_value <= adj_pval_thresh, ]
    old_iterated_final_gene_set <- iterated_final_gene_set
    iterated_final_gene_set <- unique(sort(iterated_feature_df$gene_id))
    stats_per_iter[[n + 1]] <- iterated_feature_df
    cluster_metrics_per_iter[[n+1]] <- metrics_summary
    optimal_k_per_iter[[n + 1]] <- optimal_k
    cluster_stability_per_iter[[n + 1]] <- optimal_stability_cluster_metrics
    wide_cluster_stability_per_iter[[n + 1]] <- wide_optimal_stability_cluster_metrics
    
    iter_log <- rbind(
      iter_log,
      data.frame(
        gene_set = "Use final gene set from approach one", iteration = n,
        before = length(old_iterated_final_gene_set),
        after = length(iterated_final_gene_set))
    )
    
    # --- check convergence / next iteration ---
    gene_set_check_1 <- length(old_iterated_final_gene_set)
    gene_set_check_2 <- length(iterated_final_gene_set)
    if (verbose)
      message(
        sprintf(
          "Use final gene set from approach one analysis
                                 Iter %d End: %d genes (%d prev)",
          n,
          gene_set_check_2,
          gene_set_check_1
        )
      )
    
    
    n <- n + 1
    if (n >= max_pipeline_iter) warning(sprintf("[analysis by using final gene set from approach one] Reached max_pipeline_iter (%d) in final iteration", max_pipeline_iter))
  
    # Check if there is non-immediate convergence to zero 
    if(gene_set_check_2 == 0 && gene_set_check_1 < length(new_final)){
      
      cat("Non-immediate convergence to zero detected: reverting to most recent non-zero iteration")
      
      # Convert back to old geneset 
      iterated_final_gene_set <- old_iterated_final_gene_set
      
      # update last row to reflect immediate convergence iter log 
      iter_log[nrow(iter_log), "after"] <- paste0(length(iterated_final_gene_set), "; non-immediate convergence reverted to most recent non-zero iteration") 
      
      
      # break the loop
      break
    }
    
    
    # Check if there is immediate convergence to zero 
    if (gene_set_check_2 == 0 && n == 1) {
      
      cat("Immediate convergence to zero detected: will retain the top 10% of genes by pval then eff_size. Consider relaxing adj_pval_thresh\n")
      
      # Rank old gene set by adj_p_value (ascending) and effect_size (descending)
      ranked_genes <- old_iterated_feature_df[order(old_iterated_feature_df$adj_p_value,
                                                    -old_iterated_feature_df$effect_size), ]
      
      # Keep top 10%
      top_n <- ceiling(0.10 * nrow(ranked_genes))
      iterated_feature_df <- ranked_genes[seq_len(top_n), ]
      iterated_final_gene_set <- iterated_feature_df$gene_id
      
      # update last row to reflect immediate convergence iter log 
      iter_log[nrow(iter_log), "after"] <- paste0(length(iterated_final_gene_set), "; immediate convergence to zero retained top 10% by pval then eff_size") 
      
      # Break the loop
      break
    }
    
    }
  
  
  # ---- Save per-set details ----
  details_per_set_analysis_by_gene_set_then_all_genes$final <- iterated_final_gene_set
  details_per_set_analysis_by_gene_set_then_all_genes$feature_importance_per_iter <- stats_per_iter
  details_per_set_analysis_by_gene_set_then_all_genes$optimal_k_per_iter <- optimal_k_per_iter
  details_per_set_analysis_by_gene_set_then_all_genes$cluster_metrics_per_iter <- cluster_metrics_per_iter
  details_per_set_analysis_by_gene_set_then_all_genes$cluster_stability_per_iter <- cluster_stability_per_iter
  details_per_set_analysis_by_gene_set_then_all_genes$wide_cluster_stability_per_iter <- wide_cluster_stability_per_iter
  
  iterated_use_final_gene_set_from_approach_one <- iterated_final_gene_set
  
  # Check if Gene set is empty
  if (length(iterated_use_final_gene_set_from_approach_one) == 0) {
    cat("\033[1;31mWARNING: Final Gene Set for Approach three empty! Consider Relaxing P-val\nApproach output will not be produced\033[0m\n")
    
    # toggle Switch
    no_app_three <- TRUE
    
  }
  
  })
    
  }
  
  # ---- Final Visualization -----
  
  # ---- Approach One (Iterated over all gene sets) ---- 
  
  
  
  if ("app_one" %in% filter_approach){
    
    # Set seed appropriately
    withr::with_seed(seed = 2024L, {
    
  if(no_app_one == FALSE){
  
  # set back to max_k
  working_max_k <- max_k
  
  {   
    # CHOOSE INPUT DATA TYPE FOR CLUSTERING
    # Set to "gene_set" for original gene expression data, or "umap" for UMAP coordinates, or "pca"
    input_data_type <- "pacmap"  # <-- MODIFY HERE: use "genes" or "umap" or "pca" or "pacmap"
    
    # Make matrix for delta
    #toggle here
    clustering_gene_set <- new_final
    clustering_gene_set <- clustering_gene_set[clustering_gene_set %in% colnames(clustering_matrix)]
    
    
    
    # Run UMAP on the selected gene set (same as in the original code)
    # will always be 2D 
    umap_matrix_visual <- umap::umap(clustering_matrix[,clustering_gene_set, drop=FALSE],
                              config=custom.config.visual)
    
    
    # Build embedder for analysis 
    script_pacmap_reducer <- do.call(pacmap$PaCMAP, pm_defaults)
    
    
    
    
    # Run PaCMAP
    pacmap_delta_embeddings <- script_pacmap_reducer$fit_transform(
      clustering_matrix[, clustering_gene_set, drop = FALSE]
    )
    
    # set rownames 
    rownames(pacmap_delta_embeddings) <- rownames(clustering_matrix)
    
  
    # 2D visualizations 
    # create reducer (also resets seed)
    pm_visualization_reducer <- do.call(pacmap$PaCMAP, pm_visualization)
    
    # Run PACMAP on selected gene set
    # Run PaCMAP
    pacmap_delta_embeddings_visual <- pm_visualization_reducer$fit_transform(
      clustering_matrix[, clustering_gene_set, drop = FALSE]
    )
    
    
    # set rownames 
    rownames(pacmap_delta_embeddings_visual) <- rownames(clustering_matrix)
    
    # Your input matrix: samples x genes
    pca_matrix <- clustering_matrix[, clustering_gene_set]
    
    pca_result <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)
    pca_scores <- pca_result$x[, 1:3, drop = FALSE]
    
    
    
    if (verbose) message(sprintf(
      "\n\n=== Creating Output for Approach One using %s algorithm and %s as input ===",
      clustering_alg,
      input_data_type
    ))
    
    # Define color palette for clusters
    cols <- hue_pal()(working_max_k)
    
    ########################################################################################
    
    
    # Create input matrix for consensus clustering based on chosen data type
    if (input_data_type == "genes") {
      # Use gene data
      final_input_matrix <- t(clustering_matrix[, clustering_gene_set])  # transpose to get genes as rows
      cat("Using gene expression data for clustering.\n")
    } else if (input_data_type == "umap") {
      # Use UMAP coordinates - specifically from umap_matrix to match your original approach
      final_input_matrix <- t(umap_matrix$layout)  # transpose to get dimensions as rows
      cat("Using UMAP coordinates for clustering.\n")
    } else if (input_data_type == "pca"){
      
      # Your input matrix: samples x genes
      final_input_matrix <- clustering_matrix[, clustering_gene_set]
      
      # Make sure your final_input_matrix is samples x features (rows=samples, cols=genes)
      # You can check with dim(final_input_matrix)
      
      # 1. Run Parallel Analysis to determine number of PCs
      # PA expects samples in rows, variables in columns
      centile_seq <- seq(95, 80, by = -5)
      max_pcs <- 3   # Default value in case paranormal fails
      pa_result <- NULL
      for (cent in centile_seq) {
        pa_result <- try(
          suppressMessages(suppressWarnings(
            paran(final_input_matrix,
                  iterations = 1000,
                  centile = cent,
                  graph = FALSE,
                  quietly = TRUE)
          )), silent = TRUE)
        
        if (!inherits(pa_result, "try-error")) {
          message(sprintf("paran succeeded with centile = %d", cent))
          max_pcs <- pa_result$Retained  # number of components retained by PA
          break
        }
      }
      if (inherits(pa_result, "try-error")) {
        warning("paran failed at all tested centile values (95 to 80). Using default max_pcs = 3 for this gene set.")
        # continue; max_pcs remains 3
      }
      
      cat("Number of PCs retained by Parallel Analysis:", max_pcs, "\n")
      
      
      # 2. PCA projection using selected PCs
      pca_result <- prcomp(final_input_matrix, center = TRUE, scale. = TRUE)
      pca_scores <- pca_result$x[, 1:max_pcs, drop = FALSE]
      
      # 3. Prepare input for downstream analysis, e.g., clustering
      # Here, PCs x Samples (rows = PCs, cols = samples)
      final_input_matrix <- as.matrix(t(pca_scores))
      colnames(final_input_matrix) <- umap_plot_df$sample_id  # optional, set sample IDs
      
    } else if (input_data_type == "pacmap"){
      
      final_input_matrix <- t(pacmap_delta_embeddings)
      colnames(final_input_matrix) <- rownames(clustering_matrix)
    } else {
      stop("Invalid input_data_type. Use 'gene_set' or 'umap' or 'pca' or 'pacmap'.")
    }
    
    
    
    # Initialize starting value for pItem (percentage of items to sample per iteration)
    pItem_val <- 0.8
    
    # Set maximum allowable pItem value (cannot exceed 1.0)
    max_pItem <- 1.0
    
    # Step size for incrementing pItem when NA values are found
    step <- 0.05
    
    # Flag to track presence of NA values in itemConsensus
    has_na <- TRUE
    
    # if in the case of working max k = 2 
    if (working_max_k == 2){
      
      # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
      while (has_na && pItem_val <= max_pItem) {
        # Run ConsensusClusterPlus with the current pItem_val
        consensus_results_app_one <- ConsensusClusterPlus(
          final_input_matrix, 
          maxK = 3, 
          reps = CCP_iter,            
          pItem = pItem_val,            # Percentage of items sampled per iteration
          pFeature = 1,                # Use all features
          clusterAlg = clustering_alg,      
          innerLinkage = "ward.D2",
          finalLinkage = "ward.D2",
          distance = "euclidean", 
          seed = 2024L,
          plot ="none",           
          verbose = FALSE,
          writeTable = FALSE
        )
        
        # Only extract out result for K = 2
        consensus_results_app_one <- consensus_results_app_one[1:2]
        
        # Calculate item consensus and cluster metrics (ICL)
        pdf(file = NULL)  # Start invisible plotting device
        icl_results <- calcICL(consensus_results_app_one, 
                               plot = "png", 
                               writeTable = TRUE)
        dev.off()          # Close device
        
        # Extract the itemConsensus dataframe from ICL results
        item_consensus_df <- icl_results[["itemConsensus"]]
        
        # Rename columns for clarity
        colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
        
        
        # Calculate item consensus and cluster metrics (ICL)
        pdf(file = NULL)  # Start invisible plotting device
        icl_results <- calcICL(consensus_results_app_one, 
                               plot = "png", 
                               writeTable = TRUE)
        dev.off()          # Close device
        
        # Extract the itemConsensus dataframe from ICL results
        item_consensus_df <- icl_results[["itemConsensus"]]
        
        # Rename columns for clarity
        colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
        
        # Check if any values in the item_consensus column are NA
        if (any(is.na(item_consensus_df$item_consensus))) {
          # Print message indicating NA detected and that pItem will be increased
          message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
          
          # Increase pItem by the predefined step to sample more items per iteration
          pItem_val <- pItem_val + step
        } else {
          # If no NA found, set flag to FALSE to exit the loop
          has_na <- FALSE
        }
        
      } 
      
    }
    else {
      
      # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
      while (has_na && pItem_val <= max_pItem) {
        
        # Run ConsensusClusterPlus with the current pItem_val
        consensus_results_app_one <- ConsensusClusterPlus(
          final_input_matrix, 
          maxK = working_max_k, 
          reps = CCP_iter,            
          pItem = pItem_val,            # Percentage of items sampled per iteration
          pFeature = 1,                # Use all features
          clusterAlg = clustering_alg,      
          innerLinkage = "ward.D2",
          finalLinkage = "ward.D2",
          distance = "euclidean", 
          seed = 2024L,
          plot ="none",           
          verbose = FALSE,
          writeTable = FALSE
        )
        
        # Calculate item consensus and cluster metrics (ICL)
        pdf(file = NULL)  # Start invisible plotting device
        icl_results <- calcICL(consensus_results_app_one, 
                               plot = "png", 
                               writeTable = TRUE)
        dev.off()          # Close device
        
        # Extract the itemConsensus dataframe from ICL results
        item_consensus_df <- icl_results[["itemConsensus"]]
        
        # Rename columns for clarity
        colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
        
        # Check if any values in the item_consensus column are NA
        if (any(is.na(item_consensus_df$item_consensus))) {
          # Print message indicating NA detected and that pItem will be increased
          message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
          
          # Increase pItem by the predefined step to sample more items per iteration
          pItem_val <- pItem_val + step
        } else {
          # If no NA found, set flag to FALSE to exit the loop
          has_na <- FALSE
        }
      }
    }
    
    # After loop finishes, check if NA values still exist (meaning max pItem was reached)
    if (has_na) {
      # Warn the user that maximum pItem was reached but NAs persist
      message("Reached max pItem but NA values still exist in item_consensus.")
      message("Reducing working_max_k to the maximum most stable K")
      
      # Find all K values with no NA item consensus
      ks_without_na <- item_consensus_df %>%
        group_by(k) %>%
        summarize(no_na = all(!is.na(item_consensus))) %>%
        filter(no_na) %>%
        pull(k)
      
      # if there is no solution skip it 
      if (length(ks_without_na) == 0) {
        
        message("No stable K values found. Skipping this data generation for this approach")
        
        # remove so data isn't produced  
        filter_approach <- filter_approach[filter_approach != "app_one"]
        break  # exits the  loop
        
        
        
      } else {
        # Get maximal most stable K  
        max_stable_k <- max(ks_without_na)
        message("maximum most stable K = ", max_stable_k)
        
      }
      
      # special condition for if maximum most stable K = 2 
      if (max_stable_k == 2){
        
        
        # let user know 
        message("since maximum most stable K = ", max_stable_k)
        message("Will have to extract k = 2 results using k = 3 analysis")
        
        # set that to working_k 
        working_max_k <- 2
        
        
        # Initialize starting value for pItem (percentage of items to sample per iteration)
        pItem_val <- 0.8
        
        # Set maximum allowable pItem value (cannot exceed 1.0)
        max_pItem <- 1.0
        
        # Step size for incrementing pItem when NA values are found
        step <- 0.05
        
        
        # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
        while (has_na && pItem_val <= max_pItem) {
          # Run ConsensusClusterPlus with the current pItem_val
          consensus_results_app_one <- ConsensusClusterPlus(
            final_input_matrix, 
            maxK = 3, 
            reps = CCP_iter,            
            pItem = pItem_val,            # Percentage of items sampled per iteration
            pFeature = 1,                # Use all features
            clusterAlg = clustering_alg,      
            innerLinkage = "ward.D2",
            finalLinkage = "ward.D2",
            distance = "euclidean", 
            seed = 2024L,
            plot ="none",           
            verbose = FALSE,
            writeTable = FALSE
          )
          
          # Only extract out result for K = 2
          consensus_results_app_one <- consensus_results_app_one[1:2]
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results_app_one, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results_app_one, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          # Check if any values in the item_consensus column are NA
          if (any(is.na(item_consensus_df$item_consensus))) {
            # Print message indicating NA detected and that pItem will be increased
            message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
            
            # Increase pItem by the predefined step to sample more items per iteration
            pItem_val <- pItem_val + step
          } else {
            # If no NA found, set flag to FALSE to exit the loop
            has_na <- FALSE
          }
        }
        
        
      } else {
        
        
        # set that to working_k 
        working_max_k <- max_stable_k
        
        # Initialize starting value for pItem (percentage of items to sample per iteration)
        pItem_val <- 0.8
        
        # Set maximum allowable pItem value (cannot exceed 1.0)
        max_pItem <- 1.0
        
        # Step size for incrementing pItem when NA values are found
        step <- 0.05
        
        
        # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
        while (has_na && pItem_val <= max_pItem) {
          # Run ConsensusClusterPlus with the current pItem_val
          consensus_results_app_one <- ConsensusClusterPlus(
            final_input_matrix, 
            maxK = working_max_k, 
            reps = CCP_iter,            
            pItem = pItem_val,            # Percentage of items sampled per iteration
            pFeature = 1,                # Use all features
            clusterAlg = clustering_alg,      
            innerLinkage = "ward.D2",
            finalLinkage = "ward.D2",
            distance = "euclidean", 
            seed = 2024L,
            plot ="none",           
            verbose = FALSE,
            writeTable = FALSE
          )
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results_app_one, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results_app_one, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          # Check if any values in the item_consensus column are NA
          if (any(is.na(item_consensus_df$item_consensus))) {
            # Print message indicating NA detected and that pItem will be increased
            message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
            
            # Increase pItem by the predefined step to sample more items per iteration
            pItem_val <- pItem_val + step
          } else {
            # If no NA found, set flag to FALSE to exit the loop
            has_na <- FALSE
          }
          
        }
        
      }
      
      
      
    }
    
    # If no NA found, proceed with downstream analysis steps
    if (!has_na) {
      # Classify stability based on item_consensus values using custom function
      item_consensus_df$stability <- sapply(item_consensus_df$item_consensus, get_stability_category)
      
      # Split the dataframe into a list of dataframes, one per cluster number k
      item_consensus_by_k <- split(item_consensus_df, item_consensus_df$k)
    }
    
    
    
    
    
    # Loop over each k from ConsensusClusterPlus results
    CCP_plot_list <- list() 
    
    
    for (k in 2:length(consensus_results_app_one)) {
      cm <- consensus_results_app_one[[k]]$consensusMatrix
      colnames(cm) <- colnames(final_input_matrix)
      rownames(cm) <- colnames(final_input_matrix)  # consensusMatrix is symmetric, matches samples
      
      cm <-  cm[consensus_results_app_one[[k]]$consensusTree$order,]
      
      
      clusters <- consensus_results_app_one[[k]]$consensusClass
      order_idx <- order(clusters)
      
      # Create annotation dataframe
      ann <- data.frame(Cluster = factor(clusters[order_idx]))
      #rownames(ann) <- names(clusters)
      
      # Step 1: Ensure Cluster is a factor with only the levels that exist
      ann$Cluster <- as.factor(ann$Cluster)
      ann$Cluster <- droplevels(ann$Cluster)  # drop any unused levels
      
      # Step 2: Define annotation colors correctly
      cluster_levels <- levels(ann$Cluster)
      
      # RColorBrewer requires at least 3 colors, so ensure enough levels
      num_clusters <- length(cluster_levels)
      palette <- RColorBrewer::brewer.pal(max(3, num_clusters), "Set2")[1:num_clusters]
      
      # Step 3: Set names to match factor levels exactly
      annotation_colors <- list(Cluster = setNames(palette, cluster_levels))
      
      
      # Save heatmap to list
      CCP_plot <- pheatmap(
        cm,
        cluster_rows = FALSE,
        cluster_cols = consensus_results_app_one[[k]]$consensusTree,
        annotation_col = ann,
        annotation_row = ann,
        show_colnames = FALSE,
        annotation_colors = annotation_colors,
        color = colorRampPalette(c("white", "blue"))(100),
        silent = TRUE,
        main = paste("Consensus Matrix for k =", k)
      )
      
      # Make name 
      CCP_plot_list_name <- paste0("CCP_k_", k)
      CCP_plot_list[[CCP_plot_list_name]] <- CCP_plot
    }
  
    # Create graphs for item consensus 
    
    item_consensus_overall <- as.data.frame(icl_results[["clusterConsensus"]])
    
    # If you want to round the labels nicely, e.g. 2 decimal places
    item_consensus_overall <- item_consensus_overall %>%
      mutate(label = round(clusterConsensus, 2))
    
    item_con_plot <- ggplot(item_consensus_overall,
                            aes(
                              x = factor(k),
                              y = clusterConsensus,
                              fill = factor(cluster)
                            )) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_text(aes(label = label), 
                position = position_dodge(width = 0.9),
                vjust = -0.3,  # slightly above the bar
                size = 6) +    # adjust text size
      labs(title = "Cluster Consensus by k", x = "k", y = "Consensus Score") +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        axis.text = element_text(colour = "black", size = 24),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.text.y = element_text(colour = "black", size = 24),
        axis.title = element_text(colour = "black", size = 30),
        axis.ticks = element_line(size = 1.5),
        panel.border = element_rect(
          colour = "black",
          fill = NA,
          size = 1
        )
      ) +
      scale_fill_manual(values = cols[1:working_max_k], name = "Cluster", drop = FALSE)
    
    
    
    item_consensus_by_k_plot_list <- list()
    for(k in names(item_consensus_by_k)){
      
      k_numeric <- as.numeric(k)
      item_con_dataset <- item_consensus_by_k[[as.character(k)]]
      
      
      item_consensus_by_k_plot <- ggplot(item_con_dataset,
                                         aes(
                                           x = sample_id,
                                           y = item_consensus,
                                           fill = factor(cluster)
                                         )) +
        geom_bar(stat = "identity") +
        labs(title = paste0("Item Consensus by Sample (k = ", k, ")"),
             x = "Sample ID",
             y = "Item Consensus",
             fill = "Cluster") +
        theme_classic() +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = 8
          ),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 14),
          legend.position = "right"
        ) +
        scale_fill_manual(values = cols[1:working_max_k], name = "Cluster", drop = FALSE)
      
      item_consensus_by_k_plot_list[[as.character(k)]] <- item_consensus_by_k_plot
      
    }
    
    
    # Specify custom results directory
    results_dir <- paste0(working_dir, "/CIM_states_results_", clustering_alg, "_", input_data_type, "_app_one")
    if(!dir.exists(results_dir)) {
      dir.create(results_dir, recursive = TRUE)
    }
    
    
    # Get Clustering gene set for this approach
    write.csv(clustering_gene_set,
              file = file.path(results_dir, paste0("clustering_gene_set.csv")),
              row.names = FALSE)
    
    
    # Save tables for each k into CSV files
    for (k_val in names(item_consensus_by_k)) {
      write.csv(item_consensus_by_k[[k_val]],
                file = file.path(results_dir, paste0("item_consensus_k", k_val, ".csv")),
                row.names = FALSE)
    }
    
    # Save Plots 
    
    ggsave(
      filename = paste0(results_dir, "/item_con_plot.png"),
      plot = item_con_plot,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    
    for (k in names(item_consensus_by_k_plot_list)){
      
      img <- item_consensus_by_k_plot_list[[as.character(k)]]
      
      ggsave(
        filename = paste0(results_dir, "/item_con_sample_k-", as.character(k), ".png"),
        plot = img,
        width = 10,
        height = 6,
        dpi = 300
      )
      
      
    }
    
    
    for (i in seq_along(CCP_plot_list)) {
      plot_name <- names(CCP_plot_list)[i]
      png_filename <- file.path(results_dir, paste0(plot_name, ".png"))
      png(png_filename, width = 1200, height = 1000, res = 150)
      grid.newpage()
      obj <- CCP_plot_list[[i]]
      if (inherits(obj, "grob")) {
        grid.draw(obj)
      } else if (inherits(obj, "Heatmap")) {
        draw(obj)
      } else if (!is.null(obj$gtable)) {
        grid.draw(obj$gtable)   # <-- fallback: draw the gtable element if it exists
      } else {
        warning(sprintf("Element %s is not drawable by grid.draw, draw, or print.", plot_name))
      }
      dev.off()
    }
    
    
    
    
    
    # Create cdf_plots 
    k_vec = 2:working_max_k
    cdf_data <- data.frame()
    PACs <- numeric(length(k_vec))
    
    for(i in seq_along(k_vec)) {
      k <- k_vec[i]
      consensus_matrix <- consensus_results_app_one[[k]][["consensusMatrix"]]
      cons_vals <- consensus_matrix[lower.tri(consensus_matrix)]
      Fn <- ecdf(cons_vals)  # standard consensus cluster PAC
      PACs[i] <- Fn(0.9) - Fn(0.1)
      h <- hist(cons_vals, breaks = 200, plot = FALSE)
      cdf_y <- cumsum(h$counts) / sum(h$counts)
      cdf_data <- bind_rows(
        cdf_data,
        data.frame(
          k = as.factor(k),
          x = h$mids,
          y = cdf_y
        )
      )
    }
    
    # Now join the correct PAC values
    pac_df <- data.frame(k = as.factor(k_vec), PAC = PACs)
    cdf_data <- left_join(cdf_data, pac_df, by = "k")
    cdf_data$line_label <- paste0("K = ", cdf_data$k, " (PAC=", formatC(cdf_data$PAC, format = "f", digits = 3), ")")
    
    
    cdf_plot <- ggplot(cdf_data, aes(x = x, y = y, color = line_label, group = k)) +
      geom_line(size = 1.5) +
      labs(
        x = "Consensus Index",
        y = "Cumulative Fraction",
        title = "Consensus CDFs across Ks",
        color = "K (PAC)"
      ) +
      geom_vline(xintercept = c(0.1, 0.9), linetype = "dashed", color = "gray") +
      scale_x_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
      scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.98, 0.02), # x, y from 0 (left/bottom) to 1 (right/top)
        legend.justification = c("right", "bottom"),
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill="white", color="black", size=0.8),
        axis.text = element_text(colour = "black", size = 24),
        axis.text.x = element_text(colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title = element_text(colour = "black", size = 30),
        axis.ticks = element_line(size = 1.5),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)
      )
    
    
    # save cdf plot
    ggsave(
      filename = paste0(results_dir, "/final_cdf_plot.png"),
      plot = cdf_plot,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    
    # Calculate PAC metric (Proportion of Ambiguous Clustering) using ecdf
    Kvec = 2:working_max_k
    x1 = 0.1; x2 = 0.9 
    PAC = rep(NA, length(Kvec)) 
    names(PAC) = paste("K=", Kvec, sep="") 
    for(i in Kvec){
      M = consensus_results_app_one[[i]]$consensusMatrix
      Fn = ecdf(M[lower.tri(M)])
      PAC[i-1] = Fn(x2) - Fn(x1)
    }
    optK = Kvec[which.min(PAC)]
    
    # Convert PAC to data frame for plotting
    pac_values <- data.frame(k = Kvec, pac = unname(PAC))
    
    # Plot PAC values (lower is better) and display in R environment (Zaoqu-Liu/IRLS GitHub)
    pac_plot <- ggplot(pac_values, aes(factor(k), pac, group=1)) +
      geom_line() +
      theme_bw(base_rect_size = 1.5) +
      geom_point(size=4, shape=21, color='darkred', fill='orange') +
      #ggtitle('Proportion of ambiguous clustering') +
      xlab('Number of clusters K') + ylab("PAC Score") +
      theme_classic() + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top",
            axis.text = element_text(colour = "black", size = (24)),
            axis.text.x = element_text(colour = "black", size = (24)),
            axis.text.y = element_text(colour = "black", size = (24)),
            axis.title = element_text(colour = "black", size = (30)),
            axis.ticks = element_line(size = 1.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    
    
    
    
    # Create item consensus table
    avg_cluster_consensus <- item_consensus_overall %>%
      group_by(k) %>%
      summarize(avg_cluster_consensus = mean(clusterConsensus))
    
    # Silhouette calculations using consensus cluster final
    
    orig_space_input <- t(final_input_matrix)  # transpose to get dims as samples x features
    rownames(orig_space_input) <- rownames(clustering_matrix)
    
    # Distances
    dist_space_for_clustering <- dist(orig_space_input)
    dist_ge_space <- dist(clustering_matrix[,clustering_gene_set, drop=FALSE])
    
    # Collect silhouettes for all k
    sil_results <- lapply(2:length(consensus_results_app_one), function(k) {
      
      clusters <- consensus_results_app_one[[k]]$consensusClass
      
      # Silhouette in gene space
      sil_space_for_clustering   <- silhouette(clusters, dist_space_for_clustering)
      avg_space_for_clustering   <- mean(sil_space_for_clustering[, 3])
      
      # Silhouette in PaCMAP space
      sil_ge_space <- silhouette(clusters, dist_ge_space)
      avg_ge_space <- mean(sil_ge_space[, 3])
      
      data.frame(
        k = k,
        avg_sil_space_for_clustering = avg_space_for_clustering,
        avg_sil_ge_space = avg_ge_space
      )
    })
    
    # Combine
    sil_df <- do.call(rbind, sil_results)
    sil_df$avg_sil_combined <- ((sil_df$avg_sil_space_for_clustering +
        sil_df$avg_sil_ge_space) /
        2)

    
    # Plot comparison
    sil_plot_both <- ggplot(sil_df, aes(x = k)) +
      geom_line(aes(y = avg_sil_space_for_clustering, color = paste0(input_data_type))) +
      geom_point(aes(y = avg_sil_space_for_clustering, color = paste0(input_data_type))) +
      geom_line(aes(y = avg_sil_ge_space, color = "Gene space")) +
      geom_point(aes(y = avg_sil_ge_space, color = "Gene space")) +
      labs(y = "Average silhouette width", color = "Space", title = "Silhouette comparison across k") +
      theme_classic(base_size = 16) + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top",
            axis.text = element_text(colour = "black", size = 24),
            axis.title = element_text(colour = "black", size = 30),
            axis.ticks = element_line(size = 1.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    
    

    sil_plot_original_space <- ggplot(sil_df, aes(x = k)) +
      geom_line(aes(y = avg_sil_space_for_clustering, color = paste0(input_data_type))) +
      geom_point(aes(y = avg_sil_space_for_clustering, color = paste0(input_data_type))) +
      labs(y = "Average silhouette width", color = "Space", title = "Silhouette comparison across k") +
      theme_classic(base_size = 16) + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top",
            axis.text = element_text(colour = "black", size = 24),
            axis.title = element_text(colour = "black", size = 30),
            axis.ticks = element_line(size = 1.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    

    sil_plot_ge_space <- ggplot(sil_df, aes(x = k)) +
      geom_line(aes(y = avg_sil_ge_space, color = "Gene space")) +
      geom_point(aes(y = avg_sil_ge_space, color = "Gene space")) +
      labs(y = "Average silhouette width", color = "Space", title = "Silhouette comparison across k") +
      theme_classic(base_size = 16) + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top",
            axis.text = element_text(colour = "black", size = 24),
            axis.title = element_text(colour = "black", size = 30),
            axis.ticks = element_line(size = 1.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    

  
    
    # save silhouette plot
    ggsave(
      filename = paste0(results_dir, "/sil_plot.png"),
      plot = sil_plot_both,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    
    # save silhouette plot
    ggsave(
      filename = paste0(results_dir, "/sil_plot_original_space.png"),
      plot = sil_plot_original_space,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    # save silhouette plot
    ggsave(
      filename = paste0(results_dir, "/sil_plot_ge_space.png"),
      plot = sil_plot_ge_space,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    
    # save pac plot
    ggsave(
      filename = paste0(results_dir, "/final_pac_plot.png"),
      plot = pac_plot,
      width = 10,
      height = 6,
      dpi = 300
    )
    

   

    
       
    # Get suggested optimal k values
    silhouette_dim_space_k <- sil_df$k[which.max(sil_df$avg_sil_space_for_clustering)]
    silhouette_ge_space_k <- sil_df$k[which.max(sil_df$avg_sil_ge_space)]
    silhouette_combined_avg_k <- sil_df$k[which.max(sil_df$avg_sil_combined)]
    pac_k <- Kvec[which.min(PAC)]  # Using optK from PAC calculation
    avg_cluster_consensus_k <- avg_cluster_consensus[which.max(avg_cluster_consensus$avg_cluster_consensus), ]$k
    
    # Print metric values
    cat("\n===== OPTIMAL CLUSTER SUGGESTIONS =====\n")
    cat("Silhouette analysis in dim space suggests k =", silhouette_dim_space_k, "\n")
    cat("Silhouette analysis in GE space suggests k =", silhouette_ge_space_k, "\n")
    cat("Silhouette analysis of combined avg (GE and red. dim) space suggests k =", silhouette_combined_avg_k, "\n")
    cat("PAC method suggests k =", pac_k, "\n")
    cat("cluster consensus method suggests k =", avg_cluster_consensus_k, "\n")
    
    
    
 # Create metrics summary dataframe
      metrics_summary <- data.frame(
        k = 2:working_max_k,
        silhouette_dim_reduce_space = NA, # renamed to indicate dim‑reduce (original) space
        silhouette_ge_space = NA, # silhouette computed in PacMAP (ge) space
        silhouette_combined_avg = NA, # silhouette computed in PacMAP (ge) space
        pac = NA,
        avg_cluster_consensus = NA
      )

      # Fill in the metrics
      for (i in 1:nrow(metrics_summary)) {
        k_val <- metrics_summary$k[i]

        # Silhouette in the original (dim‑reduce) space
        idx <- which(sil_df$k == k_val)
        if (length(idx) > 0) {
          metrics_summary$silhouette_dim_reduce_space[
            i
          ] <- sil_df$avg_sil_space_for_clustering[idx]

          metrics_summary$silhouette_ge_space[i] <- sil_df$avg_sil_ge_space[idx]
          metrics_summary$silhouette_combined_avg[i] <- sil_df$avg_sil_combined[
            idx
          ]
        }

        # PAC
        k_idx <- which(Kvec == k_val)
        if (length(k_idx) > 0) {
          metrics_summary$pac[i] <- PAC[k_idx]
        }

        # avg_cluster_consensus
        cc_idx <- which(avg_cluster_consensus$k == k_val)
        if (length(cc_idx) > 0) {
          metrics_summary$avg_cluster_consensus[
            i
          ] <- avg_cluster_consensus$avg_cluster_consensus[cc_idx]
        }
      }

      # Alternatively,  Optimal K can be selected using rank based compromise
      # The rank-based compromise approach assigns a rank to each k for both silhouette (descending) and PAC (ascending), then sums the ranks to identify the k with the best overall trade-off between high silhouette and low PAC. The k with the lowest combined rank is selected as the optimal choice.
      # https://dl.acm.org/doi/10.1145/371920.372165

      # Assign ranks (lowest rank = best)
      metrics_summary$rank_silhouette_dim_reduce_space <- rank(
        -metrics_summary$silhouette_dim_reduce_space,
        ties.method = "min"
      ) # higher silhouette = rank 1 (original space)
      metrics_summary$rank_silhouette_ge_space <- rank(
        -metrics_summary$silhouette_ge_space,
        ties.method = "min"
      ) # higher silhouette = rank 1 (ge space)

      metrics_summary$rank_silhouette_combined_avg <- rank(
        -metrics_summary$silhouette_combined_avg,
        ties.method = "min"
      ) # higher silhouette = rank 1 (averaged space)
      metrics_summary$rank_pac <- rank(metrics_summary$pac, ties.method = "min") # lower PAC = rank 1
      metrics_summary$rank_item_cluster_consensus <- rank(
        -metrics_summary$avg_cluster_consensus,
        ties.method = "min"
      ) # higher consensus = rank 1

      # Build the column names dynamically.
      rank_columns <- c(
        paste0("rank_", clustering_metrics) # e.g. "rank_pac", "rank_item_cluster_consensus"
      )

      # Compute the overall rank as the row sum across selected rank columns
      metrics_summary$overall_rank <- rowSums(
        metrics_summary[, rank_columns, drop = FALSE],
        na.rm = TRUE
      )

      # Find all rows with the minimum overall rank
      best_rows <- metrics_summary[
        metrics_summary$overall_rank == min(metrics_summary$overall_rank),
      ]

      # Use the last metric in clustering_metrics as the tiebreaker
      tie_metric <- clustering_metrics[length(clustering_metrics)]

      if (verbose) {
        cat(
          "Metrics being used for ranking:\n  ",
          paste(clustering_metrics, collapse = ", "),
          "\n"
        )
        cat(
          "Tie‑breaker metric (last in clustering_metrics vector):",
          tie_metric,
          "\n"
        )
      }

      # Determine if the tiebreaking metric should be maximized or minimized
      # (Assuming: max for 'silhouette' and 'item_cluster_consensus', min for 'pac')
      if (tie_metric == "silhouette_dim_reduce_space") {
        idx <- which.max(best_rows$silhouette_dim_reduce_space)
      } else if (tie_metric == "silhouette_ge") {
        idx <- which.max(best_rows$silhouette_ge_space)
      } else if (tie_metric == "item_cluster_consensus") {
        idx <- which.max(best_rows$avg_cluster_consensus)
      } else if (tie_metric == "pac") {
        idx <- which.min(best_rows$pac)
      } else {
        stop("Unknown tiebreaker metric.")
      }

      best_k_row <- best_rows[idx, ]
      best_k <- best_k_row$k


# add cluster metrics used for evaluation and tie
metrics_summary$metrics_used <- paste(clustering_metrics, collapse = ", ")
metrics_summary$metrics_used_for_tie <- tie_metric
    
    
    # Print metrics summary
    cat("\n===== CLUSTER METRICS SUMMARY =====\n")
    metrics_summary
    cat("\n")
    
    # Write metrics to file
    write.csv(metrics_summary, file.path(results_dir, "final_cluster_metrics_summary.csv"), row.names = FALSE)
    
    # Based on examining the plots and metrics, set your optimal k manually here:
    # CHANGE THIS VALUE after reviewing the plots:
    optimal_k <- best_k  # <-- MODIFY THIS BASED ON YOUR EXAMINATION IF NEEDED 
    
    cat("Using k =", optimal_k, "as the optimal number of clusters\n")
    
    # Extract cluster assignments for optimal k
    consensus_clusters <- consensus_results_app_one[[optimal_k]]$consensusClass
    
    
    ####### FINAL VISUALIZATION #######
    
    

    
    
    # Get stability/reliability metrics 
    optimal_stability_cluster_metrics <- item_consensus_by_k[[as.character(optimal_k)]]
    
    # Start from your original data frame
    wide_optimal_stability_cluster_metrics <- optimal_stability_cluster_metrics %>%
      # Only keep relevant columns
      dplyr::select(sample_id, cluster, item_consensus) %>%
      # Pivot wider: each cluster makes a new column
      pivot_wider(
        names_from = cluster,
        values_from = item_consensus,
        names_prefix = "item_consensus_cluster"
      )
    
    
    # Write metrics to file
    write.csv(optimal_stability_cluster_metrics,
              file.path(results_dir, "final_cluster_stability_metrics.csv"), row.names = FALSE)
    
    write.csv(wide_optimal_stability_cluster_metrics,
              file.path(results_dir, "final_wide_cluster_stability_metrics.csv"), row.names = FALSE)
    
    
    # Create data frame for UMAP visualization
    final_df_plot_app_one <- as.data.frame(umap_matrix_visual$layout) %>%
      rownames_to_column("sample_id") %>% 
      dplyr::rename(UMAP1 = V1, UMAP2 = V2) %>%
      dplyr::inner_join(as.data.frame(pacmap_delta_embeddings_visual) %>% 
                          rownames_to_column("sample_id") %>% 
                          dplyr::rename(PACMAP1 = V1, PACMAP2 = V2)) %>%
      dplyr::inner_join(as.data.frame(pca_scores) %>% 
                          rownames_to_column("sample_id")) %>%
      dplyr::mutate(base_id = gsub("^delta_([^_]+)_.*", "\\1", sample_id)) %>%
      dplyr::mutate(cluster_assignments = consensus_clusters) %>%
      dplyr::inner_join(as.data.frame(clustering_matrix) %>%
                          rownames_to_column("sample_id"),
                        by = "sample_id")
    
    
    # Your input matrix: samples x genes
    pca_fviz_input <- final_df_plot_app_one[, c("PC1", "PC2", "sample_id")] %>%
      tibble::column_to_rownames("sample_id")
    
    
    fviz_cluster_plot <- fviz_cluster(
      object = list(data = pca_fviz_input, cluster = as.integer(consensus_clusters)),
      palette = cols[1:optimal_k],
      geom = "text",
      ellipse.type = "convex",
      ggtheme = theme_bw(),
      ggtitle("PCA plot on Clustered Data")
    )
    
    fviz_cluster_plot
    ggsave(file.path(results_dir, paste0("final_fviz_2D_PCA_k-", optimal_k, ".png")), 
           fviz_cluster_plot, width = 10, height = 8, dpi = 300)
    
    
    # 2. PCA projection using selected PCs
    pca_scores_3D <- pca_result$x[, 1:3, drop = FALSE]
    
    # 3. Prepare input for downstream analysis, e.g., clustering
    # Here, samples x PCs (rows = samples, cols = PCs)
    pca_fviz_input_3D <- as.matrix(pca_scores_3D)
    rownames(pca_fviz_input_3D) <- final_df_plot_app_one$sample_id  # optional, set sample IDs
    
    # If pca_fviz_input_3D is samples x PCs, select first 3 PCs
    
    plot_df_3D_PCA <- data.frame(
      PC1 = pca_fviz_input_3D[, 1],
      PC2 = pca_fviz_input_3D[, 2],
      PC3 = pca_fviz_input_3D[, 3],
      cluster = as.factor(as.integer(consensus_clusters)),
      sample_id = rownames(pca_fviz_input_3D)
    )
    
    
    
    # 3D scatter plot
    fig_PCA_3D <- plot_ly(
      data = plot_df_3D_PCA,
      x = ~PC1, y = ~PC2, z = ~PC3,
      color = ~cluster,
      text = ~sample_id,
      type = "scatter3d",
      mode = "markers",
      colors = cols
    ) %>%
      plotly::layout(
        title = "3D PCA Cluster Visualization",
        scene = list(
          xaxis = list(title = "PC1"),
          yaxis = list(title = "PC2"),
          zaxis = list(title = "PC3")
        )
      )
    
    fig_PCA_3D
    # Load library if not already loaded
    library(htmlwidgets)
    
    # Save as an interactive HTML file
    saveWidget(
      widget = fig_PCA_3D,
      file = paste0(results_dir, "/final_3D_PCA_k-", optimal_k, ".html"),
      selfcontained = TRUE
    )
    
    
    # Use letters as shapes 
    shapes <- unlist(lapply(letters, utf8ToInt))
    
    # Plot using GGPlot Color by cluster
    final_cluster_plot_umap <- ggplot(
      final_df_plot_app_one,
      aes(
        x = UMAP1,
        y = UMAP2,
        label = base_id ,
        color = as.factor(cluster_assignments) 
      )
    ) +
      scale_color_manual(values = cols, name = "Cluster") +
      ggforce::geom_mark_hull(
        aes(group = as.factor(cluster_assignments),
            fill = as.factor(cluster_assignments), label = NULL),
        concavity = 2, expand = unit(2, "mm"), alpha = 0.2
      ) +
      geom_label_repel(size = 4, fontface = "bold", label.size = 0.3, fill = "white", 
                       alpha = 1, box.padding = 0.25, max.overlaps = 100) +
      labs(x = "UMAP 1", y = "UMAP 2",
           fill = "Cluster",
           color = "Cluster",
           title = paste0(
             "Final Clustering Using ", input_data_type, " input and ",
             clustering_alg, " CCP Alg (k = ", optimal_k, ")"
           )) +
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
    
    
    ggsave(file.path(results_dir, paste0("final_clustering_shown_on_umap_k-", optimal_k, ".png")), 
           final_cluster_plot_umap, width = 10, height = 8, dpi = 300)
    
    
    
    # Plot using GGPlot Color by cluster
    final_cluster_plot_pacmap <- ggplot(
      final_df_plot_app_one,
      aes(
        x = PACMAP1,
        y = PACMAP2,
        label = base_id,  # Your patient/sample label column
        color = as.factor(cluster_assignments)
      )
    ) +
      scale_color_manual(values = cols, name = "Cluster") +
      ggforce::geom_mark_hull(
        aes(group = as.factor(cluster_assignments),
            fill = as.factor(cluster_assignments), label = NULL),
        concavity = 2, expand = unit(2, "mm"), alpha = 0.2
      ) +
      geom_label_repel(size = 4, fontface = "bold", label.size = 0.3, fill = "white", 
                       alpha = 1, box.padding = 0.25, max.overlaps = 100) +
      labs(x = "PacMAP 1", y = "PacMAP 2",
           fill = "Cluster",
           color = "Cluster",
           title = paste0(
             "Final Clustering Using ", input_data_type, " input and ",
             clustering_alg, " CCP Alg (k = ", optimal_k, ")"
           )) +
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
    
    
    
    ggsave(file.path(results_dir, paste0("final_clustering_shown_on_pacmap_k-", optimal_k, ".png")), 
           final_cluster_plot_pacmap, width = 10, height = 8, dpi = 300)
    
    
  }
  
  }
  
  if (no_app_one == TRUE && 
      "app_one" %in% filter_approach) {
    cat(
      "\033[1;31mWARNING: Final Gene Set for Approach one empty! Consider Relaxing P-val
        Approach output will not be produced\033[0m\n"
    )
    
    # set string
    no_app_one_string <- "Approach one final gene set empty; Consider Relaxing P-val"
    
    
    # update iter_log
    iter_log <- rbind(
      iter_log,
      data.frame(
        gene_set = "app_one final get set", iteration = NA,
        before = NA,
        after = no_app_one_string
      )
    )
    
    
    # set df to null
    final_df_plot_app_one <- NULL
  }
  
  
  # Extra capturing of output 
  {
  # Write initial matrix to file
  write.csv(clustering_matrix,
            file.path(results_dir, "initial_clustering_mat.csv"), row.names = TRUE)
  
  # Write iter_log  to file
  write.csv(iter_log,
            file.path(results_dir, "iter_log_all_approaches.csv"), row.names = FALSE)
  
  # Convert gene_sets to to data frame
  df_gene_sets_for_file <- data.frame(
    list_name = names(all_gene_sets),
    items = sapply(all_gene_sets, function(x) paste(x, collapse = ","))
  )
  
  # Write gene_sets to CSV 
  write.csv(df_gene_sets_for_file,  
            file.path(results_dir, "all_gene_sets.csv"), row.names = FALSE)
  
  # write final_clustered df to file
  write.csv(final_df_plot_app_one,  
            file.path(results_dir, "clustered_samples_app_one.csv"), row.names = FALSE)
  
  # write pacmap settings 
  write.csv(pm_settings,  
            file.path(results_dir, "pacmap_settings.csv"), row.names = TRUE)
  
         # write removed genes setting 
    write.csv(removed_genes,  
              file.path(results_dir, "removed_genes.csv"), row.names = TRUE)
  }
  
      })
  }
  
  
  # ---- Approach Two (Iterated over all unique genes ) ---- 
  
  
  if ("app_two" %in% filter_approach){
  
    
    # Set seed appropriately
    withr::with_seed(seed = 2024L, {
      
    # set back to max_k
  working_max_k <- max_k
  
  if(no_app_two == FALSE){
  
  {   
    # CHOOSE INPUT DATA TYPE FOR CLUSTERING
    # Set to "gene_set" for original gene expression data, or "umap" for UMAP coordinates, or "pca"
    input_data_type <- "pacmap"  # <-- MODIFY HERE: use "genes" or "umap" or "pca" or "pacmap"

    
    # Make matrix for delta
    #toggle here
    clustering_gene_set <- iterated_over_all_genes
    clustering_gene_set <- clustering_gene_set[clustering_gene_set %in% colnames(clustering_matrix)]
   
    
    # Run UMAP on the selected gene set (same as in the original code)
    # will always be 2D 
    umap_matrix_visual <- umap::umap(clustering_matrix[,clustering_gene_set, drop=FALSE],
                                                          config=custom.config.visual)
    
    
    # Build embedder for analysis 
    script_pacmap_reducer <- do.call(pacmap$PaCMAP, pm_defaults)
    
 
    # Run PaCMAP
    pacmap_delta_embeddings <- script_pacmap_reducer$fit_transform(
      clustering_matrix[, clustering_gene_set, drop = FALSE]
    )
    
    # set rownames 
    rownames(pacmap_delta_embeddings) <- rownames(clustering_matrix)
    
    
    # 2D visualizations 
    # create reducer (also resets seed)
    pm_visualization_reducer <- do.call(pacmap$PaCMAP, pm_visualization)
    
    # Run PACMAP on selected gene set
    # Run PaCMAP
    pacmap_delta_embeddings_visual <- pm_visualization_reducer$fit_transform(
      clustering_matrix[, clustering_gene_set, drop = FALSE]
    )
  
    # set rownames 
    rownames(pacmap_delta_embeddings_visual) <- rownames(clustering_matrix)
    
    
    # Your input matrix: samples x genes
    pca_matrix <- clustering_matrix[, clustering_gene_set]
    
    pca_result <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)
    pca_scores <- pca_result$x[, 1:3, drop = FALSE]
    
    
    
    if (verbose) message(sprintf(
      "\n\n=== Creating Output for Approach Two using %s algorithm and %s as input ===",
      clustering_alg,
      input_data_type
    ))
    
    
    
    
    
    # Define color palette for clusters
    cols <- hue_pal()(working_max_k)
    
    ########################################################################################
    
    
    # Create input matrix for consensus clustering based on chosen data type
    if (input_data_type == "genes") {
      # Use gene data
      final_input_matrix <- t(clustering_matrix[, clustering_gene_set])  # transpose to get genes as rows
      cat("Using gene expression data for clustering.\n")
    } else if (input_data_type == "umap") {
      # Use UMAP coordinates - specifically from umap_matrix to match your original approach
      final_input_matrix <- t(umap_matrix$layout)  # transpose to get dimensions as rows
      cat("Using UMAP coordinates for clustering.\n")
    } else if (input_data_type == "pca"){
      
      # Your input matrix: samples x genes
      final_input_matrix <- clustering_matrix[, clustering_gene_set]
      
      # Make sure your final_input_matrix is samples x features (rows=samples, cols=genes)
      # You can check with dim(final_input_matrix)
      
      # 1. Run Parallel Analysis to determine number of PCs
      # PA expects samples in rows, variables in columns
      centile_seq <- seq(95, 80, by = -5)
      max_pcs <- 3   # Default value in case paranormal fails
      pa_result <- NULL
      for (cent in centile_seq) {
        pa_result <- try(
          suppressMessages(suppressWarnings(
            paran(final_input_matrix,
                  iterations = 1000,
                  centile = cent,
                  graph = FALSE,
                  quietly = TRUE)
          )), silent = TRUE)
        
        if (!inherits(pa_result, "try-error")) {
          message(sprintf("paran succeeded with centile = %d", cent))
          max_pcs <- pa_result$Retained  # number of components retained by PA
          break
        }
      }
      if (inherits(pa_result, "try-error")) {
        warning("paran failed at all tested centile values (95 to 80). Using default max_pcs = 3 for this gene set.")
        # continue; max_pcs remains 3
      }
      
      cat("Number of PCs retained by Parallel Analysis:", max_pcs, "\n")
      
      
      # 2. PCA projection using selected PCs
      pca_result <- prcomp(final_input_matrix, center = TRUE, scale. = TRUE)
      pca_scores <- pca_result$x[, 1:max_pcs, drop = FALSE]
      
      # 3. Prepare input for downstream analysis, e.g., clustering
      # Here, PCs x Samples (rows = PCs, cols = samples)
      final_input_matrix <- as.matrix(t(pca_scores))
      colnames(final_input_matrix) <- umap_plot_df$sample_id  # optional, set sample IDs
      
    } else if (input_data_type == "pacmap"){
      
      final_input_matrix <- t(pacmap_delta_embeddings)
      colnames(final_input_matrix) <- rownames(clustering_matrix)
      cat("Using PacMAP coordinates for clustering.\n")
    } else {
      stop("Invalid input_data_type. Use 'gene_set' or 'umap' or 'pca' or 'pacmap'.")
    }
    
    
    # Initialize starting value for pItem (percentage of items to sample per iteration)
    pItem_val <- 0.8
    
    # Set maximum allowable pItem value (cannot exceed 1.0)
    max_pItem <- 1.0
    
    # Step size for incrementing pItem when NA values are found
    step <- 0.05
    
    # Flag to track presence of NA values in itemConsensus
    has_na <- TRUE
    
    # if in the case of working max k = 2 
    if (working_max_k == 2){
      
      # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
      while (has_na && pItem_val <= max_pItem) {
        # Run ConsensusClusterPlus with the current pItem_val
        consensus_results_app_two <- ConsensusClusterPlus(
          final_input_matrix, 
          maxK = 3, 
          reps = CCP_iter,            
          pItem = pItem_val,            # Percentage of items sampled per iteration
          pFeature = 1,                # Use all features
          clusterAlg = clustering_alg,      
          innerLinkage = "ward.D2",
          finalLinkage = "ward.D2",
          distance = "euclidean", 
          seed = 2024L,
          plot ="none",           
          verbose = FALSE,
          writeTable = FALSE
        )
        
        # Only extract out result for K = 2
        consensus_results_app_two <- consensus_results_app_two[1:2]
        
        # Calculate item consensus and cluster metrics (ICL)
        pdf(file = NULL)  # Start invisible plotting device
        icl_results <- calcICL(consensus_results_app_two, 
                               plot = "png", 
                               writeTable = TRUE)
        dev.off()          # Close device
        
        # Extract the itemConsensus dataframe from ICL results
        item_consensus_df <- icl_results[["itemConsensus"]]
        
        # Rename columns for clarity
        colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
        
        
        # Calculate item consensus and cluster metrics (ICL)
        pdf(file = NULL)  # Start invisible plotting device
        icl_results <- calcICL(consensus_results_app_two, 
                               plot = "png", 
                               writeTable = TRUE)
        dev.off()          # Close device
        
        # Extract the itemConsensus dataframe from ICL results
        item_consensus_df <- icl_results[["itemConsensus"]]
        
        # Rename columns for clarity
        colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
        
        # Check if any values in the item_consensus column are NA
        if (any(is.na(item_consensus_df$item_consensus))) {
          # Print message indicating NA detected and that pItem will be increased
          message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
          
          # Increase pItem by the predefined step to sample more items per iteration
          pItem_val <- pItem_val + step
        } else {
          # If no NA found, set flag to FALSE to exit the loop
          has_na <- FALSE
        }
        
      } 
      
    }
    else {
      
      # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
      while (has_na && pItem_val <= max_pItem) {
        
        # Run ConsensusClusterPlus with the current pItem_val
        consensus_results_app_two <- ConsensusClusterPlus(
          final_input_matrix, 
          maxK = working_max_k, 
          reps = CCP_iter,            
          pItem = pItem_val,            # Percentage of items sampled per iteration
          pFeature = 1,                # Use all features
          clusterAlg = clustering_alg,      
          innerLinkage = "ward.D2",
          finalLinkage = "ward.D2",
          distance = "euclidean", 
          seed = 2024L,
          plot ="none",           
          verbose = FALSE,
          writeTable = FALSE
        )
        
        # Calculate item consensus and cluster metrics (ICL)
        pdf(file = NULL)  # Start invisible plotting device
        icl_results <- calcICL(consensus_results_app_two, 
                               plot = "png", 
                               writeTable = TRUE)
        dev.off()          # Close device
        
        # Extract the itemConsensus dataframe from ICL results
        item_consensus_df <- icl_results[["itemConsensus"]]
        
        # Rename columns for clarity
        colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
        
        # Check if any values in the item_consensus column are NA
        if (any(is.na(item_consensus_df$item_consensus))) {
          # Print message indicating NA detected and that pItem will be increased
          message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
          
          # Increase pItem by the predefined step to sample more items per iteration
          pItem_val <- pItem_val + step
        } else {
          # If no NA found, set flag to FALSE to exit the loop
          has_na <- FALSE
        }
      }
    }
    
    # After loop finishes, check if NA values still exist (meaning max pItem was reached)
    if (has_na) {
      # Warn the user that maximum pItem was reached but NAs persist
      message("Reached max pItem but NA values still exist in item_consensus.")
      message("Reducing working_max_k to the maximum most stable K")
      
      # Find all K values with no NA item consensus
      ks_without_na <- item_consensus_df %>%
        group_by(k) %>%
        summarize(no_na = all(!is.na(item_consensus))) %>%
        filter(no_na) %>%
        pull(k)
      
      # if there is no solution skip it 
      if (length(ks_without_na) == 0) {
        
        message("No stable K values found. Skipping this data generation for this approach")
        
        # remove so data isn't produced  
        filter_approach <- filter_approach[filter_approach != "app_two"]
        break  # exits the  loop
        
        
        
      } else {
        # Get maximal most stable K  
        max_stable_k <- max(ks_without_na)
        message("maximum most stable K = ", max_stable_k)
        
      }
      
      # special condition for if maximum most stable K = 2 
      if (max_stable_k == 2){
        
        
        # let user know 
        message("since maximum most stable K = ", max_stable_k)
        message("Will have to extract k = 2 results using k = 3 analysis")
        
        # set that to working_k 
        working_max_k <- 2
        
        
        # Initialize starting value for pItem (percentage of items to sample per iteration)
        pItem_val <- 0.8
        
        # Set maximum allowable pItem value (cannot exceed 1.0)
        max_pItem <- 1.0
        
        # Step size for incrementing pItem when NA values are found
        step <- 0.05
        
        
        # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
        while (has_na && pItem_val <= max_pItem) {
          # Run ConsensusClusterPlus with the current pItem_val
          consensus_results_app_two <- ConsensusClusterPlus(
            final_input_matrix, 
            maxK = 3, 
            reps = CCP_iter,            
            pItem = pItem_val,            # Percentage of items sampled per iteration
            pFeature = 1,                # Use all features
            clusterAlg = clustering_alg,      
            innerLinkage = "ward.D2",
            finalLinkage = "ward.D2",
            distance = "euclidean", 
            seed = 2024L,
            plot ="none",           
            verbose = FALSE,
            writeTable = FALSE
          )
          
          # Only extract out result for K = 2
          consensus_results_app_two <- consensus_results_app_two[1:2]
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results_app_two, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results_app_two, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          # Check if any values in the item_consensus column are NA
          if (any(is.na(item_consensus_df$item_consensus))) {
            # Print message indicating NA detected and that pItem will be increased
            message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
            
            # Increase pItem by the predefined step to sample more items per iteration
            pItem_val <- pItem_val + step
          } else {
            # If no NA found, set flag to FALSE to exit the loop
            has_na <- FALSE
          }
        }
        
        
      } else {
        
        
        # set that to working_k 
        working_max_k <- max_stable_k
        
        # Initialize starting value for pItem (percentage of items to sample per iteration)
        pItem_val <- 0.8
        
        # Set maximum allowable pItem value (cannot exceed 1.0)
        max_pItem <- 1.0
        
        # Step size for incrementing pItem when NA values are found
        step <- 0.05
        
        
        # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
        while (has_na && pItem_val <= max_pItem) {
          # Run ConsensusClusterPlus with the current pItem_val
          consensus_results_app_two <- ConsensusClusterPlus(
            final_input_matrix, 
            maxK = working_max_k, 
            reps = CCP_iter,            
            pItem = pItem_val,            # Percentage of items sampled per iteration
            pFeature = 1,                # Use all features
            clusterAlg = clustering_alg,      
            innerLinkage = "ward.D2",
            finalLinkage = "ward.D2",
            distance = "euclidean", 
            seed = 2024L,
            plot ="none",           
            verbose = FALSE,
            writeTable = FALSE
          )
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results_app_two, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results_app_two, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          # Check if any values in the item_consensus column are NA
          if (any(is.na(item_consensus_df$item_consensus))) {
            # Print message indicating NA detected and that pItem will be increased
            message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
            
            # Increase pItem by the predefined step to sample more items per iteration
            pItem_val <- pItem_val + step
          } else {
            # If no NA found, set flag to FALSE to exit the loop
            has_na <- FALSE
          }
          
        }
        
      }
      
      
      
    }
    
    # If no NA found, proceed with downstream analysis steps
    if (!has_na) {
      # Classify stability based on item_consensus values using custom function
      item_consensus_df$stability <- sapply(item_consensus_df$item_consensus, get_stability_category)
      
      # Split the dataframe into a list of dataframes, one per cluster number k
      item_consensus_by_k <- split(item_consensus_df, item_consensus_df$k)
    }
    
    
    
    # Loop over each k from ConsensusClusterPlus results
    CCP_plot_list <- list() 
    
    
    for (k in 2:length(consensus_results_app_two)) {
      cm <- consensus_results_app_two[[k]]$consensusMatrix
      colnames(cm) <- colnames(final_input_matrix)
      rownames(cm) <- colnames(final_input_matrix)  # consensusMatrix is symmetric, matches samples
      
      cm <-  cm[consensus_results_app_two[[k]]$consensusTree$order,]
      
      
      clusters <- consensus_results_app_two[[k]]$consensusClass
      order_idx <- order(clusters)
      
      # Create annotation dataframe
      ann <- data.frame(Cluster = factor(clusters[order_idx]))
      #rownames(ann) <- names(clusters)
      
      # Step 1: Ensure Cluster is a factor with only the levels that exist
      ann$Cluster <- as.factor(ann$Cluster)
      ann$Cluster <- droplevels(ann$Cluster)  # drop any unused levels
      
      # Step 2: Define annotation colors correctly
      cluster_levels <- levels(ann$Cluster)
      
      # RColorBrewer requires at least 3 colors, so ensure enough levels
      num_clusters <- length(cluster_levels)
      palette <- RColorBrewer::brewer.pal(max(3, num_clusters), "Set2")[1:num_clusters]
      
      # Step 3: Set names to match factor levels exactly
      annotation_colors <- list(Cluster = setNames(palette, cluster_levels))
      
      
      # Save heatmap to list
      CCP_plot <- pheatmap(
        cm,
        cluster_rows = FALSE,
        cluster_cols = consensus_results_app_two[[k]]$consensusTree,
        annotation_col = ann,
        annotation_row = ann,
        show_colnames = FALSE,
        annotation_colors = annotation_colors,
        color = colorRampPalette(c("white", "blue"))(100),
        silent = TRUE,
        main = paste("Consensus Matrix for k =", k)
      )
      
      # Make name 
      CCP_plot_list_name <- paste0("CCP_k_", k)
      CCP_plot_list[[CCP_plot_list_name]] <- CCP_plot
    }
    
    
    # Create graphs for item consensus and table 
    item_consensus_overall <- as.data.frame(icl_results[["clusterConsensus"]])
    
    
    
    
    # If you want to round the labels nicely, e.g. 2 decimal places
    item_consensus_overall <- item_consensus_overall %>%
      mutate(label = round(clusterConsensus, 2))
    
    item_con_plot <- ggplot(item_consensus_overall,
                            aes(
                              x = factor(k),
                              y = clusterConsensus,
                              fill = factor(cluster)
                            )) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_text(aes(label = label), 
                position = position_dodge(width = 0.9),
                vjust = -0.3,  # slightly above the bar
                size = 6) +    # adjust text size
      labs(title = "Cluster Consensus by k", x = "k", y = "Consensus Score") +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        axis.text = element_text(colour = "black", size = 24),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.text.y = element_text(colour = "black", size = 24),
        axis.title = element_text(colour = "black", size = 30),
        axis.ticks = element_line(size = 1.5),
        panel.border = element_rect(
          colour = "black",
          fill = NA,
          size = 1
        )
      ) +
      scale_fill_manual(values = cols[1:working_max_k], name = "Cluster", drop = FALSE)
    
    
    item_consensus_by_k_plot_list <- list()
    for(k in names(item_consensus_by_k)){
      
      k_numeric <- as.numeric(k)
      item_con_dataset <- item_consensus_by_k[[as.character(k)]]
      
      
      item_consensus_by_k_plot <- ggplot(item_con_dataset,
                                         aes(
                                           x = sample_id,
                                           y = item_consensus,
                                           fill = factor(cluster)
                                         )) +
        geom_bar(stat = "identity") +
        labs(title = paste0("Item Consensus by Sample (k = ", k, ")"),
             x = "Sample ID",
             y = "Item Consensus",
             fill = "Cluster") +
        theme_classic() +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = 8
          ),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 14),
          legend.position = "right"
        ) +
        scale_fill_manual(values = cols[1:working_max_k], name = "Cluster", drop = FALSE)
      
      item_consensus_by_k_plot_list[[as.character(k)]] <- item_consensus_by_k_plot
      
    }
    
    
    # Specify custom results directory
    results_dir <- paste0(working_dir, "/CIM_states_results_", clustering_alg, "_", input_data_type,"_app_two")
    if(!dir.exists(results_dir)) {
      dir.create(results_dir, recursive = TRUE)
    }
    
    
    # Get Clustering gene set for this approach
    write.csv(clustering_gene_set,
              file = file.path(results_dir, paste0("clustering_gene_set.csv")),
              row.names = FALSE)
    
    
    
    # Save tables for each k into CSV files
    for (k_val in names(item_consensus_by_k)) {
      write.csv(item_consensus_by_k[[k_val]],
                file = file.path(results_dir, paste0("item_consensus_k", k_val, ".csv")),
                row.names = FALSE)
    }
    
    # Save Plots 
    
    ggsave(
      filename = paste0(results_dir, "/item_con_plot.png"),
      plot = item_con_plot,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    
    for (k in names(item_consensus_by_k_plot_list)){
      
      img <- item_consensus_by_k_plot_list[[as.character(k)]]
      
      ggsave(
        filename = paste0(results_dir, "/item_con_sample_k-", as.character(k), ".png"),
        plot = img,
        width = 10,
        height = 6,
        dpi = 300
      )
      
      
    }
    
    for (i in seq_along(CCP_plot_list)) {
      plot_name <- names(CCP_plot_list)[i]
      png_filename <- file.path(results_dir, paste0(plot_name, ".png"))
      png(png_filename, width = 1200, height = 1000, res = 150)
      grid.newpage()
      obj <- CCP_plot_list[[i]]
      if (inherits(obj, "grob")) {
        grid.draw(obj)
      } else if (inherits(obj, "Heatmap")) {
        draw(obj)
      } else if (!is.null(obj$gtable)) {
        grid.draw(obj$gtable)   # <-- fallback: draw the gtable element if it exists
      } else {
        warning(sprintf("Element %s is not drawable by grid.draw, draw, or print.", plot_name))
      }
      dev.off()
    }
    
    
    
    # Create cdf_plots
    k_vec = 2:working_max_k
    cdf_data <- data.frame()
    PACs <- numeric(length(k_vec))
    
    for(i in seq_along(k_vec)) {
      k <- k_vec[i]
      consensus_matrix <- consensus_results_app_two[[k]][["consensusMatrix"]]
      cons_vals <- consensus_matrix[lower.tri(consensus_matrix)]
      Fn <- ecdf(cons_vals)  # standard consensus cluster PAC
      PACs[i] <- Fn(0.9) - Fn(0.1)
      h <- hist(cons_vals, breaks = 200, plot = FALSE)
      cdf_y <- cumsum(h$counts) / sum(h$counts)
      cdf_data <- bind_rows(
        cdf_data,
        data.frame(
          k = as.factor(k),
          x = h$mids,
          y = cdf_y
        )
      )
    }
    
    # Now join the correct PAC values
    pac_df <- data.frame(k = as.factor(k_vec), PAC = PACs)
    cdf_data <- left_join(cdf_data, pac_df, by = "k")
    cdf_data$line_label <- paste0("K = ", cdf_data$k, " (PAC=", formatC(cdf_data$PAC, format = "f", digits = 3), ")")
    
    
    cdf_plot <- ggplot(cdf_data, aes(x = x, y = y, color = line_label, group = k)) +
      geom_line(size = 1.5) +
      labs(
        x = "Consensus Index",
        y = "Cumulative Fraction",
        title = "Consensus CDFs across Ks",
        color = "K (PAC)"
      ) +
      geom_vline(xintercept = c(0.1, 0.9), linetype = "dashed", color = "gray") +
      scale_x_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
      scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.98, 0.02), # x, y from 0 (left/bottom) to 1 (right/top)
        legend.justification = c("right", "bottom"),
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill="white", color="black", size=0.8),
        axis.text = element_text(colour = "black", size = 24),
        axis.text.x = element_text(colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title = element_text(colour = "black", size = 30),
        axis.ticks = element_line(size = 1.5),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)
      )
    
    
    # save cdf plot
    ggsave(
      filename = paste0(results_dir, "/final_cdf_plot.png"),
      plot = cdf_plot,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    
    # Calculate PAC metric (Proportion of Ambiguous Clustering) using ecdf
    Kvec = 2:working_max_k
    x1 = 0.1; x2 = 0.9 
    PAC = rep(NA, length(Kvec)) 
    names(PAC) = paste("K=", Kvec, sep="") 
    for(i in Kvec){
      M = consensus_results_app_two[[i]]$consensusMatrix
      Fn = ecdf(M[lower.tri(M)])
      PAC[i-1] = Fn(x2) - Fn(x1)
    }
    optK = Kvec[which.min(PAC)]
    
    # Convert PAC to data frame for plotting
    pac_values <- data.frame(k = Kvec, pac = unname(PAC))
    
    # Plot PAC values (lower is better) and display in R environment (Zaoqu-Liu/IRLS GitHub)
    pac_plot <- ggplot(pac_values, aes(factor(k), pac, group=1)) +
      geom_line() +
      theme_bw(base_rect_size = 1.5) +
      geom_point(size=4, shape=21, color='darkred', fill='orange') +
      #ggtitle('Proportion of ambiguous clustering') +
      xlab('Number of clusters K') + ylab("PAC Score") +
      theme_classic() + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top",
            axis.text = element_text(colour = "black", size = (24)),
            axis.text.x = element_text(colour = "black", size = (24)),
            axis.text.y = element_text(colour = "black", size = (24)),
            axis.title = element_text(colour = "black", size = (30)),
            axis.ticks = element_line(size = 1.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    
    
    
    # Create item consensus table
    avg_cluster_consensus <- item_consensus_overall %>%
      group_by(k) %>%
      summarize(avg_cluster_consensus = mean(clusterConsensus))
    
    # Silhouette calculations using consensus cluster final
    
    orig_space_input <- t(final_input_matrix)  # transpose to get dims as samples x features
    rownames(orig_space_input) <- rownames(clustering_matrix)
    
    # Distances
    dist_space_for_clustering <- dist(orig_space_input)
    dist_ge_space <- dist(clustering_matrix[,clustering_gene_set, drop=FALSE])
    
    # Collect silhouettes for all k
    sil_results <- lapply(2:length(consensus_results_app_two), function(k) {
      
      clusters <- consensus_results_app_two[[k]]$consensusClass
      
      # Silhouette in gene space
      sil_space_for_clustering   <- silhouette(clusters, dist_space_for_clustering)
      avg_space_for_clustering   <- mean(sil_space_for_clustering[, 3])
      
      # Silhouette in PaCMAP space
      sil_ge_space <- silhouette(clusters, dist_ge_space)
      avg_ge_space <- mean(sil_ge_space[, 3])
      
      data.frame(
        k = k,
        avg_sil_space_for_clustering = avg_space_for_clustering,
        avg_sil_ge_space = avg_ge_space
      )
    })
    
    # Combine
    sil_df <- do.call(rbind, sil_results)
    sil_df$avg_sil_combined <- ((sil_df$avg_sil_space_for_clustering +
        sil_df$avg_sil_ge_space) /
        2)

    
    # Plot comparison
    sil_plot_both <- ggplot(sil_df, aes(x = k)) +
      geom_line(aes(y = avg_sil_space_for_clustering, color = paste0(input_data_type))) +
      geom_point(aes(y = avg_sil_space_for_clustering, color = paste0(input_data_type))) +
      geom_line(aes(y = avg_sil_ge_space, color = "Gene space")) +
      geom_point(aes(y = avg_sil_ge_space, color = "Gene space")) +
      labs(y = "Average silhouette width", color = "Space", title = "Silhouette comparison across k") +
      theme_classic(base_size = 16) + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top",
            axis.text = element_text(colour = "black", size = 24),
            axis.title = element_text(colour = "black", size = 30),
            axis.ticks = element_line(size = 1.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    
    

    sil_plot_original_space <- ggplot(sil_df, aes(x = k)) +
      geom_line(aes(y = avg_sil_space_for_clustering, color = paste0(input_data_type))) +
      geom_point(aes(y = avg_sil_space_for_clustering, color = paste0(input_data_type))) +
      labs(y = "Average silhouette width", color = "Space", title = "Silhouette comparison across k") +
      theme_classic(base_size = 16) + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top",
            axis.text = element_text(colour = "black", size = 24),
            axis.title = element_text(colour = "black", size = 30),
            axis.ticks = element_line(size = 1.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    

    sil_plot_ge_space <- ggplot(sil_df, aes(x = k)) +
      geom_line(aes(y = avg_sil_ge_space, color = "Gene space")) +
      geom_point(aes(y = avg_sil_ge_space, color = "Gene space")) +
      labs(y = "Average silhouette width", color = "Space", title = "Silhouette comparison across k") +
      theme_classic(base_size = 16) + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top",
            axis.text = element_text(colour = "black", size = 24),
            axis.title = element_text(colour = "black", size = 30),
            axis.ticks = element_line(size = 1.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    

  
    
    # save silhouette plot
    ggsave(
      filename = paste0(results_dir, "/sil_plot.png"),
      plot = sil_plot_both,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    
    # save silhouette plot
    ggsave(
      filename = paste0(results_dir, "/sil_plot_original_space.png"),
      plot = sil_plot_original_space,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    # save silhouette plot
    ggsave(
      filename = paste0(results_dir, "/sil_plot_ge_space.png"),
      plot = sil_plot_ge_space,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    
    # save pac plot
    ggsave(
      filename = paste0(results_dir, "/final_pac_plot.png"),
      plot = pac_plot,
      width = 10,
      height = 6,
      dpi = 300
    )
    

   

    
       
    # Get suggested optimal k values
    silhouette_dim_space_k <- sil_df$k[which.max(sil_df$avg_sil_space_for_clustering)]
    silhouette_ge_space_k <- sil_df$k[which.max(sil_df$avg_sil_ge_space)]
    silhouette_combined_avg_k <- sil_df$k[which.max(sil_df$avg_sil_combined)]
    pac_k <- Kvec[which.min(PAC)]  # Using optK from PAC calculation
    avg_cluster_consensus_k <- avg_cluster_consensus[which.max(avg_cluster_consensus$avg_cluster_consensus), ]$k
    
    # Print metric values
    cat("\n===== OPTIMAL CLUSTER SUGGESTIONS =====\n")
    cat("Silhouette analysis in dim space suggests k =", silhouette_dim_space_k, "\n")
    cat("Silhouette analysis in GE space suggests k =", silhouette_ge_space_k, "\n")
    cat("Silhouette analysis of combined avg (GE and red. dim) space suggests k =", silhouette_combined_avg_k, "\n")
    cat("PAC method suggests k =", pac_k, "\n")
    cat("cluster consensus method suggests k =", avg_cluster_consensus_k, "\n")
    
    
    # Create metrics summary dataframe
      metrics_summary <- data.frame(
        k = 2:working_max_k,
        silhouette_dim_reduce_space = NA, # renamed to indicate dim‑reduce (original) space
        silhouette_ge_space = NA, # silhouette computed in PacMAP (ge) space
        silhouette_combined_avg = NA, # silhouette computed in PacMAP (ge) space
        pac = NA,
        avg_cluster_consensus = NA
      )

      # Fill in the metrics
      for (i in 1:nrow(metrics_summary)) {
        k_val <- metrics_summary$k[i]

        # Silhouette in the original (dim‑reduce) space
        idx <- which(sil_df$k == k_val)
        if (length(idx) > 0) {
          metrics_summary$silhouette_dim_reduce_space[
            i
          ] <- sil_df$avg_sil_space_for_clustering[idx]

          metrics_summary$silhouette_ge_space[i] <- sil_df$avg_sil_ge_space[idx]
          metrics_summary$silhouette_combined_avg[i] <- sil_df$avg_sil_combined[
            idx
          ]
        }

        # PAC
        k_idx <- which(Kvec == k_val)
        if (length(k_idx) > 0) {
          metrics_summary$pac[i] <- PAC[k_idx]
        }

        # avg_cluster_consensus
        cc_idx <- which(avg_cluster_consensus$k == k_val)
        if (length(cc_idx) > 0) {
          metrics_summary$avg_cluster_consensus[
            i
          ] <- avg_cluster_consensus$avg_cluster_consensus[cc_idx]
        }
      }

      # Alternatively,  Optimal K can be selected using rank based compromise
      # The rank-based compromise approach assigns a rank to each k for both silhouette (descending) and PAC (ascending), then sums the ranks to identify the k with the best overall trade-off between high silhouette and low PAC. The k with the lowest combined rank is selected as the optimal choice.
      # https://dl.acm.org/doi/10.1145/371920.372165

      # Assign ranks (lowest rank = best)
      metrics_summary$rank_silhouette_dim_reduce_space <- rank(
        -metrics_summary$silhouette_dim_reduce_space,
        ties.method = "min"
      ) # higher silhouette = rank 1 (original space)
      metrics_summary$rank_silhouette_ge_space <- rank(
        -metrics_summary$silhouette_ge_space,
        ties.method = "min"
      ) # higher silhouette = rank 1 (ge space)

      metrics_summary$rank_silhouette_combined_avg <- rank(
        -metrics_summary$silhouette_combined_avg,
        ties.method = "min"
      ) # higher silhouette = rank 1 (averaged space)
      metrics_summary$rank_pac <- rank(metrics_summary$pac, ties.method = "min") # lower PAC = rank 1
      metrics_summary$rank_item_cluster_consensus <- rank(
        -metrics_summary$avg_cluster_consensus,
        ties.method = "min"
      ) # higher consensus = rank 1

      # Build the column names dynamically.
      rank_columns <- c(
        paste0("rank_", clustering_metrics) # e.g. "rank_pac", "rank_item_cluster_consensus"
      )

      # Compute the overall rank as the row sum across selected rank columns
      metrics_summary$overall_rank <- rowSums(
        metrics_summary[, rank_columns, drop = FALSE],
        na.rm = TRUE
      )

      # Find all rows with the minimum overall rank
      best_rows <- metrics_summary[
        metrics_summary$overall_rank == min(metrics_summary$overall_rank),
      ]

      # Use the last metric in clustering_metrics as the tiebreaker
      tie_metric <- clustering_metrics[length(clustering_metrics)]

      if (verbose) {
        cat(
          "Metrics being used for ranking:\n  ",
          paste(clustering_metrics, collapse = ", "),
          "\n"
        )
        cat(
          "Tie‑breaker metric (last in clustering_metrics vector):",
          tie_metric,
          "\n"
        )
      }

      # Determine if the tiebreaking metric should be maximized or minimized
      # (Assuming: max for 'silhouette' and 'item_cluster_consensus', min for 'pac')
      if (tie_metric == "silhouette_dim_reduce_space") {
        idx <- which.max(best_rows$silhouette_dim_reduce_space)
      } else if (tie_metric == "silhouette_ge") {
        idx <- which.max(best_rows$silhouette_ge_space)
      } else if (tie_metric == "item_cluster_consensus") {
        idx <- which.max(best_rows$avg_cluster_consensus)
      } else if (tie_metric == "pac") {
        idx <- which.min(best_rows$pac)
      } else {
        stop("Unknown tiebreaker metric.")
      }

      best_k_row <- best_rows[idx, ]
      best_k <- best_k_row$k

      # add cluster metrics used for evaluation and tie
      metrics_summary$metrics_used <- paste(clustering_metrics, collapse = ", ")
      metrics_summary$metrics_used_for_tie <- tie_metric

    # Print metrics summary
    cat("\n===== CLUSTER METRICS SUMMARY =====\n")
    metrics_summary
    cat("\n")
    
    # Write metrics to file
    write.csv(metrics_summary, file.path(results_dir, "final_cluster_metrics_summary.csv"), row.names = FALSE)
    
    # Based on examining the plots and metrics, set your optimal k manually here:
    # CHANGE THIS VALUE after reviewing the plots:
    optimal_k <- best_k  # <-- MODIFY THIS BASED ON YOUR EXAMINATION IF NEEDED 
    
    cat("Using k =", optimal_k, "as the optimal number of clusters\n")
    
    # Extract cluster assignments for optimal k
    consensus_clusters <- consensus_results_app_two[[optimal_k]]$consensusClass
    
    
    ####### FINAL VISUALIZATION #######
    
    
    
    # Get stability/reliability metrics 
    optimal_stability_cluster_metrics <- item_consensus_by_k[[as.character(optimal_k)]]
    
    # Start from your original data frame
    wide_optimal_stability_cluster_metrics <- optimal_stability_cluster_metrics %>%
      # Only keep relevant columns
      dplyr::select(sample_id, cluster, item_consensus) %>%
      # Pivot wider: each cluster makes a new column
      pivot_wider(
        names_from = cluster,
        values_from = item_consensus,
        names_prefix = "item_consensus_cluster"
      )
    
    
    # Write metrics to file
    write.csv(optimal_stability_cluster_metrics,
              file.path(results_dir, "final_cluster_stability_metrics.csv"), row.names = FALSE)
    
    write.csv(wide_optimal_stability_cluster_metrics,
              file.path(results_dir, "final_wide_cluster_stability_metrics.csv"), row.names = FALSE)
    
    # Create data frame for UMAP visualization
    final_df_plot_app_two <- as.data.frame(umap_matrix_visual$layout) %>%
      rownames_to_column("sample_id") %>% 
      dplyr::rename(UMAP1 = V1, UMAP2 = V2) %>%
      dplyr::inner_join(as.data.frame(pacmap_delta_embeddings_visual) %>% 
                          rownames_to_column("sample_id") %>% 
                          dplyr::rename(PACMAP1 = V1, PACMAP2 = V2)) %>%
      dplyr::inner_join(as.data.frame(pca_scores) %>% 
                          rownames_to_column("sample_id")) %>%
      dplyr::mutate(base_id = gsub("^delta_([^_]+)_.*", "\\1", sample_id)) %>%
      dplyr::mutate(cluster_assignments = consensus_clusters) %>%
      dplyr::inner_join(as.data.frame(clustering_matrix) %>%
                          rownames_to_column("sample_id"),
                        by = "sample_id")
    
    
    # Your input matrix: samples x genes
    pca_fviz_input <- final_df_plot_app_two[, c("PC1", "PC2", "sample_id")] %>%
      tibble::column_to_rownames("sample_id")
    
    
    fviz_cluster_plot <- fviz_cluster(
      object = list(data = pca_fviz_input, cluster = as.integer(consensus_clusters)),
      palette = cols[1:optimal_k],
      geom = "text",
      ellipse.type = "convex",
      ggtheme = theme_bw(),
      ggtitle("PCA plot on Clustered Data")
    )
    
    fviz_cluster_plot
    ggsave(file.path(results_dir, paste0("final_fviz_2D_PCA_k-", optimal_k, ".png")), 
           fviz_cluster_plot, width = 10, height = 8, dpi = 300)
    
    
    # 2. PCA projection using selected PCs
    pca_scores_3D <- pca_result$x[, 1:3, drop = FALSE]
    
    # 3. Prepare input for downstream analysis, e.g., clustering
    # Here, samples x PCs (rows = samples, cols = PCs)
    pca_fviz_input_3D <- as.matrix(pca_scores_3D)
    rownames(pca_fviz_input_3D) <- final_df_plot_app_two$sample_id  # optional, set sample IDs
    
    
    
    # If pca_fviz_input_3D is samples x PCs, select first 3 PCs
    
    plot_df_3D_PCA <- data.frame(
      PC1 = pca_fviz_input_3D[, 1],
      PC2 = pca_fviz_input_3D[, 2],
      PC3 = pca_fviz_input_3D[, 3],
      cluster = as.factor(as.integer(consensus_clusters)),
      sample_id = rownames(pca_fviz_input_3D)
    )
    
    
    
    # 3D scatter plot
    fig_PCA_3D <- plot_ly(
      data = plot_df_3D_PCA,
      x = ~PC1, y = ~PC2, z = ~PC3,
      color = ~cluster,
      text = ~sample_id,
      type = "scatter3d",
      mode = "markers",
      colors = cols
    ) %>%
      plotly::layout(
        title = "3D PCA Cluster Visualization",
        scene = list(
          xaxis = list(title = "PC1"),
          yaxis = list(title = "PC2"),
          zaxis = list(title = "PC3")
        )
      )
    
    fig_PCA_3D
    # Load library if not already loaded
    library(htmlwidgets)
    
    # Save as an interactive HTML file
    saveWidget(
      widget = fig_PCA_3D,
      file = paste0(results_dir, "/final_3D_PCA_k-", optimal_k, ".html"),
      selfcontained = TRUE
    )
    
    
    # Use letters as shapes 
    shapes <- unlist(lapply(letters, utf8ToInt))
    
    # Plot using GGPlot Color by cluster
    final_cluster_plot_umap <- ggplot(
      final_df_plot_app_two,
      aes(
        x = UMAP1,
        y = UMAP2,
        label = base_id ,
        color = as.factor(cluster_assignments) 
      )
    ) +
      scale_color_manual(values = cols, name = "Cluster") +
      ggforce::geom_mark_hull(
        aes(group = as.factor(cluster_assignments),
            fill = as.factor(cluster_assignments), label = NULL),
        concavity = 2, expand = unit(2, "mm"), alpha = 0.2
      ) +
      geom_label_repel(size = 4, fontface = "bold", label.size = 0.3, fill = "white", 
                       alpha = 1, box.padding = 0.25, max.overlaps = 100) +
      labs(x = "UMAP 1", y = "UMAP 2",
           fill = "Cluster",
           color = "Cluster",
           title = paste0(
             "Final Clustering Using ", input_data_type, " input and ",
             clustering_alg, " CCP Alg (k = ", optimal_k, ")"
           )) +
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
    
    ggsave(file.path(results_dir, paste0("final_clustering_shown_on_umap_k-", optimal_k, ".png")), 
           final_cluster_plot_umap, width = 10, height = 8, dpi = 300)
    
    
    # Plot using GGPlot Color by cluster
    final_cluster_plot_pacmap <- ggplot(
      final_df_plot_app_two,
      aes(
        x = PACMAP1,
        y = PACMAP2,
        label = base_id,  # Your patient/sample label column
        color = as.factor(cluster_assignments)
      )
    ) +
      scale_color_manual(values = cols, name = "Cluster") +
      ggforce::geom_mark_hull(
        aes(group = as.factor(cluster_assignments),
            fill = as.factor(cluster_assignments), label = NULL),
        concavity = 2, expand = unit(2, "mm"), alpha = 0.2
      ) +
      geom_label_repel(size = 4, fontface = "bold", label.size = 0.3, fill = "white", 
                       alpha = 1, box.padding = 0.25, max.overlaps = 100) +
      labs(x = "PacMAP 1", y = "PacMAP 2",
           fill = "Cluster",
           color = "Cluster",
           title = paste0(
             "Final Clustering Using ", input_data_type, " input and ",
             clustering_alg, " CCP Alg (k = ", optimal_k, ")"
           )) +
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
    
    
    ggsave(file.path(results_dir, paste0("final_clustering_shown_on_pacmap_k-", optimal_k, ".png")), 
           final_cluster_plot_pacmap, width = 10, height = 8, dpi = 300)
    
  }
  
  }
  
  if (no_app_two == TRUE) {
    cat(
      "\033[1;31mWARNING: Final Gene Set for Approach two empty! Consider Relaxing P-val
        Approach output will not be produced\033[0m\n"
    )
    
    # set string
    no_app_two_string <- "Approach two final gene set empty; Consider Relaxing P-val"
    
    
    # update iter_log
    iter_log <- rbind(
      iter_log,
      data.frame(
        gene_set = "app_two final get set", iteration = NA,
        before = NA,
        after = no_app_two_string
      )
    )
    
    
    # set df to null
    final_df_plot_app_two <- NULL
  }
  
  # Extra capturing of output 
  {
    # Write initial matrix to file
    write.csv(clustering_matrix,
              file.path(results_dir, "initial_clustering_mat.csv"), row.names = TRUE)
    
    # Write iter_log  to file
    write.csv(iter_log,
              file.path(results_dir, "iter_log_all_approaches.csv"), row.names = FALSE)
    
    # Convert gene_sets to to data frame
    df_gene_sets_for_file <- data.frame(
      list_name = names(all_gene_sets),
      items = sapply(all_gene_sets, function(x) paste(x, collapse = ","))
    )
    
    # Write gene_sets to CSV 
    write.csv(df_gene_sets_for_file,  
              file.path(results_dir, "all_gene_sets.csv"), row.names = FALSE)
    
    # write final_clustered df to file
    write.csv(final_df_plot_app_two,  
              file.path(results_dir, "clustered_samples_app_two.csv"), row.names = FALSE)
    
    # write pacmap settings 
    write.csv(pm_settings,  
              file.path(results_dir, "pacmap_settings.csv"), row.names = TRUE)
    
    # write removed genes setting 
    write.csv(removed_genes,  
              file.path(results_dir, "removed_genes.csv"), row.names = TRUE)
  }
  
  })
    
  }
  
  # ---- Approach Three (Iterated over all unique genes ) ---- 
  
  
  if ("app_three" %in% filter_approach){
  
    # Set seed appropriately
    withr::with_seed(seed = 2024L, {
    
  # set back to max_k
  working_max_k <- max_k
  
  if(no_app_three == FALSE){
  
  {   
    # CHOOSE INPUT DATA TYPE FOR CLUSTERING
    # Set to "gene_set" for original gene expression data, or "umap" for UMAP coordinates, or "pca"
    input_data_type <- "pacmap"  # <-- MODIFY HERE: use "genes" or "umap" or "pca" or "pacmap"
    
    # Make matrix for delta
    #toggle here
    clustering_gene_set <- iterated_use_final_gene_set_from_approach_one
    clustering_gene_set <- clustering_gene_set[clustering_gene_set %in% colnames(clustering_matrix)]
    
    # Run UMAP on the selected gene set (same as in the original code)
    # will always be 2D 
    umap_matrix_visual <- umap::umap(clustering_matrix[,clustering_gene_set, drop=FALSE],
                              config=custom.config)
    
    
    # Build embedder for analysis 
    script_pacmap_reducer <- do.call(pacmap$PaCMAP, pm_defaults)
    
    
    # Run PaCMAP
    pacmap_delta_embeddings <- script_pacmap_reducer$fit_transform(
      clustering_matrix[, clustering_gene_set, drop = FALSE]
    )
    
    # set rownames 
    rownames(pacmap_delta_embeddings) <- rownames(clustering_matrix)
    
    
    # 2D visualizations 
    # create reducer (also resets seed)
    pm_visualization_reducer <- do.call(pacmap$PaCMAP, pm_visualization)
    
    # Run PACMAP on selected gene set
    # Run PaCMAP
    pacmap_delta_embeddings_visual <- pm_visualization_reducer$fit_transform(
      clustering_matrix[, clustering_gene_set, drop = FALSE]
    )
    
    # set rownames 
    rownames(pacmap_delta_embeddings_visual) <- rownames(clustering_matrix)
    
    
    
    # Your input matrix: samples x genes
    pca_matrix <- clustering_matrix[, clustering_gene_set]
    
    pca_result <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)
    pca_scores <- pca_result$x[, 1:3, drop = FALSE]
    
    
    
    
    if (verbose) message(sprintf(
      "\n\n=== Creating Output for Approach Three using %s algorithm and %s as input ===",
      clustering_alg,
      input_data_type
    ))
    
    
    
    
    # Define color palette for clusters
    cols <- hue_pal()(working_max_k)
    
    ########################################################################################
    
    
    # Create input matrix for consensus clustering based on chosen data type
    if (input_data_type == "genes") {
      # Use gene data
      final_input_matrix <- t(clustering_matrix[, clustering_gene_set])  # transpose to get genes as rows
      cat("Using gene expression data for clustering.\n")
    } else if (input_data_type == "umap") {
      # Use UMAP coordinates - specifically from umap_matrix to match your original approach
      final_input_matrix <- t(umap_matrix$layout)  # transpose to get dimensions as rows
      cat("Using UMAP coordinates for clustering.\n")
    } else if (input_data_type == "pca"){
      
      # Your input matrix: samples x genes
      final_input_matrix <- clustering_matrix[, clustering_gene_set]
      
      # Make sure your final_input_matrix is samples x features (rows=samples, cols=genes)
      # You can check with dim(final_input_matrix)
      
      # 1. Run Parallel Analysis to determine number of PCs
      # PA expects samples in rows, variables in columns
      centile_seq <- seq(95, 80, by = -5)
      max_pcs <- 3   # Default value in case paranormal fails
      pa_result <- NULL
      for (cent in centile_seq) {
        pa_result <- try(
          suppressMessages(suppressWarnings(
            paran(final_input_matrix,
                  iterations = 1000,
                  centile = cent,
                  graph = FALSE,
                  quietly = TRUE)
          )), silent = TRUE)
        
        if (!inherits(pa_result, "try-error")) {
          message(sprintf("paran succeeded with centile = %d", cent))
          max_pcs <- pa_result$Retained  # number of components retained by PA
          break
        }
      }
      if (inherits(pa_result, "try-error")) {
        warning("paran failed at all tested centile values (95 to 80). Using default max_pcs = 3 for this gene set.")
        # continue; max_pcs remains 3
      }
      
      cat("Number of PCs retained by Parallel Analysis:", max_pcs, "\n")
      
      # 2. PCA projection using selected PCs
      pca_result <- prcomp(final_input_matrix, center = TRUE, scale. = TRUE)
      pca_scores <- pca_result$x[, 1:max_pcs, drop = FALSE]
      
      # 3. Prepare input for downstream analysis, e.g., clustering
      # Here, PCs x Samples (rows = PCs, cols = samples)
      final_input_matrix <- as.matrix(t(pca_scores))
      colnames(final_input_matrix) <- umap_plot_df$sample_id  # optional, set sample IDs
      
    } else if (input_data_type == "pacmap"){
      
      final_input_matrix <- t(pacmap_delta_embeddings)
      colnames(final_input_matrix) <- rownames(clustering_matrix)
    } else {
      stop("Invalid input_data_type. Use 'gene_set' or 'umap' or 'pca' or 'pacmap'.")
    }
    
    
    
    # Initialize starting value for pItem (percentage of items to sample per iteration)
    pItem_val <- 0.8
    
    # Set maximum allowable pItem value (cannot exceed 1.0)
    max_pItem <- 1.0
    
    # Step size for incrementing pItem when NA values are found
    step <- 0.05
    
    # Flag to track presence of NA values in itemConsensus
    has_na <- TRUE
    
    # if in the case of working max k = 2 
    if (working_max_k == 2){
      
      # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
      while (has_na && pItem_val <= max_pItem) {
        # Run ConsensusClusterPlus with the current pItem_val
        consensus_results_app_three <- ConsensusClusterPlus(
          final_input_matrix, 
          maxK = 3, 
          reps = CCP_iter,            
          pItem = pItem_val,            # Percentage of items sampled per iteration
          pFeature = 1,                # Use all features
          clusterAlg = clustering_alg,      
          innerLinkage = "ward.D2",
          finalLinkage = "ward.D2",
          distance = "euclidean", 
          seed = 2024L,
          plot ="none",           
          verbose = FALSE,
          writeTable = FALSE
        )
        
        # Only extract out result for K = 2
        consensus_results_app_three <- consensus_results_app_three[1:2]
        
        # Calculate item consensus and cluster metrics (ICL)
        pdf(file = NULL)  # Start invisible plotting device
        icl_results <- calcICL(consensus_results_app_three, 
                               plot = "png", 
                               writeTable = TRUE)
        dev.off()          # Close device
        
        # Extract the itemConsensus dataframe from ICL results
        item_consensus_df <- icl_results[["itemConsensus"]]
        
        # Rename columns for clarity
        colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
        
        
        # Calculate item consensus and cluster metrics (ICL)
        pdf(file = NULL)  # Start invisible plotting device
        icl_results <- calcICL(consensus_results_app_three, 
                               plot = "png", 
                               writeTable = TRUE)
        dev.off()          # Close device
        
        # Extract the itemConsensus dataframe from ICL results
        item_consensus_df <- icl_results[["itemConsensus"]]
        
        # Rename columns for clarity
        colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
        
        # Check if any values in the item_consensus column are NA
        if (any(is.na(item_consensus_df$item_consensus))) {
          # Print message indicating NA detected and that pItem will be increased
          message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
          
          # Increase pItem by the predefined step to sample more items per iteration
          pItem_val <- pItem_val + step
        } else {
          # If no NA found, set flag to FALSE to exit the loop
          has_na <- FALSE
        }
        
      } 
      
    }
    else {
      
      # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
      while (has_na && pItem_val <= max_pItem) {
        
        # Run ConsensusClusterPlus with the current pItem_val
        consensus_results_app_three <- ConsensusClusterPlus(
          final_input_matrix, 
          maxK = working_max_k, 
          reps = CCP_iter,            
          pItem = pItem_val,            # Percentage of items sampled per iteration
          pFeature = 1,                # Use all features
          clusterAlg = clustering_alg,      
          innerLinkage = "ward.D2",
          finalLinkage = "ward.D2",
          distance = "euclidean", 
          seed = 2024L,
          plot ="none",           
          verbose = FALSE,
          writeTable = FALSE
        )
        
        # Calculate item consensus and cluster metrics (ICL)
        pdf(file = NULL)  # Start invisible plotting device
        icl_results <- calcICL(consensus_results_app_three, 
                               plot = "png", 
                               writeTable = TRUE)
        dev.off()          # Close device
        
        # Extract the itemConsensus dataframe from ICL results
        item_consensus_df <- icl_results[["itemConsensus"]]
        
        # Rename columns for clarity
        colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
        
        # Check if any values in the item_consensus column are NA
        if (any(is.na(item_consensus_df$item_consensus))) {
          # Print message indicating NA detected and that pItem will be increased
          message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
          
          # Increase pItem by the predefined step to sample more items per iteration
          pItem_val <- pItem_val + step
        } else {
          # If no NA found, set flag to FALSE to exit the loop
          has_na <- FALSE
        }
      }
    }
    
    # After loop finishes, check if NA values still exist (meaning max pItem was reached)
    if (has_na) {
      # Warn the user that maximum pItem was reached but NAs persist
      message("Reached max pItem but NA values still exist in item_consensus.")
      message("Reducing working_max_k to the maximum most stable K")
      
      # Find all K values with no NA item consensus
      ks_without_na <- item_consensus_df %>%
        group_by(k) %>%
        summarize(no_na = all(!is.na(item_consensus))) %>%
        filter(no_na) %>%
        pull(k)
      
      # if there is no solution skip it 
      if (length(ks_without_na) == 0) {
        
        message("No stable K values found. Skipping this data generation for this approach")
        
        # remove so data isn't produced  
        filter_approach <- filter_approach[filter_approach != "app_three"]
        break  # exits the  loop
        
        
        
      } else {
        # Get maximal most stable K  
        max_stable_k <- max(ks_without_na)
        message("maximum most stable K = ", max_stable_k)
        
      }
      
      # special condition for if maximum most stable K = 2 
      if (max_stable_k == 2){
        
        
        # let user know 
        message("since maximum most stable K = ", max_stable_k)
        message("Will have to extract k = 2 results using k = 3 analysis")
        
        # set that to working_k 
        working_max_k <- 2
        
        
        # Initialize starting value for pItem (percentage of items to sample per iteration)
        pItem_val <- 0.8
        
        # Set maximum allowable pItem value (cannot exceed 1.0)
        max_pItem <- 1.0
        
        # Step size for incrementing pItem when NA values are found
        step <- 0.05
        
        
        # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
        while (has_na && pItem_val <= max_pItem) {
          # Run ConsensusClusterPlus with the current pItem_val
          consensus_results_app_three <- ConsensusClusterPlus(
            final_input_matrix, 
            maxK = 3, 
            reps = CCP_iter,            
            pItem = pItem_val,            # Percentage of items sampled per iteration
            pFeature = 1,                # Use all features
            clusterAlg = clustering_alg,      
            innerLinkage = "ward.D2",
            finalLinkage = "ward.D2",
            distance = "euclidean", 
            seed = 2024L,
            plot ="none",           
            verbose = FALSE,
            writeTable = FALSE
          )
          
          # Only extract out result for K = 2
          consensus_results_app_three <- consensus_results_app_three[1:2]
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results_app_three, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results_app_three, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          # Check if any values in the item_consensus column are NA
          if (any(is.na(item_consensus_df$item_consensus))) {
            # Print message indicating NA detected and that pItem will be increased
            message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
            
            # Increase pItem by the predefined step to sample more items per iteration
            pItem_val <- pItem_val + step
          } else {
            # If no NA found, set flag to FALSE to exit the loop
            has_na <- FALSE
          }
        }
        
        
      } else {
        
        
        # set that to working_k 
        working_max_k <- max_stable_k
        
        # Initialize starting value for pItem (percentage of items to sample per iteration)
        pItem_val <- 0.8
        
        # Set maximum allowable pItem value (cannot exceed 1.0)
        max_pItem <- 1.0
        
        # Step size for incrementing pItem when NA values are found
        step <- 0.05
        
        
        # Loop to run consensus clustering repeatedly, increasing pItem if NAs found
        while (has_na && pItem_val <= max_pItem) {
          # Run ConsensusClusterPlus with the current pItem_val
          consensus_results_app_three <- ConsensusClusterPlus(
            final_input_matrix, 
            maxK = working_max_k, 
            reps = CCP_iter,            
            pItem = pItem_val,            # Percentage of items sampled per iteration
            pFeature = 1,                # Use all features
            clusterAlg = clustering_alg,      
            innerLinkage = "ward.D2",
            finalLinkage = "ward.D2",
            distance = "euclidean", 
            seed = 2024L,
            plot ="none",           
            verbose = FALSE,
            writeTable = FALSE
          )
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results_app_three, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          
          # Calculate item consensus and cluster metrics (ICL)
          pdf(file = NULL)  # Start invisible plotting device
          icl_results <- calcICL(consensus_results_app_three, 
                                 plot = "png", 
                                 writeTable = TRUE)
          dev.off()          # Close device
          
          # Extract the itemConsensus dataframe from ICL results
          item_consensus_df <- icl_results[["itemConsensus"]]
          
          # Rename columns for clarity
          colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
          
          # Check if any values in the item_consensus column are NA
          if (any(is.na(item_consensus_df$item_consensus))) {
            # Print message indicating NA detected and that pItem will be increased
            message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
            
            # Increase pItem by the predefined step to sample more items per iteration
            pItem_val <- pItem_val + step
          } else {
            # If no NA found, set flag to FALSE to exit the loop
            has_na <- FALSE
          }
          
        }
        
      }
      
      
      
    }
    
    # If no NA found, proceed with downstream analysis steps
    if (!has_na) {
      # Classify stability based on item_consensus values using custom function
      item_consensus_df$stability <- sapply(item_consensus_df$item_consensus, get_stability_category)
      
      # Split the dataframe into a list of dataframes, one per cluster number k
      item_consensus_by_k <- split(item_consensus_df, item_consensus_df$k)
    }
    
    
    
    # Loop over each k from ConsensusClusterPlus results
    CCP_plot_list <- list() 
    
    
    for (k in 2:length(consensus_results_app_three)) {
      cm <- consensus_results_app_three[[k]]$consensusMatrix
      colnames(cm) <- colnames(final_input_matrix)
      rownames(cm) <- colnames(final_input_matrix)  # consensusMatrix is symmetric, matches samples
      
      cm <-  cm[consensus_results_app_three[[k]]$consensusTree$order,]
      
      
      clusters <- consensus_results_app_three[[k]]$consensusClass
      order_idx <- order(clusters)
      
      # Create annotation dataframe
      ann <- data.frame(Cluster = factor(clusters[order_idx]))
      #rownames(ann) <- names(clusters)
      
      # Step 1: Ensure Cluster is a factor with only the levels that exist
      ann$Cluster <- as.factor(ann$Cluster)
      ann$Cluster <- droplevels(ann$Cluster)  # drop any unused levels
      
      # Step 2: Define annotation colors correctly
      cluster_levels <- levels(ann$Cluster)
      
      # RColorBrewer requires at least 3 colors, so ensure enough levels
      num_clusters <- length(cluster_levels)
      palette <- RColorBrewer::brewer.pal(max(3, num_clusters), "Set2")[1:num_clusters]
      
      # Step 3: Set names to match factor levels exactly
      annotation_colors <- list(Cluster = setNames(palette, cluster_levels))
      
      
      # Save heatmap to list
      CCP_plot <- pheatmap(
        cm,
        cluster_rows = FALSE,
        cluster_cols = consensus_results_app_three[[k]]$consensusTree,
        annotation_col = ann,
        annotation_row = ann,
        show_colnames = FALSE,
        annotation_colors = annotation_colors,
        color = colorRampPalette(c("white", "blue"))(100),
        silent = TRUE,
        main = paste("Consensus Matrix for k =", k)
      )
      
      # Make name 
      CCP_plot_list_name <- paste0("CCP_k_", k)
      CCP_plot_list[[CCP_plot_list_name]] <- CCP_plot
    }
    
    # Calculate ICL and write tables for all k values - no dependency on optimal_k
    #it only saves Consensus-Cluster table
    pdf(file = NULL)  # Start invisible plotting device
    icl_results <- calcICL(consensus_results_app_three,
                           plot = "png",        # trick it into "thinking" it's plotting
                           writeTable = FALSE)
    dev.off()          # Close device
    
    
    
    ##Get itemConsensus score for each sample given a specific K##
    
    #Load the itemConsensus table from icl_results
    item_consensus_df <- as.data.frame(icl_results[["itemConsensus"]])
    
    # Rename columns
    colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")
    
    #function to classify stability#

    
    #stability classification
    item_consensus_df$stability <- sapply(item_consensus_df$item_consensus, get_stability_category)
    
    #Split into list of dataframes, one per k
    item_consensus_by_k <- split(item_consensus_df, item_consensus_df$k)
    
    
    # Create graphs for item consensus 
    
    item_consensus_overall <- as.data.frame(icl_results[["clusterConsensus"]])
    
    # If you want to round the labels nicely, e.g. 2 decimal places
    item_consensus_overall <- item_consensus_overall %>%
      mutate(label = round(clusterConsensus, 2))
    
    item_con_plot <- ggplot(item_consensus_overall,
                            aes(
                              x = factor(k),
                              y = clusterConsensus,
                              fill = factor(cluster)
                            )) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_text(aes(label = label), 
                position = position_dodge(width = 0.9),
                vjust = -0.3,  # slightly above the bar
                size = 6) +    # adjust text size
      labs(title = "Cluster Consensus by k", x = "k", y = "Consensus Score") +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        axis.text = element_text(colour = "black", size = 24),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.text.y = element_text(colour = "black", size = 24),
        axis.title = element_text(colour = "black", size = 30),
        axis.ticks = element_line(size = 1.5),
        panel.border = element_rect(
          colour = "black",
          fill = NA,
          size = 1
        )
      ) +
      scale_fill_manual(values = cols[1:working_max_k], name = "Cluster", drop = FALSE)
    
    
    
    item_consensus_by_k_plot_list <- list()
    for(k in names(item_consensus_by_k)){
      
      k_numeric <- as.numeric(k)
      item_con_dataset <- item_consensus_by_k[[as.character(k)]]
      
      
      item_consensus_by_k_plot <- ggplot(item_con_dataset,
                                         aes(
                                           x = sample_id,
                                           y = item_consensus,
                                           fill = factor(cluster)
                                         )) +
        geom_bar(stat = "identity") +
        labs(title = paste0("Item Consensus by Sample (k = ", k, ")"),
             x = "Sample ID",
             y = "Item Consensus",
             fill = "Cluster") +
        theme_classic() +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = 8
          ),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 14),
          legend.position = "right"
        ) +
        scale_fill_manual(values = cols[1:working_max_k], name = "Cluster", drop = FALSE)
      
      item_consensus_by_k_plot_list[[as.character(k)]] <- item_consensus_by_k_plot
      
    }
    
    
    # Specify custom results directory
    results_dir <- paste0(working_dir, "/CIM_states_results_", clustering_alg, "_", input_data_type,"_app_three")
    if(!dir.exists(results_dir)) {
      dir.create(results_dir, recursive = TRUE)
    }
    
    # Get Clustering gene set for this approach
    write.csv(clustering_gene_set,
              file = file.path(results_dir, paste0("clustering_gene_set.csv")),
              row.names = FALSE)
    
    
    
    # Save tables for each k into CSV files
    for (k_val in names(item_consensus_by_k)) {
      write.csv(item_consensus_by_k[[k_val]],
                file = file.path(results_dir, paste0("item_consensus_k", k_val, ".csv")),
                row.names = FALSE)
    }
    
    # Save Plots 
    
    ggsave(
      filename = paste0(results_dir, "/item_con_plot.png"),
      plot = item_con_plot,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    
    for (k in names(item_consensus_by_k_plot_list)){
      
      img <- item_consensus_by_k_plot_list[[as.character(k)]]
      
      ggsave(
        filename = paste0(results_dir, "/item_con_sample_k-", as.character(k), ".png"),
        plot = img,
        width = 10,
        height = 6,
        dpi = 300
      )
      
      
    }
    
    
    for (i in seq_along(CCP_plot_list)) {
      plot_name <- names(CCP_plot_list)[i]
      png_filename <- file.path(results_dir, paste0(plot_name, ".png"))
      png(png_filename, width = 1200, height = 1000, res = 150)
      grid.newpage()
      obj <- CCP_plot_list[[i]]
      if (inherits(obj, "grob")) {
        grid.draw(obj)
      } else if (inherits(obj, "Heatmap")) {
        draw(obj)
      } else if (!is.null(obj$gtable)) {
        grid.draw(obj$gtable)   # <-- fallback: draw the gtable element if it exists
      } else {
        warning(sprintf("Element %s is not drawable by grid.draw, draw, or print.", plot_name))
      }
      dev.off()
    }
    
    
    # Create cdf_plots
    k_vec = 2:working_max_k
    cdf_data <- data.frame()
    PACs <- numeric(length(k_vec))
    
    for(i in seq_along(k_vec)) {
      k <- k_vec[i]
      consensus_matrix <- consensus_results_app_three[[k]][["consensusMatrix"]]
      cons_vals <- consensus_matrix[lower.tri(consensus_matrix)]
      Fn <- ecdf(cons_vals)  # standard consensus cluster PAC
      PACs[i] <- Fn(0.9) - Fn(0.1)
      h <- hist(cons_vals, breaks = 200, plot = FALSE)
      cdf_y <- cumsum(h$counts) / sum(h$counts)
      cdf_data <- bind_rows(
        cdf_data,
        data.frame(
          k = as.factor(k),
          x = h$mids,
          y = cdf_y
        )
      )
    }
    
    # Now join the correct PAC values
    pac_df <- data.frame(k = as.factor(k_vec), PAC = PACs)
    cdf_data <- left_join(cdf_data, pac_df, by = "k")
    cdf_data$line_label <- paste0("K = ", cdf_data$k, " (PAC=", formatC(cdf_data$PAC, format = "f", digits = 3), ")")
    
    
    cdf_plot <- ggplot(cdf_data, aes(x = x, y = y, color = line_label, group = k)) +
      geom_line(size = 1.5) +
      labs(
        x = "Consensus Index",
        y = "Cumulative Fraction",
        title = "Consensus CDFs across Ks",
        color = "K (PAC)"
      ) +
      geom_vline(xintercept = c(0.1, 0.9), linetype = "dashed", color = "gray") +
      scale_x_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
      scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.98, 0.02), # x, y from 0 (left/bottom) to 1 (right/top)
        legend.justification = c("right", "bottom"),
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill="white", color="black", size=0.8),
        axis.text = element_text(colour = "black", size = 24),
        axis.text.x = element_text(colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title = element_text(colour = "black", size = 30),
        axis.ticks = element_line(size = 1.5),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)
      )
    
    
    # save cdf plot
    ggsave(
      filename = paste0(results_dir, "/final_cdf_plot.png"),
      plot = cdf_plot,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    
    # Calculate PAC metric (Proportion of Ambiguous Clustering) using ecdf
    Kvec = 2:working_max_k
    x1 = 0.1; x2 = 0.9 
    PAC = rep(NA, length(Kvec)) 
    names(PAC) = paste("K=", Kvec, sep="") 
    for(i in Kvec){
      M = consensus_results_app_three[[i]]$consensusMatrix
      Fn = ecdf(M[lower.tri(M)])
      PAC[i-1] = Fn(x2) - Fn(x1)
    }
    optK = Kvec[which.min(PAC)]
    
    # Convert PAC to data frame for plotting
    pac_values <- data.frame(k = Kvec, pac = unname(PAC))
    
    # Plot PAC values (lower is better) and display in R environment (Zaoqu-Liu/IRLS GitHub)
    pac_plot <- ggplot(pac_values, aes(factor(k), pac, group=1)) +
      geom_line() +
      theme_bw(base_rect_size = 1.5) +
      geom_point(size=4, shape=21, color='darkred', fill='orange') +
      #ggtitle('Proportion of ambiguous clustering') +
      xlab('Number of clusters K') + ylab("PAC Score") +
      theme_classic() + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top",
            axis.text = element_text(colour = "black", size = (24)),
            axis.text.x = element_text(colour = "black", size = (24)),
            axis.text.y = element_text(colour = "black", size = (24)),
            axis.title = element_text(colour = "black", size = (30)),
            axis.ticks = element_line(size = 1.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    
    
    # Create item consensus table
    avg_cluster_consensus <- item_consensus_overall %>%
      group_by(k) %>%
      summarize(avg_cluster_consensus = mean(clusterConsensus))
    
    # Silhouette calculations using consensus cluster final
    
    orig_space_input <- t(final_input_matrix)  # transpose to get dims as samples x features
    rownames(orig_space_input) <- rownames(clustering_matrix)
    
    # Distances
    dist_space_for_clustering <- dist(orig_space_input)
    dist_ge_space <- dist(clustering_matrix[,clustering_gene_set, drop=FALSE])
    
    # Collect silhouettes for all k
    sil_results <- lapply(2:length(consensus_results_app_three), function(k) {
      
      clusters <- consensus_results_app_three[[k]]$consensusClass
      
      # Silhouette in gene space
      sil_space_for_clustering   <- silhouette(clusters, dist_space_for_clustering)
      avg_space_for_clustering   <- mean(sil_space_for_clustering[, 3])
      
      # Silhouette in PaCMAP space
      sil_ge_space <- silhouette(clusters, dist_ge_space)
      avg_ge_space <- mean(sil_ge_space[, 3])
      
      data.frame(
        k = k,
        avg_sil_space_for_clustering = avg_space_for_clustering,
        avg_sil_ge_space = avg_ge_space
      )
    })
    
     # Combine
    sil_df <- do.call(rbind, sil_results)
    sil_df$avg_sil_combined <- ((sil_df$avg_sil_space_for_clustering +
        sil_df$avg_sil_ge_space) /
        2)

    
    # Plot comparison
    sil_plot_both <- ggplot(sil_df, aes(x = k)) +
      geom_line(aes(y = avg_sil_space_for_clustering, color = paste0(input_data_type))) +
      geom_point(aes(y = avg_sil_space_for_clustering, color = paste0(input_data_type))) +
      geom_line(aes(y = avg_sil_ge_space, color = "Gene space")) +
      geom_point(aes(y = avg_sil_ge_space, color = "Gene space")) +
      labs(y = "Average silhouette width", color = "Space", title = "Silhouette comparison across k") +
      theme_classic(base_size = 16) + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top",
            axis.text = element_text(colour = "black", size = 24),
            axis.title = element_text(colour = "black", size = 30),
            axis.ticks = element_line(size = 1.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    
    

    sil_plot_original_space <- ggplot(sil_df, aes(x = k)) +
      geom_line(aes(y = avg_sil_space_for_clustering, color = paste0(input_data_type))) +
      geom_point(aes(y = avg_sil_space_for_clustering, color = paste0(input_data_type))) +
      labs(y = "Average silhouette width", color = "Space", title = "Silhouette comparison across k") +
      theme_classic(base_size = 16) + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top",
            axis.text = element_text(colour = "black", size = 24),
            axis.title = element_text(colour = "black", size = 30),
            axis.ticks = element_line(size = 1.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    

    sil_plot_ge_space <- ggplot(sil_df, aes(x = k)) +
      geom_line(aes(y = avg_sil_ge_space, color = "Gene space")) +
      geom_point(aes(y = avg_sil_ge_space, color = "Gene space")) +
      labs(y = "Average silhouette width", color = "Space", title = "Silhouette comparison across k") +
      theme_classic(base_size = 16) + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "top",
            axis.text = element_text(colour = "black", size = 24),
            axis.title = element_text(colour = "black", size = 30),
            axis.ticks = element_line(size = 1.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    

  
    
    # save silhouette plot
    ggsave(
      filename = paste0(results_dir, "/sil_plot.png"),
      plot = sil_plot_both,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    
    # save silhouette plot
    ggsave(
      filename = paste0(results_dir, "/sil_plot_original_space.png"),
      plot = sil_plot_original_space,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    # save silhouette plot
    ggsave(
      filename = paste0(results_dir, "/sil_plot_ge_space.png"),
      plot = sil_plot_ge_space,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    
    # save pac plot
    ggsave(
      filename = paste0(results_dir, "/final_pac_plot.png"),
      plot = pac_plot,
      width = 10,
      height = 6,
      dpi = 300
    )
    

   

    
       
    # Get suggested optimal k values
    silhouette_dim_space_k <- sil_df$k[which.max(sil_df$avg_sil_space_for_clustering)]
    silhouette_ge_space_k <- sil_df$k[which.max(sil_df$avg_sil_ge_space)]
    silhouette_combined_avg_k <- sil_df$k[which.max(sil_df$avg_sil_combined)]
    pac_k <- Kvec[which.min(PAC)]  # Using optK from PAC calculation
    avg_cluster_consensus_k <- avg_cluster_consensus[which.max(avg_cluster_consensus$avg_cluster_consensus), ]$k
    
    # Print metric values
    cat("\n===== OPTIMAL CLUSTER SUGGESTIONS =====\n")
    cat("Silhouette analysis in dim space suggests k =", silhouette_dim_space_k, "\n")
    cat("Silhouette analysis in GE space suggests k =", silhouette_ge_space_k, "\n")
    cat("Silhouette analysis of combined avg (GE and red. dim) space suggests k =", silhouette_combined_avg_k, "\n")
    cat("PAC method suggests k =", pac_k, "\n")
    cat("cluster consensus method suggests k =", avg_cluster_consensus_k, "\n")
    
    
  # Create metrics summary dataframe
      metrics_summary <- data.frame(
        k = 2:working_max_k,
        silhouette_dim_reduce_space = NA, # renamed to indicate dim‑reduce (original) space
        silhouette_ge_space = NA, # silhouette computed in PacMAP (ge) space
        silhouette_combined_avg = NA, # silhouette computed in PacMAP (ge) space
        pac = NA,
        avg_cluster_consensus = NA
      )

      # Fill in the metrics
      for (i in 1:nrow(metrics_summary)) {
        k_val <- metrics_summary$k[i]

        # Silhouette in the original (dim‑reduce) space
        idx <- which(sil_df$k == k_val)
        if (length(idx) > 0) {
          metrics_summary$silhouette_dim_reduce_space[
            i
          ] <- sil_df$avg_sil_space_for_clustering[idx]

          metrics_summary$silhouette_ge_space[i] <- sil_df$avg_sil_ge_space[idx]
          metrics_summary$silhouette_combined_avg[i] <- sil_df$avg_sil_combined[
            idx
          ]
        }

        # PAC
        k_idx <- which(Kvec == k_val)
        if (length(k_idx) > 0) {
          metrics_summary$pac[i] <- PAC[k_idx]
        }

        # avg_cluster_consensus
        cc_idx <- which(avg_cluster_consensus$k == k_val)
        if (length(cc_idx) > 0) {
          metrics_summary$avg_cluster_consensus[
            i
          ] <- avg_cluster_consensus$avg_cluster_consensus[cc_idx]
        }
      }

      # Alternatively,  Optimal K can be selected using rank based compromise
      # The rank-based compromise approach assigns a rank to each k for both silhouette (descending) and PAC (ascending), then sums the ranks to identify the k with the best overall trade-off between high silhouette and low PAC. The k with the lowest combined rank is selected as the optimal choice.
      # https://dl.acm.org/doi/10.1145/371920.372165

      # Assign ranks (lowest rank = best)
      metrics_summary$rank_silhouette_dim_reduce_space <- rank(
        -metrics_summary$silhouette_dim_reduce_space,
        ties.method = "min"
      ) # higher silhouette = rank 1 (original space)
      metrics_summary$rank_silhouette_ge_space <- rank(
        -metrics_summary$silhouette_ge_space,
        ties.method = "min"
      ) # higher silhouette = rank 1 (ge space)

      metrics_summary$rank_silhouette_combined_avg <- rank(
        -metrics_summary$silhouette_combined_avg,
        ties.method = "min"
      ) # higher silhouette = rank 1 (averaged space)
      metrics_summary$rank_pac <- rank(metrics_summary$pac, ties.method = "min") # lower PAC = rank 1
      metrics_summary$rank_item_cluster_consensus <- rank(
        -metrics_summary$avg_cluster_consensus,
        ties.method = "min"
      ) # higher consensus = rank 1

      # Build the column names dynamically.
      rank_columns <- c(
        paste0("rank_", clustering_metrics) # e.g. "rank_pac", "rank_item_cluster_consensus"
      )

      # Compute the overall rank as the row sum across selected rank columns
      metrics_summary$overall_rank <- rowSums(
        metrics_summary[, rank_columns, drop = FALSE],
        na.rm = TRUE
      )

      # Find all rows with the minimum overall rank
      best_rows <- metrics_summary[
        metrics_summary$overall_rank == min(metrics_summary$overall_rank),
      ]

      # Use the last metric in clustering_metrics as the tiebreaker
      tie_metric <- clustering_metrics[length(clustering_metrics)]

      if (verbose) {
        cat(
          "Metrics being used for ranking:\n  ",
          paste(clustering_metrics, collapse = ", "),
          "\n"
        )
        cat(
          "Tie‑breaker metric (last in clustering_metrics vector):",
          tie_metric,
          "\n"
        )
      }

      # Determine if the tiebreaking metric should be maximized or minimized
      # (Assuming: max for 'silhouette' and 'item_cluster_consensus', min for 'pac')
      if (tie_metric == "silhouette_dim_reduce_space") {
        idx <- which.max(best_rows$silhouette_dim_reduce_space)
      } else if (tie_metric == "silhouette_ge") {
        idx <- which.max(best_rows$silhouette_ge_space)
      } else if (tie_metric == "item_cluster_consensus") {
        idx <- which.max(best_rows$avg_cluster_consensus)
      } else if (tie_metric == "pac") {
        idx <- which.min(best_rows$pac)
      } else {
        stop("Unknown tiebreaker metric.")
      }

      best_k_row <- best_rows[idx, ]
      best_k <- best_k_row$k


# add cluster metrics used for evaluation and tie
metrics_summary$metrics_used <- paste(clustering_metrics, collapse = ", ")
metrics_summary$metrics_used_for_tie <- tie_metric
    

    
    # Print metrics summary
    cat("\n===== CLUSTER METRICS SUMMARY =====\n")
    metrics_summary
    cat("\n")
    
    # Write metrics to file
    write.csv(metrics_summary, file.path(results_dir, "final_cluster_metrics_summary.csv"), row.names = FALSE)
    
    # Based on examining the plots and metrics, set your optimal k manually here:
    # CHANGE THIS VALUE after reviewing the plots:
    optimal_k <- best_k  # <-- MODIFY THIS BASED ON YOUR EXAMINATION IF NEEDED 
    
    cat("Using k =", optimal_k, "as the optimal number of clusters\n")
    
    # Extract cluster assignments for optimal k
    consensus_clusters <- consensus_results_app_three[[optimal_k]]$consensusClass
    
    
    
    
    ####### FINAL VISUALIZATION #######
    
    

    
    # Get stability/reliability metrics 
    optimal_stability_cluster_metrics <- item_consensus_by_k[[as.character(optimal_k)]]
    
    # Start from your original data frame
    wide_optimal_stability_cluster_metrics <- optimal_stability_cluster_metrics %>%
      # Only keep relevant columns
      dplyr::select(sample_id, cluster, item_consensus) %>%
      # Pivot wider: each cluster makes a new column
      pivot_wider(
        names_from = cluster,
        values_from = item_consensus,
        names_prefix = "item_consensus_cluster"
      )
    
    
    # Write metrics to file
    write.csv(optimal_stability_cluster_metrics,
              file.path(results_dir, "final_cluster_stability_metrics.csv"), row.names = FALSE)
    
    write.csv(wide_optimal_stability_cluster_metrics,
              file.path(results_dir, "final_wide_cluster_stability_metrics.csv"), row.names = FALSE)
    
    # Create data frame for UMAP visualization
    final_df_plot_app_three <- as.data.frame(umap_matrix_visual$layout) %>%
      rownames_to_column("sample_id") %>% 
      dplyr::rename(UMAP1 = V1, UMAP2 = V2) %>%
      dplyr::inner_join(as.data.frame(pacmap_delta_embeddings_visual) %>% 
                          rownames_to_column("sample_id") %>% 
                          dplyr::rename(PACMAP1 = V1, PACMAP2 = V2)) %>%
      dplyr::inner_join(as.data.frame(pca_scores) %>% 
                          rownames_to_column("sample_id")) %>%
      dplyr::mutate(base_id = gsub("^delta_([^_]+)_.*", "\\1", sample_id)) %>%
      dplyr::mutate(cluster_assignments = consensus_clusters) %>%
      dplyr::inner_join(as.data.frame(clustering_matrix) %>%
                          rownames_to_column("sample_id"),
                        by = "sample_id")
    
    
    # Your input matrix: samples x genes
    pca_fviz_input <- final_df_plot_app_three[, c("PC1", "PC2", "sample_id")] %>%
      tibble::column_to_rownames("sample_id")
    
    
    
    fviz_cluster_plot <- fviz_cluster(
      object = list(data = pca_fviz_input, cluster = as.integer(consensus_clusters)),
      palette = cols[1:optimal_k],
      geom = "text",
      ellipse.type = "convex",
      ggtheme = theme_bw(),
      ggtitle("PCA plot on Clustered Data")
    )
    
    fviz_cluster_plot
    ggsave(file.path(results_dir, paste0("final_fviz_2D_PCA_k-", optimal_k, ".png")), 
           fviz_cluster_plot, width = 10, height = 8, dpi = 300)
    
    
    # 2. PCA projection using selected PCs
    pca_scores_3D <- pca_result$x[, 1:3, drop = FALSE]
    
    # 3. Prepare input for downstream analysis, e.g., clustering
    # Here, samples x PCs (rows = samples, cols = PCs)
    pca_fviz_input_3D <- as.matrix(pca_scores_3D)
    rownames(pca_fviz_input_3D) <- final_df_plot_app_three$sample_id  # optional, set sample IDs
    
    # If pca_fviz_input_3D is samples x PCs, select first 3 PCs
    
    plot_df_3D_PCA <- data.frame(
      PC1 = pca_fviz_input_3D[, 1],
      PC2 = pca_fviz_input_3D[, 2],
      PC3 = pca_fviz_input_3D[, 3],
      cluster = as.factor(as.integer(consensus_clusters)),
      sample_id = rownames(pca_fviz_input_3D)
    )
    
    
    
    
    
    # 3D scatter plot
    fig_PCA_3D <- plot_ly(
      data = plot_df_3D_PCA,
      x = ~PC1, y = ~PC2, z = ~PC3,
      color = ~cluster,
      text = ~sample_id,
      type = "scatter3d",
      mode = "markers",
      colors = cols
    ) %>%
      plotly::layout(
        title = "3D PCA Cluster Visualization",
        scene = list(
          xaxis = list(title = "PC1"),
          yaxis = list(title = "PC2"),
          zaxis = list(title = "PC3")
        )
      )
    
    fig_PCA_3D
    # Load library if not already loaded
    library(htmlwidgets)
    
    # Save as an interactive HTML file
    saveWidget(
      widget = fig_PCA_3D,
      file = paste0(results_dir, "/final_3D_PCA_k-", optimal_k, ".html"),
      selfcontained = TRUE
    )
    
    
    # Use letters as shapes 
    shapes <- unlist(lapply(letters, utf8ToInt))
    
    # Plot using GGPlot Color by cluster
    final_cluster_plot_umap <- ggplot(
      final_df_plot_app_three,
      aes(
        x = UMAP1,
        y = UMAP2,
        label = base_id ,
        color = as.factor(cluster_assignments) 
      )
    ) +
      scale_color_manual(values = cols, name = "Cluster") +
      ggforce::geom_mark_hull(
        aes(group = as.factor(cluster_assignments),
            fill = as.factor(cluster_assignments), label = NULL),
        concavity = 2, expand = unit(2, "mm"), alpha = 0.2
      ) +
      geom_label_repel(size = 4, fontface = "bold", label.size = 0.3, fill = "white", 
                       alpha = 1, box.padding = 0.25, max.overlaps = 100) +
      labs(x = "UMAP 1", y = "UMAP 2",
           fill = "Cluster",
           color = "Cluster",
           title = paste0(
             "Final Clustering Using ", input_data_type, " input and ",
             clustering_alg, " CCP Alg (k = ", optimal_k, ")"
           )) +
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
    
    
    ggsave(file.path(results_dir, paste0("final_clustering_shown_on_umap_k-", optimal_k, ".png")), 
           final_cluster_plot_umap, width = 10, height = 8, dpi = 300)
    
    # Plot using GGPlot Color by cluster
    final_cluster_plot_pacmap <- ggplot(
      final_df_plot_app_three,
      aes(
        x = PACMAP1,
        y = PACMAP2,
        label = base_id,  # Your patient/sample label column
        color = as.factor(cluster_assignments)
      )
    ) +
      scale_color_manual(values = cols, name = "Cluster") +
      ggforce::geom_mark_hull(
        aes(group = as.factor(cluster_assignments),
            fill = as.factor(cluster_assignments), label = NULL),
        concavity = 2, expand = unit(2, "mm"), alpha = 0.2
      ) +
      geom_label_repel(size = 4, fontface = "bold", label.size = 0.3, fill = "white", 
                       alpha = 1, box.padding = 0.25, max.overlaps = 100) +
      labs(x = "PacMAP 1", y = "PacMAP 2",
           fill = "Cluster",
           color = "Cluster",
           title = paste0(
             "Final Clustering Using ", input_data_type, " input and ",
             clustering_alg, " CCP Alg (k = ", optimal_k, ")"
           )) +
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
    
    
    
    ggsave(file.path(results_dir, paste0("final_clustering_shown_on_pacmap_k-", optimal_k, ".png")), 
           final_cluster_plot_pacmap, width = 10, height = 8, dpi = 300)
    
  }
  
  }
  
  if (no_app_three == TRUE) {
    cat(
      "\033[1;31mWARNING: Final Gene Set for Approach three empty! Consider Relaxing P-val
        Approach output will not be produced\033[0m\n"
    )
    
    # set string
    no_app_three_string <- "Approach three final gene set empty; Consider Relaxing P-val"
    
    # update iter_log
    iter_log <- rbind(
      iter_log,
      data.frame(
        gene_set = "app_three final get set", iteration = NA,
        before = NA,
        after = no_app_three_string
        )
    )
    
    # set df to null
    final_df_plot_app_three <- NULL
  }
  
  # Extra capturing of output 
  {
    # Write initial matrix to file
    write.csv(clustering_matrix,
              file.path(results_dir, "initial_clustering_mat.csv"), row.names = TRUE)
    
    # Write iter_log  to file
    write.csv(iter_log,
              file.path(results_dir, "iter_log_all_approaches.csv"), row.names = FALSE)
    
    # Convert gene_sets to to data frame
    df_gene_sets_for_file <- data.frame(
      list_name = names(all_gene_sets),
      items = sapply(all_gene_sets, function(x) paste(x, collapse = ","))
    )
    
    # Write gene_sets to CSV 
    write.csv(df_gene_sets_for_file,  
              file.path(results_dir, "all_gene_sets.csv"), row.names = FALSE)
    
    # write final_clustered df to file
    write.csv(final_df_plot_app_three,  
              file.path(results_dir, "clustered_samples_app_three.csv"), row.names = FALSE)
    
    # write pacmap settings 
    write.csv(pm_settings,  
              file.path(results_dir, "pacmap_settings.csv"), row.names = TRUE)
    
     # write removed genes setting 
    write.csv(removed_genes,  
              file.path(results_dir, "removed_genes.csv"), row.names = TRUE)
  }
 
  })
    
  }
  
  
  # ---- Return result list ----
  return(list(
    iterated_by_gene_sets = if ("app_one" %in% filter_approach) new_final else "NA",
    iterated_over_all_genes = if ("app_two" %in% filter_approach) iterated_over_all_genes else "NA",
    iterated_use_final_gene_set_from_approach_one = if ("app_three" %in% filter_approach) iterated_use_final_gene_set_from_approach_one else "NA", 
    details_per_set_analysis_by_gene_sets = if ("app_one" %in% filter_approach) details_per_set_analysis_by_gene_sets else "NA",
    details_per_set_analysis_by_all_genes = if ("app_two" %in% filter_approach) details_per_set_analysis_by_all_genes else "NA" ,
    details_per_set_analysis_by_gene_set_then_all_genes = if ("app_three" %in% filter_approach) details_per_set_analysis_by_gene_set_then_all_genes else "NA",
    iter_log = iter_log,
    pacmap_settings = pm_defaults,
    removed_genes = removed_genes,
    final_df_app_one = if ("app_one" %in% filter_approach) final_df_plot_app_one else "NA",
    final_df_app_two = if ("app_two" %in% filter_approach) final_df_plot_app_two else "NA",
    final_df_app_three = if ("app_three" %in% filter_approach) final_df_plot_app_three else "NA"
  ))
  
  
}
