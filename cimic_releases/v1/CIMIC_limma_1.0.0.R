#' CIMIC pipeline (LIMMA variant of CIMIC_streamlined.R)
#' -------------------------------------------------------------------------
#' Identical to CIMIC_streamlined.R EXCEPT for the feature-selection SWAP POINT:
#' features are ranked by a single moderated linear model fit across ALL genes
#' at once (limma) rather than gene-by-gene Wilcoxon/Kruskal.
#'
#'   design = model.matrix(~ cluster)  (cluster = consensus assignment)
#'   fit    = eBayes(lmFit(t(expr), design), robust = TRUE, trend = TRUE)
#'   2 groups   -> moderated t  : p = P.Value,  effect_size = |t|,  adj = BH adj.P.Val
#'   >=3 groups -> moderated F  : p = P.Value,  effect_size = F,    adj = BH adj.P.Val
#'
#' limma's adj.P.Val (BH) IS the FDR used by the shared downstream filter
#' (p<=0.05 & adj<=thresh). Everything else (convergence, clustering, plotting)
#' is unchanged and shared with the streamlined version.
#'
#' NOTE ON OUTPUT (verified): topTable(coef = cluster terms, sort.by="none")
#' returns one row per gene with columns t/logFC (2-grp) or F (multi-grp) plus
#' P.Value and adj.P.Val; rownames are the gene ids, order preserved.
#'
#' WHAT CHANGED vs CIMIC_Release_1.0.0.R
#'   * The original file was ~7,200 lines built from SIX near-identical copies
#'     of one routine (consensus-cluster -> metrics -> rank/tiebreak ->
#'     per-gene test -> converge -> save/plot), differing only in the seed gene
#'     universe and object/label names. Those copies are now collapsed into a
#'     handful of documented helper functions. ALL logic and approach retained.
#'   * The per-gene statistical test is isolated in ONE function,
#'     `cimic_select_features()` (see "SWAP POINT"), so alternative test
#'     strategies (adaptive / limma) are drop-in replacements.
#'   * Documented bug fixes carry an inline "# BUGFIX:" tag (see list at bottom).
#'
#' Approaches (unchanged meaning):
#'   app_one   -> refine each gene set independently, then UNION survivors
#'   app_two   -> refine the pooled union of all sets at once
#'   app_three -> refine app_one's union again over that combined set
#' -------------------------------------------------------------------------
#' 
#' 
#'  NOTE: CCP/UMAP seeds intentionally fixed; only PaCMAP stochasticity is evaluated, see ARI robustness analysis
#' in original Oncoimmunology manuscript.
#' 
#' 
#' 




#### Dependency Management ####
needed_pkgs <- c(
  "umap", "cluster", "factoextra", "dplyr", "stringr", "magrittr", "tibble",
  "rstatix", "coin", "RColorBrewer", "ggforce", "grid", "gridExtra", "corrplot",
  "concaveman", "paran", "reticulate", "utils", "plotly", "scales", "withr",
  "ggrepel", "ggplot2", "tidyr", "msigdbr", "jsonlite", "pheatmap", "htmlwidgets"
)
bioc_pkgs <- c("ConsensusClusterPlus", "ComplexHeatmap", "limma")

new_pkgs <- needed_pkgs[!(needed_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)
new_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[, "Package"])]
if (length(new_bioc)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(new_bioc, ask = FALSE)
}
suppressPackageStartupMessages(lapply(c(needed_pkgs, bioc_pkgs), require, character.only = TRUE))


# =========================================================================
# Helper functions
# =========================================================================

#' Classify a consensus/stability score into a qualitative AHP category.
#' (ConsensusClusterPlus has no predefined threshold; scale from
#' DOI: 10.13140/RG.2.2.28745.83045)
get_stability_category <- function(score) {
  if (is.na(score)) return(NA)
  else if (score < 0.50) return("Very Low (<50%)")
  else if (score < 0.625) return("Low (50-62.5%)")
  else if (score < 0.75) return("Moderate (62.5-75%)")
  else if (score < 0.875) return("High (75-87.5%)")
  else return("Very High (>=87.5%)")
}


# ---- PaCMAP / conda setup (unchanged behaviour) ----
if (!requireNamespace("reticulate", quietly = TRUE)) install.packages("reticulate")
library(reticulate)

setup_conda <- function(env_name = "pacmap_env") {
  if (exists("pacmap")) {
    message("pacmap already present; skipping conda setup.")
    return(import("pacmap"))
  }
  message("Trying to set up PaCMAP...")

  cache_file <- file.path(Sys.getenv("HOME"), ".local", "share", "r-reticulate", "conda_list.json")
  if (file.exists(cache_file)) file.remove(cache_file)

  cp <- Sys.getenv("RETICULATE_CONDA", unset = NA)
  if (!(!is.na(cp) && file.exists(cp))) {
    cp <- tryCatch(reticulate::conda_binary(), error = function(e) NULL)
  }

  if (is.null(cp) || !file.exists(cp)) {
    max_tries <- 3
    install_fallback <- FALSE
    for (i in seq_len(max_tries)) {
      cat("\n=== Conda not detected ===\n")
      cat("Please type the **full path** to your conda executable (WITH NO QUOTES), e.g.\n")
      cat("  C:/Users/you/AppData/Local/anaconda3/Scripts/conda.exe\n")
      cat("  if copy/paste for path isn't working try ctrl+shift+v (Windows users)\n")
      cat("Or type **5** to let the script install Miniconda automatically.\n")
      message("\n[Info] To avoid this prompt, set RETICULATE_CONDA once, e.g.:")
      message('Sys.setenv(RETICULATE_CONDA = "C:/Users/user/AppData/Local/anaconda3/Scripts/conda.exe")')
      user_input <- readline(prompt = "conda path (or 5): ")
      if (trimws(user_input) == "5" || nzchar(user_input) == FALSE) {
        install_fallback <- TRUE
        break
      }
      if (file.exists(user_input)) {
        cp <- normalizePath(user_input, winslash = "/", mustWork = TRUE)
        Sys.setenv(RETICULATE_CONDA = cp)
        message("User provided file exists and contains valid conda...")
        break
      } else {
        message("The file you entered does not exist: ", user_input)
        if (i < max_tries) message("Please try again (attempt ", i + 1, " of ", max_tries, ").")
      }
    }
    if (install_fallback) {
      message("No valid Conda path provided - installing Miniconda locally ...")
      install_miniconda()
      cp <- reticulate::conda_binary()
      if (is.null(cp) || !file.exists(cp)) stop("Failed to install or locate Miniconda. Install it manually.")
      message("Miniconda installed at: ", cp)
    }
  }

  envs <- conda_list(conda = cp)
  if (!(env_name %in% envs$name)) {
    message(sprintf("Creating conda environment '%s' and installing pacmap.", env_name))
    conda_create(envname = env_name, conda = cp)
    conda_install(envname = env_name, packages = "pacmap", pip = TRUE, conda = cp)
  } else {
    out <- system2(cp, c("list", "-n", env_name, "--json"), stdout = TRUE)
    pkgs <- jsonlite::fromJSON(paste(out, collapse = "\n"))$name
    if (!any(pkgs == "pacmap")) {
      message("Installing pacmap into existing environment ...")
      conda_install(envname = env_name, packages = "pacmap", pip = TRUE, conda = cp)
    }
    message("Pacmap set up.")
  }

  use_condaenv(env_name, conda = cp, required = TRUE)
  import("pacmap")
}


# -------------------------------------------------------------------------
# Consensus clustering with adaptive pItem (collapses ~380 duplicated lines)
# -------------------------------------------------------------------------
# Runs ConsensusClusterPlus, bumping pItem by `step` until item-consensus has
# no NA (or pItem exceeds 1). `extract_k2 = TRUE` reproduces the original
# "run maxK = 3 but keep only k = 2" special case. calcICL is called ONCE
# (BUGFIX #4: the original called it twice back-to-back with identical args).
.ccp_pItem_loop <- function(input_matrix, maxK, extract_k2, CCP_iter, clustering_alg) {
  pItem_val <- 0.8; max_pItem <- 1.0; step <- 0.05; has_na <- TRUE
  consensus_results <- NULL; icl_results <- NULL; item_consensus_df <- NULL
  while (has_na && pItem_val <= max_pItem) {
    consensus_results <- ConsensusClusterPlus(
      input_matrix, maxK = maxK, reps = CCP_iter, pItem = pItem_val,
      pFeature = 1, clusterAlg = clustering_alg, innerLinkage = "ward.D2",
      finalLinkage = "ward.D2", distance = "euclidean", seed = 2024L,
      plot = "none", verbose = FALSE, writeTable = FALSE
    )
    if (extract_k2) consensus_results <- consensus_results[1:2]

    # IMPORTANT: calcICL ALWAYS draws (its par()/plot loop are unconditional); the `plot`
    # arg only chooses the TARGET device. With plot=NULL it opens no device and draws to
    # the CURRENT one -> in interactive R (RStudio/VSCode) that's the plot viewer, which
    # accumulates a plot every consensus run and crashes the session. The old
    # plot="png"+writeTable=TRUE instead leaked a PNG device each call and littered the
    # cwd. We only need the itemConsensus/clusterConsensus tables, so we redirect calcICL's
    # drawing into a null device (no viewer plots, no files, no leak).
    grDevices::pdf(file = nullfile())
    icl_results <- tryCatch(
      calcICL(consensus_results, plot = NULL, writeTable = FALSE),
      finally = grDevices::dev.off()
    )

    item_consensus_df <- icl_results[["itemConsensus"]]
    colnames(item_consensus_df) <- c("k", "cluster", "sample_id", "item_consensus")

    if (any(is.na(item_consensus_df$item_consensus))) {
      message(sprintf("NA found in item_consensus at pItem=%.2f, increasing pItem...", pItem_val))
      pItem_val <- pItem_val + step
    } else {
      has_na <- FALSE
    }
  }
  list(consensus_results = consensus_results, icl_results = icl_results,
       item_consensus_df = item_consensus_df, has_na = has_na)
}

#' Full adaptive consensus routine. Returns the CCP results plus the possibly
#' reduced `working_max_k`, and an `outcome`:
#'   "ok"   -> usable clustering produced
#'   "skip" -> no stable K at any pItem (caller decides what to do)
.consensus_adaptive <- function(input_matrix, working_max_k, CCP_iter, clustering_alg) {
  # First attempt (k = 2 special case runs maxK = 3 and keeps k = 2)
  if (working_max_k == 2) {
    res <- .ccp_pItem_loop(input_matrix, maxK = 3, extract_k2 = TRUE, CCP_iter, clustering_alg)
  } else {
    res <- .ccp_pItem_loop(input_matrix, maxK = working_max_k, extract_k2 = FALSE, CCP_iter, clustering_alg)
  }

  if (res$has_na) {
    message("Reached max pItem but NA values still exist in item_consensus.")
    message("Reducing working_max_k to the maximum most stable K")
    ks_without_na <- res$item_consensus_df %>%
      group_by(k) %>% summarize(no_na = all(!is.na(item_consensus))) %>%
      filter(no_na) %>% pull(k)

    if (length(ks_without_na) == 0) {
      return(list(outcome = "skip", working_max_k = working_max_k))
    }
    max_stable_k <- max(ks_without_na)
    message("maximum most stable K = ", max_stable_k)

    if (max_stable_k == 2) {
      message("since maximum most stable K = 2 will extract k = 2 results using k = 3 analysis")
      working_max_k <- 2
      res <- .ccp_pItem_loop(input_matrix, maxK = 3, extract_k2 = TRUE, CCP_iter, clustering_alg)
    } else {
      working_max_k <- max_stable_k
      res <- .ccp_pItem_loop(input_matrix, maxK = working_max_k, extract_k2 = FALSE, CCP_iter, clustering_alg)
    }
  }

  item_consensus_df <- res$item_consensus_df
  if (!res$has_na) {
    item_consensus_df$stability <- sapply(item_consensus_df$item_consensus, get_stability_category)
  }
  list(outcome = "ok", working_max_k = working_max_k,
       consensus_results = res$consensus_results, icl_results = res$icl_results,
       item_consensus_df = item_consensus_df)
}


# -------------------------------------------------------------------------
# PAC + silhouette + rank-based k selection (collapses the metrics block)
# -------------------------------------------------------------------------
#' @param samples_by_dim  dim-reduce embedding, samples x dims (for dim-space silhouette)
#' @param gene_set        genes used (for gene-expression-space silhouette)
.compute_k_metrics <- function(consensus_results, icl_results, item_consensus_df,
                               working_max_k, samples_by_dim, clustering_matrix,
                               gene_set, clustering_metrics, verbose) {
  # ---- PAC ----
  Kvec <- 2:working_max_k
  PAC <- rep(NA, length(Kvec)); names(PAC) <- paste0("K=", Kvec)
  for (i in Kvec) {
    M <- consensus_results[[i]]$consensusMatrix
    Fn <- ecdf(M[lower.tri(M)])
    PAC[i - 1] <- Fn(0.9) - Fn(0.1)
  }
  pac_values <- data.frame(k = Kvec, pac = unname(PAC))

  # ---- Silhouette (dim-reduce space and gene-expression space) ----
  dist_space_for_clustering <- dist(samples_by_dim)
  dist_ge_space <- dist(clustering_matrix[, gene_set, drop = FALSE])
  sil_results <- lapply(2:length(consensus_results), function(k) {
    clusters <- consensus_results[[k]]$consensusClass
    data.frame(
      k = k,
      avg_sil_space_for_clustering = mean(silhouette(clusters, dist_space_for_clustering)[, 3]),
      avg_sil_ge_space = mean(silhouette(clusters, dist_ge_space)[, 3])
    )
  })
  sil_df <- do.call(rbind, sil_results)
  sil_df$avg_sil_combined <- (sil_df$avg_sil_space_for_clustering + sil_df$avg_sil_ge_space) / 2

  # ---- Cluster consensus ----
  cluster_consensus_df <- as.data.frame(icl_results[["clusterConsensus"]])
  avg_cluster_consensus <- cluster_consensus_df %>%
    group_by(k) %>% summarize(avg_cluster_consensus = mean(clusterConsensus))

  # ---- Metrics summary ----
  metrics_summary <- data.frame(
    k = 2:working_max_k, silhouette_dim_reduce_space = NA, silhouette_ge_space = NA,
    silhouette_combined_avg = NA, pac = NA, avg_cluster_consensus = NA
  )
  for (i in seq_len(nrow(metrics_summary))) {
    k_val <- metrics_summary$k[i]
    idx <- which(sil_df$k == k_val)
    if (length(idx) > 0) {
      metrics_summary$silhouette_dim_reduce_space[i] <- sil_df$avg_sil_space_for_clustering[idx]
      metrics_summary$silhouette_ge_space[i] <- sil_df$avg_sil_ge_space[idx]
      metrics_summary$silhouette_combined_avg[i] <- sil_df$avg_sil_combined[idx]
    }
    k_idx <- which(Kvec == k_val)
    if (length(k_idx) > 0) metrics_summary$pac[i] <- PAC[k_idx]
    cc_idx <- which(avg_cluster_consensus$k == k_val)
    if (length(cc_idx) > 0) metrics_summary$avg_cluster_consensus[i] <- avg_cluster_consensus$avg_cluster_consensus[cc_idx]
  }

  # ---- Rank-based compromise (lowest summed rank wins) ----
  metrics_summary$rank_silhouette_dim_reduce_space <- rank(-metrics_summary$silhouette_dim_reduce_space, ties.method = "min")
  metrics_summary$rank_silhouette_ge_space <- rank(-metrics_summary$silhouette_ge_space, ties.method = "min")
  metrics_summary$rank_silhouette_combined_avg <- rank(-metrics_summary$silhouette_combined_avg, ties.method = "min")
  metrics_summary$rank_pac <- rank(metrics_summary$pac, ties.method = "min")
  metrics_summary$rank_item_cluster_consensus <- rank(-metrics_summary$avg_cluster_consensus, ties.method = "min")

  rank_columns <- paste0("rank_", clustering_metrics)
  metrics_summary$overall_rank <- rowSums(metrics_summary[, rank_columns, drop = FALSE], na.rm = TRUE)
  best_rows <- metrics_summary[metrics_summary$overall_rank == min(metrics_summary$overall_rank), ]
  tie_metric <- clustering_metrics[length(clustering_metrics)]

  if (verbose) {
    cat("Metrics being used for ranking:\n  ", paste(clustering_metrics, collapse = ", "), "\n")
    cat("Tie-breaker metric (last in clustering_metrics vector):", tie_metric, "\n")
  }

  # BUGFIX #6: tie_metric "silhouette_ge" -> "silhouette_ge_space" (the valid name).
  idx <- if (tie_metric == "silhouette_dim_reduce_space") which.max(best_rows$silhouette_dim_reduce_space)
    else if (tie_metric %in% c("silhouette_ge_space", "silhouette_ge")) which.max(best_rows$silhouette_ge_space)
    else if (tie_metric == "silhouette_combined_avg") which.max(best_rows$silhouette_combined_avg)
    else if (tie_metric == "item_cluster_consensus") which.max(best_rows$avg_cluster_consensus)
    else if (tie_metric == "pac") which.min(best_rows$pac)
    else stop("Unknown tiebreaker metric.")

  best_k <- best_rows[idx, ]$k
  metrics_summary$metrics_used <- paste(clustering_metrics, collapse = ", ")
  metrics_summary$metrics_used_for_tie <- tie_metric

  # ---- Stability metrics for the chosen k ----
  item_consensus_by_k <- split(item_consensus_df, item_consensus_df$k)
  optimal_stability <- item_consensus_by_k[[as.character(best_k)]]
  wide_optimal_stability <- optimal_stability %>%
    dplyr::select(sample_id, cluster, item_consensus) %>%
    pivot_wider(names_from = cluster, values_from = item_consensus, names_prefix = "item_consensus_cluster")

  list(
    optimal_k = best_k,
    consensus_clusters = consensus_results[[best_k]]$consensusClass,
    metrics_summary = metrics_summary, sil_df = sil_df, pac_values = pac_values,
    PAC = PAC, Kvec = Kvec, avg_cluster_consensus = avg_cluster_consensus,
    item_consensus_by_k = item_consensus_by_k,
    optimal_stability = optimal_stability, wide_optimal_stability = wide_optimal_stability
  )
}


# =========================================================================
# SWAP POINT: per-gene feature test
# =========================================================================
#' LIMMA CIMIC feature test: one moderated linear model across ALL genes per
#' call. 2 groups -> moderated t (effect = |t|); >=3 groups -> moderated F
#' (effect = F). `adj_p_value` is limma's BH FDR. Returns the FULL per-gene
#' stats table (unfiltered); the shared caller applies p<=0.05 & adj<=thresh.
#'
#' @param expr_mat  samples x genes numeric matrix (rows aligned to `clusters`)
#' @param clusters  cluster assignment per sample (aligned to rows of expr_mat)
#' @return data.frame(gene_id, effect_size, p_value, adj_p_value, test)
cimic_select_features <- function(expr_mat, clusters, adj_pval_thresh) {
  genes <- colnames(expr_mat)
  cluster <- as.factor(clusters)
  na_df <- data.frame(gene_id = genes, effect_size = NA_real_, p_value = NA_real_,
                      adj_p_value = NA_real_, test = "NA", stringsAsFactors = FALSE)
  if (nlevels(cluster) < 2) return(na_df)

  tryCatch({
    design <- model.matrix(~ cluster)
    # robust = TRUE: robust empirical-Bayes prior (Phipson et al. 2016) so a few
    # hypervariable genes don't distort the shared variance prior.
    # trend  = TRUE: allow an intensity-dependent (mean-variance) prior trend.
    # NOTE: on these delta data both barely change selection vs plain eBayes and
    # do NOT make limma track Wilcoxon (use CIMIC_limma_ranked.R for that).
    # Included for completeness / best-practice standard-limma; set either to
    # FALSE to disable.
    fit <- limma::eBayes(limma::lmFit(t(expr_mat), design), robust = TRUE, trend = TRUE)   # limma wants genes x samples
    coefs <- grep("^cluster", colnames(design), value = TRUE)
    tt <- limma::topTable(fit, coef = coefs, number = Inf, sort.by = "none")
    tt <- tt[genes, , drop = FALSE]                           # enforce original gene order
    if (nlevels(cluster) == 2) {
      data.frame(gene_id = genes, effect_size = abs(tt$t), p_value = tt$P.Value,
                 adj_p_value = tt$adj.P.Val, test = "limma_2grp(mod_t)", stringsAsFactors = FALSE)
    } else {
      data.frame(gene_id = genes, effect_size = tt$F, p_value = tt$P.Value,
                 adj_p_value = tt$adj.P.Val, test = "limma_multi(mod_F)", stringsAsFactors = FALSE)
    }
  }, error = function(e) { message("limma feature test failed: ", conditionMessage(e)); na_df })
}


# =========================================================================
# One refinement run (the shared while-loop) for a single seed gene set
# =========================================================================
#' @param universe_size  size of the ORIGINAL gene universe for this run; used
#'   by the non-immediate-convergence-to-zero revert check. (BUGFIX #1: app_one
#'   originally passed `length(gene_set_i)` = 1, so its revert branch never fired.)
#' @param select_fn       feature-test function (see SWAP POINT).
#' @return list(final, details, iter_log, converged_empty)
.refine_gene_set <- function(seed_genes, label, universe_size, clustering_matrix, pacmap,
                             pm_defaults, max_k, CCP_iter, clustering_alg, clustering_metrics,
                             adj_pval_thresh, max_pipeline_iter, select_fn, verbose) {
  working_max_k <- max_k
  iterated_final_gene_set <- seed_genes
  gene_set_check_1 <- 0; gene_set_check_2 <- length(iterated_final_gene_set); n <- 0

  stats_per_iter <- list(); cluster_metrics_per_iter <- list(); optimal_k_per_iter <- list()
  cluster_stability_per_iter <- list(); wide_cluster_stability_per_iter <- list()
  iter_log <- data.frame(gene_set = character(), iteration = integer(),
                         before = integer(), after = character(), stringsAsFactors = FALSE)

  save_details <- function(final = iterated_final_gene_set) {
    list(final = final, feature_importance_per_iter = stats_per_iter,
         optimal_k_per_iter = optimal_k_per_iter, cluster_metrics_per_iter = cluster_metrics_per_iter,
         cluster_stability_per_iter = cluster_stability_per_iter,
         wide_cluster_stability_per_iter = wide_cluster_stability_per_iter)
  }
  add_log <- function(before, after) {
    iter_log[nrow(iter_log) + 1, ] <<- list(label, n, as.integer(before), as.character(after))
  }

  while (gene_set_check_1 != gene_set_check_2 && n < max_pipeline_iter) {

     working_max_k <- max_k # Reset working_max_k to max_k at the start of each iteration to avoid it being reduced in previous iterations.

    if (verbose) message(sprintf("%s Iteration %d Start: %d genes (%d previously)",
                                 label, n, gene_set_check_2, gene_set_check_1))

    # --- PaCMAP embed (guarded) ---
    script_pacmap_reducer <- do.call(pacmap$PaCMAP, pm_defaults)
    safe_pacmap <- tryCatch(
      script_pacmap_reducer$fit_transform(clustering_matrix[, iterated_final_gene_set, drop = FALSE]),
      error = function(e) {
        message(sprintf("PaCMAP failed for %s iteration %d. Keeping most iterated gene set (if iter > 0).", label, n))
        NULL
      }
    )
    if (is.null(safe_pacmap)) {
      return(list(final = iterated_final_gene_set, details = save_details(),
                  iter_log = iter_log, converged_empty = FALSE))
    }
    rownames(safe_pacmap) <- rownames(clustering_matrix)

    input_matrix <- t(safe_pacmap)
    colnames(input_matrix) <- rownames(clustering_matrix)

    # --- Consensus clustering (adaptive pItem) ---
    cc <- .consensus_adaptive(input_matrix, working_max_k, CCP_iter, clustering_alg)
    if (cc$outcome == "skip") {
      message("No stable K values found. Keeping most iterated gene set (if iter > 0).")
      if (n >= 1) {
        add_log(length(iterated_final_gene_set),
                paste0(length(iterated_final_gene_set),
                       "; No stable clustering produced at higher iterations. Took most iterated geneset."))
      } else {
        add_log(length(iterated_final_gene_set),
                "0 ; No stable clustering produced at first iteration. Skipped geneset.")
      }
      break
    }
    working_max_k <- cc$working_max_k

    # --- Metrics + optimal k ---
    mk <- .compute_k_metrics(cc$consensus_results, cc$icl_results, cc$item_consensus_df,
                             working_max_k, safe_pacmap, clustering_matrix,
                             iterated_final_gene_set, clustering_metrics, verbose)

    # --- Feature selection (SWAP POINT) ---
    expr_mat <- clustering_matrix[, iterated_final_gene_set, drop = FALSE]
    feat_df_full <- select_fn(expr_mat, mk$consensus_clusters, adj_pval_thresh)

    # --- FDR threshold filter (shared by all test strategies) ---
    feat_df <- feat_df_full[!is.na(feat_df_full$p_value) &
      feat_df_full$adj_p_value <= adj_pval_thresh, ]


    old_iterated_final_gene_set <- iterated_final_gene_set
    iterated_final_gene_set <- unique(sort(feat_df$gene_id))

    stats_per_iter[[n + 1]] <- feat_df
    cluster_metrics_per_iter[[n + 1]] <- mk$metrics_summary
    optimal_k_per_iter[[n + 1]] <- mk$optimal_k
    cluster_stability_per_iter[[n + 1]] <- mk$optimal_stability
    wide_cluster_stability_per_iter[[n + 1]] <- mk$wide_optimal_stability

    add_log(length(old_iterated_final_gene_set), length(iterated_final_gene_set))

    gene_set_check_1 <- length(old_iterated_final_gene_set)
    gene_set_check_2 <- length(iterated_final_gene_set)
    if (verbose) message(sprintf("%s Iteration %d End: %d genes (%d previously)",
                                 label, n, gene_set_check_2, gene_set_check_1))
    n <- n + 1
    if (n >= max_pipeline_iter) warning(sprintf("Reached max_pipeline_iter (%d) for %s", max_pipeline_iter, label))

    # non-immediate convergence to zero -> revert to last non-zero set
    if (gene_set_check_2 == 0 && gene_set_check_1 < universe_size) {
      cat("Non-immediate convergence to zero detected: reverting to most recent non-zero iteration\n")
      iterated_final_gene_set <- old_iterated_final_gene_set
      iter_log[nrow(iter_log), "after"] <- paste0(length(iterated_final_gene_set),
        "; non-immediate convergence reverted to most recent non-zero iteration")
      break
    }
    # immediate convergence to zero -> retain top 10% by adj_p then effect size
    if (gene_set_check_2 == 0 && n == 1) {
      cat("Immediate convergence to zero detected: retaining top 10% of genes by pval then eff_size. Consider relaxing adj_pval_thresh\n")
      ranked <- feat_df_full[order(feat_df_full$adj_p_value, -feat_df_full$effect_size), ]
      top_n <- ceiling(0.10 * nrow(ranked))
      feat_df <- ranked[seq_len(top_n), ]
      iterated_final_gene_set <- feat_df$gene_id
      iter_log[nrow(iter_log), "after"] <- paste0(length(iterated_final_gene_set),
        "; immediate convergence to zero retained top 10% by pval then eff_size")
      break
    }
  }

  list(final = iterated_final_gene_set, details = save_details(),
       iter_log = iter_log, converged_empty = (length(iterated_final_gene_set) == 0))
}


# =========================================================================
# Phase B: final consensus clustering + visualization + saving for one approach
# =========================================================================
#' @return list(final_df_plot, keep) ; keep = FALSE means "no stable K, drop approach"
.finalize_and_visualize <- function(clustering_gene_set, approach_suffix, clustering_matrix,
                                    pacmap, pm_defaults, pm_visualization, custom.config.visual,
                                    max_k, CCP_iter, clustering_alg, clustering_metrics,
                                    working_dir, verbose) {
  input_data_type <- "pacmap"   # (retained; "genes"/"umap"/"pca" branches removed as dead code, BUGFIX #3)
  working_max_k <- max_k
  clustering_gene_set <- clustering_gene_set[clustering_gene_set %in% colnames(clustering_matrix)]

  # Embeddings for visualization
  umap_matrix_visual <- umap::umap(clustering_matrix[, clustering_gene_set, drop = FALSE], config = custom.config.visual)
  script_pacmap_reducer <- do.call(pacmap$PaCMAP, pm_defaults)
  pacmap_delta_embeddings <- script_pacmap_reducer$fit_transform(clustering_matrix[, clustering_gene_set, drop = FALSE])
  rownames(pacmap_delta_embeddings) <- rownames(clustering_matrix)
  pm_visualization_reducer <- do.call(pacmap$PaCMAP, pm_visualization)
  pacmap_delta_embeddings_visual <- pm_visualization_reducer$fit_transform(clustering_matrix[, clustering_gene_set, drop = FALSE])
  rownames(pacmap_delta_embeddings_visual) <- rownames(clustering_matrix)
  pca_result <- prcomp(clustering_matrix[, clustering_gene_set], center = TRUE, scale. = TRUE)
  pca_scores <- pca_result$x[, 1:3, drop = FALSE]

  if (verbose) message(sprintf("\n=== Creating Output for %s using %s algorithm and %s input ===",
                               approach_suffix, clustering_alg, input_data_type))
  cols <- hue_pal()(working_max_k)

  final_input_matrix <- t(pacmap_delta_embeddings)
  colnames(final_input_matrix) <- rownames(clustering_matrix)

  # Consensus clustering (adaptive pItem)
  cc <- .consensus_adaptive(final_input_matrix, working_max_k, CCP_iter, clustering_alg)
  if (cc$outcome == "skip") {
    message("No stable K values found. Skipping data generation for ", approach_suffix)
    return(list(final_df_plot = NULL, keep = FALSE))
  }
  working_max_k <- cc$working_max_k
  consensus_results <- cc$consensus_results
  icl_results <- cc$icl_results
  item_consensus_df <- cc$item_consensus_df
  item_consensus_by_k <- split(item_consensus_df, item_consensus_df$k)

  # ---- Consensus-matrix heatmaps ----
  CCP_plot_list <- list()
  for (k in 2:length(consensus_results)) {
    cm <- consensus_results[[k]]$consensusMatrix
    colnames(cm) <- colnames(final_input_matrix); rownames(cm) <- colnames(final_input_matrix)
    cm <- cm[consensus_results[[k]]$consensusTree$order, ]
    clusters <- consensus_results[[k]]$consensusClass
    ann <- data.frame(Cluster = droplevels(as.factor(clusters[order(clusters)])))
    cluster_levels <- levels(ann$Cluster); num_clusters <- length(cluster_levels)
    palette <- RColorBrewer::brewer.pal(max(3, num_clusters), "Set2")[1:num_clusters]
    annotation_colors <- list(Cluster = setNames(palette, cluster_levels))
    CCP_plot_list[[paste0("CCP_k_", k)]] <- pheatmap(
      cm, cluster_rows = FALSE, cluster_cols = consensus_results[[k]]$consensusTree,
      annotation_col = ann, annotation_row = ann, show_colnames = FALSE,
      annotation_colors = annotation_colors, color = colorRampPalette(c("white", "blue"))(100),
      silent = TRUE, main = paste("Consensus Matrix for k =", k)
    )
  }

  # ---- Cluster-consensus bar plot ----
  item_consensus_overall <- as.data.frame(icl_results[["clusterConsensus"]]) %>%
    mutate(label = round(clusterConsensus, 2))
  item_con_plot <- ggplot(item_consensus_overall, aes(factor(k), clusterConsensus, fill = factor(cluster))) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label), position = position_dodge(width = 0.9), vjust = -0.3, size = 6) +
    labs(title = "Cluster Consensus by k", x = "k", y = "Consensus Score") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "top",
          axis.text = element_text(colour = "black", size = 24),
          axis.title = element_text(colour = "black", size = 30),
          axis.ticks = element_line(size = 1.5),
          panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
    scale_fill_manual(values = cols[1:working_max_k], name = "Cluster", drop = FALSE)

  item_consensus_by_k_plot_list <- list()
  for (k in names(item_consensus_by_k)) {
    item_consensus_by_k_plot_list[[k]] <- ggplot(
      item_consensus_by_k[[k]], aes(sample_id, item_consensus, fill = factor(cluster))) +
      geom_bar(stat = "identity") +
      labs(title = paste0("Item Consensus by Sample (k = ", k, ")"), x = "Sample ID", y = "Item Consensus", fill = "Cluster") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(size = 10), axis.title = element_text(size = 14),
            legend.position = "right") +
      scale_fill_manual(values = cols[1:working_max_k], name = "Cluster", drop = FALSE)
  }

  # ---- Output directory + tables ----
  results_dir <- paste0(working_dir, "/CIM_states_results_", clustering_alg, "_", input_data_type, "_", approach_suffix)
if (dir.exists(results_dir)) unlink(results_dir, recursive = TRUE)   # ADD: clear stale outputs from prior runs
dir.create(results_dir, recursive = TRUE)
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  write.csv(clustering_gene_set, file.path(results_dir, "clustering_gene_set.csv"), row.names = FALSE)
  for (k_val in names(item_consensus_by_k)) {
    write.csv(item_consensus_by_k[[k_val]], file.path(results_dir, paste0("item_consensus_k", k_val, ".csv")), row.names = FALSE)
  }
  ggsave(file.path(results_dir, "item_con_plot.png"), item_con_plot, width = 10, height = 6, dpi = 300)
  for (k in names(item_consensus_by_k_plot_list)) {
    ggsave(file.path(results_dir, paste0("item_con_sample_k-", k, ".png")),
           item_consensus_by_k_plot_list[[k]], width = 10, height = 6, dpi = 300)
  }
  for (i in seq_along(CCP_plot_list)) {
    plot_name <- names(CCP_plot_list)[i]
    tryCatch({
      png(file.path(results_dir, paste0(plot_name, ".png")), width = 1200, height = 1000, res = 150)
      grid.newpage(); obj <- CCP_plot_list[[i]]
      if (inherits(obj, "grob")) grid.draw(obj)
      else if (inherits(obj, "Heatmap")) draw(obj)
      else if (!is.null(obj$gtable)) grid.draw(obj$gtable)
      else warning(sprintf("Element %s is not drawable.", plot_name))
    },
    error = function(e) message("Heatmap ", plot_name, " skipped: ", conditionMessage(e)),
    finally = if (!is.null(dev.list())) dev.off())   # always release the png device (no leak)
  }

  # ---- CDF plot ----
  k_vec <- 2:working_max_k; cdf_data <- data.frame(); PACs <- numeric(length(k_vec))
  for (i in seq_along(k_vec)) {
    consensus_matrix <- consensus_results[[k_vec[i]]][["consensusMatrix"]]
    cons_vals <- consensus_matrix[lower.tri(consensus_matrix)]
    Fn <- ecdf(cons_vals); PACs[i] <- Fn(0.9) - Fn(0.1)
    h <- hist(cons_vals, breaks = 200, plot = FALSE)
    cdf_data <- bind_rows(cdf_data, data.frame(k = as.factor(k_vec[i]), x = h$mids, y = cumsum(h$counts) / sum(h$counts)))
  }
  cdf_data <- left_join(cdf_data, data.frame(k = as.factor(k_vec), PAC = PACs), by = "k")
  cdf_data$line_label <- paste0("K = ", cdf_data$k, " (PAC=", formatC(cdf_data$PAC, format = "f", digits = 3), ")")
  cdf_plot <- ggplot(cdf_data, aes(x, y, color = line_label, group = k)) +
    geom_line(size = 1.5) +
    labs(x = "Consensus Index", y = "Cumulative Fraction", title = "Consensus CDFs across Ks", color = "K (PAC)") +
    geom_vline(xintercept = c(0.1, 0.9), linetype = "dashed", color = "gray") +
    scale_x_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.98, 0.02),
          legend.justification = c("right", "bottom"), legend.text = element_text(size = 16),
          legend.background = element_rect(fill = "white", color = "black", size = 0.8),
          axis.text = element_text(colour = "black", size = 16),
          axis.title = element_text(colour = "black", size = 30), axis.ticks = element_line(size = 1.5),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  ggsave(file.path(results_dir, "final_cdf_plot.png"), cdf_plot, width = 10, height = 6, dpi = 300)

  # ---- PAC plot ----
  Kvec <- 2:working_max_k; PAC <- rep(NA, length(Kvec)); names(PAC) <- paste0("K=", Kvec)
  for (i in Kvec) { M <- consensus_results[[i]]$consensusMatrix; Fn <- ecdf(M[lower.tri(M)]); PAC[i - 1] <- Fn(0.9) - Fn(0.1) }
  pac_values <- data.frame(k = Kvec, pac = unname(PAC))
  pac_plot <- ggplot(pac_values, aes(factor(k), pac, group = 1)) +
    geom_line() + theme_bw(base_rect_size = 1.5) +
    geom_point(size = 4, shape = 21, color = "darkred", fill = "orange") +
    xlab("Number of clusters K") + ylab("PAC Score") + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "top",
          axis.text = element_text(colour = "black", size = 24),
          axis.title = element_text(colour = "black", size = 30), axis.ticks = element_line(size = 1.5),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  ggsave(file.path(results_dir, "final_pac_plot.png"), pac_plot, width = 10, height = 6, dpi = 300)

  # ---- Silhouette (both spaces) ----
  orig_space_input <- t(final_input_matrix); rownames(orig_space_input) <- rownames(clustering_matrix)
  dist_space_for_clustering <- dist(orig_space_input)
  dist_ge_space <- dist(clustering_matrix[, clustering_gene_set, drop = FALSE])
  sil_df <- do.call(rbind, lapply(2:length(consensus_results), function(k) {
    clusters <- consensus_results[[k]]$consensusClass
    data.frame(k = k,
               avg_sil_space_for_clustering = mean(silhouette(clusters, dist_space_for_clustering)[, 3]),
               avg_sil_ge_space = mean(silhouette(clusters, dist_ge_space)[, 3]))
  }))
  sil_df$avg_sil_combined <- (sil_df$avg_sil_space_for_clustering + sil_df$avg_sil_ge_space) / 2
  sil_theme <- theme_classic(base_size = 16) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "top",
          axis.text = element_text(colour = "black", size = 24),
          axis.title = element_text(colour = "black", size = 30), axis.ticks = element_line(size = 1.5),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  sil_plot_both <- ggplot(sil_df, aes(x = k)) +
    geom_line(aes(y = avg_sil_space_for_clustering, color = input_data_type)) +
    geom_point(aes(y = avg_sil_space_for_clustering, color = input_data_type)) +
    geom_line(aes(y = avg_sil_ge_space, color = "Gene space")) +
    geom_point(aes(y = avg_sil_ge_space, color = "Gene space")) +
    labs(y = "Average silhouette width", color = "Space", title = "Silhouette comparison across k") + sil_theme
  ggsave(file.path(results_dir, "sil_plot.png"), sil_plot_both, width = 10, height = 6, dpi = 300)
  ggsave(file.path(results_dir, "sil_plot_original_space.png"),
         ggplot(sil_df, aes(x = k)) +
           geom_line(aes(y = avg_sil_space_for_clustering, color = input_data_type)) +
           geom_point(aes(y = avg_sil_space_for_clustering, color = input_data_type)) +
           labs(y = "Average silhouette width", color = "Space", title = "Silhouette comparison across k") + sil_theme,
         width = 10, height = 6, dpi = 300)
  ggsave(file.path(results_dir, "sil_plot_ge_space.png"),
         ggplot(sil_df, aes(x = k)) +
           geom_line(aes(y = avg_sil_ge_space, color = "Gene space")) +
           geom_point(aes(y = avg_sil_ge_space, color = "Gene space")) +
           labs(y = "Average silhouette width", color = "Space", title = "Silhouette comparison across k") + sil_theme,
         width = 10, height = 6, dpi = 300)

  # ---- Metrics summary + rank-based optimal k (BUGFIX #6 applied) ----
  avg_cluster_consensus <- item_consensus_overall %>%
    group_by(k) %>% summarize(avg_cluster_consensus = mean(clusterConsensus))
  metrics_summary <- data.frame(
    k = 2:working_max_k, silhouette_dim_reduce_space = NA, silhouette_ge_space = NA,
    silhouette_combined_avg = NA, pac = NA, avg_cluster_consensus = NA)
  for (i in seq_len(nrow(metrics_summary))) {
    k_val <- metrics_summary$k[i]
    idx <- which(sil_df$k == k_val)
    if (length(idx) > 0) {
      metrics_summary$silhouette_dim_reduce_space[i] <- sil_df$avg_sil_space_for_clustering[idx]
      metrics_summary$silhouette_ge_space[i] <- sil_df$avg_sil_ge_space[idx]
      metrics_summary$silhouette_combined_avg[i] <- sil_df$avg_sil_combined[idx]
    }
    k_idx <- which(Kvec == k_val); if (length(k_idx) > 0) metrics_summary$pac[i] <- PAC[k_idx]
    cc_idx <- which(avg_cluster_consensus$k == k_val)
    if (length(cc_idx) > 0) metrics_summary$avg_cluster_consensus[i] <- avg_cluster_consensus$avg_cluster_consensus[cc_idx]
  }
  metrics_summary$rank_silhouette_dim_reduce_space <- rank(-metrics_summary$silhouette_dim_reduce_space, ties.method = "min")
  metrics_summary$rank_silhouette_ge_space <- rank(-metrics_summary$silhouette_ge_space, ties.method = "min")
  metrics_summary$rank_silhouette_combined_avg <- rank(-metrics_summary$silhouette_combined_avg, ties.method = "min")
  metrics_summary$rank_pac <- rank(metrics_summary$pac, ties.method = "min")
  metrics_summary$rank_item_cluster_consensus <- rank(-metrics_summary$avg_cluster_consensus, ties.method = "min")
  rank_columns <- paste0("rank_", clustering_metrics)
  metrics_summary$overall_rank <- rowSums(metrics_summary[, rank_columns, drop = FALSE], na.rm = TRUE)
  best_rows <- metrics_summary[metrics_summary$overall_rank == min(metrics_summary$overall_rank), ]
  tie_metric <- clustering_metrics[length(clustering_metrics)]
  idx <- if (tie_metric == "silhouette_dim_reduce_space") which.max(best_rows$silhouette_dim_reduce_space)
    else if (tie_metric %in% c("silhouette_ge_space", "silhouette_ge")) which.max(best_rows$silhouette_ge_space)
    else if (tie_metric == "silhouette_combined_avg") which.max(best_rows$silhouette_combined_avg)
    else if (tie_metric == "item_cluster_consensus") which.max(best_rows$avg_cluster_consensus)
    else if (tie_metric == "pac") which.min(best_rows$pac)
    else stop("Unknown tiebreaker metric.")
  optimal_k <- best_rows[idx, ]$k
  metrics_summary$metrics_used <- paste(clustering_metrics, collapse = ", ")
  metrics_summary$metrics_used_for_tie <- tie_metric
  write.csv(metrics_summary, file.path(results_dir, "final_cluster_metrics_summary.csv"), row.names = FALSE)
  cat("Using k =", optimal_k, "as the optimal number of clusters\n")

  consensus_clusters <- consensus_results[[optimal_k]]$consensusClass
  optimal_stability_cluster_metrics <- item_consensus_by_k[[as.character(optimal_k)]]
  wide_optimal_stability_cluster_metrics <- optimal_stability_cluster_metrics %>%
    dplyr::select(sample_id, cluster, item_consensus) %>%
    pivot_wider(names_from = cluster, values_from = item_consensus, names_prefix = "item_consensus_cluster")
  write.csv(optimal_stability_cluster_metrics, file.path(results_dir, "final_cluster_stability_metrics.csv"), row.names = FALSE)
  write.csv(wide_optimal_stability_cluster_metrics, file.path(results_dir, "final_wide_cluster_stability_metrics.csv"), row.names = FALSE)

  # ---- Final combined embedding data frame ----
  final_df_plot <- as.data.frame(umap_matrix_visual$layout) %>%
    rownames_to_column("sample_id") %>% dplyr::rename(UMAP1 = V1, UMAP2 = V2) %>%
    dplyr::inner_join(as.data.frame(pacmap_delta_embeddings_visual) %>%
                        rownames_to_column("sample_id") %>% dplyr::rename(PACMAP1 = V1, PACMAP2 = V2), by = "sample_id") %>%
    dplyr::inner_join(as.data.frame(pca_scores) %>% rownames_to_column("sample_id"), by = "sample_id") %>%
    dplyr::mutate(base_id = gsub("^delta_([^_]+)_.*", "\\1", sample_id)) %>%
    dplyr::mutate(cluster_assignments = consensus_clusters) %>%
    dplyr::inner_join(as.data.frame(clustering_matrix) %>% rownames_to_column("sample_id"), by = "sample_id")

  # ---- fviz 2D PCA ----
  pca_fviz_input <- final_df_plot[, c("PC1", "PC2", "sample_id")] %>% tibble::column_to_rownames("sample_id")
  fviz_cluster_plot <- fviz_cluster(
    object = list(data = pca_fviz_input, cluster = as.integer(consensus_clusters)),
    palette = cols[1:optimal_k], geom = "text", ellipse.type = "convex",
    ggtheme = theme_bw(), ggtitle("PCA plot on Clustered Data"))
  ggsave(file.path(results_dir, paste0("final_fviz_2D_PCA_k-", optimal_k, ".png")), fviz_cluster_plot, width = 10, height = 8, dpi = 300)

  # ---- 3D PCA (plotly) ----
  pca_fviz_input_3D <- as.matrix(pca_result$x[, 1:3, drop = FALSE]); rownames(pca_fviz_input_3D) <- final_df_plot$sample_id
  plot_df_3D_PCA <- data.frame(PC1 = pca_fviz_input_3D[, 1], PC2 = pca_fviz_input_3D[, 2], PC3 = pca_fviz_input_3D[, 3],
                               cluster = as.factor(as.integer(consensus_clusters)), sample_id = rownames(pca_fviz_input_3D))
  fig_PCA_3D <- plot_ly(plot_df_3D_PCA, x = ~PC1, y = ~PC2, z = ~PC3, color = ~cluster, text = ~sample_id,
                        type = "scatter3d", mode = "markers", colors = cols) %>%
    plotly::layout(title = "3D PCA Cluster Visualization",
                   scene = list(xaxis = list(title = "PC1"), yaxis = list(title = "PC2"), zaxis = list(title = "PC3")))
  # selfcontained = TRUE needs pandoc; fall back to a deps-folder html if pandoc is unavailable
  .widget_file <- paste0(results_dir, "/final_3D_PCA_k-", optimal_k, ".html")
  # Fully non-fatal: the 3D widget is optional. Neither attempt may abort Phase B
  # (a throw here previously skipped the clustered_samples save, e.g. app_three).
  .w_ok <- tryCatch({ saveWidget(fig_PCA_3D, file = .widget_file, selfcontained = TRUE); TRUE },
                    error = function(e) FALSE)
  if (!.w_ok) .w_ok <- tryCatch({ saveWidget(fig_PCA_3D, file = .widget_file, selfcontained = FALSE); TRUE },
                                 error = function(e) FALSE)
  if (!.w_ok) message("3D PCA widget skipped (pandoc unavailable); continuing without it.")

  # ---- UMAP / PaCMAP hull plots ----
  hull_theme <- theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "top",
          axis.text = element_text(colour = "black", size = 24),
          axis.title = element_text(colour = "black", size = 30), axis.ticks = element_line(size = 1.5),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  make_hull_plot <- function(xvar, yvar, xlab, ylab) {
    ggplot(final_df_plot, aes(x = .data[[xvar]], y = .data[[yvar]], label = base_id, color = as.factor(cluster_assignments))) +
      scale_color_manual(values = cols, name = "Cluster") +
      ggforce::geom_mark_hull(aes(group = as.factor(cluster_assignments), fill = as.factor(cluster_assignments), label = NULL),
                              concavity = 2, expand = unit(2, "mm"), alpha = 0.2) +
      geom_label_repel(size = 4, fontface = "bold", label.size = 0.3, fill = "white", alpha = 1,
                       box.padding = 0.25, max.overlaps = 100) +
      labs(x = xlab, y = ylab, fill = "Cluster", color = "Cluster",
           title = paste0("Final Clustering Using ", input_data_type, " input and ", clustering_alg,
                          " CCP Alg (k = ", optimal_k, ")")) + hull_theme
  }
  ggsave(file.path(results_dir, paste0("final_clustering_shown_on_umap_k-", optimal_k, ".png")),
         make_hull_plot("UMAP1", "UMAP2", "UMAP 1", "UMAP 2"), width = 10, height = 8, dpi = 300)
  ggsave(file.path(results_dir, paste0("final_clustering_shown_on_pacmap_k-", optimal_k, ".png")),
         make_hull_plot("PACMAP1", "PACMAP2", "PacMAP 1", "PacMAP 2"), width = 10, height = 8, dpi = 300)

  list(final_df_plot = final_df_plot, keep = TRUE, results_dir = results_dir)
}


# =========================================================================
# Main function
# =========================================================================
CIMIC <- function(
    clustering_matrix,                       # rows: samples x cols: genes
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
    pacmap_args = NULL,
    pacmap_guardrails = TRUE,
    select_fn = cimic_select_features,       # SWAP POINT (classic / adaptive / limma)
    verbose = TRUE,
    working_dir) {

  set.seed(2024L)
  graphics.off()   # start from a clean graphics-device state (guards against a corrupted
                   # "invalid graphics state" carried over from a previous/aborted run)
  pacmap <- setup_conda(env_name = "pacmap_env")

  # ---- Resolve gene sets ----
  if (!is.null(all_gene_sets)) message("Using user provided genesets...")
  if (is.null(all_gene_sets)) {
    message("Using default CIM manuscript genesets...")
    msig_df <- msigdbr::msigdbr(species = "Homo sapiens") %>% as.data.frame()
    all_gene_set_names <- c(
      "GOBP_ADAPTIVE_IMMUNE_RESPONSE", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
      "GOBP_B_CELL_ACTIVATION", "GOBP_CELLULAR_RESPONSE_TO_STRESS",
      "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
      "GOBP_INFLAMMASOME_MEDIATED_SIGNALING_PATHWAY", "GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS",
      "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
      "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
      "GOBP_LEUKOCYTE_CHEMOTAXIS_INVOLVED_IN_INFLAMMATORY_RESPONSE",
      "GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
      "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS", "GOBP_T_CELL_ACTIVATION",
      "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
      "REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS", "REACTOME_CELLULAR_RESPONSE_TO_MITOCHONDRIAL_STRESS",
      "WP_MOLECULAR_PATHWAY_FOR_OXIDATIVE_STRESS", "WP_OXIDATIVE_STRESS_RESPONSE"
    )
    all_gene_sets <- setNames(lapply(all_gene_set_names, function(gs) {
      msig_df[msig_df$gs_name == gs, , drop = FALSE]$gene_symbol
    }), all_gene_set_names)
  }

  # ---- Optional immune/TCR gene removal ----
  immune_genes_removed_all <- character(0)
  if (remove_immune_variable_genes) {
    removed_immune_genes <- c()
    remove_prefixes <- c(
      "IGHV", "IGHD", "IGHJ", "IGHC", "IGHM", "IGHA", "IGHE", "IGHG",
      "IGLV", "IGLC", "IGLJ", "IGKV", "IGKC", "IGKJ",
      "TRAV", "TRBV", "TRDV", "TRGV", "TRAJ", "TRBJ", "TRGJ", "TRDJ",
      "TRAC", "TRBC", "TRGC", "TRDC"
    )
    pattern <- paste0("^(", paste(remove_prefixes, collapse = "|"), ")")
    all_gene_sets <- lapply(all_gene_sets, function(gene_vec) {
      removed_immune_genes <<- c(removed_immune_genes, gene_vec[grepl(pattern, gene_vec, perl = TRUE)])
      gene_vec[!grepl(pattern, gene_vec, perl = TRUE)]
    })
    immune_genes_removed_all <- unique(removed_immune_genes)
  }

  # ---- UMAP seed ----
  custom.config <- umap.defaults
  custom.config$random_state <- 2024L

  # ---- Intersect with matrix, drop zero-variance genes ----
  genes <- colnames(clustering_matrix)
  all_gene_sets <- lapply(all_gene_sets, function(gs) intersect(gs, genes))
  removed_by_membership <- unique(unlist(lapply(all_gene_sets, function(gs) setdiff(gs, genes))))
  zero_var_genes <- colnames(clustering_matrix)[apply(clustering_matrix, 2, sd) == 0]
  removed_genes_list <- lapply(all_gene_sets, function(gene_vec) intersect(gene_vec, zero_var_genes))
  all_gene_sets <- lapply(all_gene_sets, function(gene_vec) setdiff(gene_vec, zero_var_genes))
  removed_genes <- sort(unique(c(unlist(removed_genes_list), unlist(removed_by_membership), immune_genes_removed_all)))

  # ---- Validate arguments ----
  if (!(clustering_alg %in% c("km", "hc"))) stop("clustering_alg must be 'km' or 'hc'")
  valid_metrics <- c("pac", "silhouette_dim_reduce_space", "silhouette_ge_space", "silhouette_combined_avg", "item_cluster_consensus")
  if (any(!clustering_metrics %in% valid_metrics)) {
    stop(paste0("clustering_metrics must be from: ", paste(valid_metrics, collapse = ", ")))
  }

  # ---- Retained gene-count checks ----
  for (gs_name in names(all_gene_sets)) {
    n_retained <- length(all_gene_sets[[gs_name]])
    if (n_retained == 0) stop(sprintf("Gene set '%s' has no genes after filtering.", gs_name), call. = FALSE)
    orig_len <- if (exists("all_gene_set_names")) length(msig_df$gene_symbol[msig_df$gs_name == gs_name]) else n_retained
    min_warn <- max(ceiling(0.10 * orig_len), 1L)
    if (n_retained < min_warn) {
      warning(sprintf("Gene set '%s' retained only %d genes (10%% of original = %d).", gs_name, n_retained, min_warn), call. = FALSE)
    }
  }

  # app_three depends on app_one
  if ("app_three" %in% filter_approach && !"app_one" %in% filter_approach) filter_approach <- c(filter_approach, "app_one")

  # ---- PaCMAP / UMAP settings ----
  np <- reticulate::import("numpy", convert = FALSE); np$random$seed <- 2024L
  sample_number <- nrow(clustering_matrix)
  n_neigh_pacmap <- if (sample_number <= 20L) as.integer(max(5L, floor(sample_number * 0.25))) else NULL
  pm_defaults <- list(n_components = pacmap_dimensions, n_neighbors = n_neigh_pacmap, random_state = seed)
  if (!is.null(pacmap_args)) {
    if (pacmap_guardrails) {
      allowed <- c("n_components", "n_neighbors", "MN_ratio", "FP_ratio", "num_iters", "apply_pca", "random_state", "PCA_dim")
      bad <- setdiff(names(pacmap_args), allowed)
      if (length(bad)) { warning("Ignoring unrecognized pacmap_args: ", paste(bad, collapse = ", ")); pacmap_args[bad] <- NULL }
    }
    pm_defaults <- utils::modifyList(pm_defaults, pacmap_args)
  }
  pm_settings <- as.data.frame(unlist(pm_defaults))
  pm_visualization <- pm_defaults; pm_visualization$n_components <- 2L
  n_neigh_umap <- if (sample_number <= 20L) as.integer(max(5L, floor(sample_number * 0.25))) else 15L
  custom.config$n_neighbors <- min(n_neigh_umap, sample_number - 1L)
  custom.config.visual <- utils::modifyList(umap.defaults, list(n_components = 2L, n_neighbors = n_neigh_umap, random_state = 2024L))

  # ---- Output containers ----
  details_by_gene_sets <- list(); details_by_all_genes <- list(); details_gene_set_then_all <- list()
  iter_log <- data.frame(gene_set = character(), iteration = integer(), before = integer(), after = character(), stringsAsFactors = FALSE)
  new_final <- character(0)
  no_app_one <- FALSE; no_app_two <- FALSE; no_app_three <- FALSE
  iterated_over_all_genes <- character(0); iterated_use_final_gene_set_from_approach_one <- character(0)

  common <- list(clustering_matrix = clustering_matrix, pacmap = pacmap, pm_defaults = pm_defaults,
                 max_k = max_k, CCP_iter = CCP_iter, clustering_alg = clustering_alg,
                 clustering_metrics = clustering_metrics, adj_pval_thresh = adj_pval_thresh,
                 max_pipeline_iter = max_pipeline_iter, select_fn = select_fn, verbose = verbose)
  refine <- function(seed_genes, label, universe_size) {
    do.call(.refine_gene_set, c(list(seed_genes = seed_genes, label = label, universe_size = universe_size), common))
  }

  # ==== FEATURE-SELECTION PHASE ====
  # ---- Approach One: refine each set independently, union survivors ----
  if ("app_one" %in% filter_approach) {
    withr::with_seed(2024L, {
      for (gene_set_i in seq_along(all_gene_sets)) {
        gene_set_name <- names(all_gene_sets)[gene_set_i]
        if (is.null(gene_set_name)) gene_set_name <- paste0("set", gene_set_i)
        if (verbose) message(sprintf("\nGene set: %s", gene_set_name))
        seed_genes <- all_gene_sets[[gene_set_i]]
        # BUGFIX #1: pass the ORIGINAL set size (was length(gene_set_i) == 1 -> revert never fired)
        res <- refine(seed_genes, gene_set_name, length(seed_genes))
        details_by_gene_sets[[gene_set_name]] <- res$details
        iter_log <- rbind(iter_log, res$iter_log)
        new_final <- sort(unique(c(new_final, res$final)))
      }
      if (length(new_final) == 0) {
        cat("\033[1;31mWARNING: Final Gene Set for Approach one empty! Consider Relaxing P-val\033[0m\n")
        no_app_one <- TRUE
      }
    })
  }

  # ---- Approach Two: refine the pooled union ----
  if ("app_two" %in% filter_approach) {
    withr::with_seed(2024L, {
      if (verbose) message("\n[Second Approach: iteration over all unique genes]")
      pooled <- unique(unlist(all_gene_sets))
      res <- refine(pooled, "Analysis by all genes", length(pooled))
      details_by_all_genes <- res$details
      iter_log <- rbind(iter_log, res$iter_log)
      iterated_over_all_genes <- res$final
      if (length(iterated_over_all_genes) == 0) {
        cat("\033[1;31mWARNING: Final Gene Set for Approach two empty! Consider Relaxing P-val\033[0m\n")
        no_app_two <- TRUE
      }
    })
  }

  # ---- Approach Three: refine app_one's union again ----
  if ("app_three" %in% filter_approach) {
    withr::with_seed(2024L, {
      if (verbose) message("\n[Third Approach: gene set result then all genes]")
      res <- refine(new_final, "Analysis gene set then all genes", length(new_final))
      details_gene_set_then_all <- res$details
      iter_log <- rbind(iter_log, res$iter_log)
      iterated_use_final_gene_set_from_approach_one <- res$final
      if (length(iterated_use_final_gene_set_from_approach_one) == 0) {
        cat("\033[1;31mWARNING: Final Gene Set for Approach three empty! Consider Relaxing P-val\033[0m\n")
        no_app_three <- TRUE
      }
    })
  }

  # ==== PHASE B: final clustering + visualization + saving ====
  final_df_plot_app_one <- NULL; final_df_plot_app_two <- NULL; final_df_plot_app_three <- NULL

  finalize <- function(gene_set, suffix) {
    # Phase-B visualization must never abort the whole run: if embedding/clustering
    # fails for a degenerate final gene set (e.g. collapsed to 1-2 genes so PaCMAP
    # cannot build a 2D map), skip this approach's plots but KEEP its gene set.
    tryCatch(
      .finalize_and_visualize(gene_set, suffix, clustering_matrix, pacmap, pm_defaults, pm_visualization,
                              custom.config.visual, max_k, CCP_iter, clustering_alg, clustering_metrics, working_dir, verbose),
      error = function(e) {
        message(sprintf("Visualization for %s skipped (%s); gene set retained, run continues.",
                        suffix, conditionMessage(e)))
        list(final_df_plot = NULL, keep = TRUE, results_dir = NULL)
      }
    )
  }
  # write the common extra-output files into a results dir (per approach, as original did)
  save_extra <- function(results_dir, final_df, suffix) {
    if (is.null(results_dir)) return(invisible(NULL))   # viz was skipped -> nothing to save
    write.csv(clustering_matrix, file.path(results_dir, "initial_clustering_mat.csv"), row.names = TRUE)
    write.csv(iter_log, file.path(results_dir, "iter_log_all_approaches.csv"), row.names = FALSE)
    write.csv(data.frame(list_name = names(all_gene_sets),
                         items = sapply(all_gene_sets, function(x) paste(x, collapse = ","))),
              file.path(results_dir, "all_gene_sets.csv"), row.names = FALSE)
    write.csv(final_df, file.path(results_dir, paste0("clustered_samples_", suffix, ".csv")), row.names = FALSE)
    write.csv(pm_settings, file.path(results_dir, "pacmap_settings.csv"), row.names = TRUE)
    write.csv(removed_genes, file.path(results_dir, "removed_genes.csv"), row.names = TRUE)
  }

  if ("app_one" %in% filter_approach) {
    withr::with_seed(2024L, {
      if (!no_app_one) {
        fv <- finalize(new_final, "app_one")
        if (!fv$keep) { filter_approach <- filter_approach[filter_approach != "app_one"] }
        else { final_df_plot_app_one <- fv$final_df_plot; save_extra(fv$results_dir, final_df_plot_app_one, "app_one") }
      }
    })
  }
  if ("app_two" %in% filter_approach) {
    withr::with_seed(2024L, {
      if (!no_app_two) {
        fv <- finalize(iterated_over_all_genes, "app_two")
        if (!fv$keep) { filter_approach <- filter_approach[filter_approach != "app_two"] }
        else { final_df_plot_app_two <- fv$final_df_plot; save_extra(fv$results_dir, final_df_plot_app_two, "app_two") }
      }
    })
  }
  if ("app_three" %in% filter_approach) {
    withr::with_seed(2024L, {
      if (!no_app_three) {
        fv <- finalize(iterated_use_final_gene_set_from_approach_one, "app_three")
        if (!fv$keep) { filter_approach <- filter_approach[filter_approach != "app_three"] }
        else { final_df_plot_app_three <- fv$final_df_plot; save_extra(fv$results_dir, final_df_plot_app_three, "app_three") }
      }
    })
  }

  # ---- Return ----
  list(
    iterated_by_gene_sets = if ("app_one" %in% filter_approach) new_final else "NA",
    iterated_over_all_genes = if ("app_two" %in% filter_approach) iterated_over_all_genes else "NA",
    iterated_use_final_gene_set_from_approach_one = if ("app_three" %in% filter_approach) iterated_use_final_gene_set_from_approach_one else "NA",
    details_per_set_analysis_by_gene_sets = if ("app_one" %in% filter_approach) details_by_gene_sets else "NA",
    details_per_set_analysis_by_all_genes = if ("app_two" %in% filter_approach) details_by_all_genes else "NA",
    details_per_set_analysis_by_gene_set_then_all_genes = if ("app_three" %in% filter_approach) details_gene_set_then_all else "NA",
    iter_log = iter_log,
    pacmap_settings = pm_settings,
    removed_genes = removed_genes,
    final_df_app_one = if ("app_one" %in% filter_approach) final_df_plot_app_one else "NA",
    final_df_app_two = if ("app_two" %in% filter_approach) final_df_plot_app_two else "NA",
    final_df_app_three = if ("app_three" %in% filter_approach) final_df_plot_app_three else "NA"
  )
}

# -------------------------------------------------------------------------
# Documented bug fixes carried into this refactor (see "# BUGFIX" tags):
#   #1  app_one revert-on-zero used length(gene_set_i) (==1); now uses the true
#       original set size, so the "non-immediate convergence" revert can fire.
#   #2  app_two/app_three PaCMAP-failure messages referenced app_one's
#       `gene_set_name`; now each run carries its own `label`.
#   #3  Phase-B "umap"/"pca" branches referenced never-created objects
#       (`umap_matrix`, `umap_plot_df`); input is fixed to "pacmap" and those
#       dead branches were removed.
#   #4  calcICL was called twice back-to-back with identical args in every
#       consensus loop; now called once.
#   #6  Tiebreak tested "silhouette_ge" but the valid name is
#       "silhouette_ge_space"; both are now accepted.
#   #7  app_two/app_three empty-set guards now handled uniformly alongside
#       filter_approach membership.
#   Also: removed the unused/misnamed `details_per_set_gene_sets_then_all_genes`
#       initializer; container naming is now consistent.
# -------------------------------------------------------------------------
