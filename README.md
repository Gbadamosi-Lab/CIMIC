# CIMIC: Chemoimmunomodulation Induction Classifier

## Overview

CIMIC is an unsupervised, stability-guided iterative computational pipeline for
stratifying cancer patients based on chemotherapy-induced chemoimmunomodulatory
(CIM) transcriptional trajectories. It was developed and validated in the context
of neoadjuvant chemotherapy in breast cancer.

The pipeline operates on within-patient delta gene expression values
(ΔGE = post-treatment log₂[TPM+1] − pre-treatment log₂[TPM+1]) across 3,179 genes
spanning 19 CIM-related pathways, capturing treatment-induced transcriptional
responses while controlling for interpatient baseline differences.

This GitHub repository serves as a location for obtaining the code for CIMIC and a repository for data and code used to make figures included in the original publication of CIMIC on Oncoimmunology.

---

## Pipeline Versions

Two main release scripts are provided:

| Script | Use Case |
|--------|----------|
| `CIMIC_limma_1.0.0.R` | **Primary release.** For patient-level or sample-level data without replicates. Uses limma's moderated linear model for feature selection. |
| `CIMIC_limma_replicate_1.0.0.R` | **Replicate-aware release.** For replicate-level data (e.g., cell lines) where treating replicates as independent would constitute pseudoreplication. Uses limma's `duplicateCorrelation` with a sample blocking factor to model within-cell-line correlation explicitly. |

> **Which to use?** Use `CIMIC_limma_1.0.0.R` for patient cohort data or data with one replicate per sample. Use
> `CIMIC_limma_replicate_1.0.0.R` for experiments where each sample
> has multiple replicate measurements. Input for the replicate version should be
> the replicate-level matrices (e.g., `tnbc_cl_<drug>_initial_clustering_mat.csv`),
> **not** pre-averaged matrices.

---

## How It Works

CIMIC proceeds in three stages:

**Stage 1 — Per-pathway iterative clustering and feature selection**

For each of the 19 CIM pathways, ΔGE matrices are:

1. Reduced using PaCMAP dimensionality reduction
2. Clustered using ConsensusClusterPlus with hierarchical clustering (adaptive
   `pItem` to resolve NA item-consensus values)
3. Evaluated for optimal *k* using a composite ranking framework incorporating:
   - Average silhouette coefficient (s_avg) in both dimensional reduction and ΔGE space
   - Proportion of Ambiguous Clustering (PAC)
   - Cluster Consensus Score (CCS)
4. Filtered to retain only genes showing significant ΔGE differences across clusters
   using a **moderated linear model (limma)**:
   - 2 groups → moderated *t*-test; effect size = |*t*|
   - ≥3 groups → moderated *F*-test; effect size = *F*
   - FDR correction via Benjamini–Hochberg; filter: *p* ≤ 0.05 & adj. *p* ≤ threshold

   In the replicate-aware version, `duplicateCorrelation` is applied first to
   estimate the within-cell-line correlation, which is then passed to `lmFit`
   as a random block effect before empirical-Bayes moderation.

Steps 2–4 are repeated iteratively until convergence (no further gene exclusion),
up to `max_pipeline_iter` iterations (default: 50).

**Stage 2 — Pathway-refined global gene set construction**

Discriminatory genes from all pathways are pooled into a pathway-refined global gene set.

**Stage 3 — Final CIMIC classification**

The same iterative procedure is applied to the global gene set. Final cluster
assignments at convergence define the CIMIC trajectories.

### Filtering Approaches

Three complementary approaches are run (all enabled by default; can be changed to save computing time.):

| Approach | Description |
|----------|-------------|
| `app_one` | Refine each gene set independently, then take the union of survivors (i.e., just clustering using stage 1 approach)|
| `app_two` | Refine the pooled union of all gene sets at once |
| `app_three` | Refine `app_one`'s union again over the combined set (depends on `app_one`; i.e., using stage 1, stage 2, and stage 3) |

### Convergence Handling

- **Normal convergence**: iteration stops when gene count is unchanged between steps.
- **Non-immediate convergence to zero**: reverts to the most recent non-zero gene set.
- **Immediate convergence to zero** (first iteration): retains the top 10% of genes
  ranked by adjusted *p*-value then effect size, with a warning to consider relaxing
  `adj_pval_thresh`.

---

## CIMIC Trajectories

In the breast cancer neoadjuvant chemotherapy cohort, CIMIC identifies two
reproducible trajectories:

- **Fun-CIM**: Functional, immune-activating trajectory
- **Dys-CIM**: Dysfunctional, stress-adaptive trajectory

---

## Requirements

### R version
R ≥ 4.5.0

### R packages (installed automatically on first run)

**CRAN:**
`umap`, `cluster`, `factoextra`, `dplyr`, `stringr`, `magrittr`, `tibble`,
`rstatix`, `coin`, `RColorBrewer`, `ggforce`, `grid`, `gridExtra`, `corrplot`,
`concaveman`, `paran`, `reticulate`, `plotly`, `scales`, `withr`, `ggrepel`,
`ggplot2`, `tidyr`, `msigdbr`, `jsonlite`, `pheatmap`, `htmlwidgets`

**Bioconductor:**
`ConsensusClusterPlus`, `ComplexHeatmap`, `limma`

### Python (for PaCMAP)
- Anaconda or Miniconda
- The pipeline will automatically set up a `pacmap_env` conda environment on first run
- Alternatively, set the `RETICULATE_CONDA` environment variable to your conda path:

```r
Sys.setenv(RETICULATE_CONDA = "path/to/conda.exe")
```

---

## Input Data Format

### Primary release (`CIMIC_limma_1.0.0.R`)
A **clustering matrix** where:
- Rows: samples (patient IDs)
- Columns: genes (gene symbols)
- Values: ΔGE (log₂[TPM+1] post − log₂[TPM+1] pre)

### Replicate-aware release (`CIMIC_limma_replicate_1.0.0.R`)
A **replicate-level clustering matrix** where:
- Rows: replicate samples with IDs in the format `delta_<CELLLINE>_<DRUG>_<rep>`
  (the cell-line blocking factor is derived automatically by stripping the trailing `_<rep>`)
- Columns: genes (gene symbols)
- Values: ΔGE values at the replicate level (do **not** pre-average across replicates)

```r
# Example structure
dim(clustering_matrix)        # samples (or replicates) x genes
rownames(clustering_matrix)   # sample / replicate IDs
colnames(clustering_matrix)   # gene symbols
```

---

## Usage

### Primary release

```r
source("CIMIC_limma_1.0.0.R")

results <- CIM_feature_selection_by_gene_set_pacmap(
  clustering_matrix  = your_delta_ge_matrix,  # samples x genes
  clustering_alg     = "hc",                  # "hc" (hierarchical) or "km" (k-means)
  max_k              = 5,                     # maximum clusters to test
  CCP_iter           = 5000,                  # ConsensusClusterPlus iterations
  adj_pval_thresh    = 0.05,                  # FDR threshold
  max_pipeline_iter  = 50,                    # max convergence iterations
  seed               = 2024L,
  clustering_metrics = c("pac", "silhouette_combined_avg", "item_cluster_consensus"),
  filter_approach    = c("app_one", "app_two", "app_three"),
  working_dir        = "path/to/your/output/directory"
)
```

### Replicate-aware release

```r
source("CIMIC_limma_replicate_1.0.0.R")

results <- CIM_feature_selection_by_gene_set_pacmap(
  clustering_matrix  = your_replicate_delta_ge_matrix,  # replicates x genes
  clustering_alg     = "hc",
  max_k              = 5,
  CCP_iter           = 5000,
  adj_pval_thresh    = 0.05,
  max_pipeline_iter  = 50,
  seed               = 2024L,
  clustering_metrics = c("pac", "silhouette_combined_avg", "item_cluster_consensus"),
  filter_approach    = c("app_one", "app_two", "app_three"),
  working_dir        = "path/to/your/output/directory"
)
```

### Return value

Both versions return a named list:

```r
results$iterated_by_gene_sets                          # Approach 1 final gene set
results$iterated_over_all_genes                        # Approach 2 final gene set
results$iterated_use_final_gene_set_from_approach_one  # Approach 3 final gene set
results$final_df_app_one                               # Clustered samples dataframe (Approach 1)
results$final_df_app_two                               # Clustered samples dataframe (Approach 2)
results$final_df_app_three                             # Clustered samples dataframe (Approach 3)
results$iter_log                                       # Iteration log (all approaches)
results$pacmap_settings                                # PaCMAP parameters used
results$removed_genes                                  # Genes removed during filtering
```

---

## Output Files

Results are saved to `working_dir/CIM_states_results_{alg}_pacmap_{approach}/` and include:

| File | Description |
|------|-------------|
| `clustered_samples_{approach}.csv` | Final cluster assignments with UMAP, PaCMAP, and PCA embeddings |
| `final_cluster_metrics_summary.csv` | PAC, silhouette, and CCS for each *k* tested, with ranks and optimal *k* selection |
| `final_cluster_stability_metrics.csv` | Per-sample item consensus scores at optimal *k* |
| `final_wide_cluster_stability_metrics.csv` | Wide-format per-sample item consensus scores |
| `clustering_gene_set.csv` | Final discriminatory gene set used for clustering |
| `iter_log_all_approaches.csv` | Gene counts before and after each iteration across all approaches |
| `all_gene_sets.csv` | Input gene sets after filtering |
| `removed_genes.csv` | Genes removed (zero variance, not in matrix, or immune variable genes) |
| `pacmap_settings.csv` | PaCMAP parameters used |
| `item_consensus_k{k}.csv` | Per-sample item consensus scores for each *k* tested |
| `final_pac_plot.png` | PAC curve across *k* |
| `final_cdf_plot.png` | Consensus CDF curves across *k* |
| `sil_plot.png` | Silhouette comparison (dim-reduce and gene-expression space) across *k* |
| `sil_plot_original_space.png` | Silhouette in dim-reduce space only |
| `sil_plot_ge_space.png` | Silhouette in gene-expression space only |
| `item_con_plot.png` | Cluster consensus bar plot across *k* |
| `item_con_sample_k-{k}.png` | Per-sample item consensus bar plot for each *k* |
| `CCP_k_{k}.png` | Consensus matrix heatmap for each *k* |
| `final_fviz_2D_PCA_k-{k}.png` | 2D PCA cluster visualization |
| `final_clustering_shown_on_pacmap_k-{k}.png` | PaCMAP visualization with cluster hulls |
| `final_clustering_shown_on_umap_k-{k}.png` | UMAP visualization with cluster hulls |
| `final_3D_PCA_k-{k}.html` | Interactive 3D PCA plot (requires pandoc; skipped gracefully if unavailable) |

---

## Contact

Mohammed Gbadamosi — mgbadamosi@ufl.edu  
University of Florida