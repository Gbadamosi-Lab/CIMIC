# CIMIC User Guide & Technical Reference

This guide contains the full technical reference for **CIMIC** (`cimic_releases/v1/CIMIC_1.0.0.R`), including mode selection, `block_metadata` configuration, full argument details, pipeline internals, worked examples, output file reference, and troubleshooting.

For a quick overview and minimal example, see the [README](../README.md).

**Contents:** [The two modes](#the-two-modes) · [Requirements](#requirements) · [Input data](#input-data) · [Replicate structure (`block_metadata`)](#replicate-structure-block_metadata) · [Full usage & arguments](#full-usage--arguments) · [How it works](#how-it-works) · [Interpreting output](#interpreting-cimic-output) · [Worked examples](#worked-examples) · [Output files](#output-files) · [Troubleshooting](#troubleshooting) · [Release notes](#release-notes--v100)

---

## The two modes

Both modes share the entire pipeline — embedding, consensus clustering, *k* selection, iterative refinement, and outputs. They differ only in how the per-gene model treats samples that are repeated measurements of the same biological unit.

| | **Independent-samples mode** | **Replicate-aware mode** |
|---|---|---|
| **Set** | `block_metadata = NULL` *(default)* | `block_metadata = <table>` |
| **Assumes** | Every row is an independent observation | Rows within a block are correlated |
| **Model** | `lmFit` → `eBayes` (moderated *t* / *F*) | `duplicateCorrelation` → `lmFit(block=, correlation=)` → `eBayes` |
| **`test` label in output** | `limma_2grp(mod_t)` / `limma_multi(mod_F)` | `limma_dupcor_2grp(mod_t)` / `limma_dupcor_multi(mod_F)` |
| **Use for** | Patient cohorts with one sample each | Cell lines in triplicate; repeated biopsies or timepoints per patient |

### Choosing a mode

The question is not "patients or cell lines" — it is **how many rows come from the same biological unit**:

| Your design | Mode | `block` is |
|---|---|---|
| One pre/post pair per patient | Independent-samples | — |
| 9 cell lines × 3 replicates | Replicate-aware | cell line |
| 20 patients × 2 timepoints each | Replicate-aware | patient |
| One sample per donor, several donors | Independent-samples | — |
| Technical replicates of the same library | Replicate-aware | library/sample |

If you are unsure, count: if no two rows share a unit, use independent-samples mode.

> [!IMPORTANT]
> Using replicate-aware mode on data with no replicates is **safe and gives the identical answer** — with one sample per block, `duplicateCorrelation` is skipped and the fit reduces exactly to unblocked limma. The reverse is not safe: using independent-samples mode on replicate data is pseudoreplication and will produce anticonservative p-values.

### Relationship to the earlier release

| Script | Status |
|---|---|
| `cimic_releases/v1/CIMIC_1.0.0.R` | **Current — use this.** Both modes in one file. |
| `cimic_releases/v1/CIMIC_limma_1.0.0.R` | *Superseded.* Earlier patient-only release, retained for provenance of published results. |

Independent-samples mode reproduces `CIMIC_limma_1.0.0.R` **exactly** — verified on the published cohorts at *k* = 2, 3 and 4: identical p-values, identical effect sizes, identical selected genes, and the same `test` label. Existing patient analyses therefore need no re-run when migrating to `CIMIC_1.0.0.R`.

---

## Requirements

- **R** — developed and tested on R 4.5.0.
- **R packages** — installed automatically on first `source()`. To pre-install manually:
  - *CRAN:* `umap`, `cluster`, `factoextra`, `dplyr`, `stringr`, `magrittr`, `tibble`, `rstatix`, `coin`, `RColorBrewer`, `ggforce`, `grid`, `gridExtra`, `corrplot`, `concaveman`, `paran`, `reticulate`, `plotly`, `scales`, `withr`, `ggrepel`, `ggplot2`, `tidyr`, `msigdbr`, `jsonlite`, `pheatmap`, `htmlwidgets`
  - *Bioconductor:* `ConsensusClusterPlus`, `ComplexHeatmap`, `limma`
- **Python (for PaCMAP)** — requires Anaconda or Miniconda. On first run the pipeline creates a `pacmap_env` conda environment and installs `pacmap` into it. If conda is not auto-detected you will be prompted for its path; to avoid the prompt, set it once:
  ```r
  Sys.setenv(RETICULATE_CONDA = "path/to/conda.exe")
  ```
- **Optional:** `pandoc`, for the self-contained interactive 3D PCA HTML. Without it, that one output is skipped and the run continues.

---

## Input data

A **clustering matrix**: numeric, samples in rows, genes in columns.

- **Rows:** sample IDs (`rownames` must be set — they are how samples are tracked)
- **Columns:** gene symbols
- **Values:** ΔGE (log₂[TPM+1] post − log₂[TPM+1] pre)

```r
clustering_matrix <- as.matrix(read.csv(input_csv, row.names = 1, check.names = FALSE))
storage.mode(clustering_matrix) <- "double"
```

For replicate-level data, supply **one row per replicate** — do **not** pre-average across replicates — and declare the structure via `block_metadata`.

Genes are dropped automatically when they have zero variance, are absent from the matrix, or (with `remove_immune_variable_genes = TRUE`, the default) are immunoglobulin/T-cell-receptor variable-region genes whose expression reflects clonal composition rather than regulation. All removals are logged to `removed_genes.csv`.

---

## Replicate structure (`block_metadata`)

This is the one argument specific to your experimental design: **which rows are repeated measurements of the same cell line or subject.**

| Your data | Pass |
|---|---|
| One sample per patient/subject | `block_metadata = NULL` *(default)* |
| Several replicates per cell line/subject | a table mapping sample → block |

**A "block" is the unit that repeated measurements share** — a cell line, a patient, a donor, a library. It is *not* specific to cell lines; for a study with two biopsies per patient, the block is the patient.

**Why it matters.** Replicates are neither averaged away (which loses *n* and the within-block variance) nor treated as independent (pseudoreplication). Instead, the within-block correlation is estimated with limma's `duplicateCorrelation` and passed to `lmFit` as a random block effect.

With `NULL`, every sample is its own block, `duplicateCorrelation` is skipped, and the fit reduces to plain unblocked limma — numerically identical to `CIMIC_limma_1.0.0.R`.

### Accepted forms

```r
# 1. data.frame (any column order; extra columns ignored)
block_metadata = data.frame(sample_id = rownames(mat), block = cell_line_per_row)

# 2. path to a two-column CSV
block_metadata = "path/to/replicate_metadata.csv"

# 3. named vector, sample id -> block
block_metadata = c(delta_BT549_EPI_R1 = "BT549", delta_BT549_EPI_R2 = "BT549")

# 4. NULL -> independent-samples mode (the default)
block_metadata = NULL

# 5. "auto" -> derive by stripping a trailing "_1" / "_R1" / "_rep1" from the IDs
block_metadata = "auto"
```

> [!CAUTION]
> `"auto"` infers a statistical model from a naming convention, so it is opt-in rather than the default. It is only as good as your IDs — it can fail in both directions:
>
> | Sample IDs | `"auto"` gives | Correct? |
> |---|---|---|
> | `delta_BT549_EPI_R1` … (3 lines × 3) | 3 blocks of 3 | ✅ |
> | `sample_1` … `sample_12` (12 independent) | **1 block of 12** | ❌ over-merged |
> | `BT549-1`, `BT549-2`, … (dash, 2 reps × 2 lines) | **4 blocks of 1** | ❌ replicates missed |
> | `BT549_A`, `BT549_B`, … (letter suffix) | **4 blocks of 1** | ❌ replicates missed |
>
> Always read the printed block table when using `"auto"`. Declaring the structure explicitly (forms 1–3) is what the published runners do.

A table looks like this — a shipped example is [tnbc_cl_epi_replicate_metadata.csv](../oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epi_replicate_metadata.csv):

| sample_id | block |
|---|---|
| `delta_BT549_EPI_R1` | BT549 |
| `delta_BT549_EPI_R2` | BT549 |
| `delta_BT549_EPI_R3` | BT549 |
| `delta_DU4475_EPI_R1` | DU4475 |

Column names are detected automatically (sample id from `sample_id`/`sample`/`id`/…, block from `block`/`base_id`/`cell_line`/`subject`/…), and a plain two-column table is read positionally. Rows are matched to matrix rownames by sample id, so row order and extra rows or columns are irrelevant. **A sample missing from the table is a hard error**, never a silently mis-blocked model.

### Verify the structure before trusting a run

Every run prints the resolved structure once, before clustering starts:

```text
# replicate-aware mode
Replicate structure (declared in block_metadata): 9 block(s) over 27 samples, block sizes 3.

# independent-samples mode
Replicate structure: none declared - all 36 samples treated as independent
(identical to unblocked limma).
```

**Read that line.** If it reports one block per sample when your samples share cell lines or subjects, the block labels are wrong and the replicate model is doing nothing. When blocking is requested and the result is all singletons, CIMIC raises an immediate warning.

### What blocking is *not* for

`duplicateCorrelation` is designed for repeated measurements of the same unit — many small blocks. Do **not** use `block_metadata` to absorb batch, site, or cohort effects (e.g., a two-cohort merged dataset): a handful of very large blocks is not the intended regime, and a fixed-effect covariate in the design is the standard way to handle batch.

---

## Full usage & arguments

A full call with every argument shown explicitly. The **only** line that differs between the two modes is `block_metadata`.

```r
source("cimic_releases/v1/CIMIC_1.0.0.R")

clustering_matrix <- as.matrix(read.csv(input_csv, row.names = 1, check.names = FALSE))
storage.mode(clustering_matrix) <- "double"

results <- CIMIC(
  clustering_matrix  = clustering_matrix,     # samples x genes
  block_metadata     = NULL,                  # NULL = independent; or data.frame/CSV/named vector
  all_gene_sets      = NULL,                  # NULL = the 19 default CIM pathways
  clustering_alg     = "hc",
  max_k              = 5,
  CCP_iter           = 5000,
  adj_pval_thresh    = 0.05,
  max_pipeline_iter  = 50,
  seed               = 2026L,
  clustering_metrics = c("pac", "silhouette_combined_avg", "item_cluster_consensus"),
  filter_approach    = c("app_one", "app_two", "app_three"),
  verbose            = TRUE,
  working_dir        = "path/to/output"
)
saveRDS(results, file.path("path/to/output", "CIMIC_result.rds"))
```

### Arguments

| Argument | Default | Description |
|---|---|---|
| `clustering_matrix` | *required* | Numeric matrix, samples x genes, with rownames |
| `clustering_alg` | *required* | `"hc"` (hierarchical) or `"km"` (k-means) |
| `seed` | *required* | PaCMAP `random_state` (see [Reproducibility](#reproducibility)) |
| `working_dir` | *required* | Output directory — **see warning below** |
| `block_metadata` | `NULL` | Replicate structure; see [Replicate structure](#replicate-structure-block_metadata) |
| `all_gene_sets` | `NULL` | Named list of gene-symbol vectors. `NULL` uses the 19 default CIM pathways (13 GO:BP, 2 Hallmark, 2 Reactome, 2 WikiPathways) fetched via `msigdbr` |
| `remove_immune_variable_genes` | `TRUE` | Drop immunoglobulin/TCR variable-region genes |
| `max_k` | `5` | Maximum *k* tested |
| `CCP_iter` | `5000` | ConsensusClusterPlus resampling iterations; lower (e.g., 500) for a quick test |
| `adj_pval_thresh` | `0.05` | FDR threshold for retaining a gene |
| `max_pipeline_iter` | `50` | Iteration cap for the refinement loop |
| `clustering_metrics` | `c("pac", "silhouette_combined_avg", "item_cluster_consensus")` | Metrics used to rank *k*; the **last** one breaks ties. Valid: `pac`, `silhouette_dim_reduce_space`, `silhouette_ge_space`, `silhouette_combined_avg`, `item_cluster_consensus` |
| `filter_approach` | all three | Any of `app_one`, `app_two`, `app_three`; drop some to save time (`app_three` requires `app_one`) |
| `pacmap_dimensions` | `2L` | PaCMAP components used for clustering |
| `pacmap_args` | `NULL` | Named list passed to PaCMAP, e.g., `list(n_neighbors = 7)` |
| `pacmap_guardrails` | `TRUE` | Warn about and drop unrecognized `pacmap_args` |
| `select_fn` | `cimic_select_features` | The per-gene test; see [Custom statistical test](#custom-statistical-test) |
| `verbose` | `TRUE` | Progress messages |

> [!WARNING]
> **`working_dir` subdirectories are deleted and recreated at the start of every run** (`CIM_states_results_{alg}_pacmap_{approach}/`). Copy anything you want to keep before rerunning.

### Reproducibility

The `seed` argument sets PaCMAP's `random_state` only. The global RNG, ConsensusClusterPlus, and UMAP seeds are fixed internally at `2024L`, so PaCMAP is the only stochastic component you control — deliberate, so embedding stability can be assessed independently.

### Public functions

Sourcing the release defines four public functions (internal helpers are prefixed with `.`):

| Function | Purpose |
|---|---|
| `CIMIC()` | The pipeline. The only function most users call. |
| `cimic_select_features()` | The per-gene test; replaceable via `select_fn` |
| `setup_conda()` | PaCMAP / conda bootstrap, called automatically |
| `get_stability_category()` | Maps a consensus score to a qualitative label |

### Custom statistical test

The per-gene test is isolated in one function so it can be replaced via `select_fn`. A replacement must accept `(expr_mat, clusters, adj_pval_thresh)` and return one row per gene, **in the same order as `colnames(expr_mat)`**, with columns `gene_id`, `effect_size`, `p_value`, `adj_p_value`, `test`. Return results *unfiltered* — CIMIC applies the FDR threshold itself. Handle the fewer-than-two-clusters case and wrap the body in `tryCatch`.

```r
my_wilcox_test <- function(expr_mat, clusters, adj_pval_thresh) {
  genes   <- colnames(expr_mat)
  cluster <- as.factor(clusters)
  na_df <- data.frame(gene_id = genes, effect_size = NA_real_, p_value = NA_real_,
                      adj_p_value = NA_real_, test = "NA", stringsAsFactors = FALSE)
  if (nlevels(cluster) < 2) return(na_df)

  p <- apply(expr_mat, 2, function(g) tryCatch(
    if (nlevels(cluster) == 2) stats::wilcox.test(g ~ cluster)$p.value
    else stats::kruskal.test(g ~ cluster)$p.value,
    error = function(e) NA_real_))
  eff <- apply(expr_mat, 2, function(g) diff(range(tapply(g, cluster, median))))

  data.frame(gene_id = genes, effect_size = as.numeric(eff), p_value = as.numeric(p),
             adj_p_value = p.adjust(p, method = "BH"),
             test = if (nlevels(cluster) == 2) "wilcox" else "kruskal",
             stringsAsFactors = FALSE)
}

results <- CIMIC(..., select_fn = my_wilcox_test)
```

Add a fourth argument named `block_metadata` if your test needs the replicate structure. CIMIC passes it only to functions that declare it.

### Return value

```r
results$iterated_by_gene_sets                          # Approach 1 final gene set
results$iterated_over_all_genes                        # Approach 2 final gene set
results$iterated_use_final_gene_set_from_approach_one  # Approach 3 final gene set
results$final_df_app_one                               # Clustered samples (Approach 1)
results$final_df_app_two                               # Clustered samples (Approach 2)
results$final_df_app_three                             # Clustered samples (Approach 3)
results$iter_log                                       # Iteration log, all approaches
results$pacmap_settings                                # PaCMAP parameters used
results$removed_genes                                  # Genes removed during filtering
```

Approaches not requested return `"NA"`.

---

## How it works

**Stage 1 — Per-pathway iterative clustering and feature selection.** For each of the 19 CIM pathways, the ΔGE matrix is:

1. Reduced with PaCMAP.
2. Clustered with ConsensusClusterPlus (hierarchical or k-means; adaptive `pItem` to resolve NA item-consensus values).
3. Evaluated for optimal *k* by a composite rank across average silhouette (in both embedding and ΔGE space), Proportion of Ambiguous Clustering (PAC), and Cluster Consensus Score (CCS).
4. Filtered to genes with significant ΔGE differences across clusters via a **moderated linear model (`limma`)**:
   - 2 clusters → moderated *t*; effect size = |*t*|
   - ≥3 clusters → moderated *F*; effect size = *F*
   - Benjamini–Hochberg FDR; retained if *p* ≤ 0.05 **and** adj. *p* ≤ `adj_pval_thresh`
   - In **replicate-aware mode**, `duplicateCorrelation` estimates within-block correlation and `lmFit` receives it as a random block effect.

Steps 2–4 repeat until convergence (gene set stops changing), capped at `max_pipeline_iter`.

**Stage 2 — Pathway-refined global gene set.** Discriminatory genes from all pathways are pooled.

**Stage 3 — Final classification.** The same iterative procedure is applied to the pooled set; cluster assignments at convergence define the CIMIC trajectories.

### Filtering approaches

| Approach | Description |
|---|---|
| `app_one` | Refine each gene set independently, then union the survivors (Stage 1 only) |
| `app_two` | Refine the pooled union of all gene sets at once |
| `app_three` | Refine `app_one`'s union again over the combined set (Stages 1–3; requires `app_one`) |

### Convergence handling

- **Normal:** stops when the gene count is unchanged between iterations.
- **Non-immediate convergence to zero:** reverts to the most recent non-zero gene set.
- **Immediate convergence to zero** (first iteration): retains the top 10% of genes by adjusted *p*-value then effect size, and advises relaxing `adj_pval_thresh`.
- **No stable *k*:** if no *k* yields NA-free item consensus at any `pItem`, the gene set is skipped (first iteration) or the most-refined set is kept (later).

### CIMIC trajectories

In the breast cancer neoadjuvant cohort, CIMIC identifies two reproducible trajectories:
- **Fun-CIM** — functional, immune-activating
- **Dys-CIM** — dysfunctional, stress-adaptive

---

## Interpreting CIMIC output

CIMIC is an unsupervised clustering pipeline: it groups samples into candidate CIM trajectories using iteratively refined, cluster-informative gene subsets. It does **not** perform confirmatory statistical testing. Its internal per-gene statistics (p-values, adjusted p-values, effect sizes) are *selection criteria* used to refine gene sets for clustering — they are not, and must not be reported as, differential expression results.

After running CIMIC, you should:
1. **Independently test** whether genes driving a cluster are differentially induced between clusters using a standard DE framework.
2. **Characterize clusters biologically** (pathway enrichment, marker genes, functional annotation).
3. **Validate against data not used in clustering** (clinical variables, outcomes, orthogonal datasets, independent cohorts).

---

## Worked examples

Three standalone scripts reproduce the manuscript analyses. Edit the paths at the top of each script before running.

| Data type | Script | Mode | `block_metadata` |
|---|---|---|---|
| Patient cohort (36 samples) | [nki_smc_combine_cimic.R](../oncoimmunology_paper/Datasets/NKI_SMC/nki_smc_combine_cimic.R) | Independent-samples | `NULL` |
| Patient cohort (19 samples) | [neo_cimic.R](../oncoimmunology_paper/Datasets/NEO/neo_cimic.R) | Independent-samples | `NULL` |
| Cell lines (9 lines × 3 replicates) | [tnbc_cl_cimic.R](../oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_cimic.R) | Replicate-aware | 2-column CSV |

All three source `CIMIC_1.0.0.R` with identical settings (`clustering_alg = "hc"`, `max_k = 5`, `CCP_iter = 5000`, `adj_pval_thresh = 0.05`, `max_pipeline_iter = 50`, `seed = 2026L`, all three approaches) and differ substantively only in `block_metadata`. Each writes `<output_dir>/CIMIC_result.rds`.

---

## Output files

Written to `working_dir/CIM_states_results_{alg}_pacmap_{approach}/`, one folder per requested approach.

**Tables**

| File | Description |
|---|---|
| `clustered_samples_{approach}.csv` | Final cluster assignments with UMAP, PaCMAP, and PCA embeddings |
| `final_cluster_metrics_summary.csv` | PAC, silhouette, CCS per *k*, with ranks and optimal-*k* selection |
| `final_cluster_stability_metrics.csv` | Per-sample item consensus at optimal *k* |
| `final_wide_cluster_stability_metrics.csv` | The same, wide format |
| `clustering_gene_set.csv` | Final discriminatory gene set used for clustering |
| `iter_log_all_approaches.csv` | Gene counts before/after each iteration, all approaches |
| `all_gene_sets.csv` | Input gene sets after filtering |
| `removed_genes.csv` | Genes removed (zero variance, absent, or immune-variable) |
| `pacmap_settings.csv` | PaCMAP parameters used |
| `item_consensus_k{k}.csv` | Per-sample item consensus for each *k* |
| `initial_clustering_mat.csv` | The input matrix, as used |

**Figures**

| File | Description |
|---|---|
| `final_pac_plot.png` | PAC curve across *k* |
| `final_cdf_plot.png` | Consensus CDF curves across *k* |
| `sil_plot.png` | Silhouette in both spaces across *k* |
| `sil_plot_original_space.png` / `sil_plot_ge_space.png` | Silhouette, one space each |
| `item_con_plot.png` | Cluster consensus bar plot across *k* |
| `item_con_sample_k-{k}.png` | Per-sample item consensus per *k* |
| `CCP_k_{k}.png` | Consensus matrix heatmap per *k* |
| `final_fviz_2D_PCA_k-{k}.png` | 2D PCA cluster visualization |
| `final_clustering_shown_on_pacmap_k-{k}.png` | PaCMAP embedding with cluster hulls |
| `final_clustering_shown_on_umap_k-{k}.png` | UMAP embedding with cluster hulls |
| `final_3D_PCA_k-{k}.html` | Interactive 3D PCA (requires pandoc; skipped if unavailable) |

---

## Troubleshooting

- **`Sample size is smaller than the total number of assigned points…`** — A harmless PaCMAP notice for ~21–35 samples when `n_neighbors` is left at its default (10). Setting `pacmap_args = list(n_neighbors = …)` silences it but changes the embedding.
- **`Every block came out a singleton`** — You asked for blocking, but the labels gave one sample per block. Check the `block` column in `block_metadata`.
- **Conda not detected** — Set `RETICULATE_CONDA` (see [Requirements](#requirements)) or let the pipeline install Miniconda when prompted.
- **Empty final gene set** — The FDR filter removed everything; relax `adj_pval_thresh`. On immediate convergence to zero, CIMIC retains the top 10% of genes automatically.

---

## Release notes — v1.0.0

- **One release, two modes.** `CIMIC_1.0.0.R` replaces the previous pair of scripts. Independent-samples mode reproduces `CIMIC_limma_1.0.0.R` exactly (verified: identical p-values, effect sizes, selected genes, and `test` label).
- **Replicate structure is declared, not inferred.** Passed explicitly via `block_metadata` (`NULL` = no replicates; `"auto"` is opt-in).
- **`duplicateCorrelation` is skipped when there is no replication**, making independent-samples mode faster and preserving the plain-limma `test` label.
- **Resolved block structure is reported once per run**, with a warning if blocking was requested but produced an all-singleton factor.
