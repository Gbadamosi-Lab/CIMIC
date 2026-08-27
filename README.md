# CIMIC: Chemoimmunomodulation Induction Classifier

CIMIC is an unsupervised, stability-guided iterative pipeline for stratifying
cancer patients by chemotherapy-induced chemoimmunomodulatory (CIM)
transcriptional trajectories. It was developed and validated in neoadjuvant
chemotherapy for breast cancer.

The pipeline operates on within-patient delta gene expression
(ΔGE = post-treatment log₂[TPM+1] − pre-treatment log₂[TPM+1]) across 3,179 genes
spanning 19 CIM-related pathways, capturing treatment-induced transcriptional
responses while controlling for interpatient baseline differences.

This repository contains the CIMIC source releases plus the data and code used to
produce the figures in the original *Oncoimmunology* publication.

**Contents:** [Quick start](#quick-start) · [Requirements](#requirements) ·
[Input data](#input-data) · [Replicate structure](#replicate-structure-block_metadata) ·
[Usage](#usage) · [How it works](#how-it-works) ·
[Interpreting output](#interpreting-cimic-output) ·
[Worked examples](#worked-examples) · [Output files](#output-files)

---

## Quick start

```r
source("cimic_releases/v1/CIMIC_1.0.0.R")   # installs missing packages on first run

results <- CIMIC(
  clustering_matrix = my_delta_ge_matrix,   # samples x genes, rownames = sample IDs
  block_metadata    = NULL,                 # NULL = no replicates (see below)
  clustering_alg    = "hc",
  seed              = 2026L,
  working_dir       = "path/to/output"
)
```

`clustering_matrix`, `clustering_alg`, `seed`, and `working_dir` are required;
everything else has a default. See [Usage](#usage) for the full argument list.

### Which release?

| Script | Status |
|---|---|
| `cimic_releases/v1/CIMIC_1.0.0.R` | **Current — use this.** Handles data with or without replicates. |
| `cimic_releases/v1/CIMIC_limma_1.0.0.R` | *Superseded.* Earlier patient-only release, retained for provenance of published results. Reproduced exactly by `CIMIC_1.0.0.R` with `block_metadata = NULL`. |

---

## Requirements

**R** — developed and tested on R 4.5.0.

**R packages** are installed automatically on first `source()`. If you would
rather control this yourself, install them before sourcing:

- *CRAN:* `umap`, `cluster`, `factoextra`, `dplyr`, `stringr`, `magrittr`,
  `tibble`, `rstatix`, `coin`, `RColorBrewer`, `ggforce`, `grid`, `gridExtra`,
  `corrplot`, `concaveman`, `paran`, `reticulate`, `plotly`, `scales`, `withr`,
  `ggrepel`, `ggplot2`, `tidyr`, `msigdbr`, `jsonlite`, `pheatmap`, `htmlwidgets`
- *Bioconductor:* `ConsensusClusterPlus`, `ComplexHeatmap`, `limma`

**Python (for PaCMAP)** — requires Anaconda or Miniconda. On first run the
pipeline creates a `pacmap_env` conda environment and installs `pacmap` into it.
If conda is not auto-detected you will be prompted for its path; to avoid the
prompt, set it once:

```r
Sys.setenv(RETICULATE_CONDA = "path/to/conda.exe")
```

**Optional:** pandoc, for the self-contained interactive 3D PCA HTML. Without it
that one output is skipped and the run continues.

---

## Input data

A **clustering matrix**: numeric, samples in rows, genes in columns.

- Rows: sample IDs (`rownames` must be set — they are how samples are tracked)
- Columns: gene symbols
- Values: ΔGE (log₂[TPM+1] post − log₂[TPM+1] pre)

```r
clustering_matrix <- as.matrix(read.csv(input_csv, row.names = 1, check.names = FALSE))
storage.mode(clustering_matrix) <- "double"
```

For replicate-level data, supply **one row per replicate** — do **not**
pre-average across replicates — and declare the structure via `block_metadata`.

Genes are dropped automatically when they have zero variance, are absent from the
matrix, or (with `remove_immune_variable_genes = TRUE`, the default) are
immunoglobulin/T-cell-receptor variable-region genes whose expression reflects
clonal composition rather than regulation. All removals are logged to
`removed_genes.csv`.

---

## Replicate structure (`block_metadata`)

This is the one argument specific to your experimental design: **which rows are
repeated measurements of the same cell line or subject.**

| Your data | Pass |
|---|---|
| One sample per patient/subject | `block_metadata = NULL` *(default)* |
| Several replicates per cell line/subject | a table mapping sample → block |

**Why it matters.** Replicates are neither averaged away (which loses *n* and the
within-block variance) nor treated as independent (pseudoreplication). Instead
the within-block correlation is estimated with limma's `duplicateCorrelation` and
passed to `lmFit` as a random block effect.

With `NULL`, every sample is its own block, `duplicateCorrelation` is skipped, and
the fit reduces to plain unblocked limma — numerically identical to
`CIMIC_limma_1.0.0.R` (same p-values, same selected genes, same `test` label).

**Accepted forms:**

```r
# data.frame (any column order; extra columns ignored)
block_metadata = data.frame(sample_id = rownames(mat), block = cell_line_per_row)

# path to a two-column CSV
block_metadata = "path/to/replicate_metadata.csv"

# named vector
block_metadata = c(delta_BT549_EPI_R1 = "BT549", delta_BT549_EPI_R2 = "BT549", ...)

# derive from a trailing "_1" / "_R1" / "_rep1" in the sample IDs (opt-in)
block_metadata = "auto"
```

A table looks like this — a shipped example is
[tnbc_cl_epi_replicate_metadata.csv](oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_epi_replicate_metadata.csv):

| sample_id | block |
|---|---|
| `delta_BT549_EPI_R1` | BT549 |
| `delta_BT549_EPI_R2` | BT549 |
| `delta_BT549_EPI_R3` | BT549 |
| `delta_DU4475_EPI_R1` | DU4475 |

Column names are detected automatically (sample id from
`sample_id`/`sample`/`id`/…, block from `block`/`base_id`/`cell_line`/`subject`/…),
and a plain two-column table is read positionally, so no configuration is needed.
Rows are matched to the matrix rownames by sample id, so row order and extra rows
or columns are irrelevant. **A sample missing from the table is a hard error**,
never a silently mis-blocked model.

Every run prints the resolved structure once before starting:

```
Replicate structure (declared in block_metadata): 9 block(s) over 27 samples, block sizes 3.
```

**Check that line.** If it reports one block per sample when your samples do
share cell lines, the block labels are wrong and the replicate model is doing
nothing. `"auto"` additionally warns if it produces an all-singleton factor.

---

## Usage

```r
source("cimic_releases/v1/CIMIC_1.0.0.R")

results <- CIMIC(
  clustering_matrix  = my_delta_ge_matrix,
  block_metadata     = NULL,
  clustering_alg     = "hc",
  max_k              = 5,
  CCP_iter           = 5000,
  adj_pval_thresh    = 0.05,
  max_pipeline_iter  = 50,
  seed               = 2026L,
  clustering_metrics = c("pac", "silhouette_combined_avg", "item_cluster_consensus"),
  filter_approach    = c("app_one", "app_two", "app_three"),
  working_dir        = "path/to/output"
)
```

### Arguments

| Argument | Default | Description |
|---|---|---|
| `clustering_matrix` | *required* | Numeric matrix, samples x genes, with rownames |
| `clustering_alg` | *required* | `"hc"` (hierarchical) or `"km"` (k-means) |
| `seed` | *required* | PaCMAP `random_state` (see [Reproducibility](#reproducibility)) |
| `working_dir` | *required* | Output directory — **see the warning below** |
| `block_metadata` | `NULL` | Replicate structure; see [above](#replicate-structure-block_metadata) |
| `all_gene_sets` | `NULL` | Named list of gene-symbol vectors. `NULL` uses the 19 default CIM pathways (13 GO:BP, 2 Hallmark, 2 Reactome, 2 WikiPathways) fetched via `msigdbr` |
| `remove_immune_variable_genes` | `TRUE` | Drop immunoglobulin/TCR variable-region genes |
| `max_k` | `5` | Maximum *k* tested |
| `CCP_iter` | `5000` | ConsensusClusterPlus resampling iterations; lower (e.g. 500) for a quick test |
| `adj_pval_thresh` | `0.05` | FDR threshold for retaining a gene |
| `max_pipeline_iter` | `50` | Iteration cap for the refinement loop |
| `clustering_metrics` | `c("pac", "silhouette_combined_avg", "item_cluster_consensus")` | Metrics used to rank *k*; the **last** one breaks ties. Valid: `pac`, `silhouette_dim_reduce_space`, `silhouette_ge_space`, `silhouette_combined_avg`, `item_cluster_consensus` |
| `filter_approach` | all three | Any of `app_one`, `app_two`, `app_three`; drop some to save time (`app_three` requires `app_one`) |
| `pacmap_dimensions` | `2L` | PaCMAP components used for clustering |
| `pacmap_args` | `NULL` | Named list passed to PaCMAP, e.g. `list(n_neighbors = 7)` |
| `pacmap_guardrails` | `TRUE` | Warn about and drop unrecognized `pacmap_args` |
| `select_fn` | `cimic_select_features` | The per-gene test; see [Custom test](#custom-statistical-test) |
| `verbose` | `TRUE` | Progress messages |

> [!WARNING]
> **`working_dir` subdirectories are deleted and recreated at the start of every
> run** (`CIM_states_results_{alg}_pacmap_{approach}/`). Copy anything you want to
> keep before rerunning.

### Reproducibility

The `seed` argument sets PaCMAP's `random_state` only. The global RNG,
ConsensusClusterPlus, and UMAP seeds are fixed internally at `2024L`, so PaCMAP
is the only stochastic component you control — this was deliberate, so that
embedding stability could be assessed independently (see the ARI robustness
analysis in the manuscript).

### Custom statistical test

The per-gene test is isolated in one function so it can be replaced without
touching the pipeline. A replacement must accept
`(expr_mat, clusters, adj_pval_thresh)` and return one row per gene, **in the same
order as `colnames(expr_mat)`**, with columns
`gene_id`, `effect_size`, `p_value`, `adj_p_value`, `test`. Return results
*unfiltered* — CIMIC applies the FDR threshold itself. Handle the
fewer-than-two-clusters case and wrap the body in `tryCatch`; the function is
called hundreds of times per run and an error aborts everything.

```r
results <- CIMIC(..., select_fn = my_own_test)
```

Add a fourth argument named `block_metadata` if your test needs the replicate
structure; CIMIC passes it only to functions that declare it, so three-argument
tests keep working unchanged.

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

Approaches not requested return the string `"NA"`.

---

## How it works

**Stage 1 — Per-pathway iterative clustering and feature selection.** For each of
the 19 CIM pathways, the ΔGE matrix is:

1. Reduced with PaCMAP
2. Clustered with ConsensusClusterPlus (hierarchical or k-means; adaptive `pItem`
   to resolve NA item-consensus values)
3. Evaluated for optimal *k* by a composite rank across average silhouette (in
   both embedding and ΔGE space), Proportion of Ambiguous Clustering (PAC), and
   Cluster Consensus Score (CCS)
4. Filtered to genes with significant ΔGE differences across clusters, via a
   **moderated linear model (limma)**:
   - 2 clusters → moderated *t*; effect size = |*t*|
   - ≥3 clusters → moderated *F*; effect size = *F*
   - Benjamini–Hochberg FDR; retained if *p* ≤ 0.05 **and** adj. *p* ≤ `adj_pval_thresh`
   - With a declared replicate structure, `duplicateCorrelation` estimates the
     within-block correlation first and `lmFit` receives it as a random block
     effect before empirical-Bayes moderation

Steps 2–4 repeat until convergence — the surviving gene set stops changing —
capped at `max_pipeline_iter`.

**Stage 2 — Pathway-refined global gene set.** Discriminatory genes from all
pathways are pooled.

**Stage 3 — Final classification.** The same iterative procedure is applied to
the pooled set; cluster assignments at convergence define the CIMIC trajectories.

### Filtering approaches

| Approach | Description |
|---|---|
| `app_one` | Refine each gene set independently, then union the survivors (Stage 1 only) |
| `app_two` | Refine the pooled union of all gene sets at once |
| `app_three` | Refine `app_one`'s union again over the combined set (Stages 1–3; requires `app_one`) |

### Convergence handling

- **Normal:** stops when the gene count is unchanged between iterations.
- **Non-immediate convergence to zero:** reverts to the most recent non-zero gene set.
- **Immediate convergence to zero** (first iteration): retains the top 10% of
  genes by adjusted *p*-value then effect size, and advises relaxing
  `adj_pval_thresh`.
- **No stable *k*:** if no *k* yields NA-free item consensus at any `pItem`, the
  gene set is skipped (first iteration) or the most-refined set is kept (later).

### CIMIC trajectories

In the breast cancer neoadjuvant cohort, CIMIC identifies two reproducible
trajectories:

- **Fun-CIM** — functional, immune-activating
- **Dys-CIM** — dysfunctional, stress-adaptive

---

## Interpreting CIMIC output

CIMIC is an unsupervised clustering pipeline: it groups samples into candidate CIM
trajectories using iteratively refined, cluster-informative gene subsets. It does
**not** perform confirmatory statistical testing. Its internal per-gene statistics
(p-values, adjusted p-values, effect sizes) are *selection criteria* used to refine
gene sets for clustering — they are not, and must not be reported as,
differential expression results.

CIMIC output is hypothesis-generating. After running it you should:

1. **Independently test** whether genes driving a cluster are differentially
   induced between clusters, using a standard DE framework appropriate to your
   data — not the pipeline's internal selection statistics.
2. **Characterize clusters biologically** (pathway enrichment, marker genes,
   functional annotation) to confirm the partition reflects a coherent biological
   process rather than an artifact of feature refinement.
3. **Validate against data not used in clustering** — clinical variables,
   outcomes, orthogonal datasets, independent cohorts. This is the strongest
   evidence that a trajectory is real and not a description of the algorithm's
   own optimization objective.

Every partition CIMIC returns should be treated as a hypothesis to be tested, not
a validated finding.

---

## Worked examples

Three standalone scripts reproduce the manuscript analyses. Each is
self-contained: open it, edit the paths at the top, run it.

| Data type | Script | `block_metadata` |
|---|---|---|
| Patient cohort (36 samples) | [nki_smc_combine_cimic.R](oncoimmunology_paper/Datasets/NKI_SMC/nki_smc_combine_cimic.R) | `NULL` |
| Patient cohort (19 samples) | [neo_cimic.R](oncoimmunology_paper/Datasets/NEO/neo_cimic.R) | `NULL` |
| Cell lines (9 lines x 3 replicates) | [tnbc_cl_cimic.R](oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_cimic.R) | 2-column CSV |

**Patient cohorts.** Combined NKI + SMC breast-cancer cohort (one pre/post pair
per patient) and the NEO cohort. One measurement per patient, so
`block_metadata = NULL` and no blocking factor is needed.

- Inputs: `nki_smc_combine_initial_clustering_mat.csv`, `neo_initial_clustering_mat.csv`
- Outputs: `oncoimmunology_paper/Results/{nki_smc,neo}_cimic_results/`

**Cell lines.** TNBC panel treated with epirubicin, 9 cell lines x 3 biological
replicates = 27 samples. Because each line contributes several replicates,
`block_metadata` points at a mapping table so `duplicateCorrelation` models the
within-line correlation instead of treating replicates as independent.

- Input: `tnbc_cl_epi_initial_clustering_mat.csv` (replicate-level, **not** pre-averaged)
- Blocking: `tnbc_cl_epi_replicate_metadata.csv`
- Output: `oncoimmunology_paper/Results/tnbc_cl_epirubicin_cimic_results/`

**What differs between them.** All three source the same release and call `CIMIC()`
with identical settings — `clustering_alg = "hc"`, `max_k = 5`,
`CCP_iter = 5000`, `adj_pval_thresh = 0.05`, `max_pipeline_iter = 50`,
`seed = 2026L`, all three approaches. The only substantive difference is
`block_metadata`. That choice, driven purely by whether the data contain
replicates, is the whole of what separates the analyses.

Each script writes the full result object to `<output_dir>/CIMIC_result.rds`
alongside the per-approach folders described below.

> [!NOTE]
> These scripts hard-code absolute paths for `script_dir`, `input_csv`,
> `output_dir`, and (where applicable) the blocking CSV. Edit them for your
> machine, or replace them with paths relative to the repository root.

---

## Output files

Written to `working_dir/CIM_states_results_{alg}_pacmap_{approach}/`, one folder
per requested approach.

**Tables**

| File | Description |
|---|---|
| `clustered_samples_{approach}.csv` | Final cluster assignments with UMAP, PaCMAP, and PCA embeddings |
| `final_cluster_metrics_summary.csv` | PAC, silhouette, CCS per *k*, with ranks and the optimal-*k* selection |
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

## Repository layout

```
cimic_releases/v1/         CIMIC source releases
oncoimmunology_paper/
  Datasets/                Input matrices, blocking metadata, and runner scripts
  Results/                 Published pipeline outputs
```

---

## Troubleshooting

**`Sample size is smaller than the total number of assigned points…`** — a PaCMAP
notice, not an error. When `n_neighbors` is left at PaCMAP's default of 10, the
requested neighbour/mid-near/far pairs exceed what a small sample can supply, so
PaCMAP scales them down and continues. It appears for datasets of roughly 21–35
samples and is harmless. Setting `pacmap_args = list(n_neighbors = …)` silences it
but changes the embedding, and therefore your results.

**`Every block came out a singleton`** — you asked for blocking but the labels
gave one sample per block, so the replicate model is doing nothing. Check the
block column in your `block_metadata`.

**Conda not detected** — set `RETICULATE_CONDA` (see [Requirements](#requirements))
or let the pipeline install Miniconda when prompted.

**Empty final gene set** — the FDR filter removed everything; relax
`adj_pval_thresh`. On immediate convergence to zero CIMIC retains the top 10% of
genes automatically and says so.

---

## Citation

If you use CIMIC, please cite the original publication in *Oncoimmunology*.

<!-- TODO: replace with the full citation and DOI once assigned. -->

## License

MIT — see [LICENSE](LICENSE).

## Contact

Mohammed Gbadamosi — mgbadamosi@ufl.edu  
University of Florida
