# CIMIC: Chemoimmunomodulation Induction Classifier

CIMIC is an unsupervised, stability-guided iterative pipeline for stratifying
cancer patients by chemotherapy-induced chemoimmunomodulatory (CIM)
transcriptional trajectories. It was developed and validated in neoadjuvant
chemotherapy for breast cancer.

The pipeline operates on within-patient delta gene expression
(ΔGE = post-treatment log₂[TPM+1] − pre-treatment log₂[TPM+1]) across 3,189 genes
spanning 19 CIM-related pathways, capturing treatment-induced transcriptional
responses while controlling for interpatient baseline differences.

This repository contains the CIMIC source releases plus the data and code used to
produce the figures in the original *Oncoimmunology* publication.

**Contents:** [Quick start](#quick-start) · [The two modes](#the-two-modes) ·
[Requirements](#requirements) · [Input data](#input-data) ·
[Replicate structure](#replicate-structure-block_metadata) · [Usage](#usage) ·
[How it works](#how-it-works) ·
[Interpreting output](#interpreting-cimic-output) ·
[Worked examples](#worked-examples) · [Output files](#output-files) ·
[Troubleshooting](#troubleshooting) · [Release notes](#release-notes--v100)

---

## Quick start

There is **one** release, `cimic_releases/v1/CIMIC_1.0.0.R`, which runs in one of
two modes. The mode is chosen entirely by the `block_metadata` argument.

```r
source("cimic_releases/v1/CIMIC_1.0.0.R")   # installs missing packages on first run

# --- Independent-samples mode: one measurement per patient/subject -------------
results <- CIMIC(
  clustering_matrix = my_delta_ge_matrix,   # samples x genes, rownames = sample IDs
  block_metadata    = NULL,                 # <- no replicate structure
  clustering_alg    = "hc",
  seed              = 2026L,
  working_dir       = "path/to/output"
)

# --- Replicate-aware mode: several replicates per cell line/subject ------------
results <- CIMIC(
  clustering_matrix = my_replicate_delta_ge_matrix,   # one row per replicate
  block_metadata    = data.frame(                     # <- declare the structure
                        sample_id = rownames(my_replicate_delta_ge_matrix),
                        block     = cell_line_per_row),
  clustering_alg    = "hc",
  seed              = 2026L,
  working_dir       = "path/to/output"
)
```

`clustering_matrix`, `clustering_alg`, `seed`, and `working_dir` are required;
everything else has a default. See [Usage](#usage) for the full argument list.

---

## The two modes

Both modes share the entire pipeline — embedding, consensus clustering, *k*
selection, iterative refinement, outputs. They differ only in how the per-gene
model treats samples that are repeated measurements of the same biological unit.

| | **Independent-samples mode** | **Replicate-aware mode** |
|---|---|---|
| **Set** | `block_metadata = NULL` *(default)* | `block_metadata = <table>` |
| **Assumes** | every row is an independent observation | rows within a block are correlated |
| **Model** | `lmFit` → `eBayes` (moderated *t* / *F*) | `duplicateCorrelation` → `lmFit(block=, correlation=)` → `eBayes` |
| **`test` label in output** | `limma_2grp(mod_t)` / `limma_multi(mod_F)` | `limma_dupcor_2grp(mod_t)` / `limma_dupcor_multi(mod_F)` |
| **Use for** | patient cohorts with one sample each | cell lines in triplicate; repeated biopsies or timepoints per patient |

### Choosing a mode

The question is not "patients or cell lines" — it is **how many rows come from
the same biological unit**:

| Your design | Mode | `block` is |
|---|---|---|
| One pre/post pair per patient | Independent-samples | — |
| 9 cell lines × 3 replicates | Replicate-aware | cell line |
| 20 patients × 2 timepoints each | Replicate-aware | patient |
| One sample per donor, several donors | Independent-samples | — |
| Technical replicates of the same library | Replicate-aware | library/sample |

If you are unsure, count: if no two rows share a unit, use independent-samples
mode.

> [!IMPORTANT]
> Using replicate-aware mode on data with no replicates is **safe and gives the
> identical answer** — with one sample per block, `duplicateCorrelation` is
> skipped and the fit reduces exactly to unblocked limma. The reverse is not
> safe: using independent-samples mode on replicate data is pseudoreplication and
> will produce anticonservative p-values.

### Relationship to the earlier release

| Script | Status |
|---|---|
| `cimic_releases/v1/CIMIC_1.0.0.R` | **Current — use this.** Both modes in one file. |
| `cimic_releases/v1/CIMIC_limma_1.0.0.R` | *Superseded.* Earlier patient-only release, retained for provenance of published results. |

Independent-samples mode reproduces `CIMIC_limma_1.0.0.R` **exactly** — verified
on the published cohorts at *k* = 2, 3 and 4: identical p-values, identical
effect sizes, identical selected genes, and the same `test` label. Existing
patient analyses therefore need no re-run when migrating to `CIMIC_1.0.0.R`.

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

**A "block" is the unit that repeated measurements share** — a cell line, a
patient, a donor, a library. It is *not* specific to cell lines; for a study with
two biopsies per patient, the block is the patient.

**Why it matters.** Replicates are neither averaged away (which loses *n* and the
within-block variance) nor treated as independent (pseudoreplication). Instead
the within-block correlation is estimated with limma's `duplicateCorrelation` and
passed to `lmFit` as a random block effect.

With `NULL`, every sample is its own block, `duplicateCorrelation` is skipped, and
the fit reduces to plain unblocked limma — numerically identical to
`CIMIC_limma_1.0.0.R` (same p-values, same selected genes, same `test` label).

**Accepted forms:**

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
> `"auto"` infers a statistical model from a naming convention, so it is opt-in
> rather than the default. It is only as good as your IDs — it can fail in both
> directions:
>
> | Sample IDs | `"auto"` gives | Correct? |
> |---|---|---|
> | `delta_BT549_EPI_R1` … (3 lines × 3) | 3 blocks of 3 | ✅ |
> | `sample_1` … `sample_12` (12 independent) | **1 block of 12** | ❌ over-merged |
> | `BT549-1`, `BT549-2`, … (dash, 2 reps × 2 lines) | **4 blocks of 1** | ❌ replicates missed |
> | `BT549_A`, `BT549_B`, … (letter suffix) | **4 blocks of 1** | ❌ replicates missed |
>
> Always read the printed block table when using `"auto"`. Declaring the structure
> explicitly (forms 1–3) is what the published runners do.

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

### Verify the structure before trusting a run

Every run prints the resolved structure once, before the (long) clustering starts:

```
# replicate-aware mode
Replicate structure (declared in block_metadata): 9 block(s) over 27 samples, block sizes 3.

# independent-samples mode
Replicate structure: none declared - all 36 samples treated as independent
(identical to unblocked limma).
```

**Read that line.** If it reports one block per sample when your samples really do
share cell lines or subjects, the block labels are wrong and the replicate model
is doing nothing. When you asked for blocking and the result is all singletons,
CIMIC raises an immediate warning; a sample present in the matrix but missing from
`block_metadata` is a hard error, so a mismatch fails in seconds rather than
part-way through a run.

### What blocking is *not* for

`duplicateCorrelation` is designed for repeated measurements of the same unit —
many small blocks. Do **not** use `block_metadata` to absorb batch, site, or
cohort effects (e.g. a two-cohort merged dataset): a handful of very large blocks
is not the intended regime, and a fixed-effect covariate in the design is the
standard way to handle batch. Blocking answers "which rows are the same thing
measured twice", not "which rows were processed together".

---

## Usage

A full call with every argument shown explicitly. The **only** line that differs
between the two modes is `block_metadata`.

### Independent-samples mode

```r
source("cimic_releases/v1/CIMIC_1.0.0.R")

clustering_matrix <- as.matrix(read.csv(input_csv, row.names = 1, check.names = FALSE))
storage.mode(clustering_matrix) <- "double"

results <- CIMIC(
  clustering_matrix  = clustering_matrix,     # samples x genes
  block_metadata     = NULL,                  # <- one measurement per subject
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

### Replicate-aware mode

```r
source("cimic_releases/v1/CIMIC_1.0.0.R")

clustering_matrix <- as.matrix(read.csv(input_csv, row.names = 1, check.names = FALSE))
storage.mode(clustering_matrix) <- "double"

# one row per replicate; block = the cell line each replicate came from
block_metadata <- data.frame(
  sample_id = rownames(clustering_matrix),
  block     = my_cell_line_per_row,
  stringsAsFactors = FALSE
)
table(block_metadata$block)                   # sanity-check before the long run

results <- CIMIC(
  clustering_matrix  = clustering_matrix,     # replicates x genes
  block_metadata     = block_metadata,        # <- the only changed line
  all_gene_sets      = NULL,
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

### Public functions

Sourcing the release defines four public functions; everything else is internal
and prefixed with a dot (`.cimic_resolve_block`, `.refine_gene_set`, …), which
keeps it out of `ls()`. You should never need to call an internal helper — if a
task seems to require one, it is a gap in the public API.

| Function | Purpose |
|---|---|
| `CIMIC()` | The pipeline. The only function most users call. |
| `cimic_select_features()` | The per-gene test; replaceable via `select_fn` (below) |
| `setup_conda()` | PaCMAP / conda bootstrap, called automatically |
| `get_stability_category()` | Maps a consensus score to a qualitative label |

### Custom statistical test

The per-gene test is isolated in one function so it can be replaced without
touching the rest of the pipeline. A replacement must accept
`(expr_mat, clusters, adj_pval_thresh)` and return one row per gene, **in the same
order as `colnames(expr_mat)`**, with columns
`gene_id`, `effect_size`, `p_value`, `adj_p_value`, `test`. Return results
*unfiltered* — CIMIC applies the FDR threshold itself. Handle the
fewer-than-two-clusters case and wrap the body in `tryCatch`; the function is
called hundreds of times per run and an unhandled error aborts everything.

```r
my_wilcox_test <- function(expr_mat, clusters, adj_pval_thresh) {
  genes   <- colnames(expr_mat)
  cluster <- as.factor(clusters)
  na_df <- data.frame(gene_id = genes, effect_size = NA_real_, p_value = NA_real_,
                      adj_p_value = NA_real_, test = "NA", stringsAsFactors = FALSE)
  if (nlevels(cluster) < 2) return(na_df)          # required: CIMIC reaches this case

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

Add a fourth argument named `block_metadata` if your test needs the replicate
structure. CIMIC passes it only to functions that declare it, so three-argument
tests like the one above keep working unchanged.

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
   - In **replicate-aware mode**, `duplicateCorrelation` estimates the
     within-block correlation first and `lmFit` receives it as a random block
     effect before empirical-Bayes moderation. In **independent-samples mode**
     that step is skipped entirely

   This step — and only this step — differs between the two modes.

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

| Data type | Script | Mode | `block_metadata` |
|---|---|---|---|
| Patient cohort (36 samples) | [nki_smc_combine_cimic.R](oncoimmunology_paper/Datasets/NKI_SMC/nki_smc_combine_cimic.R) | Independent-samples | `NULL` |
| Patient cohort (19 samples) | [neo_cimic.R](oncoimmunology_paper/Datasets/NEO/neo_cimic.R) | Independent-samples | `NULL` |
| Cell lines (9 lines × 3 replicates) | [tnbc_cl_cimic.R](oncoimmunology_paper/Datasets/TNBC_CL_Epirubicin/tnbc_cl_cimic.R) | Replicate-aware | 2-column CSV |

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

**What differs between them.** All three source the same release, `CIMIC_1.0.0.R`,
and call the same entry point, `CIMIC()`, with identical settings —
`clustering_alg = "hc"`, `max_k = 5`, `CCP_iter = 5000`,
`adj_pval_thresh = 0.05`, `max_pipeline_iter = 50`, `seed = 2026L`, all three
approaches. The only substantive difference is `block_metadata`, i.e. the mode.
That choice, driven purely by whether the data contain replicates, is the whole
of what separates the analyses.

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

## Release notes — v1.0.0

- **One release, two modes.** `CIMIC_1.0.0.R` replaces the previous pair of
  scripts. Independent-samples mode reproduces `CIMIC_limma_1.0.0.R` exactly
  (verified: identical p-values, effect sizes, selected genes and `test` label),
  so migrating a patient analysis changes nothing.
- **Replicate structure is now declared, not inferred.** It is passed as the
  explicit `block_metadata` argument to `CIMIC()`. `NULL` means "no replicates",
  and deriving blocks from sample-ID suffixes is opt-in via `"auto"`.
- **`duplicateCorrelation` is skipped when there is no replication**, since limma
  would set the intrablock correlation to zero anyway. This makes
  independent-samples mode faster and keeps the plain-limma `test` label.
- **The resolved block structure is reported once per run**, with a warning if
  blocking was requested but produced an all-singleton factor.

> [!IMPORTANT]
> **If you produced replicate-level results with a pre-1.0.0 script, re-run them.**
> Earlier versions derived the block by stripping a trailing `_<digits>` from the
> sample ID, which does not match the `_R1` / `_R2` / `_R3` convention used by the
> replicate matrices in this repository. Every block came out a singleton, so the
> replicate model silently reduced to an unblocked fit. Patient-level results are
> unaffected.

---

## Citation

If you use CIMIC, please cite the original publication in *Oncoimmunology*.

<!-- TODO: replace with the full citation and DOI once assigned. -->

## License

MIT — see [LICENSE](LICENSE).

## Contact

Mohammed Gbadamosi — mgbadamosi@ufl.edu  
University of Florida
