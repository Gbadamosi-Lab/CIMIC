# CIMIC: Chemoimmunomodulation Induction Classifier

CIMIC is an unsupervised, stability-guided iterative pipeline for stratifying cancer samples by chemotherapy-induced chemoimmunomodulatory (CIM) transcriptional trajectories. It operates on within-sample delta gene expression ($\Delta\text{GE} = \text{post} - \text{pre}$ $\log_2[\text{TPM}+1]$) across 3,189 genes spanning 19 CIM-related pathways, identifying two reproducible trajectories in breast cancer: **Fun-CIM** (immune-activating) and **Dys-CIM** (stress-adaptive).

This repository contains the CIMIC release (`cimic_releases/v1/CIMIC_1.0.0.R`) and the datasets and scripts for the *Oncoimmunology* manuscript.

> **Looking for full documentation?** See the **[User Guide](docs/USER_GUIDE.md)** for complete arguments, `block_metadata` details, pipeline internals, worked examples, output file reference, troubleshooting, and release notes.

---

## Quick start

```r
source("cimic_releases/v1/CIMIC_1.0.0.R")   # installs missing R packages on first run

# 1. Independent-samples mode (one measurement per patient/subject)
results <- CIMIC(
  clustering_matrix = my_delta_ge_matrix,   # samples x genes, rownames = sample IDs
  block_metadata    = NULL,                 # no replicate structure
  clustering_alg    = "hc",
  seed              = 2026L,
  working_dir       = "path/to/output"
)

# 2. Replicate-aware mode (e.g., cell lines in triplicate)
results <- CIMIC(
  clustering_matrix = my_replicate_delta_ge_matrix,   # one row per replicate
  block_metadata    = data.frame(
    sample_id = rownames(my_replicate_delta_ge_matrix),
    block     = cell_line_per_row
  ),
  clustering_alg    = "hc",
  seed              = 2026L,
  working_dir       = "path/to/output"
)
```

### Choosing a mode

| Design | Mode | `block_metadata` |
|---|---|---|
| One pre/post pair per patient | Independent-samples | `NULL` *(default)* |
| Repeated measurements (e.g., cell line replicates) | Replicate-aware | Data frame, CSV path, or named vector mapping `sample_id` → `block` |

> [!WARNING]
> `working_dir` subdirectories are deleted and recreated at the start of every run. Copy anything you want to keep before rerunning.

---

## Requirements

- **R** (tested on 4.5.0) — CRAN and Bioconductor packages are installed automatically on first `source()`.
- **Python (Conda/Miniconda)** — used for PaCMAP embedding (a `pacmap_env` environment is created automatically on first run; set `Sys.setenv(RETICULATE_CONDA = "path/to/conda.exe")` if not auto-detected).
- **Optional:** `pandoc` for interactive 3D PCA HTML output.

---

## Input & output at a glance

- **Input:** A numeric $\Delta\text{GE}$ matrix (samples in rows with `rownames`, gene symbols in columns). Do **not** pre-average replicates; pass them as separate rows with `block_metadata`.
- **Output:** Written to `working_dir/CIM_states_results_{alg}_pacmap_{approach}/`, including cluster assignments (`clustered_samples_*.csv`), stability metrics, refined gene sets, and UMAP/PaCMAP/PCA plots.
- **Note:** CIMIC is **hypothesis-generating** — its internal `limma` statistics are feature-selection criteria for clustering, not confirmatory differential expression results.

---

## Documentation

Extensive notes and technical details are in **[docs/USER_GUIDE.md](docs/USER_GUIDE.md)**:

- [Mode selection & `CIMIC_limma_1.0.0.R` equivalence](docs/USER_GUIDE.md#the-two-modes)
- [Replicate structure (`block_metadata`) & verification](docs/USER_GUIDE.md#replicate-structure-block_metadata)
- [Full argument list, custom statistical tests (`select_fn`), & return values](docs/USER_GUIDE.md#full-usage--arguments)
- [How the 3-stage iterative pipeline works](docs/USER_GUIDE.md#how-it-works)
- [Worked manuscript examples (NKI/SMC, NEO, TNBC cell lines)](docs/USER_GUIDE.md#worked-examples)
- [Output files, troubleshooting, & v1.0.0 release notes](docs/USER_GUIDE.md#output-files)

---

## Repository layout

```text
cimic_releases/v1/     CIMIC source releases (use CIMIC_1.0.0.R)
docs/                  Extended user guide and technical reference
oncoimmunology_paper/  Manuscript datasets, runner scripts, and published results
```

---

## Citation, License & Contact

- **Citation:** If you use CIMIC, please cite the original publication in *Oncoimmunology*.
- **License:** MIT — see [LICENSE](LICENSE).
- **Contact:** Mohammed Gbadamosi — mgbadamosi@ufl.edu (University of Florida)
