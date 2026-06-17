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

## How It Works

CIMIC proceeds in three stages:

**Stage 1 — Per-pathway iterative clustering and feature selection**  
For each of the 19 CIM pathways, ΔGE matrices are:
1. Reduced using PaCMAP dimensionality reduction (v0.8.0)
2. Clustered using ConsensusClusterPlus (v1.73.0) with hierarchical clustering
3. Evaluated for optimal *k* using a composite ranking framework incorporating:
   - Average silhouette coefficient (s_avg) in both dimensional reduction and ΔGE space
   - Proportion of Ambiguous Clustering (PAC)
   - Cluster Consensus Score (CCS)
4. Filtered to retain only genes showing significant ΔGE differences across clusters 
   (Wilcoxon rank-sum, FDR < 0.05)

Steps 2–4 are repeated iteratively until convergence (no further gene exclusion), 
typically within 3–5 iterations (maximum: 30).

**Stage 2 — Pathway-refined global gene set construction**  
Discriminatory genes from all pathways are pooled into a pathway-refined global gene set.

**Stage 3 — Final CIMIC classification**  
The same iterative procedure is applied to the global gene set. Final cluster 
assignments at convergence define the CIMIC trajectories.

## CIMIC Trajectories

In the breast cancer neoadjuvant chemotherapy cohort, CIMIC identifies two 
reproducible trajectories:
- **Fun-CIM**: Functional, immune-activating trajectory
- **Dys-CIM**: Dysfunctional, stress-adaptive trajectory

## Requirements

### R version
R ≥ 4.5.0

### R packages (installed automatically on first run)
**CRAN:**
`umap`, `cluster`, `factoextra`, `dplyr`, `stringr`, `magrittr`, `tibble`,
`rstatix`, `coin`, `RColorBrewer`, `ggforce`, `grid`, `gridExtra`, `corrplot`,
`concaveman`, `paran`, `reticulate`, `plotly`, `scales`, `withr`, `ggrepel`,
`ggplot2`, `tidyr`, `msigdbr`, `jsonlite`, `htmlwidgets`

**Bioconductor:**
`ConsensusClusterPlus`, `ComplexHeatmap`

### Python (for PaCMAP)
- Anaconda or Miniconda
- The pipeline will automatically set up a `pacmap_env` conda environment on first run
- Alternatively, set `RETICULATE_CONDA` environment variable to your conda path

## Input Data Format

The primary input is a **clustering matrix**:
- Rows: samples
- Columns: genes (gene symbols)
- Values: ΔGE (log₂[TPM+1] post − log₂[TPM+1] pre)

```r
# Example structure
dim(clustering_matrix)   # samples x genes
rownames(clustering_matrix)  # sample IDs
colnames(clustering_matrix)  # gene symbols
