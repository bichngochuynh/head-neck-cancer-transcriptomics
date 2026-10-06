# Characterising Immune-hot and Immune-cold Tumour Microenvironment in Head and Neck Squamous Cell Carcinoma: Integrating TCGA Bulk RNA-Sequencing Immune Deconvolution with Single-Cell Transcriptomic Profiling

![R](https://img.shields.io/badge/R-4.6-276DC3?logo=r&logoColor=white)
![Seurat](https://img.shields.io/badge/Seurat-scRNA--seq-6A5ACD)
![DESeq2](https://img.shields.io/badge/Bioconductor-DESeq2-1F8DD6)
![renv](https://img.shields.io/badge/reproducible-renv-2E8B57)

*Analysis code for my Master's thesis.* This repository contains the R
pipeline only. No data, results or figures are included.

## Abstract

Head and neck squamous cell carcinoma (HNSCC) is a diverse cancer with large differences in its tumour immune microenvironment. HPV-positive tumours often have better clinical outcomes than HPV-negative tumours, but HPV status alone does not fully explain immune variation between patients. This study therefore focused first on immune phenotype, rather than HPV status alone. Using TCGA-HNSC bulk RNA-seq data, immune infiltration was estimated with MCP-counter, ESTIMATE and xCell. Tumours were classified into immune-hot and immune-cold groups. In total, 139 tumours were classified as immune-hot and 381 as immune-cold. Immune-hot tumours showed stronger immune-related signals, including T-cell, B-cell, plasma-cell and antigen-presentation programmes. Immune-cold tumours showed weaker immune activation and relatively stronger metabolic, ribosomal and structural gene-expression programmes.

HPV-positive tumours were significantly enriched in the immune-hot group, but HPV status and immune phenotype were not identical. Many immune-hot tumours were HPV-negative, showing that immune phenotype provides information beyond HPV status. Survival analysis showed that HPV-positive status was associated with better overall survival, while immune-hot/cold status alone was not significantly associated with overall survival.

Single-cell RNA-seq data from GSE181919 were then used to interpret the bulk findings at cell-type level. HPV-associated gene-expression differences were mainly linked to epithelial cells, while B/plasma cells showed the strongest immune-cell-associated signal. Overall, this study shows that HPV status and immune phenotype are related but distinct features of HNSCC biology. Classifying tumours by immune-hot and immune-cold states may help describe tumour immune heterogeneity and guide future studies, but further validation is needed before clinical use.

## Workflow

```mermaid
flowchart TB
    subgraph BULK["Stage 1 · Bulk RNA-seq · TCGA-HNSC"]
        direction LR
        A[GDC STAR counts<br/>TCGAbiolinks] --> B[Immune deconvolution<br/>MCP-counter · xCell · ESTIMATE]
        B --> C[k-means k=2<br/>Immune-hot / cold]
        C --> D[DESeq2 + apeglm<br/>GSEA / ORA · clusterProfiler]
        C --> E[HPV association<br/>Kaplan–Meier · Cox PH]
    end
    subgraph SC["Stage 2 · Single-cell RNA-seq · GEO GSE181919"]
        direction LR
        F[UMI counts] --> G[Per-tissue MAD QC<br/>SCTransform]
        G --> H[PCA · Louvain · UMAP]
        H --> I[SingleR + marker panels<br/>cell-type annotation]
        I --> J[Composition by HPV<br/>pseudobulk DE · ORA / GSEA]
    end
    BULK ==> SC
```

## Repository structure

```
.
├── data/                       # README with public accessions (no data committed)
├── bulk_seq_analysis/          # Stage 1: bulk RNA-seq, run in numeric order
│   ├── 01_download_HNSC_TCGAbiolinks.R          # GDC download (STAR counts)
│   ├── 02_import_HNSC_clinical_2018.R           # clinical / HPV / survival import
│   ├── 03_run_immune_deconvolution.R            # MCP-counter, xCell, ESTIMATE
│   ├── 04_plot_immune_classification.R          # k-means immune-hot / cold
│   ├── 05_plot_immune_scores_by_tool.R          # validation across tools
│   ├── 06_run_de_pathway_immune_hot_vs_cold.R   # DESeq2 + GO/KEGG ORA and GSEA
│   ├── 07_visualise_de_immune_hot_vs_cold.R     # volcano, heatmaps
│   ├── 08_correlate_genes_with_immunescore.R    # gene–ImmuneScore correlation
│   └── 09_clinical_association_immune_hot_vs_cold.R  # Fisher, KM, Cox PH
├── scRNA_analysis/             # Stage 2: single-cell RNA-seq
│   └── 01_run_scRNAseq_pipeline.R   # QC → SCTransform → clustering → annotation →
│                                    # HPV composition → pseudobulk DE → ORA / GSEA
└── renv.lock                   # pinned package versions
```

## Methods implemented

### Stage 1: Bulk RNA-seq (TCGA-HNSC)

| Step | Implementation | Rationale |
|---|---|---|
| Immune deconvolution | MCP-counter, xCell and ESTIMATE (immunedeconv) | Three independent estimates of immune infiltration |
| Immune groups | k-means on MCP-counter immune scores; ESTIMATE only names clusters; xCell held out | Keeps an independent method for validation |
| Differential expression | DESeq2 Wald test with apeglm LFC shrinkage | Shrinkage stabilises fold changes for low-count genes |
| Pathways | clusterProfiler ORA and GSEA (GO BP, KEGG) | GSEA catches coordinated sub-threshold shifts that ORA misses |
| Survival | Kaplan–Meier, log-rank, uni- and multivariable Cox PH | Tests whether the immune phenotype is prognostic independently of HPV |

### Stage 2: Single-cell RNA-seq

| Step | Implementation | Rationale |
|---|---|---|
| QC | 3 × MAD on log10 counts / genes / mito%, computed **per tissue type** | Library complexity differs by tissue, so one global cutoff would over- or under-filter |
| Normalisation | SCTransform v2 (glmGamPoi), regressing mito% | Removes sequencing-depth and mitochondrial technical variance |
| Clustering | PCA → shared-nearest-neighbour graph → Louvain → UMAP | PCs chosen from the elbow plot; checkpoints saved after expensive steps |
| Annotation | SingleR (celldex reference) + canonical marker panels | Combines reference-based labels with marker validation |
| Differential expression | Pseudobulk DESeq2 (counts summed per patient per cell type), then ORA / GSEA | Avoids pseudo-replication from treating cells of one patient as independent |

## Data

All input data are public and are **not** included in this repository. See
[`data/README.md`](data/README.md) for accessions and download instructions
(TCGA-HNSC via GDC/cBioPortal; GEO
[GSE181919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181919)).

## Running the pipeline

```r
# 1. Clone, open the .Rproj in RStudio, then restore the pinned packages
renv::restore()

# 2. Bulk RNA-seq track (downloads TCGA data on first run)
source("bulk_seq_analysis/01_download_HNSC_TCGAbiolinks.R")
# ... run 02 → 09 in order

# 3. Single-cell track (place the GSE181919 *.gz files in the project root first)
source("scRNA_analysis/01_run_scRNAseq_pipeline.R")
```

- All paths are **relative to the project root**. No machine-specific
  absolute paths are used.
- Package versions are pinned with [`renv`](https://rstudio.github.io/renv/).
- Each script creates its own output folders under `outputs/`, which is not
  tracked in the repository.
