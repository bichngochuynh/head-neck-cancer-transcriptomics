## Working directory must be the project root (paths below are root-relative);
## this script lives one folder down in scRNA_analysis/, hence the extra dirname()
tryCatch(
  setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path))),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    f    <- sub("--file=", "", args[grep("--file=", args)])
    if (length(f) > 0) setwd(dirname(dirname(normalizePath(f))))
  }
)

# Load libraries
library(Seurat)
library(data.table)
library(Matrix)
library(dplyr)
library(tidyr)
library(R.utils)
library(ggplot2)
if (!requireNamespace("MAST", quietly = TRUE)) {
  BiocManager::install("MAST", update = FALSE, ask = FALSE)
}
library(MAST)
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2", update = FALSE, ask = FALSE)
}
library(DESeq2)

# Pipeline output tree, organised by analysis stage
scrnaseq_out <- "outputs/scRNAseq_outputs"
for (d in c("checkpoints", "01_qc", "02_clustering", "03_marker_validation",
            "04_singler_annotation", "05_final_celltype_annotation",
            "06_composition_by_HPV", "07_exploratory_celltype_DE",
            "08_tcell_subtypes", "09_differential_expression",
            "10_pathway_enrichment", "de_tables")) {
  dir.create(file.path(scrnaseq_out, d), recursive = TRUE, showWarnings = FALSE)
}

# Checkpoint covering data load through SCTransform/clustering/UMAP — the
# most expensive and crash-prone steps. Lets any failure in FindAllMarkers/
# SingleR/annotation below resume without redoing all of this.
post_clustering_checkpoint <- "outputs/scRNAseq_outputs/checkpoints/seurat_post_clustering.rds"

if (file.exists(post_clustering_checkpoint)) {
  
  cat("Checkpoint found — loading post-clustering seurat_obj\n")
  seurat_obj <- readRDS(post_clustering_checkpoint)
  
} else {

## STEP 1 — Load raw counts and metadata
## Why: the GEO deposit ships counts and per-cell metadata as two separate
## flat files. They must be parsed, barcode-matched, and combined into a
## single Seurat object before any QC or analysis is possible.

# 1. Read raw file: First row = barcodes
raw_lines <- readLines("GSE181919_UMI_counts.txt.gz", n = 1)
barcodes <- strsplit(raw_lines, "\t")[[1]]

# Fix barcode format: change ".1" or ".2" to "-1", "-2".
barcodes <- gsub("\\.(\\d+)$", "-\\1", barcodes)

# 2. Read actual expression data (skip first row)
umi_counts_df <- fread("GSE181919_UMI_counts.txt.gz", skip = 1, header = FALSE)

# Extract gene names and convert to matrix
gene_names <- umi_counts_df[[1]]
umi_counts_mat <- as.matrix(umi_counts_df[, -1, with = FALSE])
rownames(umi_counts_mat) <- gene_names
colnames(umi_counts_mat) <- barcodes

# Convert to sparse matrix
umi_counts_sparse <- as(umi_counts_mat, "dgCMatrix")

# Free the dense intermediates (multi-GB) now that the sparse matrix exists —
# leaving them in .GlobalEnv inflates SCTransform's future-globals export later
rm(umi_counts_df, umi_counts_mat)
gc()

# 3. Load metadata
barcode_metadata <- fread("GSE181919_Barcode_metadata.txt.gz", header = FALSE, skip = 1)
colnames(barcode_metadata) <- c("barcodes", "patientid", "sampleid", "gender", "age",
                                "tissuetype", "subsite", "HPVstatus", "celltype")
barcode_metadata <- as.data.frame(barcode_metadata)

# Match barcode suffix format in metadata
barcode_metadata$barcodes <- gsub("\\.(\\d+)$", "-\\1", barcode_metadata$barcodes)
rownames(barcode_metadata) <- barcode_metadata$barcodes

# 4. Intersect barcodes
shared_cells <- intersect(colnames(umi_counts_sparse), rownames(barcode_metadata))
length(shared_cells)

# 5. Subset both datasets
umi_counts_sparse <- umi_counts_sparse[, shared_cells]
barcode_metadata <- barcode_metadata[shared_cells, ]

# 6. Verify everything matches
dim(barcode_metadata)
length(shared_cells)
head(rownames(barcode_metadata))
all(colnames(umi_counts_sparse) == rownames(barcode_metadata))  # Should be TRUE

# These should all be TRUE:
ncol(umi_counts_sparse) == nrow(barcode_metadata)
length(shared_cells) == ncol(umi_counts_sparse)
all(colnames(umi_counts_sparse) == rownames(barcode_metadata))

# Check dimensions
cat("Count matrix:", nrow(umi_counts_sparse), "genes x", ncol(umi_counts_sparse), "cells\n")
cat("Metadata:", nrow(barcode_metadata), "cells\n")
cat("Shared cells:", length(shared_cells), "\n")

# Keep genes that are expressed in at least 3 cells
## Why: genes detected in only 1-2 cells out of ~50k can't support any
## statistical test and just bloat memory/compute for every downstream step.
genes_to_keep <- rowSums(umi_counts_sparse > 0) >= 3
# Subset the sparse matrix to only include those genes
umi_counts_sparse <- umi_counts_sparse[genes_to_keep, ]

# 7. Create Seurat object
seurat_obj <- CreateSeuratObject(counts = umi_counts_sparse, meta.data = barcode_metadata)

# Final check
table(seurat_obj@meta.data$HPVstatus)

## STEP 2 — Quality control (remove low-quality / dying / doublet-like cells)
## Why: raw droplet-based scRNA-seq always contains empty droplets, dying
## cells (high mito%), and library-prep artifacts (too few/too many genes
## detected). Leaving them in would distort normalisation, clustering, and
## every downstream comparison. Thresholds are computed per tissue type
## (3 x MAD on log10 scale) because library size/complexity differs systematically by tissue, 
## so a single global cutoff would over- or under-filter some tissues.

# 1. Calculate mitochondrial gene percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-") # identify dying cells

# View QC metrics
head(seurat_obj@meta.data)

# 2. Visualise raw QC metrics, split by tissue, before filtering
p_qc_raw <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "tissuetype", ncol = 3, pt.size = 0)
ggsave("outputs/scRNAseq_outputs/01_qc/01_QC_metrics_raw.pdf", p_qc_raw, width = 12, height = 5)

seurat_before_qc <- seurat_obj

# 3. Calculate per-tissue MAD-based thresholds (3 x MAD on log10 scale)
qc_metrics <- seurat_obj@meta.data %>%
  mutate(
    log10_nCount   = log10(nCount_RNA),
    log10_nFeature = log10(nFeature_RNA),
    log10_mt       = log10(percent.mt + 1)   # +1 so 0% mito cells don't give -Inf
  )

mad_thresholds <- qc_metrics %>%
  group_by(tissuetype) %>%
  summarise(
    nCount_med   = median(log10_nCount),   nCount_mad   = mad(log10_nCount),
    nFeature_med = median(log10_nFeature), nFeature_mad = mad(log10_nFeature),
    mt_med       = median(log10_mt),       mt_mad       = mad(log10_mt),
    .groups = "drop"
  ) %>%
  mutate(
    nCount_lower   = nCount_med   - 3 * nCount_mad,
    nCount_upper   = nCount_med   + 3 * nCount_mad,
    nFeature_lower = nFeature_med - 3 * nFeature_mad,
    nFeature_upper = nFeature_med + 3 * nFeature_mad,
    # Mito%: only high values indicate dying cells, so upper bound only
    mt_upper       = mt_med + 3 * mt_mad
  )

cat("Per-tissue MAD-based QC thresholds (log10 scale):\n")
print(mad_thresholds)

# 4. Flag outlier cells using their own tissue's thresholds
qc_metrics <- qc_metrics %>%
  left_join(mad_thresholds, by = "tissuetype") %>%
  mutate(
    outlier_nCount   = log10_nCount   < nCount_lower   | log10_nCount   > nCount_upper,
    outlier_nFeature = log10_nFeature < nFeature_lower | log10_nFeature > nFeature_upper,
    outlier_mt       = log10_mt > mt_upper,
    discard          = outlier_nCount | outlier_nFeature | outlier_mt
  )

seurat_obj$discard <- qc_metrics$discard

cat("Outliers flagged per tissue:\n")
qc_metrics %>%
  group_by(tissuetype) %>%
  summarise(
    n_total     = n(),
    n_discard   = sum(discard),
    pct_discard = round(100 * n_discard / n_total, 2),
    .groups = "drop"
  ) %>%
  print()

# 5. Visualise flagged vs. kept cells before removing them
p_qc_flagged <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "tissuetype", split.by = "discard", ncol = 3, pt.size = 0)
ggsave("outputs/scRNAseq_outputs/01_qc/02_QC_metrics_flagged_outliers.pdf", p_qc_flagged, width = 12, height = 5)

# 6. Apply the filter
seurat_obj <- subset(seurat_obj, subset = discard == FALSE)

cells_before <- ncol(seurat_before_qc)
cells_after <- ncol(seurat_obj)

cat("Cells before QC:", cells_before, "\n")
cat("Cells after QC:", cells_after, "\n")
cat("Cells removed:", cells_before - cells_after, "\n")
cat("Percentage removed:", round((cells_before - cells_after) / cells_before * 100, 2), "%\n")

# 7. Confirm filtering
p_qc_filtered <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "tissuetype", ncol = 3, pt.size = 0)
ggsave("outputs/scRNAseq_outputs/01_qc/03_QC_metrics_post_filter.pdf", p_qc_filtered, width = 12, height = 5)

## STEP 3 — Cell cycle scoring
## Why: proliferating cells (tumour/malignant cells) can cluster by cell
## cycle phase rather than by cell identity, which would confound clustering
## and annotation. Scoring S/G2M phase here lets us check for that effect and,
## if needed, regress it out later.

# Cell cycle assessment
seurat_phase <- NormalizeData(seurat_obj)

seurat_phase <- CellCycleScoring(
  seurat_phase,
  s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes
)

seurat_phase <- FindVariableFeatures(seurat_phase)
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)

p_cc_pca <- DimPlot(
  seurat_phase,
  reduction = "pca",
  group.by = "Phase"
)
ggsave("outputs/scRNAseq_outputs/01_qc/04_cellcycle_PCA.pdf", p_cc_pca, width = 7, height = 6)
#Store scores
seurat_obj$S.Score <- seurat_phase$S.Score
seurat_obj$G2M.Score <- seurat_phase$G2M.Score
seurat_obj$Phase <- seurat_phase$Phase

## STEP 4 — SCTransform normalisation
## Why: raw UMI counts can't be compared across cells directly because
## sequencing depth (library size) varies cell-to-cell. SCTransform models
## and removes this technical variance (and the mito% effect, via
## vars.to.regress) using a regularised negative-binomial regression, which
## is the standard normalisation Seurat clustering/UMAP expect.

# SCTransform
options(future.globals.maxSize = 4 * 1024^3)

seurat_obj <- SCTransform(
  seurat_obj,
  vars.to.regress = "percent.mt",
  do.correct.umi = FALSE,
  method = "glmGamPoi",
  vst.flavor = "v2",
  verbose = TRUE
)

## Note:  a manual decision point (look at the 04_cellcycle_PCA plot, 
## decide not to regress out cell cycle).

## STEP 5 — Dimensionality reduction and clustering
## Why: with ~20k genes per cell, distances in full gene-expression space are
## dominated by noise. PCA compresses the data to a handful of components
## that capture the dominant biological variance; FindNeighbors/FindClusters
## then groups cells into clusters using those PCs (graph-based Louvain
## clustering) instead of raw expression. UMAP is purely for 2D visualisation
## of that same neighbour graph — it doesn't affect clustering itself.

# PCA
seurat_obj <- RunPCA(seurat_obj)
p_elbow <- ElbowPlot(seurat_obj)
ggsave("outputs/scRNAseq_outputs/02_clustering/05_ElbowPlot.pdf", p_elbow, width = 7, height = 5)
#Select PCs: 1:20

# Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
p_umap_clusters <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
ggsave("outputs/scRNAseq_outputs/02_clustering/06_UMAP_clusters.pdf", p_umap_clusters, width = 8, height = 7)

saveRDS(seurat_obj, post_clustering_checkpoint)
cat("Saved checkpoint:", post_clustering_checkpoint, "\n")

}

post_annotation_checkpoint <- "outputs/scRNAseq_outputs/checkpoints/seurat_post_annotation.rds"

if (file.exists(post_annotation_checkpoint)) {

  cat("Checkpoint found — loading post-annotation seurat_obj (skips marker validation/SingleR)\n")
  seurat_obj <- readRDS(post_annotation_checkpoint)

} else {

DefaultAssay(seurat_obj) <- "SCT"

## STEP 6 — Marker gene validation per cluster
## Why: clustering only groups cells by similarity — it doesn't say what
## those cells are. FindAllMarkers identifies genes specifically
## upregulated in each cluster, which is the evidence needed to assign a
## biological cell-type label to each cluster in the next steps.

#CELL ANNOTATION
#MARKER GENE VALIDATION PER CLUSTER

library(Seurat)
library(dplyr)
library(ggplot2)


DefaultAssay(seurat_obj) <- "SCT"
Idents(seurat_obj) <- "seurat_clusters"

n_clusters <- length(unique(seurat_obj$seurat_clusters))

if (file.exists("outputs/scRNAseq_outputs/03_marker_validation/markers_all_clusters.csv")) {
  
  cat("Checkpoint found — loading markers_all_clusters.csv\n")
  
  markers_all <- read.csv(
    "outputs/scRNAseq_outputs/03_marker_validation/markers_all_clusters.csv",
    stringsAsFactors = FALSE
  )
  
} else {
  
  markers_all <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    verbose = FALSE
  )
  
  write.csv(
    markers_all,
    "outputs/scRNAseq_outputs/03_marker_validation/markers_all_clusters.csv",
    row.names = FALSE
  )
  
  cat("Saved: outputs/scRNAseq_outputs/03_marker_validation/markers_all_clusters.csv\n")
}

# Top markers per cluster
top5 <- markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup()

write.csv(
  top5,
  "outputs/scRNAseq_outputs/03_marker_validation/markers_top5_per_cluster.csv",
  row.names = FALSE
)

cat("\nTop 5 markers per cluster:\n")

top5 %>%
  group_by(cluster) %>%
  summarise(
    markers = paste(gene, collapse = ", ")
  ) %>%
  print(n = Inf)

# DotPlot of top markers

top3_genes <- markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3) %>%
  pull(gene) %>%
  unique()

pdf(
  "outputs/scRNAseq_outputs/03_marker_validation/S3a_dotplot_top3_markers.pdf",
  width = max(14, length(top3_genes) * 0.35 + 4),
  height = max(6, n_clusters * 0.4 + 3)
)

DotPlot(
  seurat_obj,
  features = top3_genes,
  group.by = "seurat_clusters"
) +
  labs(title = "Top 3 marker genes per cluster") +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 7
    )
  )

dev.off()

cat("Saved: S3a_dotplot_top3_markers.pdf\n")

# Heatmap

top5_genes <- unique(top5$gene)

set.seed(42)

cells_use <- seurat_obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  group_by(seurat_clusters) %>%
  group_modify(~ slice_sample(.x, n = min(200, nrow(.x)))) %>%
  ungroup() %>%
  pull(cell)

pdf(
  "outputs/scRNAseq_outputs/03_marker_validation/S3b_heatmap_top5_markers.pdf",
  width = 14,
  height = max(8, length(top5_genes) * 0.18 + 3)
)

DoHeatmap(
  seurat_obj[, cells_use],
  features = top5_genes,
  group.by = "seurat_clusters"
)

dev.off()

cat("Saved: S3b_heatmap_top5_markers.pdf\n")

# Canonical marker expression

canonical_markers <- c(
  "CD3D",   # T cells
  "CD79A",  # B / plasma cells
  "MS4A1",  # B cells
  "MZB1",   # Plasma cells
  "CD68",   # Macrophages
  "CLEC9A", # Dendritic cells
  "DCN",    # Fibroblasts
  "PECAM1", # Endothelial cells
  "EPCAM",  # Epithelial / malignant
  "TPSAB1", # Mast cells
  "MYH11"   # Myocytes
)

canonical_present <- canonical_markers[
  canonical_markers %in% rownames(seurat_obj)
]

pdf(
  "outputs/scRNAseq_outputs/03_marker_validation/S3c_UMAP_canonical_markers.pdf",
  width = 16,
  height = ceiling(length(canonical_present) / 4) * 4
)

FeaturePlot(
  seurat_obj,
  features = canonical_present,
  ncol = 4,
  pt.size = 0.2,
  order = TRUE
)

dev.off()

png(
  "outputs/scRNAseq_outputs/03_marker_validation/S3c_UMAP_canonical_markers.png",
  width = 16,
  height = ceiling(length(canonical_present) / 4) * 4,
  units = "in",
  res = 300
)

FeaturePlot(
  seurat_obj,
  features = canonical_present,
  ncol = 4,
  pt.size = 0.2,
  order = TRUE
)

dev.off()

cat("Saved: S3c_UMAP_canonical_markers.pdf / .png\n\n")

## STEP 7 — SingleR cross-check
## Why: manual marker-based annotation is somewhat subjective. SingleR
## independently predicts a cell type per cluster by correlating its
## expression profile against a labelled reference atlas (HPCA), giving a
## second, automated opinion to validate (or flag disagreement with) the
## marker-based calls before finalising cell-type labels.

#SINGLER CROSS-CHECK

library(SingleR)
library(celldex)
library(SingleCellExperiment)

if (
  file.exists("outputs/scRNAseq_outputs/checkpoints/singler_results.rds") &
  file.exists("outputs/scRNAseq_outputs/04_singler_annotation/singler_cluster_labels.csv")
) {
  
  cat("Checkpoint found — loading SingleR results\n")
  
  singler_results <- readRDS(
    "outputs/scRNAseq_outputs/checkpoints/singler_results.rds"
  )
  
  cluster_singler_map <- read.csv(
    "outputs/scRNAseq_outputs/04_singler_annotation/singler_cluster_labels.csv",
    stringsAsFactors = FALSE
  )
  
} else {
  
  cat("Loading HumanPrimaryCellAtlasData reference...\n")
  
  ref_hpca <- HumanPrimaryCellAtlasData()
  
  DefaultAssay(seurat_obj) <- "SCT"
  
  sce <- as.SingleCellExperiment(seurat_obj)
  
  clusters_vec <- seurat_obj$seurat_clusters
  
  cat("Running SingleR...\n")
  
  singler_results <- SingleR(
    test = sce,
    ref = ref_hpca,
    labels = ref_hpca$label.main,
    clusters = clusters_vec
  )
  
  cluster_singler_map <- data.frame(
    seurat_clusters = rownames(singler_results),
    singler_label = singler_results$labels,
    singler_score = apply(
      singler_results$scores,
      1,
      max
    ),
    stringsAsFactors = FALSE
  )
  
  saveRDS(
    singler_results,
    "outputs/scRNAseq_outputs/checkpoints/singler_results.rds"
  )
  
  write.csv(
    cluster_singler_map,
    "outputs/scRNAseq_outputs/04_singler_annotation/singler_cluster_labels.csv",
    row.names = FALSE
  )
}

seurat_obj$singler_label <- cluster_singler_map$singler_label[
  match(
    as.character(seurat_obj$seurat_clusters),
    cluster_singler_map$seurat_clusters
  )
]

print(cluster_singler_map)

# UMAP by SingleR

pdf(
  "outputs/scRNAseq_outputs/04_singler_annotation/S4a_UMAP_SingleR_labels.pdf",
  width = 10,
  height = 8
)

print(
  DimPlot(
    seurat_obj,
    group.by = "singler_label",
    label = TRUE,
    repel = TRUE
  )
)

dev.off()

cat("Saved: S4a_UMAP_SingleR_labels.pdf\n")

# SingleR score heatmap


pdf(
  "outputs/scRNAseq_outputs/04_singler_annotation/S4b_SingleR_score_heatmap.pdf",
  width = 10,
  height = 8
)

plotScoreHeatmap(singler_results)

dev.off()

cat("Saved: S4b_SingleR_score_heatmap.pdf\n\n")

# COMPARE MARKERS VS SINGLER

cluster_summary <- data.frame(
  Cluster = cluster_singler_map$seurat_clusters,
  SingleR = cluster_singler_map$singler_label,
  Score = round(cluster_singler_map$singler_score, 3)
)

print(cluster_summary)

write.csv(
  cluster_summary,
  "outputs/scRNAseq_outputs/04_singler_annotation/cluster_SingleR_summary.csv",
  row.names = FALSE
)

cat("Saved: cluster_SingleR_summary.csv\n")

cat("\nUse:\n")
cat("1. Top marker genes\n")
cat("2. DotPlot\n")
cat("3. Heatmap\n")
cat("4. Canonical marker FeaturePlots\n")
cat("5. SingleR predictions\n")
cat("to manually assign final cell identities.\n\n")

# Checkpoint before manual cluster annotation — lets PART 2 resume without
# rerunning QC/SCTransform/clustering/FindAllMarkers/SingleR
saveRDS(seurat_obj, post_annotation_checkpoint)
cat("Saved checkpoint:", post_annotation_checkpoint, "\n")
cat("n_clusters =", n_clusters, "\n")

}

## STEP 8 — Final cell-type annotation
## Why: this is where the marker genes (Step 6) and SingleR predictions
## (Step 7) get synthesised into one human-decided label per cluster. It has
## to be a manual mapping (cluster number -> cell type) because no automated
## method can be fully trusted alone for a tumour microenvironment with
## atypical cell states (e.g. malignant cells, CAFs) absent from reference
## atlases.

#FINAL CELL TYPE ANNOTATION

# EDIT THIS SECTION AFTER REVIEWING MARKERS + SINGLER
## Note: the genes listed per cluster below are that cluster's actual top-5
## FindAllMarkers hits (Step 6), not a literature panel chosen in advance —
## only the cell-type label is the literature/SingleR-informed call.

cluster_annotation <- c(
  "T.cells",          # 0  - GNLY, GZMH, KRT86, FASLG, NKG7 (cytotoxic)
  "Fibroblasts",       # 1  - PLA2G2A, SFRP1, GDF10
  "T.cells",          # 2  - FOXP3, TNFRSF4, IL2RA, LAIR2 (Treg)
  "T.cells",          # 3  - IL7R, CD69
  "B_Plasma.cells",   # 4  - MS4A1, BANK1, CD22, VPREB3
  "Fibroblasts",       # 5  - COL11A1, COL10A1, POSTN (CAF)
  "Fibroblasts",       # 6  - PRG4, PCOLCE2, ADAMTS5
  "Endothelial",       # 7  - CCL14, SOX17, SELE, ACKR1
  "Dendritic.cells",  # 8  - FCN1, CFP, CD1C, CPVL (myeloid/cDC)
  "Macrophages",       # 9  - C1QC, C1QB, C1QA, CCL18
  "T.cells",          # 10 - GZMK, GZMM, CD8A, TRAT1 (CD8)
  "B_Plasma.cells",   # 11 - IGHGP, IGHG1-4 (IgG plasma)
  "Epithelial",        # 12 - RHOXF1-AS1, MAGEA1, KRT19 (malignant)
  "Epithelial",        # 13 - LGALS7B, KLK5, COL4A6
  "Fibroblasts",       # 14 - SCN7A, MYOC, SMOC2 (Schwann-like marker overlap; called CAF)
  "Dendritic.cells",  # 15 - CLEC4C, LILRA4, PTCRA (pDC)
  "Myocytes",          # 16 - RERGL, MYH11, PLN
  "B_Plasma.cells",   # 17 - IGHA2, IGHA1, IGHM, JCHAIN (IgA plasma)
  "Epithelial",        # 18 - DUOXA2, CEACAM7, FLG
  "Epithelial",        # 19 - BPIFB1, PIGR, MSLN, CAPS
  "Endothelial",       # 20 - CCL21, MMRN1, LYVE1 (lymphatic)
  "Myocytes",          # 21 - MYF5, PAX7, MYBPHL
  "Myocytes",          # 22 - LMOD2, CSRP3, MYH2
  "Mast.cells"         # 23 - TPSAB1, TPSB2, CPA3, CTSG
)

names(cluster_annotation) <- levels(seurat_obj)

seurat_obj$celltype_final <- unname(cluster_annotation[
  as.character(seurat_obj$seurat_clusters)
])

table(seurat_obj$celltype_final)

pdf(
  "outputs/scRNAseq_outputs/05_final_celltype_annotation/S6_final_celltype_annotation.pdf",
  width = 10,
  height = 8
)

print(
  DimPlot(
    seurat_obj,
    group.by = "celltype_final",
    label = TRUE,
    repel = TRUE
  )
)

dev.off()

cat("Saved: S6_final_celltype_annotation.pdf\n")

## STEP 9 — Validate the final annotation with marker violin plots
## Why: a sanity check that each assigned cell type actually expresses its
## defining canonical markers (e.g. CD3D for T cells) and not some other
## cell type's markers, catching mislabelled clusters before they propagate
## into every downstream comparison.

# VIOLIN PLOTS — MARKER GENES BY CLUSTER-DERIVED CELLTYPE (celltype_final)

dir.create("outputs/scRNAseq_outputs/05_final_celltype_annotation/violin_celltype_final", recursive = TRUE, showWarnings = FALSE)

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(seurat_obj) <- "celltype_final"

marker_sets_final <- list(
  T.cells         = c("CD3D", "CD3E", "CD8A", "CD4"),
  B_Plasma.cells  = c("CD79A", "MS4A1", "IGHG1", "MZB1"),
  Macrophages     = c("CD68", "CD14", "LYZ"),
  Dendritic.cells = c("FCER1A", "CLEC10A", "CD1C"),
  Fibroblasts     = c("COL1A1", "COL1A2", "FAP"),
  Endothelial     = c("PECAM1", "VWF", "CDH5"),
  Mast.cells      = c("TPSAB1", "CPA3", "KIT"),
  Epithelial      = c("EPCAM", "KRT5", "KRT14", "KRT17"),
  Myocytes        = c("ACTA1", "MYH11")
)

for (ct_name in names(marker_sets_final)) {
  genes <- intersect(marker_sets_final[[ct_name]], rownames(seurat_obj))
  if (length(genes) == 0) { cat("  No genes found for", ct_name, "\n"); next }
  tryCatch({
    p <- VlnPlot(seurat_obj, features = genes, ncol = length(genes), pt.size = 0) &
      theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
            axis.title.x = element_blank())
    ggsave(
      file.path("outputs/scRNAseq_outputs/05_final_celltype_annotation/violin_celltype_final",
                 paste0("violin_", gsub("\\.", "_", ct_name), ".pdf")),
      p, width = max(5, length(genes) * 3.2), height = 5
    )
    cat("  Saved violin:", ct_name, "\n")
  }, error = function(e) cat("  Warning in violin", ct_name, ":", conditionMessage(e), "\n"))
}

top1_final <- sapply(marker_sets_final, function(g) intersect(g, rownames(seurat_obj))[1])
top1_final <- top1_final[!is.na(top1_final)]
tryCatch({
  p <- VlnPlot(seurat_obj, features = unique(top1_final), ncol = 5, pt.size = 0) &
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
          axis.title.x = element_blank())
  ggsave("outputs/scRNAseq_outputs/05_final_celltype_annotation/violin_celltype_final/violin_combined_top_markers.pdf",
         p, width = 18, height = 10)
  cat("  Saved combined top-marker violin (celltype_final)\n")
}, error = function(e) cat("  Warning in combined violin:", conditionMessage(e), "\n"))

# Restore state expected by the rest of the script
DefaultAssay(seurat_obj) <- "SCT"
Idents(seurat_obj) <- "seurat_clusters"

## STEP 10 — Cell-type composition by HPV status

# ALL CELL TYPES — PROPORTION PLOTS. Here celltype_final is the 24-cluster
# annotation collapsed to 9 types (Malignant folded into Epithelial).

dir.create("outputs/scRNAseq_outputs/06_composition_by_HPV/proportions_all_celltypes", recursive = TRUE, showWarnings = FALSE)

ct_final_levels <- sort(unique(seurat_obj$celltype_final))
ct_final_pal <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(ct_final_levels)),
  ct_final_levels
)
hpv_pal <- c("HPV+" = "#E41A1C", "HPV-" = "#377EB8")

prop_df_all <- seurat_obj@meta.data %>%
  group_by(HPVstatus, celltype_final) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(HPVstatus) %>%
  mutate(prop = n / sum(n))

write.csv(prop_df_all,
          "outputs/scRNAseq_outputs/06_composition_by_HPV/proportions_all_celltypes/celltype_final_proportions_by_HPVstatus.csv",
          row.names = FALSE)

# Stacked bar (composition, normalised to 1) — all 9 cell types
p_prop_all_stacked <- ggplot(prop_df_all, aes(x = HPVstatus, y = prop, fill = celltype_final)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = ct_final_pal) +
  labs(title = "Cell Type Proportions by HPV Status (9 cell types)",
       x = "HPV Status", y = "Proportion", fill = "Cell Type") +
  theme_minimal()
ggsave("outputs/scRNAseq_outputs/06_composition_by_HPV/proportions_all_celltypes/09b_all_celltype_proportions_stacked.pdf",
       p_prop_all_stacked, width = 7, height = 6)

# Grouped bar — each of the 9 cell types side-by-side across HPV groups
p_prop_all_grouped <- ggplot(prop_df_all, aes(x = reorder(celltype_final, -prop), y = prop, fill = HPVstatus)) +
  geom_bar(stat = "identity", position = "dodge", colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = hpv_pal) +
  labs(title = "Cell Type Proportions: HPV+ vs HPV- (9 cell types)",
       x = NULL, y = "Proportion", fill = "HPV Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("outputs/scRNAseq_outputs/06_composition_by_HPV/proportions_all_celltypes/10b_all_celltype_proportions_grouped.pdf",
       p_prop_all_grouped, width = 10, height = 6)

cat("All-celltype (9 types) proportion plots saved.\n")
cat("B_Plasma.cells proportion by HPV status (celltype_final, all cells):\n")
print(prop_df_all %>% filter(celltype_final == "B_Plasma.cells"))

# Patient-level statistical testing of composition by HPV status
## Why: prop_df_all above is descriptive only (per-HPV-group proportions with
## no test of significance). Cell-level pooled tests would treat thousands of
## correlated cells from a handful of patients as independent observations
## (pseudoreplication), so the primary test here uses patient as the unit of
## replication (Wilcoxon on per-patient proportions), with the pooled-cell
## chi-squared/Fisher tests kept only as a descriptive complement.
prop_stats_dir <- "outputs/scRNAseq_outputs/06_composition_by_HPV/proportions_all_celltypes"

md <- seurat_obj@meta.data %>%
  select(patientid, HPVstatus, celltype_final) %>%
  filter(!is.na(celltype_final), !is.na(HPVstatus), !is.na(patientid))

# One patient should carry a single HPV status; this is a sanity check, not a
# correction step — if it fails, the metadata itself is inconsistent.
patient_hpv <- md %>% distinct(patientid, HPVstatus)
stopifnot(!any(duplicated(patient_hpv$patientid)))

cat("Patients per HPV group:\n")
print(table(patient_hpv$HPVstatus))

# Patient-level proportions (replicate unit = patient, avoids cell-level
# pseudoreplication). This is what the combined box/violin figure and the
# Wilcoxon tests below are built on.
patient_counts <- md %>%
  dplyr::count(patientid, HPVstatus, celltype_final, name = "n") %>%
  complete(patientid, celltype_final, fill = list(n = 0)) %>%
  group_by(patientid) %>%
  fill(HPVstatus, .direction = "downup") %>%
  ungroup() %>%
  group_by(patientid) %>%
  mutate(total = sum(n), prop = n / total) %>%
  ungroup()

write.csv(patient_counts,
          file.path(prop_stats_dir, "celltype_final_proportions_by_patient.csv"),
          row.names = FALSE)

# Statistical testing per cell type:
#  a) Wilcoxon rank-sum test on patient-level proportions -> answers "which
#     cell types shift", respecting patients as the unit of replication.
#  b) Two-proportion test (chi-squared / Fisher) on pooled cell counts
#     (cell type vs all others, HPV+ vs HPV-) -> classic composition test,
#     but pools all cells from all patients into one number per group.
#  c) Overall chi-squared test on the full celltype_final x HPVstatus
#     contingency table of pooled cell counts -> single omnibus p-value.
stats_rows <- list()
pooled_counts <- md %>% dplyr::count(HPVstatus, celltype_final, name = "n")

for (ct in ct_final_levels) {

  df_ct <- patient_counts %>% filter(celltype_final == ct)

  pos_vals <- df_ct$prop[df_ct$HPVstatus == "HPV+"]
  neg_vals <- df_ct$prop[df_ct$HPVstatus == "HPV-"]

  wilcox_p <- tryCatch(
    wilcox.test(pos_vals, neg_vals, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )

  n_pos_cells <- pooled_counts$n[pooled_counts$HPVstatus == "HPV+" & pooled_counts$celltype_final == ct]
  n_neg_cells <- pooled_counts$n[pooled_counts$HPVstatus == "HPV-" & pooled_counts$celltype_final == ct]
  n_pos_cells <- if (length(n_pos_cells) == 0) 0 else n_pos_cells
  n_neg_cells <- if (length(n_neg_cells) == 0) 0 else n_neg_cells

  total_pos_cells <- sum(md$HPVstatus == "HPV+")
  total_neg_cells <- sum(md$HPVstatus == "HPV-")

  tab2x2 <- matrix(
    c(n_pos_cells, total_pos_cells - n_pos_cells,
      n_neg_cells, total_neg_cells - n_neg_cells),
    nrow = 2,
    dimnames = list(c("celltype", "other"), c("HPV+", "HPV-"))
  )

  expected_ok <- tryCatch(all(suppressWarnings(chisq.test(tab2x2)$expected) >= 5), error = function(e) FALSE)
  prop_test_method <- if (isTRUE(expected_ok)) "chi-squared" else "fisher.exact"
  prop_test_p <- tryCatch({
    if (isTRUE(expected_ok)) {
      suppressWarnings(prop.test(c(n_pos_cells, n_neg_cells), c(total_pos_cells, total_neg_cells))$p.value)
    } else {
      fisher.test(tab2x2)$p.value
    }
  }, error = function(e) NA_real_)

  stats_rows[[ct]] <- data.frame(
    celltype           = ct,
    n_patients_HPVpos  = length(pos_vals),
    n_patients_HPVneg  = length(neg_vals),
    mean_prop_HPVpos   = mean(pos_vals),
    mean_prop_HPVneg   = mean(neg_vals),
    median_prop_HPVpos = median(pos_vals),
    median_prop_HPVneg = median(neg_vals),
    n_cells_HPVpos     = n_pos_cells,
    n_cells_HPVneg     = n_neg_cells,
    pooled_prop_HPVpos = n_pos_cells / total_pos_cells,
    pooled_prop_HPVneg = n_neg_cells / total_neg_cells,
    wilcox_p_patientlevel = wilcox_p,
    proptest_method_pooled = prop_test_method,
    proptest_p_pooled      = prop_test_p,
    stringsAsFactors = FALSE
  )
}

stats_df <- bind_rows(stats_rows)
stats_df$wilcox_p_adj_BH   <- p.adjust(stats_df$wilcox_p_patientlevel, method = "BH")
stats_df$proptest_p_adj_BH <- p.adjust(stats_df$proptest_p_pooled, method = "BH")
stats_df <- stats_df %>% arrange(wilcox_p_adj_BH)

write.csv(stats_df, file.path(prop_stats_dir, "celltype_HPV_proportion_test_summary.csv"), row.names = FALSE)

cat("\nPer-cell-type statistical summary (sorted by Wilcoxon BH-adjusted p):\n")
print(stats_df, row.names = FALSE)

# Omnibus chi-squared test — full celltype_final x HPVstatus contingency
# table, pooled cells
ct_table_all <- table(md$celltype_final, md$HPVstatus)
chi_all <- chisq.test(ct_table_all)

cat("\nOmnibus chi-squared test (", length(ct_final_levels), "cell types x HPV status, pooled cells):\n")
print(ct_table_all)
print(chi_all)
cat("Any expected count < 5 (chi-squared assumption check):",
    any(chi_all$expected < 5), "\n")

omnibus_summary <- data.frame(
  test = paste0("omnibus_chisq_", length(ct_final_levels), "celltypes_x_HPVstatus"),
  statistic = unname(chi_all$statistic),
  df = unname(chi_all$parameter),
  p_value = chi_all$p.value,
  any_expected_lt5 = any(chi_all$expected < 5)
)
write.csv(omnibus_summary, file.path(prop_stats_dir, "celltype_HPV_omnibus_chisq.csv"), row.names = FALSE)

# Combined faceted figure — all cell types in one panel for quick view
combined_df <- patient_counts %>%
  left_join(stats_df %>% select(celltype_final = celltype, wilcox_p_adj_BH), by = "celltype_final") %>%
  mutate(facet_label = sprintf("%s\n(BH p = %s)", celltype_final,
                                ifelse(is.na(wilcox_p_adj_BH), "NA", formatC(wilcox_p_adj_BH, format = "g", digits = 2))))

p_combined <- ggplot(combined_df, aes(x = HPVstatus, y = prop, fill = HPVstatus)) +
  geom_violin(alpha = 0.4, trim = FALSE, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.08, size = 1.2, shape = 21, colour = "black", stroke = 0.25) +
  scale_fill_manual(values = hpv_pal) +
  facet_wrap(~ facet_label, scales = "free_y", nrow = 3) +
  labs(title = "Cell Type Proportions by HPV Status (patient-level, Wilcoxon BH-adjusted)",
       x = "HPV Status", y = "Proportion of cells (per patient)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(prop_stats_dir, "boxviolin_all_celltypes_combined_HPVstatus.pdf"), p_combined, width = 12, height = 11)
ggsave(file.path(prop_stats_dir, "boxviolin_all_celltypes_combined_HPVstatus.png"), p_combined, width = 12, height = 11, dpi = 300, bg = "white")

cat("\nComposition statistics written to:", prop_stats_dir, "\n")
cat(" - celltype_final_proportions_by_patient.csv  (per-patient proportions, raw data)\n")
cat(" - boxviolin_all_celltypes_combined_HPVstatus.pdf/.png (combined facet figure, all cell types)\n")
cat(" - celltype_HPV_proportion_test_summary.csv    (Wilcoxon + pooled proportion test per cell type)\n")
cat(" - celltype_HPV_omnibus_chisq.csv               (omnibus chi-squared across all cell types)\n")

# IMMUNE CELLS (5-type subset, kept for the original immune-focused comparison)
# Exact labels verified against table(seurat_obj$celltype) from the GEO metadata
immune_celltype_labels <- c("T.cells", "B_Plasma.cells", "Macrophages",
                            "Dendritic.cells", "Mast.cells")

# Confirm all labels exist and show their counts before subsetting
cat("Immune celltype counts in full dataset:\n")
print(table(seurat_obj$celltype)[immune_celltype_labels])

# Subset using the original celltype column
immune_obj <- subset(seurat_obj, subset = celltype %in% immune_celltype_labels)

# Confirm subset composition
cat("\nCells in immune_obj by celltype:\n")
print(table(immune_obj$celltype))

# Stratify by HPV status
cat("\nimmune_obj cells by HPV status:\n")
print(table(immune_obj$HPVstatus))

# Plot cell type proportions by HPV status
library(dplyr)
library(ggplot2)

# Tabulate proportions
prop_df <- immune_obj@meta.data %>%
  group_by(HPVstatus, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(HPVstatus) %>%
  mutate(prop = n / sum(n))

# Now print the resulting data frame
print(prop_df)

# UMAP Coloured by HPV Status
p_immune_umap_hpv <- DimPlot(immune_obj, reduction = "umap", group.by = "HPVstatus", label = TRUE) +
  ggtitle("UMAP of Immune Cells Coloured by HPV Status")
ggsave("outputs/scRNAseq_outputs/06_composition_by_HPV/07_immune_UMAP_HPVstatus.pdf", p_immune_umap_hpv, width = 8, height = 7)

# Split by HPV status, coloured by celltype
p_immune_umap_split <- DimPlot(immune_obj, reduction = "umap", split.by = "HPVstatus", group.by = "celltype", label = TRUE) +
  ggtitle("UMAP Split by HPV Status")
ggsave("outputs/scRNAseq_outputs/06_composition_by_HPV/08_immune_UMAP_split_HPVstatus.pdf", p_immune_umap_split, width = 14, height = 7)

# Stacked bar: absolute proportions per HPV group
p_prop_stacked <- ggplot(prop_df, aes(x = HPVstatus, y = prop, fill = celltype)) +
  geom_bar(stat = "identity") +
  labs(title = "Immune Cell Proportions by HPV Status", y = "Proportion", fill = "Immune Cell Type") +
  theme_minimal()
ggsave("outputs/scRNAseq_outputs/06_composition_by_HPV/09_immune_proportions_stacked.pdf", p_prop_stacked, width = 7, height = 6)

# Stacked bar: normalised to 1 (composition view)
p_prop_composition <- ggplot(prop_df, aes(x = HPVstatus, y = prop, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Proportion of Immune Cell Types") +
  xlab("HPV Status") +
  theme_minimal() +
  labs(fill = "Immune Cell Type")
ggsave("outputs/scRNAseq_outputs/06_composition_by_HPV/10_immune_proportions_composition.pdf", p_prop_composition, width = 7, height = 6)

# Grouped bar: each cell type side-by-side across HPV groups
p_prop_grouped <- ggplot(prop_df, aes(x = celltype, y = prop, fill = HPVstatus)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Immune Cell Proportions by Cell Type and HPV Status",
       x = "Cell Type",
       y = "Proportion",
       fill = "HPV Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("outputs/scRNAseq_outputs/06_composition_by_HPV/11_immune_proportions_grouped.pdf", p_prop_grouped, width = 9, height = 6)

## STEP 11 — Cross-check available cell-type labels before subsetting
## Why: before subsetting to specific cell types for HPV comparisons below,
## confirm exactly which label strings exist in both the original GEO
## `celltype` column and the SingleR predictions — subsetting on a typo'd or
## mismatched label string would silently produce an empty/wrong subset.

#Validation of Cell Types
# 4.1 First: Explore Available Cell Types

# STEP 1: Examine what cell types SingleR predicted
available_types <- sort(unique(seurat_obj$singler_label))
cat("Available cell types in SingleR annotation:\n")
print(available_types)

# Distribution of SingleR-predicted labels
cat("\nSingleR label counts:\n")
print(table(seurat_obj$singler_label))

# STEP 2: Look for immune-related patterns in SingleR labels
# More precise pattern to avoid false matches
immune_pattern <- "^T\\.|^B_|^Macrophages|^Dendritic|^Mast|NK|Mono|Neutro"
potential_immune <- available_types[grepl(immune_pattern, available_types, ignore.case = TRUE)]
cat("\nImmune cell types found in SingleR labels:\n")
print(potential_immune)

# Cross-tab: original celltype annotation vs SingleR predicted labels
cat("\nCross-tabulation: original celltype (rows) vs SingleR_label (columns):\n")
print(table(seurat_obj$celltype, seurat_obj$singler_label))

# =============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS: HPV+ vs HPV-
# =============================================================================

library(ggplot2)
library(ggrepel)
library(dplyr)

## STEP 12 — Systematic differential expression: HPV+ vs HPV- (pseudobulk DESeq2)
## sum raw counts into one pseudobulk profile per patient (per cell type) and run DESeq2 
## with design = ~ HPVstatus on those profiles. Patient_ID can't be added to the design: 
## it's collinear with HPVstatus, since every patient has only ever had one HPV status.

DefaultAssay(seurat_obj) <- "RNA"
raw_counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
meta       <- seurat_obj@meta.data

min_cells_per_sample   <- 10  # min cells contributing to a pseudobulk profile
min_patients_per_group <- 3   # min patients per HPV group to attempt DESeq2

aggregate_pseudobulk <- function(cell_ids_by_patient) {
  mat <- sapply(cell_ids_by_patient, function(ids) Matrix::rowSums(raw_counts[, ids, drop = FALSE]))
  colnames(mat) <- names(cell_ids_by_patient)
  mat
}

run_pseudobulk_deseq2 <- function(count_mat, sample_meta) {
  # DESeqDataSet runs factor levels through make.names(), which collapses
  # both "HPV+" and "HPV-" to "HPV." (it strips +/- to "."), so the levels
  # must be recoded to names that survive make.names() unchanged.
  sample_meta$HPVstatus <- factor(ifelse(sample_meta$HPVstatus == "HPV+", "HPVpos", "HPVneg"),
                                   levels = c("HPVneg", "HPVpos"))
  count_mat <- count_mat[rowSums(count_mat) >= 10, , drop = FALSE]

  dds <- DESeqDataSetFromMatrix(countData = round(count_mat),
                                 colData   = sample_meta,
                                 design    = ~ HPVstatus)
  dds <- DESeq(dds)
  contrast <- c("HPVstatus", "HPVpos", "HPVneg")
  res <- results(dds, contrast = contrast, alpha = 0.05)
  res_stat <- res$stat  ## Wald stat, saved before lfcShrink (which drops it) -- used as the GSEA ranking metric below
  ## apeglm requires coef= (not contrast=) matching resultsNames(dds); since
  ## HPVneg is the reference level, that coefficient is algebraically the same contrast.
  coef_name <- resultsNames(dds)[grepl("HPVstatus", resultsNames(dds))]
  res <- lfcShrink(dds, coef = coef_name, res = res, type = "apeglm")

  res_df <- as.data.frame(res)
  res_df$gene       <- rownames(res_df)
  res_df$avg_log2FC <- res_df$log2FoldChange
  res_df$p_val      <- res_df$pvalue
  res_df$p_val_adj  <- res_df$padj
  res_df$stat       <- res_stat
  res_df
}

patients_for_cells <- function(cell_meta) {
  ids_by_patient <- split(rownames(cell_meta), cell_meta$patientid)
  ids_by_patient[sapply(ids_by_patient, length) >= min_cells_per_sample]
}

patient_metadata_for <- function(cell_meta, pb_colnames) {
  pm <- cell_meta %>%
    distinct(patientid, HPVstatus) %>%
    filter(patientid %in% pb_colnames)
  rownames(pm) <- pm$patientid
  pm[pb_colnames, ]
}

# 1. All-cell DE: HPV+ vs HPV- (pseudobulk per patient, all cell types pooled)
de_all_checkpoint <- "outputs/scRNAseq_outputs/checkpoints/de_all.rds"
if (file.exists(de_all_checkpoint)) {
  cat("Checkpoint found — loading de_all\n")
  de_all <- readRDS(de_all_checkpoint)
} else {
  ids_by_patient  <- patients_for_cells(meta)
  pb_all          <- aggregate_pseudobulk(ids_by_patient)
  patient_meta    <- patient_metadata_for(meta, colnames(pb_all))

  de_all <- run_pseudobulk_deseq2(pb_all, patient_meta)

  saveRDS(de_all, de_all_checkpoint)
  cat("Saved checkpoint:", de_all_checkpoint, "\n")
}

head(de_all[order(de_all$avg_log2FC, decreasing = TRUE),  ], 20)
head(de_all[order(de_all$avg_log2FC, decreasing = FALSE), ], 20)

# Volcano plot: all cells
de_all$sig <- ifelse(de_all$p_val_adj < 0.05 & abs(de_all$avg_log2FC) > 0.5,
                     ifelse(de_all$avg_log2FC > 0, "HPV+ up", "HPV- up"), "NS")

top_labels <- de_all %>%
  filter(sig != "NS") %>%
  slice_max(order_by = abs(avg_log2FC), n = 20)

p_de_all_volcano <- ggplot(de_all, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sig)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_text_repel(data = top_labels, aes(label = gene), size = 3, max.overlaps = 20) +
  scale_colour_manual(values = c("HPV+ up" = "#E41A1C", "HPV- up" = "#377EB8", "NS" = "grey70")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  labs(title = "DE: HPV+ vs HPV- (all cells)", x = "log2 Fold Change", y = "-log10 adj. p-value") +
  theme_classic()
ggsave("outputs/scRNAseq_outputs/09_differential_expression/23_DE_all_cells_volcano.pdf", p_de_all_volcano, width = 9, height = 7)

# 2. Cell-type-specific DE (pseudobulk per patient, within each cell type)
## Uses celltype_final (this pipeline's own SingleR+marker-validated
## reclustering annotation) for consistency with STEP 10/11 composition
## analysis and the 06/07 bulk-signature-integration scripts -- NOT the
## original GEO-provided "celltype" column, which uses a different label set
## (e.g. a separate Malignant.cells category folded into Epithelial here).
Idents(seurat_obj) <- "celltype_final"

de_by_celltype_checkpoint <- "outputs/scRNAseq_outputs/checkpoints/de_by_celltype.rds"
if (file.exists(de_by_celltype_checkpoint)) {
  cat("Checkpoint found — loading de_combined_hpv\n")
  de_combined_hpv <- readRDS(de_by_celltype_checkpoint)
} else {
  cell_types_de  <- unique(meta$celltype_final)
  de_by_celltype <- list()

  for (ct in cell_types_de) {
    meta_ct        <- meta[meta$celltype_final == ct, ]
    ids_by_patient <- patients_for_cells(meta_ct)
    if (length(ids_by_patient) == 0) next

    pb_ct        <- aggregate_pseudobulk(ids_by_patient)
    patient_meta <- patient_metadata_for(meta_ct, colnames(pb_ct))

    n_per_group <- table(patient_meta$HPVstatus)
    if (length(n_per_group) < 2 || any(n_per_group < min_patients_per_group)) {
      cat("Skipping", ct, "- fewer than", min_patients_per_group, "patients in one HPV group\n")
      next
    }

    de_result            <- run_pseudobulk_deseq2(pb_ct, patient_meta)
    de_result$celltype   <- ct
    de_by_celltype[[ct]] <- de_result
  }

  de_combined_hpv <- bind_rows(de_by_celltype)
  saveRDS(de_combined_hpv, de_by_celltype_checkpoint)
  cat("Saved checkpoint:", de_by_celltype_checkpoint, "\n")
}

de_sig <- de_combined_hpv %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  arrange(celltype, desc(abs(avg_log2FC)))
print(de_sig)

# Volcano plots per cell type
for (ct in unique(de_combined_hpv$celltype)) {
  df     <- de_combined_hpv %>% filter(celltype == ct)
  df$sig <- ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) > 0.5,
                   ifelse(df$avg_log2FC > 0, "HPV+ up", "HPV- up"), "NS")
  top_l  <- df %>% filter(sig != "NS") %>% slice_max(order_by = abs(avg_log2FC), n = 15)

  p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = sig)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_text_repel(data = top_l, aes(label = gene), size = 3, max.overlaps = 15) +
    scale_colour_manual(values = c("HPV+ up" = "#E41A1C", "HPV- up" = "#377EB8", "NS" = "grey70")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    labs(title = paste("DE:", ct, "— HPV+ vs HPV-"),
         x = "log2 Fold Change", y = "-log10 adj. p-value") +
    theme_classic()
  print(p)
  ct_safe <- gsub("[^A-Za-z0-9]+", "_", ct)
  ggsave(paste0("outputs/scRNAseq_outputs/09_differential_expression/24_DE_volcano_", ct_safe, ".pdf"), p, width = 9, height = 7)
}

# Key visualisations

# Bar plot: significant DEGs per cell type
deg_counts <- de_combined_hpv %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  mutate(direction = ifelse(avg_log2FC > 0, "HPV+ up", "HPV- up")) %>%
  dplyr::count(celltype, direction)

write.csv(deg_counts, "outputs/scRNAseq_outputs/de_tables/DEG_counts_per_celltype.csv", row.names = FALSE)
cat("Saved: outputs/scRNAseq_outputs/de_tables/DEG_counts_per_celltype.csv\n")

p_deg_counts <- ggplot(deg_counts, aes(x = reorder(celltype, n), y = n, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("HPV+ up" = "#E41A1C", "HPV- up" = "#377EB8")) +
  coord_flip() +
  labs(title = "Number of Significant DEGs per Cell Type (HPV+ vs HPV-)",
       x = "Cell Type", y = "Number of DEGs") +
  theme_classic()
ggsave("outputs/scRNAseq_outputs/09_differential_expression/25_DEG_counts_per_celltype.pdf", p_deg_counts, width = 8, height = 6)

# Heatmap: top 40 DE genes, all cells split by HPV status
Idents(seurat_obj) <- "HPVstatus"
top_de_genes <- de_all %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(order_by = abs(avg_log2FC), n = 40) %>%
  pull(gene)

seurat_obj <- ScaleData(seurat_obj, features = top_de_genes, assay = "RNA")
p_heatmap_all <- DoHeatmap(seurat_obj, features = top_de_genes, group.by = "HPVstatus",
          group.colors = c("HPV+" = "#E41A1C", "HPV-" = "#377EB8")) +
  ggtitle("Top 40 DE Genes: HPV+ vs HPV-")
ggsave("outputs/scRNAseq_outputs/09_differential_expression/26_heatmap_top40_DEgenes_allcells.pdf", p_heatmap_all, width = 10, height = 12)

# Heatmap: top 5 DE genes per cell type
Idents(seurat_obj) <- "celltype_final"
top5_per_ct <- de_combined_hpv %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  group_by(celltype) %>%
  slice_max(order_by = abs(avg_log2FC), n = 5) %>%
  pull(gene) %>%
  unique()

seurat_obj <- ScaleData(seurat_obj, features = union(top_de_genes, top5_per_ct), assay = "RNA")
p_heatmap_ct <- DoHeatmap(seurat_obj, features = top5_per_ct, group.by = "celltype_final", size = 3) +
  ggtitle("Top DE Genes per Cell Type")
ggsave("outputs/scRNAseq_outputs/09_differential_expression/27_heatmap_top5_DEgenes_per_celltype.pdf", p_heatmap_ct, width = 10, height = 12)

# VlnPlot: top 6 DE genes, split by HPV status
Idents(seurat_obj) <- "HPVstatus"
top6_genes <- de_all %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(order_by = abs(avg_log2FC), n = 6) %>%
  pull(gene)

p_vln_top6 <- VlnPlot(seurat_obj, features = top6_genes, ncol = 3,
        cols = c("HPV+" = "#E41A1C", "HPV-" = "#377EB8")) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("outputs/scRNAseq_outputs/09_differential_expression/28_VlnPlot_top6_DEgenes.pdf", p_vln_top6, width = 12, height = 8)

# FeaturePlot: top 4 DE genes on UMAP, split by HPV status
top4_genes <- de_all %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(order_by = abs(avg_log2FC), n = 4) %>%
  pull(gene)

p_feature_top4 <- FeaturePlot(seurat_obj, features = top4_genes, split.by = "HPVstatus",
            ncol = 4, pt.size = 0.3) &
  theme(legend.position = "right")
ggsave("outputs/scRNAseq_outputs/09_differential_expression/29_FeaturePlot_top4_DEgenes.pdf", p_feature_top4, width = 14, height = 8)

# DotPlot: top DE genes across cell types and HPV status
Idents(seurat_obj) <- "celltype_final"
seurat_obj$celltype_hpv <- paste(seurat_obj$celltype_final, seurat_obj$HPVstatus, sep = "_")
Idents(seurat_obj) <- "celltype_hpv"

top10_de <- de_all %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(order_by = abs(avg_log2FC), n = 10) %>%
  pull(gene)

p_dotplot_de <- DotPlot(seurat_obj, features = top10_de, dot.scale = 5) +
  RotatedAxis() +
  ggtitle("Top DE Genes Expression Across Cell Types and HPV Status")
ggsave("outputs/scRNAseq_outputs/09_differential_expression/30_DotPlot_top10_DEgenes.pdf", p_dotplot_de, width = 12, height = 8)

Idents(seurat_obj) <- "celltype_final"

# =============================================================================
# PATHWAY ENRICHMENT ANALYSIS
# =============================================================================
## STEP 13 — Pathway enrichment on the pseudobulk DESeq2 results
## Why: GO/KEGG enrichment run systematically — for the pooled all-cell DE
## result and for every cell type that passed the pseudobulk DESeq2 testing
## threshold — to ask whether the statistically robust, patient-level DE
## genes converge on coherent biological pathways.

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

symbol_to_entrez <- function(genes) {
  res <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = org.Hs.eg.db, drop = TRUE)
  res$ENTREZID
}

# 1. All-cell GO and KEGG enrichment

hpv_pos_genes <- de_all %>% filter(p_val_adj < 0.05, avg_log2FC >  0.5) %>% pull(gene)
hpv_neg_genes <- de_all %>% filter(p_val_adj < 0.05, avg_log2FC < -0.5) %>% pull(gene)

go_kegg_allcells_checkpoint <- "outputs/scRNAseq_outputs/checkpoints/go_kegg_allcells.rds"
if (file.exists(go_kegg_allcells_checkpoint)) {
  cat("Checkpoint found — loading go_pos/go_neg/kegg_pos/kegg_neg\n")
  ck       <- readRDS(go_kegg_allcells_checkpoint)
  go_pos   <- ck$go_pos
  go_neg   <- ck$go_neg
  kegg_pos <- ck$kegg_pos
  kegg_neg <- ck$kegg_neg
} else {
  go_pos <- enrichGO(
    gene          = hpv_pos_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )

  go_neg <- enrichGO(
    gene          = hpv_neg_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )

  kegg_pos <- enrichKEGG(
    gene          = symbol_to_entrez(hpv_pos_genes),
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )

  kegg_neg <- enrichKEGG(
    gene          = symbol_to_entrez(hpv_neg_genes),
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )

  saveRDS(list(go_pos = go_pos, go_neg = go_neg, kegg_pos = kegg_pos, kegg_neg = kegg_neg),
          go_kegg_allcells_checkpoint)
  cat("Saved checkpoint:", go_kegg_allcells_checkpoint, "\n")
}

if (!is.null(go_pos) && nrow(go_pos) > 0) {
  p <- dotplot(go_pos, showCategory = 20) + ggtitle("GO BP: HPV+ Upregulated (all cells)")
  ggsave("outputs/scRNAseq_outputs/10_pathway_enrichment/31_GO_BP_HPVpos_allcells.pdf", p, width = 9, height = 8)
}
if (!is.null(go_neg) && nrow(go_neg) > 0) {
  p <- dotplot(go_neg, showCategory = 20) + ggtitle("GO BP: HPV- Upregulated (all cells)")
  ggsave("outputs/scRNAseq_outputs/10_pathway_enrichment/32_GO_BP_HPVneg_allcells.pdf", p, width = 9, height = 8)
}
if (!is.null(kegg_pos) && nrow(kegg_pos) > 0) {
  p <- dotplot(kegg_pos, showCategory = 20) + ggtitle("KEGG: HPV+ Upregulated (all cells)")
  ggsave("outputs/scRNAseq_outputs/10_pathway_enrichment/33_KEGG_HPVpos_allcells.pdf", p, width = 9, height = 8)
}
if (!is.null(kegg_neg) && nrow(kegg_neg) > 0) {
  p <- dotplot(kegg_neg, showCategory = 20) + ggtitle("KEGG: HPV- Upregulated (all cells)")
  ggsave("outputs/scRNAseq_outputs/10_pathway_enrichment/34_KEGG_HPVneg_allcells.pdf", p, width = 9, height = 8)
}

if (!is.null(go_pos) && nrow(go_pos) > 0) {
  go_pos_sim <- pairwise_termsim(go_pos)
  p_emap <- emapplot(go_pos_sim, showCategory = 30) + ggtitle("GO BP Enrichment Map: HPV+ Upregulated")
  print(p_emap)
  ggsave("outputs/scRNAseq_outputs/10_pathway_enrichment/35_GO_BP_emapplot_HPVpos.pdf", p_emap, width = 10, height = 9)
}

# 2. Cell-type-specific GO and KEGG

go_kegg_by_celltype_checkpoint <- "outputs/scRNAseq_outputs/checkpoints/go_kegg_by_celltype.rds"
if (file.exists(go_kegg_by_celltype_checkpoint)) {
  cat("Checkpoint found — loading go_by_celltype/kegg_by_celltype\n")
  ck               <- readRDS(go_kegg_by_celltype_checkpoint)
  go_by_celltype   <- ck$go_by_celltype
  kegg_by_celltype <- ck$kegg_by_celltype
} else {
  go_by_celltype   <- list()
  kegg_by_celltype <- list()

  for (ct in unique(de_combined_hpv$celltype)) {
    df        <- de_combined_hpv %>% filter(celltype == ct)
    pos_genes <- df %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>% pull(gene)

    if (length(pos_genes) >= 5) {
      go_res <- enrichGO(gene = pos_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                         ont = "BP", pAdjustMethod = "BH",
                         pvalueCutoff = 0.05, readable = TRUE)
      go_by_celltype[[ct]] <- go_res

      entrez <- symbol_to_entrez(pos_genes)
      if (length(entrez) >= 5) {
        kegg_res <- enrichKEGG(gene = entrez, organism = "hsa",
                               pAdjustMethod = "BH", pvalueCutoff = 0.05)
        kegg_by_celltype[[ct]] <- kegg_res
      }
    }
  }

  saveRDS(list(go_by_celltype = go_by_celltype, kegg_by_celltype = kegg_by_celltype),
          go_kegg_by_celltype_checkpoint)
  cat("Saved checkpoint:", go_kegg_by_celltype_checkpoint, "\n")
}

# Plots are always (re)generated from the cached or freshly computed results
for (ct in names(go_by_celltype)) {
  go_res <- go_by_celltype[[ct]]
  if (!is.null(go_res) && nrow(go_res) > 0) {
    ct_safe <- gsub("[^A-Za-z0-9]+", "_", ct)
    p_go_ct <- dotplot(go_res, showCategory = 15) +
            ggtitle(paste("GO BP:", ct, "— HPV+ upregulated"))
    print(p_go_ct)
    ggsave(paste0("outputs/scRNAseq_outputs/10_pathway_enrichment/36_GO_BP_", ct_safe, ".pdf"), p_go_ct, width = 9, height = 8)
  }
}
for (ct in names(kegg_by_celltype)) {
  kegg_res <- kegg_by_celltype[[ct]]
  if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
    ct_safe <- gsub("[^A-Za-z0-9]+", "_", ct)
    p_kegg_ct <- dotplot(kegg_res, showCategory = 15) +
            ggtitle(paste("KEGG:", ct, "— HPV+ upregulated"))
    print(p_kegg_ct)
    ggsave(paste0("outputs/scRNAseq_outputs/10_pathway_enrichment/37_KEGG_", ct_safe, ".pdf"), p_kegg_ct, width = 9, height = 8)
  }
}

# 3. Comparative dotplot: top pathways across cell types

gene_list_by_ct <- lapply(unique(de_combined_hpv$celltype), function(ct) {
  de_combined_hpv %>%
    filter(celltype == ct, p_val_adj < 0.05, avg_log2FC > 0.5) %>%
    pull(gene)
})
names(gene_list_by_ct) <- unique(de_combined_hpv$celltype)
gene_list_by_ct <- Filter(function(x) length(x) >= 5, gene_list_by_ct)

cc_go_checkpoint <- "outputs/scRNAseq_outputs/checkpoints/cc_go.rds"
if (length(gene_list_by_ct) > 0) {
  if (file.exists(cc_go_checkpoint)) {
    cat("Checkpoint found — loading cc_go\n")
    cc_go <- readRDS(cc_go_checkpoint)
  } else {
    cc_go <- compareCluster(
      geneCluster   = gene_list_by_ct,
      fun           = "enrichGO",
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    )
    saveRDS(cc_go, cc_go_checkpoint)
    cat("Saved checkpoint:", cc_go_checkpoint, "\n")
  }

  if (!is.null(cc_go) && nrow(cc_go) > 0) {
    p_cc_go <- dotplot(cc_go, showCategory = 5) +
      RotatedAxis() +
      ggtitle("Top GO BP Pathways per Cell Type (HPV+ upregulated)")
    ggsave("outputs/scRNAseq_outputs/10_pathway_enrichment/38_comparative_GO_BP_per_celltype.pdf", p_cc_go, width = 12, height = 9)
  }
} else {
  cat("Skipping comparative GO enrichment: no cell type had >=5 significant upregulated genes.\n")
}

# =============================================================================
# GSEA (GO BP + KEGG) on the full ranked all-cell pseudobulk gene list
# =============================================================================
## STEP 13b — Complements the ORA above: ORA only sees genes that pass the
## padj/logFC cutoff, tested per direction separately; GSEA ranks every gene
## by the DESeq2 Wald stat and finds pathways skewed toward one end in a
## single run, picking up coordinated sub-threshold shifts ORA misses.
## Mirrors bulk_seq_analysis/06_run_de_pathway_immune_hot_vs_cold.R's GSEA
## block; run here on the all-cell pseudobulk comparison (de_all) only, not
## per cell type -- add a per-celltype loop the same way if that's wanted.
cat("\nRunning GSEA (GO BP + KEGG) on the all-cell pseudobulk ranked gene list...\n")

## One rank per gene symbol: keep the strongest |stat| when a symbol maps to
## multiple rows, then sort descending (gseGO/gseKEGG require unique,
## sorted-decreasing named vectors).
gsea_ranks <- de_all %>% filter(!is.na(stat), !is.na(gene), gene != "")
gsea_ranks <- gsea_ranks[order(-abs(gsea_ranks$stat)), ]
gsea_ranks <- gsea_ranks[!duplicated(gsea_ranks$gene), ]
gsea_ranks <- gsea_ranks[order(-gsea_ranks$stat), ]
gene_rank_go <- setNames(gsea_ranks$stat, gsea_ranks$gene)
cat("Ranked gene list for GSEA:", length(gene_rank_go), "unique gene symbols\n")

gsea_allcells_checkpoint <- "outputs/scRNAseq_outputs/checkpoints/gsea_allcells.rds"
if (file.exists(gsea_allcells_checkpoint)) {
  cat("Checkpoint found — loading gsea_go/gsea_kegg\n")
  ck        <- readRDS(gsea_allcells_checkpoint)
  gsea_go   <- ck$gsea_go
  gsea_kegg <- ck$gsea_kegg
} else {
  set.seed(42)  ## GSEA permutation p-values are stochastic; fix seed for reproducibility

  gsea_go <- gseGO(geneList = gene_rank_go, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                    minGSSize = 10, maxGSSize = 500, eps = 0)

  ## KEGG needs Entrez IDs -- map, then re-collapse/re-sort the same way as
  ## above since bitr's SYMBOL->ENTREZID map is not 1:1.
  entrez_map    <- bitr(gsea_ranks$gene, fromType = "SYMBOL", toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db, drop = TRUE)
  rank_kegg_df  <- merge(gsea_ranks, entrez_map, by.x = "gene", by.y = "SYMBOL")
  rank_kegg_df  <- rank_kegg_df[order(-abs(rank_kegg_df$stat)), ]
  rank_kegg_df  <- rank_kegg_df[!duplicated(rank_kegg_df$ENTREZID), ]
  rank_kegg_df  <- rank_kegg_df[order(-rank_kegg_df$stat), ]
  gene_rank_kegg <- setNames(rank_kegg_df$stat, rank_kegg_df$ENTREZID)

  gsea_kegg <- gseKEGG(geneList = gene_rank_kegg, organism = "hsa",
                        pAdjustMethod = "BH", pvalueCutoff = 0.05,
                        minGSSize = 10, maxGSSize = 500, eps = 0)

  saveRDS(list(gsea_go = gsea_go, gsea_kegg = gsea_kegg), gsea_allcells_checkpoint)
  cat("Saved checkpoint:", gsea_allcells_checkpoint, "\n")
}

## split=".sign" + facet puts HPV+ (activated) and HPV- (suppressed)
## enrichments in separate panels, same idea as the ORA hot/cold dotplot
## pairs above but from one GSEA run instead of two ORA runs.
if (!is.null(gsea_go) && nrow(gsea_go) > 0) {
  p <- dotplot(gsea_go, showCategory = 15, split = ".sign") +
    facet_grid(. ~ .sign, labeller = as_labeller(c(activated = "Higher in HPV+",
                                                     suppressed = "Higher in HPV-"))) +
    ggtitle("GSEA GO BP: HPV+ vs HPV- (all cells, pseudobulk)")
  ggsave("outputs/scRNAseq_outputs/10_pathway_enrichment/39_GSEA_GO_BP_allcells.pdf", p, width = 13, height = 9)
}
if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
  p <- dotplot(gsea_kegg, showCategory = 15, split = ".sign") +
    facet_grid(. ~ .sign, labeller = as_labeller(c(activated = "Higher in HPV+",
                                                     suppressed = "Higher in HPV-"))) +
    ggtitle("GSEA KEGG: HPV+ vs HPV- (all cells, pseudobulk)")
  ggsave("outputs/scRNAseq_outputs/10_pathway_enrichment/40_GSEA_KEGG_allcells.pdf", p, width = 13, height = 9)
}

## Classic GSEA running-enrichment-score plot for the single strongest hit at
## each end of the ranking (most HPV+-skewed, most HPV--skewed).
if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
  kegg_ord    <- gsea_kegg@result[order(-gsea_kegg@result$NES), ]
  top_hpvpos_id <- kegg_ord$ID[1]
  top_hpvneg_id <- kegg_ord$ID[nrow(kegg_ord)]

  p_hpvpos <- gseaplot2(gsea_kegg, geneSetID = top_hpvpos_id,
                         title = paste0("Top HPV+ pathway: ", kegg_ord$Description[1]))
  ggsave("outputs/scRNAseq_outputs/10_pathway_enrichment/41_GSEA_KEGG_top_HPVpos_gseaplot.pdf", p_hpvpos, width = 8, height = 7)

  p_hpvneg <- gseaplot2(gsea_kegg, geneSetID = top_hpvneg_id,
                         title = paste0("Top HPV- pathway: ", kegg_ord$Description[nrow(kegg_ord)]))
  ggsave("outputs/scRNAseq_outputs/10_pathway_enrichment/42_GSEA_KEGG_top_HPVneg_gseaplot.pdf", p_hpvneg, width = 8, height = 7)
}

write.csv(as.data.frame(gsea_go),   file = "outputs/scRNAseq_outputs/de_tables/GSEA_GO_BP_allcells.csv",   row.names = FALSE)
write.csv(as.data.frame(gsea_kegg), file = "outputs/scRNAseq_outputs/de_tables/GSEA_KEGG_allcells.csv",    row.names = FALSE)

cat("\n=== GSEA summary ===\n")
cat("GO BP gene sets (padj<0.05):", ifelse(is.null(gsea_go),   0, nrow(gsea_go)),   "\n")
cat("KEGG gene sets (padj<0.05):",  ifelse(is.null(gsea_kegg), 0, nrow(gsea_kegg)), "\n")
cat("GSEA results saved to outputs/scRNAseq_outputs/de_tables/\n")

# =============================================================================
# SAVE RESULTS
# =============================================================================
## STEP 14 — Export final tables
## Why: persist the DE and pathway enrichment results as plain CSVs, so the
## thesis write-up and any downstream analysis can consume them without
## having to re-run this entire pipeline (or even have R installed).

dir.create("outputs/scRNAseq_outputs/de_tables", recursive = TRUE, showWarnings = FALSE)

write.csv(de_all,          file = "outputs/scRNAseq_outputs/de_tables/DE_all_cells_HPVpos_vs_HPVneg.csv",   row.names = FALSE)
write.csv(de_combined_hpv, file = "outputs/scRNAseq_outputs/de_tables/DE_by_celltype_HPVpos_vs_HPVneg.csv", row.names = FALSE)
cat("DE results saved.\n")

write.csv(as.data.frame(go_pos),   file = "outputs/scRNAseq_outputs/de_tables/GO_BP_HPVpos_upregulated.csv",  row.names = FALSE)
write.csv(as.data.frame(go_neg),   file = "outputs/scRNAseq_outputs/de_tables/GO_BP_HPVneg_upregulated.csv",  row.names = FALSE)
write.csv(as.data.frame(kegg_pos), file = "outputs/scRNAseq_outputs/de_tables/KEGG_HPVpos_upregulated.csv",   row.names = FALSE)
write.csv(as.data.frame(kegg_neg), file = "outputs/scRNAseq_outputs/de_tables/KEGG_HPVneg_upregulated.csv",   row.names = FALSE)
cat("Pathway enrichment results saved.\n")