# Download TCGA-HNSC RNA-seq + clinical data (GDC, via TCGAbiolinks)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) BiocManager::install("TCGAbiolinks")
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) BiocManager::install("SummarizedExperiment")

library(TCGAbiolinks)
library(SummarizedExperiment)

dir.create("GDC_HNSC", showWarnings = FALSE)

# Workaround: TCGAbiolinks 2.40.0 GDCdownload.by.chunk bug
## GDCdownload.aux() doesn't return a value on success, so the caller's
## `if (ret == 1) break` check crashes on NULL. Errors already raise via
## stop(), so it's safe to call once per chunk and let them propagate.
patched_chunk <- function(server = "https://api.gdc.cancer.gov/data/",
                           manifest, name = "TCGAbiolinks_download",
                           path = ".", step = 1) {
  for (idx in 0:ceiling(nrow(manifest) / step - 1)) {
    end <- ifelse(((idx + 1) * step) > nrow(manifest), nrow(manifest), ((idx + 1) * step))
    manifest.aux <- manifest[((idx * step) + 1):end, ]
    size <- TCGAbiolinks:::humanReadableByteCount(sum(as.numeric(manifest.aux$size)))
    if (nrow(manifest.aux) > 1) {
      name.aux <- gsub("\\.tar", paste0("_", idx, ".tar"), name)
    } else {
      name.aux <- manifest.aux$filename
    }
    message(paste0("Downloading chunk ", idx + 1, " of ", ceiling(nrow(manifest) / step),
                    " (", nrow(manifest.aux), " files, size = ", size, ") ", "as ", name.aux))
    TCGAbiolinks:::GDCdownload.aux(server, manifest.aux, name.aux, path)
  }
}
assignInNamespace("GDCdownload.by.chunk", patched_chunk, ns = "TCGAbiolinks")

# 1. RNA-seq expression (STAR-Counts, GDC harmonised, hg38)
query_expr <- GDCquery(
  project       = "TCGA-HNSC",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query_expr, directory = "GDC_HNSC", files.per.chunk = 10)
hnsc_expr <- GDCprepare(query_expr, directory = "GDC_HNSC")

saveRDS(hnsc_expr, "GDC_HNSC/HNSC_RNAseq_SummarizedExperiment.rds")

## Raw counts and gene metadata, pulled out for convenience
counts <- assay(hnsc_expr, "unstranded")
gene_info <- as.data.frame(rowData(hnsc_expr))
write.table(
  cbind(gene_id = rownames(counts), counts),
  "GDC_HNSC/HNSC_RNAseq_counts.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 2. GDC clinical data
clinical_gdc <- GDCquery_clinic(project = "TCGA-HNSC", type = "clinical")
write.table(clinical_gdc, "GDC_HNSC/HNSC_clinical_GDC.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nDone. Outputs in ./GDC_HNSC:\n",
    "- HNSC_RNAseq_SummarizedExperiment.rds (full SE: counts, gene/sample metadata)\n",
    "- HNSC_RNAseq_counts.tsv (raw counts matrix)\n",
    "- HNSC_clinical_GDC.tsv (standard GDC clinical fields)\n")
cat("HPV status: use the cBioPortal PanCancer Atlas 2018 'Subtype' column",
    "(see import_HNSC_clinical_2018.R) instead of TCGAbiolinks' curated subtype table.\n")
