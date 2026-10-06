# Immune infiltration estimation (ESTIMATE, MCP-counter, xCell) — TCGA-HNSC bulk RNA-seq (TPM)

# Install packages if missing
if (!requireNamespace("remotes",    quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("estimate",   quietly = TRUE))
  install.packages("estimate", repos = "http://r-forge.r-project.org")
if (!requireNamespace("MCPcounter", quietly = TRUE))
  remotes::install_github("ebecht/MCPcounter", subdir = "Source")
if (!requireNamespace("xCell",      quietly = TRUE))
  remotes::install_github("dviraran/xCell")

library(SummarizedExperiment)
library(estimate)
library(MCPcounter)
library(xCell)

dir.create("outputs/tcga_bulk_seq_outputs", recursive = TRUE, showWarnings = FALSE)

# Load TPM matrix, map to gene symbols
se  <- readRDS("GDC_HNSC/HNSC_RNAseq_SummarizedExperiment.rds")
tpm <- assay(se, "tpm_unstrand")

## Map Ensembl IDs -> gene symbols; drop rows with no symbol
symbols <- rowData(se)$gene_name
keep     <- !is.na(symbols) & symbols != ""
tpm      <- tpm[keep, ]
rownames(tpm) <- symbols[keep]

## Remove duplicate gene symbols (keep highest mean expression)
dup <- duplicated(rownames(tpm))
if (any(dup)) {
  means    <- rowMeans(tpm)
  tpm      <- tpm[order(rownames(tpm), -means), ]
  tpm      <- tpm[!duplicated(rownames(tpm)), ]
}

## Keep tumour samples only (sample type "01" at barcode positions 14-15)
## Avoids duplicate patient IDs after trimming (tumour -01A vs normal -11A)
sample_type <- substr(colnames(tpm), 14, 15)
tpm <- tpm[, sample_type == "01"]
colnames(tpm) <- substr(colnames(tpm), 1, 12)

cat("TPM matrix:", nrow(tpm), "genes x", ncol(tpm), "samples (tumour only)\n")

# 1. ESTIMATE
## ESTIMATE uses a file-based GCT interface
gct_filtered <- "outputs/tcga_bulk_seq_outputs/tpm_filtered.gct"
gct_scores   <- "outputs/tcga_bulk_seq_outputs/estimate_scores.gct"

## Write GCT (version 1.2 format) — single write to avoid mixed line endings
write_gct <- function(mat, file) {
  tmp <- tempfile()
  df  <- data.frame(Name = rownames(mat), Description = rownames(mat),
                    mat, check.names = FALSE, stringsAsFactors = FALSE)
  write.table(df, tmp, sep = "\t", quote = FALSE, row.names = FALSE)
  body <- readLines(tmp)
  unlink(tmp)
  writeLines(c("#1.2", paste(nrow(mat), ncol(mat), sep = "\t"), body), file)
}

## Pre-filter to ESTIMATE's common gene list (bypasses filterCommonGenes,
## avoids writing a 59k-gene GCT that causes read.table issues on Windows)
data("common_genes", package = "estimate")
tpm_est <- tpm[rownames(tpm) %in% common_genes$GeneSymbol, ]
cat("ESTIMATE: matched", nrow(tpm_est), "of", nrow(common_genes), "common genes\n")
write_gct(tpm_est, gct_filtered)
estimateScore(gct_filtered, gct_scores, platform = "illumina")

## Parse GCT output into a tidy data frame
parse_gct <- function(file) {
  lines  <- readLines(file)
  header <- strsplit(lines[3], "\t")[[1]]
  mat    <- do.call(rbind, lapply(lines[-(1:3)], function(l) strsplit(l, "\t")[[1]]))
  df     <- as.data.frame(t(mat[, -(1:2)]), stringsAsFactors = FALSE)
  colnames(df) <- mat[, 1]
  df$sample     <- header[-(1:2)]
  score_cols    <- setdiff(names(df), "sample")
  df[, score_cols] <- lapply(df[, score_cols], as.numeric)
  df
}

estimate_scores <- parse_gct(gct_scores)
estimate_scores$sample <- gsub("\\.", "-", estimate_scores$sample)
write.csv(estimate_scores, "outputs/tcga_bulk_seq_outputs/ESTIMATE_scores.csv", row.names = FALSE)
cat("ESTIMATE done:", nrow(estimate_scores), "samples\n")

# 2. MCP-counter
mcp_scores <- MCPcounter.estimate(
  expression    = tpm,
  featuresType  = "HUGO_symbols"
)
mcp_df <- as.data.frame(t(mcp_scores))
mcp_df$sample <- rownames(mcp_df)
write.csv(mcp_df, "outputs/tcga_bulk_seq_outputs/MCPcounter_scores.csv", row.names = FALSE)
cat("MCP-counter done:", nrow(mcp_df), "samples,", ncol(mcp_df) - 1, "cell types\n")

# 3. xCell
xcell_scores <- xCellAnalysis(tpm)
xcell_df     <- as.data.frame(t(xcell_scores))
xcell_df$sample <- rownames(xcell_df)
write.csv(xcell_df, "outputs/tcga_bulk_seq_outputs/xCell_scores.csv", row.names = FALSE)
cat("xCell done:", nrow(xcell_df), "samples,", ncol(xcell_df) - 1, "cell types\n")

# Save combined summary
## Add method prefix to avoid duplicate column names (MCP-counter and xCell share some cell types)
names(mcp_df)   <- ifelse(names(mcp_df)   == "sample", "sample", paste0("MCP.",   names(mcp_df)))
names(xcell_df) <- ifelse(names(xcell_df) == "sample", "sample", paste0("xCell.", names(xcell_df)))
combined <- Reduce(function(a, b) merge(a, b, by = "sample"), list(estimate_scores, mcp_df, xcell_df))
write.csv(combined, "outputs/tcga_bulk_seq_outputs/all_deconvolution_scores.csv", row.names = FALSE)

cat("\nOutputs saved to outputs/tcga_bulk_seq_outputs/\n")
cat("  ESTIMATE_scores.csv     - ImmuneScore, StromalScore, ESTIMATEScore, TumorPurity\n")
cat("  MCPcounter_scores.csv   - 10 immune/stromal cell populations\n")
cat("  xCell_scores.csv        - 64 cell types + ImmuneScore, StromaScore\n")
cat("  all_deconvolution_scores.csv - merged\n")
