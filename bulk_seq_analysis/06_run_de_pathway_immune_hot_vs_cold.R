# Differential expression + pathway analysis: Immune-hot vs Immune-cold (TCGA-HNSC)
## DESeq2 on raw counts, using ImmuneGroup labels from
## 04_plot_immune_classification.R (MCP-counter clustering, ESTIMATE-named)

library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

dir.create("outputs/tcga_bulk_seq_outputs", recursive = TRUE, showWarnings = FALSE)

# 1. Load raw counts and gene symbol map
se     <- readRDS("GDC_HNSC/HNSC_RNAseq_SummarizedExperiment.rds")
counts <- assay(se, "unstranded")
gene_symbols <- rowData(se)$gene_name

## Tumor samples only, trimmed to 12-char patient barcode (matches immune_group_labels.csv)
sample_type <- substr(colnames(counts), 14, 15)
counts <- counts[, sample_type == "01"]
colnames(counts) <- substr(colnames(counts), 1, 12)

## Drop any duplicate patient barcodes (keep first aliquot) to guarantee unique columns
dup_sample <- duplicated(colnames(counts))
if (any(dup_sample)) {
  cat("Dropping", sum(dup_sample), "duplicate tumour aliquot(s) for patients already represented\n")
  counts <- counts[, !dup_sample]
}

cat("Raw counts matrix:", nrow(counts), "genes x", ncol(counts), "tumour samples\n")

# 2. Attach ImmuneGroup labels
labels <- read.csv("outputs/tcga_bulk_seq_outputs/immune_group_labels.csv")
common <- intersect(colnames(counts), labels$sample)
cat("Samples with both counts and ImmuneGroup label:", length(common), "of", ncol(counts), "\n")

counts    <- counts[, common]
col_data  <- labels[match(common, labels$sample), c("sample", "ImmuneGroup")]
rownames(col_data) <- col_data$sample
col_data$ImmuneGroup <- factor(col_data$ImmuneGroup, levels = c("Immune-cold", "Immune-hot"))
stopifnot(identical(colnames(counts), rownames(col_data)))

cat("\nGroup sizes going into DESeq2:\n")
print(table(col_data$ImmuneGroup))

# 3. Filter low-count genes and run DESeq2
keep   <- rowSums(counts) >= 10
counts <- counts[keep, ]
gene_symbols <- gene_symbols[keep]
cat("\nGenes retained after low-count filter (rowSums >= 10):", nrow(counts), "\n")

dds <- DESeqDataSetFromMatrix(countData = round(counts),
                               colData   = col_data,
                               design    = ~ ImmuneGroup)

checkpoint <- "outputs/tcga_bulk_seq_outputs/de_immune_hot_vs_cold_DESeq2_checkpoint.rds"
if (file.exists(checkpoint)) {
  cat("Checkpoint found — loading DESeq2 fit\n")
  dds <- readRDS(checkpoint)
} else {
  cat("Running DESeq2 (this may take a few minutes for", ncol(counts), "samples)...\n")
  dds <- DESeq(dds)
  saveRDS(dds, checkpoint)
  cat("Saved checkpoint:", checkpoint, "\n")
}

# 4. Extract results: Immune-hot vs Immune-cold
## apeglm (Zhu, Ibrahim & Love 2018) gives a more accurate posterior LFC
## estimator than type="normal" and is the DESeq2-recommended default for
## shrinkage. It requires the `res` passed in to be built via `name=` (not
## `contrast=`) matching its `coef=` exactly, or lfcShrink errors out even
## when the two are algebraically identical (as here, since "Immune-cold"
## is the reference level so name and contrast agree).
coef_name <- resultsNames(dds)[grepl("ImmuneGroup", resultsNames(dds))]
res <- results(dds, name = coef_name, alpha = 0.05)
res_stat <- res$stat  ## Wald stat, saved before lfcShrink (which drops it) -- used as the GSEA ranking metric below
res <- lfcShrink(dds, coef = coef_name, res = res, type = "apeglm")

res_df <- as.data.frame(res)
res_df$ensembl_id <- rownames(res_df)
res_df$gene       <- gene_symbols[match(res_df$ensembl_id, rownames(counts))]
res_df$stat       <- res_stat
res_df <- res_df[, c("ensembl_id", "gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
res_df <- res_df[order(res_df$padj), ]

## Positive log2FoldChange = higher in Immune-hot; negative = higher in Immune-cold
sig_hot_up  <- res_df$padj < 0.05 & !is.na(res_df$padj) & res_df$log2FoldChange >  0.5
sig_cold_up <- res_df$padj < 0.05 & !is.na(res_df$padj) & res_df$log2FoldChange < -0.5
res_df$sig <- "NS"
res_df$sig[sig_hot_up]  <- "Higher in Immune-hot"
res_df$sig[sig_cold_up] <- "Higher in Immune-cold"

cat("\n=== DE summary: Immune-hot vs Immune-cold ===\n")
cat("Genes tested:", nrow(res_df), "\n")
cat("Significant (padj<0.05):", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("Higher in Immune-hot  (padj<0.05 & log2FC>0.5):", sum(sig_hot_up), "\n")
cat("Higher in Immune-cold (padj<0.05 & log2FC<-0.5):", sum(sig_cold_up), "\n")

cat("\nTop 15 genes higher in Immune-hot:\n")
print(head(res_df[order(-res_df$log2FoldChange), c("gene", "log2FoldChange", "padj")], 15))
cat("\nTop 15 genes higher in Immune-cold:\n")
print(head(res_df[order(res_df$log2FoldChange), c("gene", "log2FoldChange", "padj")], 15))

write.csv(res_df, "outputs/tcga_bulk_seq_outputs/DE_ImmuneHot_vs_ImmuneCold.csv", row.names = FALSE)
cat("\nFull DE table saved: outputs/tcga_bulk_seq_outputs/DE_ImmuneHot_vs_ImmuneCold.csv\n")

## Note: volcano plot is generated by 07_visualise_de_immune_hot_vs_cold.R
## (padj-ranked gene labels), run that script next to produce it.

# 5. Pathway enrichment (GO BP + KEGG) on the DE genes
## clusterProfiler conventions: org.Hs.eg.db, SYMBOL input for GO,
## SYMBOL->ENTREZID via bitr for KEGG, BH-adjusted, readable GO terms.
symbol_to_entrez <- function(genes) {
  res <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = org.Hs.eg.db, drop = TRUE)
  res$ENTREZID
}

cat("\nRunning GO/KEGG enrichment on Immune-hot vs Immune-cold DE genes...\n")
hot_genes  <- res_df %>% filter(sig == "Higher in Immune-hot")  %>% pull(gene)
cold_genes <- res_df %>% filter(sig == "Higher in Immune-cold") %>% pull(gene)
cat("Hot-up genes:", length(hot_genes), " | Cold-up genes:", length(cold_genes), "\n")

go_kegg_checkpoint <- "outputs/tcga_bulk_seq_outputs/go_kegg_ImmuneHotCold_checkpoint.rds"
if (file.exists(go_kegg_checkpoint)) {
  cat("Checkpoint found — loading go_hot/go_cold/kegg_hot/kegg_cold\n")
  ck       <- readRDS(go_kegg_checkpoint)
  go_hot   <- ck$go_hot
  go_cold  <- ck$go_cold
  kegg_hot <- ck$kegg_hot
  kegg_cold<- ck$kegg_cold
} else {
  go_hot <- enrichGO(gene = hot_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                     ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2, readable = TRUE)

  go_cold <- enrichGO(gene = cold_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                      ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2, readable = TRUE)

  kegg_hot <- enrichKEGG(gene = symbol_to_entrez(hot_genes), organism = "hsa",
                         pAdjustMethod = "BH", pvalueCutoff = 0.05)

  kegg_cold <- enrichKEGG(gene = symbol_to_entrez(cold_genes), organism = "hsa",
                          pAdjustMethod = "BH", pvalueCutoff = 0.05)

  saveRDS(list(go_hot = go_hot, go_cold = go_cold, kegg_hot = kegg_hot, kegg_cold = kegg_cold),
          go_kegg_checkpoint)
  cat("Saved checkpoint:", go_kegg_checkpoint, "\n")
}

## showCategory=15 (not 20) + taller canvas + wrapped labels -- GO BP term names
## are long and at height=8/20 categories they overlapped illegibly
wrap_labels <- function(width = 45) function(x) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

if (!is.null(go_hot) && nrow(go_hot) > 0) {
  p <- dotplot(go_hot, showCategory = 15) + ggtitle("GO BP: Higher in Immune-hot") +
    scale_y_discrete(labels = wrap_labels())
  ggsave("outputs/tcga_bulk_seq_outputs/GO_BP_ImmuneHot_upregulated.pdf", p, width = 9, height = 9)
  ggsave("outputs/tcga_bulk_seq_outputs/GO_BP_ImmuneHot_upregulated.png", p, width = 9, height = 9, dpi = 300)
}
if (!is.null(go_cold) && nrow(go_cold) > 0) {
  p <- dotplot(go_cold, showCategory = 15) + ggtitle("GO BP: Higher in Immune-cold") +
    scale_y_discrete(labels = wrap_labels())
  ggsave("outputs/tcga_bulk_seq_outputs/GO_BP_ImmuneCold_upregulated.pdf", p, width = 9, height = 9)
  ggsave("outputs/tcga_bulk_seq_outputs/GO_BP_ImmuneCold_upregulated.png", p, width = 9, height = 9, dpi = 300)
}
if (!is.null(kegg_hot) && nrow(kegg_hot) > 0) {
  p <- dotplot(kegg_hot, showCategory = 15) + ggtitle("KEGG: Higher in Immune-hot") +
    scale_y_discrete(labels = wrap_labels())
  ggsave("outputs/tcga_bulk_seq_outputs/KEGG_ImmuneHot_upregulated.pdf", p, width = 9, height = 9)
  ggsave("outputs/tcga_bulk_seq_outputs/KEGG_ImmuneHot_upregulated.png", p, width = 9, height = 9, dpi = 300)
}
if (!is.null(kegg_cold) && nrow(kegg_cold) > 0) {
  p <- dotplot(kegg_cold, showCategory = 15) + ggtitle("KEGG: Higher in Immune-cold") +
    scale_y_discrete(labels = wrap_labels())
  ggsave("outputs/tcga_bulk_seq_outputs/KEGG_ImmuneCold_upregulated.pdf", p, width = 9, height = 9)
  ggsave("outputs/tcga_bulk_seq_outputs/KEGG_ImmuneCold_upregulated.png", p, width = 9, height = 9, dpi = 300)
}

if (!is.null(go_hot) && nrow(go_hot) > 0) {
  go_hot_sim <- pairwise_termsim(go_hot)
  p_emap <- emapplot(go_hot_sim, showCategory = 30) + ggtitle("GO BP Enrichment Map: Higher in Immune-hot")
  ggsave("outputs/tcga_bulk_seq_outputs/GO_BP_emapplot_ImmuneHot.pdf", p_emap, width = 10, height = 9)
  ggsave("outputs/tcga_bulk_seq_outputs/GO_BP_emapplot_ImmuneHot.png", p_emap, width = 10, height = 9, dpi = 300)
}
if (!is.null(go_cold) && nrow(go_cold) > 0) {
  go_cold_sim <- pairwise_termsim(go_cold)
  p_emap <- emapplot(go_cold_sim, showCategory = 30) + ggtitle("GO BP Enrichment Map: Higher in Immune-cold")
  ggsave("outputs/tcga_bulk_seq_outputs/GO_BP_emapplot_ImmuneCold.pdf", p_emap, width = 10, height = 9)
  ggsave("outputs/tcga_bulk_seq_outputs/GO_BP_emapplot_ImmuneCold.png", p_emap, width = 10, height = 9, dpi = 300)
}

write.csv(as.data.frame(go_hot),   file = "outputs/tcga_bulk_seq_outputs/GO_BP_ImmuneHot_upregulated.csv",  row.names = FALSE)
write.csv(as.data.frame(go_cold),  file = "outputs/tcga_bulk_seq_outputs/GO_BP_ImmuneCold_upregulated.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_hot), file = "outputs/tcga_bulk_seq_outputs/KEGG_ImmuneHot_upregulated.csv",   row.names = FALSE)
write.csv(as.data.frame(kegg_cold),file = "outputs/tcga_bulk_seq_outputs/KEGG_ImmuneCold_upregulated.csv",  row.names = FALSE)

cat("\n=== Pathway enrichment summary ===\n")
cat("GO BP terms (Immune-hot):", ifelse(is.null(go_hot), 0, nrow(go_hot)), "\n")
cat("GO BP terms (Immune-cold):", ifelse(is.null(go_cold), 0, nrow(go_cold)), "\n")
cat("KEGG pathways (Immune-hot):", ifelse(is.null(kegg_hot), 0, nrow(kegg_hot)), "\n")
cat("KEGG pathways (Immune-cold):", ifelse(is.null(kegg_cold), 0, nrow(kegg_cold)), "\n")
cat("\nPathway enrichment results saved to outputs/tcga_bulk_seq_outputs/\n")

# 6. GSEA (GO BP + KEGG) on the full ranked gene list
## Unlike the ORA above (hard padj/LFC cutoff, hot-up/cold-up tested
## separately), GSEA ranks every gene by the DESeq2 Wald stat and finds
## pathways skewed toward one end -- picks up coordinated shifts too subtle
## for the ORA cutoff, both directions in a single run via the enrichment
## score's sign. Ranking metric: Wald stat (weights log2FC by its SE).
## Positive NES = enriched toward Immune-hot; negative = Immune-cold.
cat("\nRunning GSEA (GO BP + KEGG) on the full ranked gene list...\n")

## One rank per gene symbol: keep the strongest |stat| when a symbol maps to
## multiple Ensembl IDs, then sort descending (gseGO/gseKEGG require unique,
## sorted-decreasing named vectors).
gsea_ranks <- res_df %>% filter(!is.na(stat), !is.na(gene), gene != "")
gsea_ranks <- gsea_ranks[order(-abs(gsea_ranks$stat)), ]
gsea_ranks <- gsea_ranks[!duplicated(gsea_ranks$gene), ]
gsea_ranks <- gsea_ranks[order(-gsea_ranks$stat), ]
gene_rank_go <- setNames(gsea_ranks$stat, gsea_ranks$gene)
cat("Ranked gene list for GSEA:", length(gene_rank_go), "unique gene symbols\n")

gsea_checkpoint <- "outputs/tcga_bulk_seq_outputs/gsea_ImmuneHotCold_checkpoint.rds"
if (file.exists(gsea_checkpoint)) {
  cat("Checkpoint found — loading gsea_go/gsea_kegg\n")
  ck        <- readRDS(gsea_checkpoint)
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

  saveRDS(list(gsea_go = gsea_go, gsea_kegg = gsea_kegg), gsea_checkpoint)
  cat("Saved checkpoint:", gsea_checkpoint, "\n")
}

## split=".sign" + facet puts Immune-hot (activated) and Immune-cold
## (suppressed) enrichments in separate panels, same idea as the hot/cold
## dotplot pairs above but from one GSEA run instead of two ORA runs.
if (!is.null(gsea_go) && nrow(gsea_go) > 0) {
  p <- dotplot(gsea_go, showCategory = 15, split = ".sign") +
    facet_grid(. ~ .sign, labeller = as_labeller(c(activated = "Higher in Immune-hot",
                                                     suppressed = "Higher in Immune-cold"))) +
    ggtitle("GSEA GO BP: Immune-hot vs Immune-cold") +
    scale_y_discrete(labels = wrap_labels())
  ggsave("outputs/tcga_bulk_seq_outputs/GSEA_GO_BP_ImmuneHotCold.pdf", p, width = 13, height = 9)
  ggsave("outputs/tcga_bulk_seq_outputs/GSEA_GO_BP_ImmuneHotCold.png", p, width = 13, height = 9, dpi = 300)
}
if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
  p <- dotplot(gsea_kegg, showCategory = 15, split = ".sign") +
    facet_grid(. ~ .sign, labeller = as_labeller(c(activated = "Higher in Immune-hot",
                                                     suppressed = "Higher in Immune-cold"))) +
    ggtitle("GSEA KEGG: Immune-hot vs Immune-cold") +
    scale_y_discrete(labels = wrap_labels())
  ggsave("outputs/tcga_bulk_seq_outputs/GSEA_KEGG_ImmuneHotCold.pdf", p, width = 13, height = 9)
  ggsave("outputs/tcga_bulk_seq_outputs/GSEA_KEGG_ImmuneHotCold.png", p, width = 13, height = 9, dpi = 300)
}

## Classic GSEA running-enrichment-score plot for the single strongest hit at
## each end of the ranking (most Immune-hot-skewed, most Immune-cold-skewed).
if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
  kegg_ord    <- gsea_kegg@result[order(-gsea_kegg@result$NES), ]
  top_hot_id  <- kegg_ord$ID[1]
  top_cold_id <- kegg_ord$ID[nrow(kegg_ord)]

  p_hot <- gseaplot2(gsea_kegg, geneSetID = top_hot_id,
                      title = paste0("Top Immune-hot pathway: ", kegg_ord$Description[1]))
  ggsave("outputs/tcga_bulk_seq_outputs/GSEA_KEGG_top_ImmuneHot_gseaplot.pdf", p_hot, width = 8, height = 7)
  ggsave("outputs/tcga_bulk_seq_outputs/GSEA_KEGG_top_ImmuneHot_gseaplot.png", p_hot, width = 8, height = 7, dpi = 300)

  p_cold <- gseaplot2(gsea_kegg, geneSetID = top_cold_id,
                       title = paste0("Top Immune-cold pathway: ", kegg_ord$Description[nrow(kegg_ord)]))
  ggsave("outputs/tcga_bulk_seq_outputs/GSEA_KEGG_top_ImmuneCold_gseaplot.pdf", p_cold, width = 8, height = 7)
  ggsave("outputs/tcga_bulk_seq_outputs/GSEA_KEGG_top_ImmuneCold_gseaplot.png", p_cold, width = 8, height = 7, dpi = 300)
}

write.csv(as.data.frame(gsea_go),   file = "outputs/tcga_bulk_seq_outputs/GSEA_GO_BP_ImmuneHotCold.csv",   row.names = FALSE)
write.csv(as.data.frame(gsea_kegg), file = "outputs/tcga_bulk_seq_outputs/GSEA_KEGG_ImmuneHotCold.csv",    row.names = FALSE)

cat("\n=== GSEA summary ===\n")
cat("GO BP gene sets (padj<0.05):",   ifelse(is.null(gsea_go),   0, nrow(gsea_go)),   "\n")
cat("KEGG gene sets (padj<0.05):",    ifelse(is.null(gsea_kegg), 0, nrow(gsea_kegg)), "\n")
cat("\nGSEA results saved to outputs/tcga_bulk_seq_outputs/\n")
