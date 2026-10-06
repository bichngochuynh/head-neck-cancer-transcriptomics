# Genes correlated with continuous ESTIMATE ImmuneScore (Spearman)
## Complements the discrete hot/cold DESeq2 test by treating the labelling
## variable as continuous, capturing the hot<->cold gradient and the
## ambiguous middle-zone samples seen in the classification scatter plot.

library(DESeq2)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(ggrepel)

dir.create("outputs/tcga_bulk_seq_outputs", recursive = TRUE, showWarnings = FALSE)

# 1. Reuse the cached DESeq2 fit (same samples as the hot vs cold DE)
dds <- readRDS("outputs/tcga_bulk_seq_outputs/de_immune_hot_vs_cold_DESeq2_checkpoint.rds")

# 2. Variance-stabilised expression (for correlation across many genes)
cat("Running variance-stabilising transformation...\n")
vsd  <- vst(dds, blind = FALSE)
expr <- assay(vsd)   # genes x samples
cat("Expression matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")

# 3. Attach continuous ESTIMATE ImmuneScore for the same samples
scores <- read.csv("outputs/tcga_bulk_seq_outputs/all_deconvolution_scores.csv", check.names = FALSE)
immune_score <- scores$ImmuneScore[match(colnames(expr), scores$sample)]
stopifnot(!any(is.na(immune_score)))

# 4. Map Ensembl IDs -> gene symbols
se <- readRDS("GDC_HNSC/HNSC_RNAseq_SummarizedExperiment.rds")
gene_symbols <- rowData(se)$gene_name[match(rownames(expr), rownames(se))]

# 5. Spearman correlation for every gene, vectorised
## Rank-transform expression and score, then Pearson-correlate the ranks
## (mathematically identical to Spearman rho) for all ~54k genes in one matrix
## operation, instead of looping cor.test() gene-by-gene (prohibitively slow).
n          <- length(immune_score)
rank_score <- rank(immune_score)
rank_expr  <- t(apply(expr, 1, rank))   # genes x samples, ranked within each gene

rho <- as.numeric(cor(t(rank_expr), rank_score))

## p-value via the standard t-approximation -- the same fallback cor.test()
## uses for Spearman when ties are present / n is large (our case: n=520, and
## RNA-seq expression routinely has ties even after vst).
t_stat <- rho * sqrt((n - 2) / (1 - rho^2))
p_val  <- 2 * pt(-abs(t_stat), df = n - 2)
padj   <- p.adjust(p_val, method = "BH")

result <- data.frame(
  ensembl_id = rownames(expr),
  gene       = gene_symbols,
  rho        = rho,
  p_value    = p_val,
  padj       = padj,
  stringsAsFactors = FALSE
)

# 6. Cross-check against the discrete Immune-hot vs Immune-cold DE test
de_discrete <- read.csv("outputs/tcga_bulk_seq_outputs/DE_ImmuneHot_vs_ImmuneCold.csv")
result$also_sig_in_discrete_DE <- result$ensembl_id %in%
  de_discrete$ensembl_id[de_discrete$sig != "NS"]

result <- result[order(-result$rho), ]
write.csv(result, "outputs/tcga_bulk_seq_outputs/GeneCorrelation_with_ImmuneScore.csv", row.names = FALSE)

# 7. Summary
sig_pos <- result %>% filter(padj < 0.05, rho > 0) %>% arrange(desc(rho))
sig_neg <- result %>% filter(padj < 0.05, rho < 0) %>% arrange(rho)

cat("\n=== Spearman correlation with ImmuneScore: summary ===\n")
cat("Genes tested:", nrow(result), "\n")
cat("Significantly POSITIVELY correlated (padj<0.05, rho>0):", nrow(sig_pos), "\n")
cat("Significantly NEGATIVELY correlated (padj<0.05, rho<0):", nrow(sig_neg), "\n")
cat("Of those, also significant in the discrete Hot-vs-Cold DESeq2 test:",
    sum(c(sig_pos$also_sig_in_discrete_DE, sig_neg$also_sig_in_discrete_DE)),
    "/", nrow(sig_pos) + nrow(sig_neg),
    "-- genes flagged by ONE test but not the other are worth a closer look\n")

cat("\nTop 15 genes positively correlated with ImmuneScore:\n")
print(head(sig_pos[, c("gene", "rho", "padj", "also_sig_in_discrete_DE")], 15))
cat("\nTop 15 genes negatively correlated with ImmuneScore:\n")
print(head(sig_neg[, c("gene", "rho", "padj", "also_sig_in_discrete_DE")], 15))

cat("\nFull table saved: outputs/tcga_bulk_seq_outputs/GeneCorrelation_with_ImmuneScore.csv\n")

# 8. Visual check — expression vs ImmuneScore for top 3 genes each direction
labels_df <- read.csv("outputs/tcga_bulk_seq_outputs/immune_group_labels.csv")
group_colors <- c("Immune-cold" = "#2166AC", "Immune-hot" = "#D6604D")

top_genes <- rbind(head(sig_pos, 3), head(sig_neg, 3))

plot_df <- lapply(seq_len(nrow(top_genes)), function(i) {
  eid <- top_genes$ensembl_id[i]
  data.frame(
    sample       = colnames(expr),
    expression   = expr[eid, ],
    ImmuneScore  = immune_score,
    gene         = top_genes$gene[i],
    rho          = top_genes$rho[i],
    padj         = top_genes$padj[i]
  )
}) %>% bind_rows() %>%
  left_join(labels_df[, c("sample", "ImmuneGroup")], by = "sample") %>%
  mutate(gene_label = sprintf("%s (rho=%.2f, padj=%.1e)", gene, rho, padj),
         gene_label = factor(gene_label, levels = unique(gene_label)))

p_corr <- ggplot(plot_df, aes(x = ImmuneScore, y = expression)) +
  geom_point(aes(colour = ImmuneGroup), alpha = 0.5, size = 1) +
  geom_smooth(method = "loess", colour = "black", linewidth = 0.6, se = TRUE) +
  scale_colour_manual(values = group_colors, name = "Immune group") +
  facet_wrap(~ gene_label, scales = "free_y", ncol = 3) +
  labs(
    title = "Top genes correlated with ESTIMATE ImmuneScore (Spearman)",
    subtitle = "Continuous view -- includes the ambiguous middle-zone samples, not just discrete hot/cold",
    x = "ESTIMATE ImmuneScore", y = "VST-normalised expression"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 9.5),
        strip.background = element_rect(fill = "grey90", colour = NA))

ggsave("outputs/tcga_bulk_seq_outputs/GeneCorrelation_top_genes_scatter.pdf", p_corr, width = 12, height = 8)
ggsave("outputs/tcga_bulk_seq_outputs/GeneCorrelation_top_genes_scatter.png", p_corr, width = 12, height = 8, dpi = 300)
cat("Scatter panel saved: outputs/tcga_bulk_seq_outputs/GeneCorrelation_top_genes_scatter.pdf/png\n")

# 9. Three-region split of the classification scatter
## Hot-only/Cold-only = inside that group's 90% ellipse, outside the other's;
## Overlap = inside both (ambiguous middle zone); Outside both = neither.
## Ellipse = 90% bivariate-normal region (Mahalanobis distance vs chi-square),
## same level as stat_ellipse in 04_plot_immune_classification.R.
scores_xy <- scores[match(colnames(expr), scores$sample), ]  # same sample order as expr/dds
scores_xy$ImmuneGroup <- labels_df$ImmuneGroup[match(scores_xy$sample, labels_df$sample)]
scores_xy$ImmuneGroup <- factor(scores_xy$ImmuneGroup, levels = c("Immune-cold", "Immune-hot"))

ellipse_membership <- function(data_xy, ref_mask, level = 0.90) {
  mu    <- colMeans(data_xy[ref_mask, ])
  Sigma <- cov(data_xy[ref_mask, ])
  d2    <- mahalanobis(data_xy, center = mu, cov = Sigma)
  d2 <= qchisq(level, df = 2)
}

xy <- scores_xy[, c("ImmuneScore", "StromalScore")]
in_hot_ellipse  <- ellipse_membership(xy, scores_xy$ImmuneGroup == "Immune-hot")
in_cold_ellipse <- ellipse_membership(xy, scores_xy$ImmuneGroup == "Immune-cold")

scores_xy$Region <- case_when(
  in_hot_ellipse  & in_cold_ellipse  ~ "Overlap",
  in_hot_ellipse  & !in_cold_ellipse ~ "Hot-only",
  !in_hot_ellipse & in_cold_ellipse  ~ "Cold-only",
  TRUE                                ~ "Outside both ellipses"
)
scores_xy$Region <- factor(scores_xy$Region,
                            levels = c("Cold-only", "Overlap", "Hot-only", "Outside both ellipses"))

cat("\nRegion sizes:\n")
print(table(scores_xy$Region))
cat("\nRegion x ImmuneGroup cross-tab:\n")
print(table(scores_xy$Region, scores_xy$ImmuneGroup))

write.csv(scores_xy[, c("sample", "ImmuneGroup", "ImmuneScore", "StromalScore", "Region")],
          "outputs/tcga_bulk_seq_outputs/scatter_region_labels.csv", row.names = FALSE)
cat("Region labels saved: outputs/tcga_bulk_seq_outputs/scatter_region_labels.csv\n")

# New scatter plot, coloured by region
region_colors <- c("Hot-only" = "#D6604D", "Cold-only" = "#2166AC",
                    "Overlap" = "#762A83", "Outside both ellipses" = "grey70")

scatter_regions <- ggplot(scores_xy, aes(x = ImmuneScore, y = StromalScore, colour = Region)) +
  geom_point(alpha = 0.65, size = 1.4) +
  stat_ellipse(data = subset(scores_xy, ImmuneGroup == "Immune-hot"),
               aes(x = ImmuneScore, y = StromalScore), level = 0.90,
               colour = "#D6604D", linewidth = 0.6, inherit.aes = FALSE) +
  stat_ellipse(data = subset(scores_xy, ImmuneGroup == "Immune-cold"),
               aes(x = ImmuneScore, y = StromalScore), level = 0.90,
               colour = "#2166AC", linewidth = 0.6, inherit.aes = FALSE) +
  scale_colour_manual(values = region_colors, name = "Region") +
  labs(
    title = "Classification Scatter Split into Three Regions",
    subtitle = paste0("Hot-only n=", sum(scores_xy$Region == "Hot-only"),
                      "  |  Overlap n=", sum(scores_xy$Region == "Overlap"),
                      "  |  Cold-only n=", sum(scores_xy$Region == "Cold-only"),
                      "  |  Outside both n=", sum(scores_xy$Region == "Outside both ellipses"),
                      "\nRegions = 90% confidence ellipse membership per group (Mahalanobis distance)"),
    x = "ESTIMATE ImmuneScore", y = "ESTIMATE StromalScore"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 9))

ggsave("outputs/tcga_bulk_seq_outputs/scatter_three_regions.pdf", scatter_regions, width = 8, height = 6.5)
ggsave("outputs/tcga_bulk_seq_outputs/scatter_three_regions.png", scatter_regions, width = 8, height = 6.5, dpi = 300)
cat("New 3-region scatter saved: outputs/tcga_bulk_seq_outputs/scatter_three_regions.pdf/png\n")

# 10. DE between regions — markers of Hot-only, Cold-only, Overlap
## Reuses the raw counts inside the cached `dds` object, refit with a
## 3-level Region design.
keep_region <- scores_xy$Region != "Outside both ellipses"
cat("\nExcluding", sum(!keep_region), "sample(s) outside both ellipses from the region DE test\n")

region_counts <- counts(dds)[, keep_region]
region_coldata <- data.frame(
  sample = scores_xy$sample[keep_region],
  Region = droplevels(scores_xy$Region[keep_region])
)
rownames(region_coldata) <- region_coldata$sample
stopifnot(identical(colnames(region_counts), rownames(region_coldata)))

dds_region <- DESeqDataSetFromMatrix(countData = region_counts,
                                      colData   = region_coldata,
                                      design    = ~ Region)
dds_region <- dds_region[rowSums(counts(dds_region)) >= 10, ]

region_checkpoint <- "outputs/tcga_bulk_seq_outputs/de_scatter_regions_DESeq2_checkpoint.rds"
if (file.exists(region_checkpoint)) {
  cat("Checkpoint found — loading region DESeq2 fit\n")
  dds_region <- readRDS(region_checkpoint)
} else {
  cat("Running DESeq2 for the 3-region design (Hot-only / Overlap / Cold-only)...\n")
  dds_region <- DESeq(dds_region)
  saveRDS(dds_region, region_checkpoint)
  cat("Saved checkpoint:", region_checkpoint, "\n")
}

get_contrast <- function(grp1, grp2) {
  res <- results(dds_region, contrast = c("Region", grp1, grp2), alpha = 0.05)
  res <- lfcShrink(dds_region, contrast = c("Region", grp1, grp2), res = res, type = "normal")
  df  <- as.data.frame(res)
  df$ensembl_id <- rownames(df)
  df$gene       <- gene_symbols[match(df$ensembl_id, rownames(expr))]
  df$contrast   <- paste0(grp1, " vs ", grp2)
  df$sig        <- with(df, ifelse(!is.na(padj) & padj < 0.05 & log2FoldChange >  0.5, paste0("Higher in ", grp1),
                             ifelse(!is.na(padj) & padj < 0.05 & log2FoldChange < -0.5, paste0("Higher in ", grp2), "NS")))
  df[order(df$padj), c("contrast", "ensembl_id", "gene", "baseMean", "log2FoldChange", "padj", "sig")]
}

de_hot_vs_overlap  <- get_contrast("Hot-only", "Overlap")
de_cold_vs_overlap <- get_contrast("Cold-only", "Overlap")
de_hot_vs_cold     <- get_contrast("Hot-only", "Cold-only")

de_regions_all <- bind_rows(de_hot_vs_overlap, de_cold_vs_overlap, de_hot_vs_cold)
write.csv(de_regions_all, "outputs/tcga_bulk_seq_outputs/DE_three_regions_scatter.csv", row.names = FALSE)

cat("\n=== DE between scatter regions: summary ===\n")
for (ct in unique(de_regions_all$contrast)) {
  sub <- de_regions_all[de_regions_all$contrast == ct, ]
  cat(sprintf("%-28s significant genes: %d\n", ct, sum(sub$sig != "NS")))
}

cat("\nTop 10 genes marking the Overlap zone vs Hot-only (higher in Overlap):\n")
print(head(de_hot_vs_overlap[de_hot_vs_overlap$sig == "Higher in Overlap",
                              c("gene", "log2FoldChange", "padj")], 10))
cat("\nTop 10 genes marking the Overlap zone vs Cold-only (higher in Overlap):\n")
print(head(de_cold_vs_overlap[de_cold_vs_overlap$sig == "Higher in Overlap",
                               c("gene", "log2FoldChange", "padj")], 10))

cat("\nFull region DE table saved: outputs/tcga_bulk_seq_outputs/DE_three_regions_scatter.csv\n")

# 11. Marker-gene heatmap
## A gene counts as a region's marker if it's significantly higher there vs
## BOTH other regions (consistent across contrasts, not just one comparison).
library(pheatmap)
library(RColorBrewer)

hot_markers <- intersect(
  de_hot_vs_overlap$ensembl_id[de_hot_vs_overlap$sig == "Higher in Hot-only"],
  de_hot_vs_cold$ensembl_id[de_hot_vs_cold$sig == "Higher in Hot-only"]
)
cold_markers <- intersect(
  de_cold_vs_overlap$ensembl_id[de_cold_vs_overlap$sig == "Higher in Cold-only"],
  de_hot_vs_cold$ensembl_id[de_hot_vs_cold$sig == "Higher in Cold-only"]
)
overlap_markers <- intersect(
  de_hot_vs_overlap$ensembl_id[de_hot_vs_overlap$sig == "Higher in Overlap"],
  de_cold_vs_overlap$ensembl_id[de_cold_vs_overlap$sig == "Higher in Overlap"]
)

cat("\nConsistent markers -- Hot-only:", length(hot_markers),
    "| Cold-only:", length(cold_markers),
    "| Overlap:", length(overlap_markers), "\n")

## Rank each region's markers by significance (this is what worked cleanly for
## Hot-only before), but FIRST drop low-expression/sparse genes (baseMean < 50)
## from the candidate pool. Sparse genes (near-zero in most samples, high in a
## few) produce inflated fold-change estimates that look dramatic but are
## driven by outliers -- ranking by raw |log2FC| alone (tried and reverted)
## surfaced exactly these genes (e.g. KRT38, TDRG1) instead of clean markers.
top_n_markers <- function(ids, contrast_df1, contrast_df2, n = 12, min_base_mean = 50) {
  bm1 <- contrast_df1$baseMean[match(ids, contrast_df1$ensembl_id)]
  bm2 <- contrast_df2$baseMean[match(ids, contrast_df2$ensembl_id)]
  ids <- ids[bm1 >= min_base_mean & bm2 >= min_base_mean]

  combined_padj <- pmax(contrast_df1$padj[match(ids, contrast_df1$ensembl_id)],
                        contrast_df2$padj[match(ids, contrast_df2$ensembl_id)])
  ids[order(combined_padj)][seq_len(min(n, length(ids)))]
}

top_hot     <- top_n_markers(hot_markers,     de_hot_vs_overlap,  de_hot_vs_cold)
top_cold    <- top_n_markers(cold_markers,    de_cold_vs_overlap, de_hot_vs_cold)
top_overlap <- top_n_markers(overlap_markers, de_hot_vs_overlap,  de_cold_vs_overlap)

marker_ids <- c(top_hot, top_overlap, top_cold)   # ordered Hot -> Overlap -> Cold
marker_gene_labels <- gene_symbols[match(marker_ids, rownames(expr))]

# Build the heatmap matrix: z-scored expression, samples grouped by Region
region_order <- order(scores_xy$Region[keep_region])
sample_order <- scores_xy$sample[keep_region][region_order]

mat <- t(scale(t(expr[marker_ids, sample_order])))
mat[mat >  3] <-  3
mat[mat < -3] <- -3
rownames(mat) <- marker_gene_labels

ann_col <- data.frame(
  Region = scores_xy$Region[keep_region][region_order],
  row.names = sample_order
)
ann_row <- data.frame(
  `Marker of` = rep(c("Hot-only", "Overlap", "Cold-only"),
                     times = c(length(top_hot), length(top_overlap), length(top_cold))),
  row.names = marker_gene_labels,
  check.names = FALSE
)

ann_colors <- list(
  Region = c("Hot-only" = "#D6604D", "Overlap" = "#762A83", "Cold-only" = "#2166AC"),
  `Marker of` = c("Hot-only" = "#D6604D", "Overlap" = "#762A83", "Cold-only" = "#2166AC")
)

for (ext in c("pdf", "png")) {
  pheatmap(mat,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    gaps_row          = cumsum(c(length(top_hot), length(top_overlap)))[1:2],
    show_colnames     = FALSE,
    annotation_col    = ann_col,
    annotation_row    = ann_row,
    annotation_colors = ann_colors,
    color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    border_color      = NA,
    fontsize_row      = 8,
    main = paste0("Marker Genes by Scatter Region -- Hot-only (n=", length(top_hot),
                  "), Overlap (n=", length(top_overlap), "), Cold-only (n=", length(top_cold), ")"),
    filename = paste0("outputs/tcga_bulk_seq_outputs/heatmap_three_regions_markers.", ext),
    width = 12, height = 7
  )
}
cat("Marker-gene heatmap saved: outputs/tcga_bulk_seq_outputs/heatmap_three_regions_markers.pdf/png\n")

# 12. Feature-style scatter — top marker gene per region
## Plotted on the same ImmuneScore/StromalScore coordinates, coloured by
## expression; checks whether high-expression samples fall in their region.
feature_gene_ids   <- c(top_hot[1], top_overlap[1], top_cold[1])
feature_gene_names <- gene_symbols[match(feature_gene_ids, rownames(expr))]
feature_region     <- c("Hot-only", "Overlap", "Cold-only")

feature_df <- lapply(seq_along(feature_gene_ids), function(i) {
  data.frame(
    sample        = scores_xy$sample,
    ImmuneScore   = scores_xy$ImmuneScore,
    StromalScore  = scores_xy$StromalScore,
    expression    = expr[feature_gene_ids[i], ],   # scores_xy is already in expr's column order
    ## Z-score within each gene so all three panels share a comparable colour
    ## range (Raw VST is on a different absolute scale per gene -- without
    ## this, a highly-expressed gene like CD48 saturates the shared legend and
    ## visually flattens genes with a lower expression range).
    z_expression  = as.numeric(scale(expr[feature_gene_ids[i], ])),
    gene          = feature_gene_names[i],
    marker_of     = feature_region[i]
  )
}) %>% bind_rows() %>%
  mutate(panel_label = sprintf("%s  (marker of %s)", gene, marker_of),
         panel_label = factor(panel_label, levels = unique(panel_label)))

p_feature <- ggplot(feature_df, aes(x = ImmuneScore, y = StromalScore, colour = z_expression)) +
  geom_point(size = 1.4, alpha = 0.8) +
  stat_ellipse(data = subset(scores_xy, ImmuneGroup == "Immune-hot"),
               aes(x = ImmuneScore, y = StromalScore), level = 0.90,
               colour = "#D6604D", linewidth = 0.5, inherit.aes = FALSE) +
  stat_ellipse(data = subset(scores_xy, ImmuneGroup == "Immune-cold"),
               aes(x = ImmuneScore, y = StromalScore), level = 0.90,
               colour = "#2166AC", linewidth = 0.5, inherit.aes = FALSE) +
  scale_colour_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", midpoint = 0,
                         limits = c(-3, 3), oob = scales::squish, name = "Expression\n(z-score)") +
  facet_wrap(~ panel_label, nrow = 1) +
  labs(
    title = "Top Marker Gene per Region, Plotted on the Classification Scatter",
    subtitle = "Colour = per-gene z-scored expression -- checks whether each marker's high-expression samples fall in its named region",
    x = "ESTIMATE ImmuneScore", y = "ESTIMATE StromalScore"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 9),
        strip.background = element_rect(fill = "grey90", colour = NA))

ggsave("outputs/tcga_bulk_seq_outputs/scatter_regions_marker_gene_expression.pdf", p_feature, width = 13, height = 5)
ggsave("outputs/tcga_bulk_seq_outputs/scatter_regions_marker_gene_expression.png", p_feature, width = 13, height = 5, dpi = 300)
cat("Marker-gene feature scatter saved: outputs/tcga_bulk_seq_outputs/scatter_regions_marker_gene_expression.pdf/png\n")
