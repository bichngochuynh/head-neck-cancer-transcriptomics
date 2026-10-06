# Immune classification: immune-hot vs immune-cold
##
## Workflow: k-means (k=2) on MCP-counter scores only; ESTIMATE and xCell play
## no role in clustering. Clusters are labelled "hot"/"cold" by mean ESTIMATE
## ImmuneScore (ImmuneScore is thus a labelling variable, not validation).
## StromalScore, ESTIMATEScore, and xCell are independent validation
## (see 05_plot_immune_scores_by_tool.R). Unsupervised, score-driven grouping.

library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(patchwork)

dir.create("outputs/tcga_bulk_seq_outputs", recursive = TRUE, showWarnings = FALSE)

# Load scores
scores <- read.csv("outputs/tcga_bulk_seq_outputs/all_deconvolution_scores.csv",
                   check.names = FALSE)

# STEP 1 — Cluster: k-means (k=2) on MCP-counter scores
clust_cols <- c(
  "MCP.T cells", "MCP.CD8 T cells", "MCP.Cytotoxic lymphocytes",
  "MCP.NK cells", "MCP.B lineage", "MCP.Monocytic lineage",
  "MCP.Myeloid dendritic cells"
)

set.seed(42)
clust_mat  <- scale(scores[, clust_cols])
km         <- kmeans(clust_mat, centers = 2, nstart = 50)
scores$MCP_CompositeScore <- rowMeans(clust_mat)   # kept for reference/ordering only

# STEP 2/3 — Label: name clusters using ESTIMATE ImmuneScore
## ESTIMATE played no part in clustering -- used only to name hot vs cold.
group_means <- tapply(scores$ImmuneScore, km$cluster, mean)
hot_cluster <- as.integer(names(which.max(group_means)))
scores$ImmuneGroup <- ifelse(km$cluster == hot_cluster, "Immune-hot", "Immune-cold")
scores$ImmuneGroup <- factor(scores$ImmuneGroup, levels = c("Immune-cold", "Immune-hot"))

cat("Immune group sizes:\n")
print(table(scores$ImmuneGroup))

group_colors <- c("Immune-cold" = "#2166AC", "Immune-hot" = "#D6604D")

# 1. Heatmap — MCP-counter clustering basis, ordered by ImmuneScore
heatmap_cols <- c(
  ## Clustering features (rows 1-7)
  "MCP.T cells", "MCP.CD8 T cells", "MCP.Cytotoxic lymphocytes",
  "MCP.NK cells", "MCP.B lineage", "MCP.Monocytic lineage",
  "MCP.Myeloid dendritic cells",
  ## Held-out MCP-counter populations, not used in clustering (rows 8-10)
  "MCP.Neutrophils", "MCP.Endothelial cells", "MCP.Fibroblasts"
)

scores_ord <- scores[order(scores$ImmuneScore), ]
mat        <- t(apply(scores_ord[, heatmap_cols], 2, scale))
colnames(mat) <- scores_ord$sample
mat[mat >  3] <-  3
mat[mat < -3] <- -3

rownames(mat) <- c(
  "T cells", "CD8+ T cells", "Cytotoxic lymphocytes",
  "NK cells", "B cells", "Monocytic lineage", "Myeloid DCs",
  "Neutrophils", "Endothelial cells", "Fibroblasts"
)

## Column annotation: immune group
ann_col <- data.frame(
  `Immune group` = scores_ord$ImmuneGroup,
  row.names      = scores_ord$sample,
  check.names    = FALSE
)

## Row annotation: whether the MCP-counter population fed the k-means clustering
ann_row <- data.frame(
  Role = c(rep("Clustering feature", 7), rep("Held-out (MCP-counter)", 3)),
  row.names = rownames(mat)
)

ann_colors <- list(
  `Immune group` = group_colors,
  Role = c("Clustering feature" = "#FF7F00", "Held-out (MCP-counter)" = "#984EA3")
)

for (ext in c("pdf", "png")) {
  pheatmap(mat,
    cluster_rows      = FALSE,          # fixed order: clustering features vs held-out
    cluster_cols      = FALSE,
    gaps_row          = 7,              # gap between clustering features and held-out rows
    show_colnames     = FALSE,
    annotation_col    = ann_col,
    annotation_row    = ann_row,
    annotation_colors = ann_colors,
    color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    border_color      = NA,
    fontsize_row      = 8.5,
    main = "MCP-counter Cell-Type Scores — Clustering Basis — TCGA-HNSC (n = 520, ordered by ESTIMATE ImmuneScore)",
    filename = paste0("outputs/tcga_bulk_seq_outputs/heatmap_immune_classification.", ext),
    width = 14, height = 5.5
  )
}
cat("Heatmap saved\n")

# 2. Scatter — ImmuneScore vs StromalScore, coloured by group
scatter <- ggplot(scores, aes(x = ImmuneScore, y = StromalScore,
                               colour = ImmuneGroup, shape = ImmuneGroup)) +
  geom_point(alpha = 0.55, size = 1.4) +
  stat_ellipse(level = 0.90, linewidth = 0.6) +
  scale_colour_manual(values = group_colors, name = "Immune group") +
  scale_shape_manual(values = c("Immune-cold" = 16, "Immune-hot" = 17),
                     name = "Immune group") +
  geom_vline(xintercept = median(scores$ImmuneScore),
             linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = median(scores$StromalScore),
             linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  labs(
    x     = "ESTIMATE ImmuneScore (labelling variable)",
    y     = "ESTIMATE StromalScore (independent)",
    title = "Immune vs Stromal Score (ESTIMATE)",
    subtitle = "Clusters from MCP-counter; ImmuneScore named them (x, circular); StromalScore validates (y, independent)\n90% confidence ellipses per group"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 9))

ggsave("outputs/tcga_bulk_seq_outputs/immune_classification_scatter.pdf",
       scatter, width = 7.5, height = 5)
ggsave("outputs/tcga_bulk_seq_outputs/immune_classification_scatter.png",
       scatter, width = 7.5, height = 5, dpi = 300)

## Save classification labels for downstream use
saveRDS(scores[, c("sample", "ImmuneGroup", "MCP_CompositeScore")],
        "outputs/tcga_bulk_seq_outputs/immune_group_labels.rds")
write.csv(scores[, c("sample", "ImmuneGroup", "MCP_CompositeScore")],
          "outputs/tcga_bulk_seq_outputs/immune_group_labels.csv", row.names = FALSE)

cat("Scatter plot saved\n")
cat("Group labels saved (immune_group_labels.csv)\n")
cat("\nGroup summary:\n")
print(table(scores$ImmuneGroup))
