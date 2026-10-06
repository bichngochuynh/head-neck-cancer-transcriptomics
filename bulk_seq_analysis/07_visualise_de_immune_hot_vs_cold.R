# Visualisations for the Immune-hot vs Immune-cold DESeq2 results
## 1. Volcano (labelled by padj, not fold-change) 2. DE summary table
## 3. Heat map (ComplexHeatmap) of top 10 up/down DEGs, row-scaled VST expression

library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)

dir.create("outputs/tcga_bulk_seq_outputs", recursive = TRUE, showWarnings = FALSE)

res_df <- read.csv("outputs/tcga_bulk_seq_outputs/DE_ImmuneHot_vs_ImmuneCold.csv")

## Loaded once here and reused below (subtitle sample counts + heat map
## expression matrix in section 3) so the group sizes shown in every figure
## always reflect the samples actually used in the DESeq2 fit.
dds <- readRDS("outputs/tcga_bulk_seq_outputs/de_immune_hot_vs_cold_DESeq2_checkpoint.rds")
sample_counts <- table(colData(dds)$ImmuneGroup)
n_total <- ncol(dds)
n_hot   <- unname(sample_counts["Immune-hot"])
n_cold  <- unname(sample_counts["Immune-cold"])

group_colors <- c("Higher in Immune-hot"  = "#D6604D",
                   "Higher in Immune-cold" = "#2166AC",
                   "NS" = "grey70")

n_hot_sig  <- sum(res_df$sig == "Higher in Immune-hot")
n_cold_sig <- sum(res_df$sig == "Higher in Immune-cold")

## Single ranking criterion (padj, smallest first) and label count shared by
## the volcano plot and the summary table below, so both show the same genes.
N_LABEL <- 15

# 1. Volcano — top genes labelled by significance (padj), not fold-change
top_labels <- res_df %>%
  filter(sig != "NS") %>%
  group_by(sig) %>%
  slice_min(order_by = padj, n = N_LABEL) %>%
  ungroup()

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), colour = sig)) +
  geom_point(alpha = 0.55, size = 1.1) +
  geom_point(data = top_labels, size = 2, shape = 21, colour = "black",
             aes(fill = sig), stroke = 0.4, show.legend = FALSE) +
  geom_text_repel(data = top_labels, aes(label = gene), size = 3.2, fontface = "bold",
                   max.overlaps = Inf, min.segment.length = 0, seed = 42,
                   segment.size = 0.3, box.padding = 0.4) +
  scale_colour_manual(values = group_colors, name = NULL) +
  scale_fill_manual(values = group_colors, guide = "none") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  labs(
    title = "DESeq2: Immune-hot vs Immune-cold (TCGA-HNSC bulk RNA-seq)",
    subtitle = paste0("n = ", n_total, "  |  Immune-hot (n=", n_hot, ") vs Immune-cold (n=", n_cold, ")  |  ",
                      n_hot_sig, " genes up in hot, ", n_cold_sig, " up in cold (padj<0.05, |log2FC|>0.5)",
                      "\nLabelled genes = top ", N_LABEL, " per direction by adjusted p-value"),
    x = "log2 Fold Change (Immune-hot / Immune-cold)",
    y = "-log10 adjusted p-value"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 9.5))

ggsave("outputs/tcga_bulk_seq_outputs/volcano_ImmuneHot_vs_ImmuneCold.pdf", p_volcano, width = 9.5, height = 8)
ggsave("outputs/tcga_bulk_seq_outputs/volcano_ImmuneHot_vs_ImmuneCold.png", p_volcano, width = 9.5, height = 8, dpi = 300)
cat("Volcano plot updated (labelled by significance, not fold-change)\n")

# 2. DE summary table — top N_LABEL significant genes per direction
## Same padj ranking/N as the volcano labels, so both show identical genes.
make_table_df <- function(direction, n = N_LABEL) {
  res_df %>%
    filter(sig == direction) %>%
    arrange(padj) %>%
    slice_head(n = n) %>%
    transmute(
      Gene       = gene,
      `log2FC`   = sprintf("%.2f", log2FoldChange),
      `Adj. p`   = ifelse(padj < 1e-300, "< 1e-300", formatC(padj, format = "e", digits = 2)),
      Direction  = ifelse(direction == "Higher in Immune-hot", "Hot", "Cold")
    )
}

table_df <- bind_rows(make_table_df("Higher in Immune-hot"),
                       make_table_df("Higher in Immune-cold"))

## Build a two-column layout (Hot | Cold) side by side, 15 rows each
hot_tbl  <- table_df %>% filter(Direction == "Hot")  %>% select(Gene, log2FC, `Adj. p`)
cold_tbl <- table_df %>% filter(Direction == "Cold") %>% select(Gene, log2FC, `Adj. p`)

base_theme <- ttheme_minimal(
  core = list(
    fg_params = list(hjust = 0, x = 0.05, fontsize = 10),
    bg_params = list(fill = c("white", "grey96"), col = NA)
  ),
  colhead = list(
    fg_params = list(fontface = "bold", fontsize = 10, col = "white"),
    bg_params = list(fill = "grey30", col = NA)
  )
)

hot_grob  <- tableGrob(hot_tbl,  rows = NULL, theme = base_theme)
cold_grob <- tableGrob(cold_tbl, rows = NULL, theme = base_theme)

## Colour the two table headers to match the volcano plot's group colours
hot_grob$grobs[hot_grob$layout$name == "colhead-bg"][[1]]$gp$fill   <- group_colors["Higher in Immune-hot"]
cold_grob$grobs[cold_grob$layout$name == "colhead-bg"][[1]]$gp$fill <- group_colors["Higher in Immune-cold"]

hot_title  <- textGrob(paste0("Higher in Immune-hot (top ", N_LABEL, " by padj, n sig = ", n_hot_sig, ")"),
                        gp = gpar(fontface = "bold", fontsize = 12, col = group_colors["Higher in Immune-hot"]))
cold_title <- textGrob(paste0("Higher in Immune-cold (top ", N_LABEL, " by padj, n sig = ", n_cold_sig, ")"),
                        gp = gpar(fontface = "bold", fontsize = 12, col = group_colors["Higher in Immune-cold"]))

main_title <- textGrob("Differential Expression: Immune-hot vs Immune-cold (TCGA-HNSC, DESeq2)",
                        gp = gpar(fontface = "bold", fontsize = 14), just = "left", x = 0.01)
subtitle   <- textGrob(paste0("Positive log2FC = higher in Immune-hot  |  n = ", n_total,
                              " (", n_hot, " hot / ", n_cold, " cold)"),
                        gp = gpar(fontsize = 10, col = "grey40"), just = "left", x = 0.01)

## Explicit per-row heights (lines for text rows, null/expanding for the tables)
## so title text never overlaps the table grobs below it
row_grobs   <- list(main_title, subtitle, hot_title, hot_grob, cold_title, cold_grob)
row_heights <- grid::unit.c(unit(2, "lines"), unit(1.5, "lines"),
                            unit(1.8, "lines"), unit(1, "null"),
                            unit(1.8, "lines"), unit(1, "null"))

pdf("outputs/tcga_bulk_seq_outputs/DE_table_ImmuneHot_vs_ImmuneCold.pdf", width = 10, height = 11)
grid.arrange(grobs = row_grobs, ncol = 1, heights = row_heights)
dev.off()

png("outputs/tcga_bulk_seq_outputs/DE_table_ImmuneHot_vs_ImmuneCold.png", width = 10, height = 11, units = "in", res = 300)
grid.arrange(grobs = row_grobs, ncol = 1, heights = row_heights)
dev.off()

cat("DE table figure saved: outputs/tcga_bulk_seq_outputs/DE_table_ImmuneHot_vs_ImmuneCold.pdf/png\n")

# 3. Heat map — top 10 up- and downregulated DEGs
## VST-normalised expression, row-scaled (z-score) for cross-sample comparability.
N_HEATMAP <- 10

top_up   <- res_df %>% filter(sig == "Higher in Immune-hot")  %>% arrange(padj) %>% slice_head(n = N_HEATMAP)
top_down <- res_df %>% filter(sig == "Higher in Immune-cold") %>% arrange(padj) %>% slice_head(n = N_HEATMAP)
heatmap_genes <- bind_rows(top_up, top_down)

## dds loaded once near the top of the script and reused here
vsd <- vst(dds, blind = FALSE)
expr_mat <- assay(vsd)[heatmap_genes$ensembl_id, ]
rownames(expr_mat) <- heatmap_genes$gene

## Row-scale (z-score) so each gene's own dynamic range drives its colour
expr_z <- t(scale(t(expr_mat)))

col_data <- as.data.frame(colData(dds))[colnames(expr_z), , drop = FALSE]
sample_order <- order(col_data$ImmuneGroup)
expr_z    <- expr_z[, sample_order]
col_data  <- col_data[sample_order, , drop = FALSE]

top_ann <- HeatmapAnnotation(
  ImmuneGroup = col_data$ImmuneGroup,
  col = list(ImmuneGroup = c("Immune-hot" = "#D6604D", "Immune-cold" = "#2166AC")),
  annotation_name_side = "left"
)

row_direction <- c(rep("Higher in Immune-hot", nrow(top_up)), rep("Higher in Immune-cold", nrow(top_down)))
row_ann <- rowAnnotation(
  Direction = row_direction,
  col = list(Direction = c("Higher in Immune-hot" = "#D6604D", "Higher in Immune-cold" = "#2166AC")),
  show_annotation_name = FALSE
)

col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#D6604D"))

ht <- Heatmap(
  expr_z,
  name = "z-score",
  col = col_fun,
  top_annotation = top_ann,
  left_annotation = row_ann,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  row_split = factor(row_direction, levels = c("Higher in Immune-hot", "Higher in Immune-cold")),
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  column_title = paste0("Top ", N_HEATMAP, " up- and downregulated DEGs: Immune-hot vs Immune-cold (TCGA-HNSC)"),
  column_title_gp = gpar(fontface = "bold", fontsize = 12),
  heatmap_legend_param = list(title = "Row z-score (VST)")
)

pdf("outputs/tcga_bulk_seq_outputs/heatmap_ImmuneHot_vs_ImmuneCold_top10.pdf", width = 10, height = 7)
draw(ht)
dev.off()

png("outputs/tcga_bulk_seq_outputs/heatmap_ImmuneHot_vs_ImmuneCold_top10.png", width = 10, height = 7, units = "in", res = 300)
draw(ht)
dev.off()

cat("Heat map saved: outputs/tcga_bulk_seq_outputs/heatmap_ImmuneHot_vs_ImmuneCold_top10.pdf/png\n")
