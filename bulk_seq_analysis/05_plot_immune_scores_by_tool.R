# Immune infiltration scores by tool — immune-hot vs immune-cold
## Groups: k-means on MCP-counter only, named "hot"/"cold" via ESTIMATE
## ImmuneScore (see 04_plot_immune_classification.R). Three figures:
##   ESTIMATE     - ImmuneScore is the labelling variable (circular); StromalScore
##                  & ESTIMATEScore are independent (partially derived, resp.)
##   MCP-counter  - 7/10 populations are the clustering features (circular,
##                  confirms clustering); 3 held out as independent checks
##   xCell        - fully independent validation

library(ggplot2)
library(patchwork)

dir.create("outputs/tcga_bulk_seq_outputs", recursive = TRUE, showWarnings = FALSE)

# Load scores and immune group labels
scores <- read.csv("outputs/tcga_bulk_seq_outputs/all_deconvolution_scores.csv",
                   check.names = FALSE)
labels <- read.csv("outputs/tcga_bulk_seq_outputs/immune_group_labels.csv")
scores <- merge(scores, labels, by = "sample")
scores$ImmuneGroup <- factor(scores$ImmuneGroup, levels = c("Immune-cold", "Immune-hot"))

group_colors <- c("Immune-cold" = "#2166AC", "Immune-hot" = "#D6604D")
n_cold <- sum(scores$ImmuneGroup == "Immune-cold")
n_hot  <- sum(scores$ImmuneGroup == "Immune-hot")
cat(sprintf("Immune-cold: %d  |  Immune-hot: %d\n", n_cold, n_hot))

# Shared violin+boxplot helper
## Returns list(plot = <ggplot>, stat = <1-row data.frame>) so p-values used
## in the plot annotations can also be exported to a stats table.
make_violin <- function(col, y_label, title, tool, role = "validation (independent)") {
  hot_vals  <- scores[[col]][scores$ImmuneGroup == "Immune-hot"]
  cold_vals <- scores[[col]][scores$ImmuneGroup == "Immune-cold"]
  wt    <- wilcox.test(hot_vals, cold_vals)
  pval  <- wt$p.value
  stars <- ifelse(pval < 0.001, "***",
           ifelse(pval < 0.01,  "**",
           ifelse(pval < 0.05,  "*", "ns")))
  p_below_floor <- pval < 2.2e-16
  ptext <- if (p_below_floor) {
    paste0(stars, "\np < 2.2e-16")
  } else {
    paste0(stars, "\np = ", formatC(pval, format = "e", digits = 1))
  }

  rng    <- range(scores[[col]], na.rm = TRUE)
  y_max  <- rng[2]
  y_step <- diff(rng) * 0.09

  plot <- ggplot(scores, aes(x = ImmuneGroup, y = .data[[col]], fill = ImmuneGroup)) +
    geom_violin(alpha = 0.55, trim = FALSE, linewidth = 0.3) +
    geom_boxplot(width = 0.13, fill = "white",
                 outlier.size = 0.3, linewidth = 0.35) +
    annotate("segment",
             x = 1, xend = 2,
             y = y_max + y_step, yend = y_max + y_step,
             linewidth = 0.4) +
    annotate("text",
             x = 1.5, y = y_max + y_step * 2.2,
             label = ptext, size = 2.6, lineheight = 0.9) +
    scale_fill_manual(values = group_colors) +
    scale_x_discrete(labels = c(paste0("Immune-cold\n(n=", n_cold, ")"),
                                paste0("Immune-hot\n(n=", n_hot, ")"))) +
    labs(x = NULL, y = y_label, title = title) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none",
          plot.title  = element_text(size = 9, face = "bold"),
          axis.text.x = element_text(size = 8),
          axis.title.y = element_text(size = 8))

  stat <- data.frame(
    tool                = tool,
    cell_type           = gsub("\n", " ", title),
    column              = col,
    role                = role,
    n_cold              = length(cold_vals),
    n_hot               = length(hot_vals),
    median_cold         = median(cold_vals, na.rm = TRUE),
    median_hot          = median(hot_vals, na.rm = TRUE),
    W_statistic         = unname(wt$statistic),
    p_value             = pval,
    p_value_display     = if (p_below_floor) "< 2.2e-16" else formatC(pval, format = "e", digits = 2),
    significance        = stars,
    stringsAsFactors = FALSE
  )

  list(plot = plot, stat = stat)
}

subtitle_theme <- theme(
  plot.title    = element_text(size = 13, face = "bold"),
  plot.subtitle = element_text(size = 9.5, colour = "grey40")
)

# FIGURE 1 — ESTIMATE (3 scores), independent validation panel
estimate_vars <- list(
  list(col = "ImmuneScore",   label = "ESTIMATE score", title = "ImmuneScore\n(labelling variable)",
       role = "labelling variable (circular)"),
  list(col = "StromalScore",  label = "ESTIMATE score", title = "StromalScore",
       role = "validation (independent)"),
  list(col = "ESTIMATEScore", label = "ESTIMATE score", title = "ESTIMATEScore\n(Immune + Stromal)",
       role = "validation (partially derived: includes ImmuneScore)")
)

res_est   <- lapply(estimate_vars, function(s) make_violin(s$col, s$label, s$title, tool = "ESTIMATE", role = s$role))
p_est     <- lapply(res_est, `[[`, "plot")
stats_est <- do.call(rbind, lapply(res_est, `[[`, "stat"))
fig_estimate <- wrap_plots(p_est, nrow = 1) +
  plot_annotation(
    title    = "ESTIMATE Scores — ImmuneScore Named the Groups; StromalScore Validates",
    subtitle = paste0("TCGA-HNSC bulk RNA-seq  |  n = 520  |  Clusters from MCP-counter; ImmuneScore assigned hot/cold labels",
                      "  |  Wilcoxon test"),
    theme = subtitle_theme
  )

ggsave("outputs/tcga_bulk_seq_outputs/fig_ESTIMATE_by_group.pdf",
       fig_estimate, width = 9, height = 5)
ggsave("outputs/tcga_bulk_seq_outputs/fig_ESTIMATE_by_group.png",
       fig_estimate, width = 9, height = 5, dpi = 300)
cat("ESTIMATE figure saved\n")

# FIGURE 2 — MCP-counter (10 cell types, 2×5 grid)
## First 7 are the clustering features (circular, confirms separation); last 3
## (Neutrophils, Endothelial cells, Fibroblasts) were held out as independent checks.
mcp_vars <- list(
  list(col = "MCP.T cells",              label = "MCP-counter score", title = "T cells",                 role = "clustering (circular)"),
  list(col = "MCP.CD8 T cells",          label = "MCP-counter score", title = "CD8+ T cells",             role = "clustering (circular)"),
  list(col = "MCP.Cytotoxic lymphocytes",label = "MCP-counter score", title = "Cytotoxic\nlymphocytes",   role = "clustering (circular)"),
  list(col = "MCP.B lineage",            label = "MCP-counter score", title = "B cells",                  role = "clustering (circular)"),
  list(col = "MCP.NK cells",             label = "MCP-counter score", title = "NK cells",                 role = "clustering (circular)"),
  list(col = "MCP.Monocytic lineage",    label = "MCP-counter score", title = "Monocytic\nlineage",       role = "clustering (circular)"),
  list(col = "MCP.Myeloid dendritic cells", label = "MCP-counter score", title = "Myeloid DCs",           role = "clustering (circular)"),
  list(col = "MCP.Neutrophils",          label = "MCP-counter score", title = "Neutrophils\n(held out)",  role = "validation (independent)"),
  list(col = "MCP.Endothelial cells",    label = "MCP-counter score", title = "Endothelial cells\n(held out)", role = "validation (independent)"),
  list(col = "MCP.Fibroblasts",          label = "MCP-counter score", title = "Fibroblasts\n(held out)",  role = "validation (independent)")
)

res_mcp   <- lapply(mcp_vars, function(s) make_violin(s$col, s$label, s$title, tool = "MCP-counter", role = s$role))
p_mcp     <- lapply(res_mcp, `[[`, "plot")
stats_mcp <- do.call(rbind, lapply(res_mcp, `[[`, "stat"))
fig_mcp <- wrap_plots(p_mcp, nrow = 2, ncol = 5) +
  plot_annotation(
    title    = "MCP-counter: Clustering Basis (7) + Held-out Populations (3)",
    subtitle = paste0("TCGA-HNSC  |  n = 520  |  Immune-hot (n=", n_hot, ") vs Immune-cold (n=", n_cold, ")",
                      "  |  first 7 panels are the k-means features — expected to separate  |  Wilcoxon test"),
    theme = subtitle_theme
  )

ggsave("outputs/tcga_bulk_seq_outputs/fig_MCPcounter_by_group.pdf",
       fig_mcp, width = 16, height = 8)
ggsave("outputs/tcga_bulk_seq_outputs/fig_MCPcounter_by_group.png",
       fig_mcp, width = 16, height = 8, dpi = 300)
cat("MCP-counter figure saved\n")

# FIGURE 3 — xCell validation (12 selected cell types, 3×4 grid)
## Selected to mirror MCP-counter cell types + unique xCell insights
xcell_vars <- list(
  ## Lymphocytes
  list(col = "xCell.CD8+ T-cells",   label = "xCell score", title = "CD8+ T cells"),
  list(col = "xCell.CD4+ T-cells",   label = "xCell score", title = "CD4+ T cells"),
  list(col = "xCell.Tregs",          label = "xCell score", title = "Tregs"),
  list(col = "xCell.NK cells",       label = "xCell score", title = "NK cells"),
  ## B cells
  list(col = "xCell.Memory B-cells", label = "xCell score", title = "Memory B cells"),
  list(col = "xCell.Plasma cells",   label = "xCell score", title = "Plasma cells"),
  ## Myeloid
  list(col = "xCell.Macrophages M1", label = "xCell score", title = "Macrophages M1"),
  list(col = "xCell.Macrophages M2", label = "xCell score", title = "Macrophages M2"),
  list(col = "xCell.Monocytes",      label = "xCell score", title = "Monocytes"),
  list(col = "xCell.DC",             label = "xCell score", title = "Dendritic cells"),
  ## Stromal
  list(col = "xCell.Fibroblasts",    label = "xCell score", title = "Fibroblasts"),
  list(col = "xCell.Endothelial cells", label = "xCell score", title = "Endothelial\ncells")
)

res_xc      <- lapply(xcell_vars, function(s) make_violin(s$col, s$label, s$title, tool = "xCell"))
p_xc        <- lapply(res_xc, `[[`, "plot")
stats_xcell <- do.call(rbind, lapply(res_xc, `[[`, "stat"))
fig_xcell <- wrap_plots(p_xc, nrow = 3, ncol = 4) +
  plot_annotation(
    title    = "xCell Enrichment Scores — Independent Validation Panel",
    subtitle = paste0("TCGA-HNSC  |  n = 520  |  Groups defined by MCP-counter clustering only  |",
                      "  xCell not used for clustering  |  Wilcoxon test"),
    theme = subtitle_theme
  )

ggsave("outputs/tcga_bulk_seq_outputs/fig_xCell_by_group.pdf",
       fig_xcell, width = 14, height = 11)
ggsave("outputs/tcga_bulk_seq_outputs/fig_xCell_by_group.png",
       fig_xcell, width = 14, height = 11, dpi = 300)
cat("xCell figure saved\n")

# STATS TABLE — Wilcoxon tests underlying the plot annotations
all_stats <- rbind(stats_est, stats_mcp, stats_xcell)
write.csv(all_stats, "outputs/tcga_bulk_seq_outputs/immune_group_stats.csv", row.names = FALSE)
cat("Stats table saved (immune_group_stats.csv):", nrow(all_stats), "rows\n")

cat("\nAll figures saved to outputs/tcga_bulk_seq_outputs/\n")
