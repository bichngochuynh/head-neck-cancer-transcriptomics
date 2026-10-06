# Clinical association: Immune-hot vs Immune-cold (TCGA-HNSC)
## 1. HPV enrichment (Fisher's exact) 2. Survival (KM + log-rank + Cox), OS + PFS endpoints
## Source: GDC_HNSC/cBioPortal/hnsc_tcga_pan_can_atlas_2018_clinical_data.tsv
## (standard GDC export has no HPV column; this cBioPortal file's "Subtype"
## field does, and also carries Overall Survival time/status)

library(dplyr)
library(ggplot2)
library(survival)

dir.create("outputs/tcga_bulk_seq_outputs", recursive = TRUE, showWarnings = FALSE)

# Load clinical data and derive HPV status
clinical <- read.delim(
  "GDC_HNSC/cBioPortal/hnsc_tcga_pan_can_atlas_2018_clinical_data.tsv",
  check.names = FALSE
)
clinical$HPVstatus <- case_when(
  clinical$Subtype == "HNSC_HPV+" ~ "HPV+",
  clinical$Subtype == "HNSC_HPV-" ~ "HPV-",
  TRUE ~ NA_character_
)
cat("HPV status distribution (all patients with a Subtype call):\n")
print(table(clinical$HPVstatus, useNA = "always"))

# Load ImmuneGroup labels and merge on 12-char patient barcode
labels <- read.csv("outputs/tcga_bulk_seq_outputs/immune_group_labels.csv")
clinical_merged <- merge(labels, clinical, by.x = "sample", by.y = "Patient ID")
cat("\nSamples with both ImmuneGroup and clinical data:", nrow(clinical_merged),
    "of", nrow(labels), "\n")

group_colors <- c("Immune-cold" = "#2166AC", "Immune-hot" = "#D6604D")
hpv_colors   <- c("HPV+" = "#F1A340", "HPV-" = "#998EC3")

# 1. HPV enrichment in Immune-hot vs Immune-cold
hpv_data <- clinical_merged %>% filter(!is.na(HPVstatus))
cat("\nSamples with a definitive HPV call:", nrow(hpv_data), "\n")

hpv_table <- table(hpv_data$ImmuneGroup, hpv_data$HPVstatus)
cat("\nContingency table (ImmuneGroup x HPVstatus):\n")
print(hpv_table)
cat("\nRow percentages:\n")
print(round(100 * prop.table(hpv_table, margin = 1), 1))

## Fisher's exact is the preferred test for a 2x2 table: exact (not an
## asymptotic approximation like chi-square) at any sample size, and it
## directly yields an odds ratio + CI rather than just a test statistic.
fisher_res <- fisher.test(hpv_table)

cat("\n=== Fisher's exact test ===\n")
cat("Odds ratio:", round(fisher_res$estimate, 3),
    " | 95% CI:", round(fisher_res$conf.int[1], 3), "-", round(fisher_res$conf.int[2], 3),
    " | p =", format.pval(fisher_res$p.value, digits = 3), "\n")

write.csv(as.data.frame.matrix(hpv_table),
          "outputs/tcga_bulk_seq_outputs/HPV_vs_ImmuneGroup_contingency_table.csv")

hpv_test_summary <- data.frame(
  test       = "Fisher's exact",
  odds_ratio = unname(fisher_res$estimate),
  ci_lower   = fisher_res$conf.int[1],
  ci_upper   = fisher_res$conf.int[2],
  p_value    = fisher_res$p.value
)
write.csv(hpv_test_summary, "outputs/tcga_bulk_seq_outputs/HPV_ImmuneGroup_test_summary.csv", row.names = FALSE)
cat("\nHPV association results saved to outputs/tcga_bulk_seq_outputs/\n")

# Stacked proportion bar chart
p_hpv <- ggplot(hpv_data, aes(x = ImmuneGroup, fill = HPVstatus)) +
  geom_bar(position = "fill", width = 0.6) +
  geom_text(stat = "count", aes(label = after_stat(count)),
            position = position_fill(vjust = 0.5), colour = "white", fontface = "bold") +
  scale_fill_manual(values = hpv_colors, name = "HPV status") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "HPV Status by Immune Group",
    subtitle = paste0("Fisher's exact p = ", format.pval(fisher_res$p.value, digits = 3),
                      "  |  OR = ", round(fisher_res$estimate, 2),
                      "  |  n = ", nrow(hpv_data)),
    x = NULL, y = "Proportion of samples"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 9.5))

ggsave("outputs/tcga_bulk_seq_outputs/HPV_by_ImmuneGroup_barplot.pdf", p_hpv, width = 6, height = 5.5)
ggsave("outputs/tcga_bulk_seq_outputs/HPV_by_ImmuneGroup_barplot.png", p_hpv, width = 6, height = 5.5, dpi = 300)
cat("HPV bar plot saved: outputs/tcga_bulk_seq_outputs/HPV_by_ImmuneGroup_barplot.pdf/png\n")

km_xmax <- 120  # cap follow-up at 10 years for all KM plots

## Reusable KM plotting helper (manual ggplot step function; survminer not installed)
km_survival_plot <- function(fit, strata_prefix, colour_values, legend_name,
                              plot_title, subtitle, file_stem, xmax = km_xmax,
                              legend_nrow = 1, width = 7.5, height = 6,
                              y_label = "Overall survival probability") {
  km_df <- data.frame(
    time     = fit$time,
    surv     = fit$surv,
    n_censor = fit$n.censor,
    strata   = rep(names(fit$strata), fit$strata)
  )
  km_df$strata <- gsub(strata_prefix, "", km_df$strata)

  start_rows <- data.frame(time = 0, surv = 1, n_censor = 0, strata = unique(km_df$strata))
  km_df <- rbind(start_rows, km_df)
  km_df <- km_df[order(km_df$strata, km_df$time), ]
  km_df <- km_df[km_df$time <= xmax, ]

  p <- ggplot(km_df, aes(x = time, y = surv, colour = strata)) +
    geom_step(linewidth = 0.7) +
    geom_point(data = subset(km_df, n_censor > 0), shape = 3, size = 2) +
    scale_colour_manual(values = colour_values, name = legend_name) +
    scale_x_continuous(limits = c(0, xmax), breaks = seq(0, xmax, by = 12)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(title = plot_title, subtitle = subtitle,
         x = "Time (months)", y = y_label) +
    guides(colour = guide_legend(nrow = legend_nrow, byrow = TRUE)) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom",
          legend.text   = element_text(size = 9),
          plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey40", size = 9))

  ggsave(paste0("outputs/tcga_bulk_seq_outputs/", file_stem, ".pdf"), p, width = width, height = height)
  ggsave(paste0("outputs/tcga_bulk_seq_outputs/", file_stem, ".png"), p, width = width, height = height, dpi = 300)
  cat("Kaplan-Meier plot saved: outputs/tcga_bulk_seq_outputs/", file_stem, ".pdf/png\n", sep = "")
  p
}

# 2. Survival analysis by Immune Group
surv_data <- clinical_merged %>%
  filter(!is.na(`Overall Survival (Months)`), !is.na(`Overall Survival Status`)) %>%
  transmute(
    sample      = sample,
    ImmuneGroup = ImmuneGroup,
    time_months = `Overall Survival (Months)`,
    event       = as.numeric(substr(`Overall Survival Status`, 1, 1))  # "0:LIVING"/"1:DECEASED" -> 0/1
  )
cat("\nSamples with usable survival data:", nrow(surv_data), "\n")
cat("Events (deaths):", sum(surv_data$event), "/ censored:", sum(surv_data$event == 0), "\n")

surv_obj <- Surv(surv_data$time_months, surv_data$event)
fit      <- survfit(surv_obj ~ ImmuneGroup, data = surv_data)
logrank  <- survdiff(surv_obj ~ ImmuneGroup, data = surv_data)
logrank_p <- pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)

cox_fit <- coxph(surv_obj ~ ImmuneGroup, data = surv_data)
cox_summary <- summary(cox_fit)

cat("\n=== Log-rank test (Immune-hot vs Immune-cold survival) ===\n")
cat("Chi-square =", round(logrank$chisq, 3), ", p =", format.pval(logrank_p, digits = 3), "\n")

cat("\n=== Cox proportional hazards (Immune-hot vs Immune-cold) ===\n")
cat("Hazard ratio:", round(cox_summary$coefficients[1, "exp(coef)"], 3),
    " | 95% CI:", round(cox_summary$conf.int[1, "lower .95"], 3), "-",
    round(cox_summary$conf.int[1, "upper .95"], 3),
    " | p =", format.pval(cox_summary$coefficients[1, "Pr(>|z|)"], digits = 3), "\n")

survival_test_summary <- data.frame(
  test = c("Log-rank", "Cox PH (Immune-hot vs Immune-cold)"),
  statistic = c(logrank$chisq, cox_summary$coefficients[1, "z"]),
  hazard_ratio = c(NA, cox_summary$coefficients[1, "exp(coef)"]),
  ci_lower = c(NA, cox_summary$conf.int[1, "lower .95"]),
  ci_upper = c(NA, cox_summary$conf.int[1, "upper .95"]),
  p_value = c(logrank_p, cox_summary$coefficients[1, "Pr(>|z|)"])
)
write.csv(survival_test_summary, "outputs/tcga_bulk_seq_outputs/Survival_ImmuneGroup_test_summary.csv", row.names = FALSE)
cat("\nSurvival test results saved: outputs/tcga_bulk_seq_outputs/Survival_ImmuneGroup_test_summary.csv\n")

km_survival_plot(
  fit, strata_prefix = "ImmuneGroup=", colour_values = group_colors, legend_name = "Immune group",
  plot_title = "Overall Survival by Immune Group",
  subtitle = paste0("Log-rank p = ", format.pval(logrank_p, digits = 3),
                    "  |  Cox HR (hot vs cold) = ", round(cox_summary$coefficients[1, "exp(coef)"], 2),
                    " (95% CI ", round(cox_summary$conf.int[1, "lower .95"], 2), "-",
                    round(cox_summary$conf.int[1, "upper .95"], 2), ")",
                    "  |  n = ", nrow(surv_data)),
  file_stem = "KM_survival_ImmuneGroup"
)

# 3. Survival analysis by HPV status
surv_data_hpv <- clinical_merged %>%
  filter(!is.na(`Overall Survival (Months)`), !is.na(`Overall Survival Status`), !is.na(HPVstatus)) %>%
  transmute(
    sample      = sample,
    HPVstatus   = HPVstatus,
    time_months = `Overall Survival (Months)`,
    event       = as.numeric(substr(`Overall Survival Status`, 1, 1))
  )
cat("\nSamples with usable survival + HPV data:", nrow(surv_data_hpv), "\n")
cat("Events (deaths):", sum(surv_data_hpv$event), "/ censored:", sum(surv_data_hpv$event == 0), "\n")

surv_obj_hpv  <- Surv(surv_data_hpv$time_months, surv_data_hpv$event)
fit_hpv       <- survfit(surv_obj_hpv ~ HPVstatus, data = surv_data_hpv)
logrank_hpv   <- survdiff(surv_obj_hpv ~ HPVstatus, data = surv_data_hpv)
logrank_p_hpv <- pchisq(logrank_hpv$chisq, df = length(logrank_hpv$n) - 1, lower.tail = FALSE)

cox_fit_hpv     <- coxph(surv_obj_hpv ~ HPVstatus, data = surv_data_hpv)
cox_summary_hpv <- summary(cox_fit_hpv)
hpv_levels      <- sort(unique(surv_data_hpv$HPVstatus))  # [1] = Cox reference level
hpv_cox_label   <- paste0("Cox HR (", hpv_levels[2], " vs ", hpv_levels[1], ")")

cat("\n=== Log-rank test (HPV+ vs HPV- survival) ===\n")
cat("Chi-square =", round(logrank_hpv$chisq, 3), ", p =", format.pval(logrank_p_hpv, digits = 3), "\n")

cat("\n===", hpv_cox_label, "===\n")
cat("Hazard ratio:", round(cox_summary_hpv$coefficients[1, "exp(coef)"], 3),
    " | 95% CI:", round(cox_summary_hpv$conf.int[1, "lower .95"], 3), "-",
    round(cox_summary_hpv$conf.int[1, "upper .95"], 3),
    " | p =", format.pval(cox_summary_hpv$coefficients[1, "Pr(>|z|)"], digits = 3), "\n")

survival_test_summary_hpv <- data.frame(
  test = c("Log-rank", paste0("Cox PH (", hpv_cox_label, ")")),
  statistic = c(logrank_hpv$chisq, cox_summary_hpv$coefficients[1, "z"]),
  hazard_ratio = c(NA, cox_summary_hpv$coefficients[1, "exp(coef)"]),
  ci_lower = c(NA, cox_summary_hpv$conf.int[1, "lower .95"]),
  ci_upper = c(NA, cox_summary_hpv$conf.int[1, "upper .95"]),
  p_value = c(logrank_p_hpv, cox_summary_hpv$coefficients[1, "Pr(>|z|)"])
)
write.csv(survival_test_summary_hpv, "outputs/tcga_bulk_seq_outputs/Survival_HPVstatus_test_summary.csv", row.names = FALSE)
cat("\nSurvival test results saved: outputs/tcga_bulk_seq_outputs/Survival_HPVstatus_test_summary.csv\n")

km_survival_plot(
  fit_hpv, strata_prefix = "HPVstatus=", colour_values = hpv_colors, legend_name = "HPV status",
  plot_title = "Overall Survival by HPV Status",
  subtitle = paste0("Log-rank p = ", format.pval(logrank_p_hpv, digits = 3),
                    "  |  ", hpv_cox_label, " = ", round(cox_summary_hpv$coefficients[1, "exp(coef)"], 2),
                    " (95% CI ", round(cox_summary_hpv$conf.int[1, "lower .95"], 2), "-",
                    round(cox_summary_hpv$conf.int[1, "upper .95"], 2), ")",
                    "  |  n = ", nrow(surv_data_hpv)),
  file_stem = "KM_survival_HPVstatus"
)

# 4. Survival analysis by Immune Group x HPV status (4-group)
surv_data_combo <- clinical_merged %>%
  filter(!is.na(`Overall Survival (Months)`), !is.na(`Overall Survival Status`), !is.na(HPVstatus)) %>%
  transmute(
    sample      = sample,
    ImmuneGroup = ImmuneGroup,
    HPVstatus   = HPVstatus,
    ComboGroup  = paste(ImmuneGroup, HPVstatus, sep = " / "),
    time_months = `Overall Survival (Months)`,
    event       = as.numeric(substr(`Overall Survival Status`, 1, 1))
  )
cat("\nSamples with usable survival + HPV + ImmuneGroup data:", nrow(surv_data_combo), "\n")
print(table(surv_data_combo$ComboGroup))

surv_obj_combo  <- Surv(surv_data_combo$time_months, surv_data_combo$event)
fit_combo       <- survfit(surv_obj_combo ~ ComboGroup, data = surv_data_combo)
logrank_combo   <- survdiff(surv_obj_combo ~ ComboGroup, data = surv_data_combo)
logrank_p_combo <- pchisq(logrank_combo$chisq, df = length(logrank_combo$n) - 1, lower.tail = FALSE)

cat("\n=== Log-rank test (4-group: ImmuneGroup x HPV status) ===\n")
cat("Chi-square =", round(logrank_combo$chisq, 3), ", df =", length(logrank_combo$n) - 1,
    ", p =", format.pval(logrank_p_combo, digits = 3), "\n")

write.csv(
  data.frame(test = "Log-rank", statistic = logrank_combo$chisq,
             df = length(logrank_combo$n) - 1, p_value = logrank_p_combo),
  "outputs/tcga_bulk_seq_outputs/Survival_ImmuneGroup_x_HPVstatus_test_summary.csv",
  row.names = FALSE
)
cat("\nSurvival test results saved: outputs/tcga_bulk_seq_outputs/Survival_ImmuneGroup_x_HPVstatus_test_summary.csv\n")

combo_colors <- c(
  "Immune-cold / HPV-" = "#2166AC",
  "Immune-cold / HPV+" = "#92C5DE",
  "Immune-hot / HPV-"  = "#B2182B",
  "Immune-hot / HPV+"  = "#F4A582"
)

km_survival_plot(
  fit_combo, strata_prefix = "ComboGroup=", colour_values = combo_colors,
  legend_name = "Immune group / HPV status",
  plot_title = "Overall Survival by Immune Group and HPV Status",
  subtitle = paste0("Log-rank p = ", format.pval(logrank_p_combo, digits = 3),
                    "  |  n = ", nrow(surv_data_combo)),
  file_stem = "KM_survival_ImmuneGroup_x_HPVstatus",
  legend_nrow = 2, width = 7.5, height = 6.5
)

## Multivariable Cox: does ImmuneGroup remain prognostic after adjusting for
## HPV status? Needed because HPV status is itself associated with ImmuneGroup
## (Fisher's test above) -- a univariate ImmuneGroup HR alone can't tell
## you whether that's an independent immune-phenotype effect or just riding
## on the HPV correlation.
cox_fit_multi_os <- coxph(surv_obj_combo ~ ImmuneGroup + HPVstatus, data = surv_data_combo)
cox_summary_multi_os <- summary(cox_fit_multi_os)

cat("\n=== Multivariable Cox (OS ~ ImmuneGroup + HPVstatus) ===\n")
print(cox_summary_multi_os$coefficients)

multivariable_cox_os <- data.frame(
  term         = rownames(cox_summary_multi_os$coefficients),
  hazard_ratio = cox_summary_multi_os$coefficients[, "exp(coef)"],
  ci_lower     = cox_summary_multi_os$conf.int[, "lower .95"],
  ci_upper     = cox_summary_multi_os$conf.int[, "upper .95"],
  p_value      = cox_summary_multi_os$coefficients[, "Pr(>|z|)"]
)
write.csv(multivariable_cox_os, "outputs/tcga_bulk_seq_outputs/Cox_multivariable_OS_ImmuneGroup_HPVstatus.csv", row.names = FALSE)
cat("\nMultivariable Cox results saved: outputs/tcga_bulk_seq_outputs/Cox_multivariable_OS_ImmuneGroup_HPVstatus.csv\n")

# =============================================================================
# PROGRESSION-FREE SURVIVAL (PFS) -- second endpoint, mirrors blocks 2-4 above
# =============================================================================
## Same cBioPortal file also carries PFS (event = progression, not death) --
## the standard secondary endpoint alongside OS for tumour-microenvironment /
## immune-classification survival analyses. "0:CENSORED"/"1:PROGRESSION" -> 0/1,
## same substr() parse as OS.

# 5. PFS by Immune Group
surv_data_pfs <- clinical_merged %>%
  filter(!is.na(`Progress Free Survival (Months)`), !is.na(`Progression Free Status`)) %>%
  transmute(
    sample      = sample,
    ImmuneGroup = ImmuneGroup,
    time_months = `Progress Free Survival (Months)`,
    event       = as.numeric(substr(`Progression Free Status`, 1, 1))
  )
cat("\nSamples with usable PFS data:", nrow(surv_data_pfs), "\n")
cat("Events (progression):", sum(surv_data_pfs$event), "/ censored:", sum(surv_data_pfs$event == 0), "\n")

surv_obj_pfs  <- Surv(surv_data_pfs$time_months, surv_data_pfs$event)
fit_pfs       <- survfit(surv_obj_pfs ~ ImmuneGroup, data = surv_data_pfs)
logrank_pfs   <- survdiff(surv_obj_pfs ~ ImmuneGroup, data = surv_data_pfs)
logrank_p_pfs <- pchisq(logrank_pfs$chisq, df = length(logrank_pfs$n) - 1, lower.tail = FALSE)

cox_fit_pfs     <- coxph(surv_obj_pfs ~ ImmuneGroup, data = surv_data_pfs)
cox_summary_pfs <- summary(cox_fit_pfs)

cat("\n=== Log-rank test (PFS, Immune-hot vs Immune-cold) ===\n")
cat("Chi-square =", round(logrank_pfs$chisq, 3), ", p =", format.pval(logrank_p_pfs, digits = 3), "\n")

cat("\n=== Cox proportional hazards (PFS, Immune-hot vs Immune-cold) ===\n")
cat("Hazard ratio:", round(cox_summary_pfs$coefficients[1, "exp(coef)"], 3),
    " | 95% CI:", round(cox_summary_pfs$conf.int[1, "lower .95"], 3), "-",
    round(cox_summary_pfs$conf.int[1, "upper .95"], 3),
    " | p =", format.pval(cox_summary_pfs$coefficients[1, "Pr(>|z|)"], digits = 3), "\n")

survival_test_summary_pfs <- data.frame(
  test = c("Log-rank", "Cox PH (Immune-hot vs Immune-cold)"),
  statistic = c(logrank_pfs$chisq, cox_summary_pfs$coefficients[1, "z"]),
  hazard_ratio = c(NA, cox_summary_pfs$coefficients[1, "exp(coef)"]),
  ci_lower = c(NA, cox_summary_pfs$conf.int[1, "lower .95"]),
  ci_upper = c(NA, cox_summary_pfs$conf.int[1, "upper .95"]),
  p_value = c(logrank_p_pfs, cox_summary_pfs$coefficients[1, "Pr(>|z|)"])
)
write.csv(survival_test_summary_pfs, "outputs/tcga_bulk_seq_outputs/PFS_ImmuneGroup_test_summary.csv", row.names = FALSE)
cat("\nSurvival test results saved: outputs/tcga_bulk_seq_outputs/PFS_ImmuneGroup_test_summary.csv\n")

km_survival_plot(
  fit_pfs, strata_prefix = "ImmuneGroup=", colour_values = group_colors, legend_name = "Immune group",
  plot_title = "Progression-Free Survival by Immune Group",
  subtitle = paste0("Log-rank p = ", format.pval(logrank_p_pfs, digits = 3),
                    "  |  Cox HR (hot vs cold) = ", round(cox_summary_pfs$coefficients[1, "exp(coef)"], 2),
                    " (95% CI ", round(cox_summary_pfs$conf.int[1, "lower .95"], 2), "-",
                    round(cox_summary_pfs$conf.int[1, "upper .95"], 2), ")",
                    "  |  n = ", nrow(surv_data_pfs)),
  file_stem = "KM_PFS_ImmuneGroup",
  y_label = "Progression-free survival probability"
)

# 6. PFS by HPV status
surv_data_hpv_pfs <- clinical_merged %>%
  filter(!is.na(`Progress Free Survival (Months)`), !is.na(`Progression Free Status`), !is.na(HPVstatus)) %>%
  transmute(
    sample      = sample,
    HPVstatus   = HPVstatus,
    time_months = `Progress Free Survival (Months)`,
    event       = as.numeric(substr(`Progression Free Status`, 1, 1))
  )
cat("\nSamples with usable PFS + HPV data:", nrow(surv_data_hpv_pfs), "\n")
cat("Events (progression):", sum(surv_data_hpv_pfs$event), "/ censored:", sum(surv_data_hpv_pfs$event == 0), "\n")

surv_obj_hpv_pfs  <- Surv(surv_data_hpv_pfs$time_months, surv_data_hpv_pfs$event)
fit_hpv_pfs       <- survfit(surv_obj_hpv_pfs ~ HPVstatus, data = surv_data_hpv_pfs)
logrank_hpv_pfs   <- survdiff(surv_obj_hpv_pfs ~ HPVstatus, data = surv_data_hpv_pfs)
logrank_p_hpv_pfs <- pchisq(logrank_hpv_pfs$chisq, df = length(logrank_hpv_pfs$n) - 1, lower.tail = FALSE)

cox_fit_hpv_pfs     <- coxph(surv_obj_hpv_pfs ~ HPVstatus, data = surv_data_hpv_pfs)
cox_summary_hpv_pfs <- summary(cox_fit_hpv_pfs)
hpv_levels_pfs      <- sort(unique(surv_data_hpv_pfs$HPVstatus))  # [1] = Cox reference level
hpv_cox_label_pfs   <- paste0("Cox HR (", hpv_levels_pfs[2], " vs ", hpv_levels_pfs[1], ")")

cat("\n=== Log-rank test (PFS, HPV+ vs HPV-) ===\n")
cat("Chi-square =", round(logrank_hpv_pfs$chisq, 3), ", p =", format.pval(logrank_p_hpv_pfs, digits = 3), "\n")

cat("\n===", hpv_cox_label_pfs, "(PFS) ===\n")
cat("Hazard ratio:", round(cox_summary_hpv_pfs$coefficients[1, "exp(coef)"], 3),
    " | 95% CI:", round(cox_summary_hpv_pfs$conf.int[1, "lower .95"], 3), "-",
    round(cox_summary_hpv_pfs$conf.int[1, "upper .95"], 3),
    " | p =", format.pval(cox_summary_hpv_pfs$coefficients[1, "Pr(>|z|)"], digits = 3), "\n")

survival_test_summary_hpv_pfs <- data.frame(
  test = c("Log-rank", paste0("Cox PH (", hpv_cox_label_pfs, ")")),
  statistic = c(logrank_hpv_pfs$chisq, cox_summary_hpv_pfs$coefficients[1, "z"]),
  hazard_ratio = c(NA, cox_summary_hpv_pfs$coefficients[1, "exp(coef)"]),
  ci_lower = c(NA, cox_summary_hpv_pfs$conf.int[1, "lower .95"]),
  ci_upper = c(NA, cox_summary_hpv_pfs$conf.int[1, "upper .95"]),
  p_value = c(logrank_p_hpv_pfs, cox_summary_hpv_pfs$coefficients[1, "Pr(>|z|)"])
)
write.csv(survival_test_summary_hpv_pfs, "outputs/tcga_bulk_seq_outputs/PFS_HPVstatus_test_summary.csv", row.names = FALSE)
cat("\nSurvival test results saved: outputs/tcga_bulk_seq_outputs/PFS_HPVstatus_test_summary.csv\n")

km_survival_plot(
  fit_hpv_pfs, strata_prefix = "HPVstatus=", colour_values = hpv_colors, legend_name = "HPV status",
  plot_title = "Progression-Free Survival by HPV Status",
  subtitle = paste0("Log-rank p = ", format.pval(logrank_p_hpv_pfs, digits = 3),
                    "  |  ", hpv_cox_label_pfs, " = ", round(cox_summary_hpv_pfs$coefficients[1, "exp(coef)"], 2),
                    " (95% CI ", round(cox_summary_hpv_pfs$conf.int[1, "lower .95"], 2), "-",
                    round(cox_summary_hpv_pfs$conf.int[1, "upper .95"], 2), ")",
                    "  |  n = ", nrow(surv_data_hpv_pfs)),
  file_stem = "KM_PFS_HPVstatus",
  y_label = "Progression-free survival probability"
)

# 7. PFS by Immune Group x HPV status (4-group)
surv_data_combo_pfs <- clinical_merged %>%
  filter(!is.na(`Progress Free Survival (Months)`), !is.na(`Progression Free Status`), !is.na(HPVstatus)) %>%
  transmute(
    sample      = sample,
    ImmuneGroup = ImmuneGroup,
    HPVstatus   = HPVstatus,
    ComboGroup  = paste(ImmuneGroup, HPVstatus, sep = " / "),
    time_months = `Progress Free Survival (Months)`,
    event       = as.numeric(substr(`Progression Free Status`, 1, 1))
  )
cat("\nSamples with usable PFS + HPV + ImmuneGroup data:", nrow(surv_data_combo_pfs), "\n")
print(table(surv_data_combo_pfs$ComboGroup))

surv_obj_combo_pfs  <- Surv(surv_data_combo_pfs$time_months, surv_data_combo_pfs$event)
fit_combo_pfs       <- survfit(surv_obj_combo_pfs ~ ComboGroup, data = surv_data_combo_pfs)
logrank_combo_pfs   <- survdiff(surv_obj_combo_pfs ~ ComboGroup, data = surv_data_combo_pfs)
logrank_p_combo_pfs <- pchisq(logrank_combo_pfs$chisq, df = length(logrank_combo_pfs$n) - 1, lower.tail = FALSE)

cat("\n=== Log-rank test (PFS, 4-group: ImmuneGroup x HPV status) ===\n")
cat("Chi-square =", round(logrank_combo_pfs$chisq, 3), ", df =", length(logrank_combo_pfs$n) - 1,
    ", p =", format.pval(logrank_p_combo_pfs, digits = 3), "\n")

write.csv(
  data.frame(test = "Log-rank", statistic = logrank_combo_pfs$chisq,
             df = length(logrank_combo_pfs$n) - 1, p_value = logrank_p_combo_pfs),
  "outputs/tcga_bulk_seq_outputs/PFS_ImmuneGroup_x_HPVstatus_test_summary.csv",
  row.names = FALSE
)
cat("\nSurvival test results saved: outputs/tcga_bulk_seq_outputs/PFS_ImmuneGroup_x_HPVstatus_test_summary.csv\n")

km_survival_plot(
  fit_combo_pfs, strata_prefix = "ComboGroup=", colour_values = combo_colors,
  legend_name = "Immune group / HPV status",
  plot_title = "Progression-Free Survival by Immune Group and HPV Status",
  subtitle = paste0("Log-rank p = ", format.pval(logrank_p_combo_pfs, digits = 3),
                    "  |  n = ", nrow(surv_data_combo_pfs)),
  file_stem = "KM_PFS_ImmuneGroup_x_HPVstatus",
  legend_nrow = 2, width = 7.5, height = 6.5,
  y_label = "Progression-free survival probability"
)

## Multivariable Cox (PFS): same independence check as the OS block above.
cox_fit_multi_pfs <- coxph(surv_obj_combo_pfs ~ ImmuneGroup + HPVstatus, data = surv_data_combo_pfs)
cox_summary_multi_pfs <- summary(cox_fit_multi_pfs)

cat("\n=== Multivariable Cox (PFS ~ ImmuneGroup + HPVstatus) ===\n")
print(cox_summary_multi_pfs$coefficients)

multivariable_cox_pfs <- data.frame(
  term         = rownames(cox_summary_multi_pfs$coefficients),
  hazard_ratio = cox_summary_multi_pfs$coefficients[, "exp(coef)"],
  ci_lower     = cox_summary_multi_pfs$conf.int[, "lower .95"],
  ci_upper     = cox_summary_multi_pfs$conf.int[, "upper .95"],
  p_value      = cox_summary_multi_pfs$coefficients[, "Pr(>|z|)"]
)
write.csv(multivariable_cox_pfs, "outputs/tcga_bulk_seq_outputs/Cox_multivariable_PFS_ImmuneGroup_HPVstatus.csv", row.names = FALSE)
cat("\nMultivariable Cox results saved: outputs/tcga_bulk_seq_outputs/Cox_multivariable_PFS_ImmuneGroup_HPVstatus.csv\n")
