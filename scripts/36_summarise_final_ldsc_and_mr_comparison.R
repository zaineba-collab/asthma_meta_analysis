#!/usr/bin/env Rscript

# Formatting, joins, validation, and plotting only. No MR or LDSC model is run.
suppressPackageStartupMessages(library(ggplot2))

root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))] |> sub("^--file=", "", x = _)), ".."))
setwd(root)

read_tsv <- function(path) read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
                                      stringsAsFactors = FALSE, na.strings = "NA")
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
yn <- function(x) ifelse(is.na(x), NA, ifelse(tolower(as.character(x)) %in% c("yes", "true"), "TRUE", "FALSE"))
fmt_num <- function(x, digits = 3) ifelse(is.na(x), "—", formatC(x, format = "f", digits = digits))
fmt_p <- function(x) vapply(x, function(value) {
  if (is.na(value)) "—" else format(value, digits = 3, scientific = value < 0.001, trim = TRUE)
}, character(1))
fmt_sig <- function(x, digits = 4) ifelse(is.na(x), "—", format(x, digits = digits, scientific = FALSE, trim = TRUE))

rg <- read_tsv("results/ldsc/qc/ldsc_rg_summary.tsv")
rg_all <- read_tsv("results/ldsc/qc/ldsc_rg_all_12_outcomes.tsv")
h2 <- read_tsv("results/ldsc/qc/ldsc_h2_summary.tsv")
h2_ready <- read_tsv("results/ldsc/qc/ldsc_h2_readiness.tsv")
mr <- read_tsv("results/mr/finngen/summary/finngen_mr_final_summary.tsv")

stopifnot(nrow(rg) == 12, nrow(rg_all) == 12, nrow(h2) == 13, nrow(mr) == 12)
stopifnot(!anyDuplicated(rg$outcome), !anyDuplicated(mr$outcome), !anyDuplicated(h2$trait))
stopifnot(all(rg$outcome == rg_all$outcome))

display_order <- c("Allergic rhinitis", "Bronchiectasis", "Eosinophilic disease", "BMI IRN",
                   "NAFLD", "CRP", "ALT", "AST", "GGT", "HDL cholesterol",
                   "LDL cholesterol", "Fasting triglycerides")
stopifnot(identical(rg$outcome, display_order))

h2_idx <- match(rg$outcome, h2$trait)
stopifnot(!anyNA(h2_idx))
stopifnot(abs(rg$outcome_univariate_h2 - h2$h2_obs[h2_idx]) < 1e-15 | is.na(rg$outcome_univariate_h2))

ldsc_table <- data.frame(
  Outcome = rg$outcome,
  univariate_h2 = h2$h2_obs[h2_idx],
  h2_SE = h2$h2_se[h2_idx],
  h2_Z = h2$h2_z[h2_idx],
  rg = rg$rg,
  rg_SE = rg$rg_se,
  rg_lower_95CI = rg$rg_lower_95CI,
  rg_upper_95CI = rg$rg_upper_95CI,
  rg_P = rg$rg_p,
  Bonferroni_adjusted_P = rg$p_bonferroni_12,
  BH_FDR_P = rg$p_fdr_bh_11,
  Bonferroni_significant = yn(rg$bonferroni12_significant),
  FDR_significant = yn(rg$fdr11_significant),
  cross_trait_intercept = rg$cross_trait_intercept,
  cross_trait_intercept_SE = rg$cross_trait_intercept_se,
  SNPs = rg$SNPs,
  h2_readiness = rg$h2_readiness_category,
  status = rg$run_status,
  check.names = FALSE
)
write_tsv(ldsc_table, "results/ldsc/summary/ldsc_rg_dissertation_table.tsv")

ldsc_formatted <- data.frame(
  Outcome = ldsc_table$Outcome,
  h2 = fmt_sig(ldsc_table$univariate_h2),
  `h2 SE` = fmt_sig(ldsc_table$h2_SE),
  `h2 Z` = fmt_num(ldsc_table$h2_Z),
  rg = fmt_num(ldsc_table$rg),
  `rg SE` = fmt_num(ldsc_table$rg_SE),
  `rg 95% CI` = ifelse(is.na(ldsc_table$rg), "—", paste0(fmt_num(ldsc_table$rg_lower_95CI), " to ", fmt_num(ldsc_table$rg_upper_95CI))),
  `rg P` = fmt_p(ldsc_table$rg_P),
  `Bonferroni P` = fmt_p(ldsc_table$Bonferroni_adjusted_P),
  `BH-FDR P` = fmt_p(ldsc_table$BH_FDR_P),
  `Multiple-testing result` = ifelse(is.na(ldsc_table$rg), "Not estimable",
                                     ifelse(ldsc_table$Bonferroni_significant == "TRUE", "Bonferroni significant", "Not Bonferroni significant")),
  `Cross-trait intercept` = fmt_num(ldsc_table$cross_trait_intercept, 4),
  SNPs = ifelse(is.na(ldsc_table$SNPs), "—", format(ldsc_table$SNPs, scientific = FALSE, big.mark = ",")),
  `h2 readiness` = ldsc_table$h2_readiness,
  Status = ldsc_table$status,
  check.names = FALSE
)
write_tsv(ldsc_formatted, "results/ldsc/summary/ldsc_rg_dissertation_table_formatted.tsv")

estimable <- !is.na(rg$rg)
nominal <- estimable & rg$rg_p < 0.05
bonf <- estimable & tolower(rg$bonferroni12_significant) == "yes"
fdr <- estimable & tolower(rg$fdr11_significant) == "yes"
key_lines <- c(
  "LDSC key result counts (recalculated from authoritative source tables)",
  sprintf("planned outcomes = %d", nrow(rg)),
  sprintf("estimable rg outcomes = %d", sum(estimable)),
  sprintf("non-estimable outcomes = %d", sum(!estimable)),
  sprintf("nominal rg P < 0.05 = %d", sum(nominal)),
  sprintf("Bonferroni significant = %d", sum(bonf)),
  sprintf("BH-FDR significant = %d", sum(fdr)),
  sprintf("significant positive rg = %d", sum(bonf & rg$rg > 0, na.rm = TRUE)),
  sprintf("significant negative rg = %d", sum(bonf & rg$rg < 0, na.rm = TRUE)),
  paste("nominal outcomes =", paste(rg$outcome[nominal], collapse = "; ")),
  paste("Bonferroni outcomes =", paste(rg$outcome[bonf], collapse = "; ")),
  paste("BH-FDR outcomes =", paste(rg$outcome[fdr], collapse = "; "))
)
dir.create("results/ldsc/summary", recursive = TRUE, showWarnings = FALSE)
writeLines(key_lines, "results/ldsc/summary/ldsc_key_results.txt")

# Forest plot: no point is plotted for the non-estimable outcome.
plot_rg <- rg
plot_rg$outcome_plot <- factor(plot_rg$outcome, levels = rev(display_order))
plot_rg$multiple_testing <- ifelse(!estimable, "Not estimable",
                                   ifelse(bonf, "Bonferroni significant", "Not Bonferroni significant"))
p_rg <- ggplot(plot_rg, aes(y = outcome_plot)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.45, linetype = 2) +
  geom_errorbar(data = subset(plot_rg, !is.na(rg)), aes(xmin = rg_lower_95CI, xmax = rg_upper_95CI),
                orientation = "y", width = 0.18, linewidth = 0.55) +
  geom_point(data = subset(plot_rg, !is.na(rg)), aes(x = rg, shape = multiple_testing), size = 2.5, fill = "white") +
  geom_text(data = subset(plot_rg, is.na(rg)), aes(x = 0, label = "Not estimable"), hjust = 0.5, size = 3.2, fontface = "italic") +
  scale_shape_manual(values = c("Bonferroni significant" = 16, "Not Bonferroni significant" = 21)) +
  scale_y_discrete(limits = rev(display_order), drop = FALSE) +
  labs(x = "Genetic correlation with asthma (rg)", y = NULL, shape = NULL,
       title = "Asthma genetic correlations", caption = "Points show rg with 95% confidence intervals; eosinophilic disease was not estimable.") +
  theme_classic(base_size = 11) + theme(legend.position = "bottom", plot.title.position = "plot")
dir.create("results/ldsc/summary/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/ldsc/summary/figures/ldsc_asthma_rg_forest.pdf", p_rg, width = 8.2, height = 6.2)
ggsave("results/ldsc/summary/figures/ldsc_asthma_rg_forest.png", p_rg, width = 8.2, height = 6.2, dpi = 300)

# Observed-scale h2 QC plot.
h2$lower <- h2$h2_obs - 1.96 * h2$h2_se
h2$upper <- h2$h2_obs + 1.96 * h2$h2_se
h2$trait_plot <- factor(h2$trait, levels = rev(h2$trait))
h2$estimate_status <- ifelse(h2$h2_obs <= 0, "Non-positive h2", "Positive h2")
p_h2 <- ggplot(h2, aes(y = trait_plot)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.45, linetype = 2) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.18, linewidth = 0.55) +
  geom_point(aes(x = h2_obs, shape = estimate_status), size = 2.5, fill = "white") +
  scale_shape_manual(values = c("Non-positive h2" = 4, "Positive h2" = 21)) +
  labs(x = "Observed-scale LDSC h2 (95% CI)", y = NULL, shape = NULL,
       title = "Univariate LDSC heritability QC", caption = "Observed-scale estimates; eosinophilic disease has a non-positive point estimate.") +
  theme_classic(base_size = 11) + theme(legend.position = "bottom", plot.title.position = "plot")
ggsave("results/ldsc/summary/figures/ldsc_h2_qc.pdf", p_h2, width = 8.2, height = 6.4)
ggsave("results/ldsc/summary/figures/ldsc_h2_qc.png", p_h2, width = 8.2, height = 6.4, dpi = 300)

# Normalize outcome names only for an explicit one-to-one join.
mr$join_outcome <- mr$outcome
mr$join_outcome[mr$join_outcome == "BMI (inverse-rank normalized)"] <- "BMI IRN"
mr$join_outcome[mr$join_outcome == "Triglycerides (fasting)"] <- "Fasting triglycerides"
stopifnot(!anyDuplicated(mr$join_outcome), setequal(mr$join_outcome, rg$outcome))
mi <- match(rg$outcome, mr$join_outcome)
stopifnot(!anyNA(mi))
mrj <- mr[mi, ]
mr_sig <- as.logical(mrj$ivw_bonferroni_significant)
ldsc_sig <- bonf
ldsc_estimable <- estimable
pattern <- ifelse(!ldsc_estimable, "LDSC_NOT_ESTIMABLE",
                  ifelse(mr_sig & ldsc_sig, "BOTH_SIGNIFICANT",
                         ifelse(mr_sig, "MR_SIGNIFICANT_LDSC_NOT_SIGNIFICANT",
                                ifelse(ldsc_sig, "LDSC_SIGNIFICANT_MR_NOT_SIGNIFICANT", "NEITHER_SIGNIFICANT"))))

comparison <- data.frame(
  Outcome = rg$outcome,
  MR_ready_SNPs = mrj$final_mr_ready_snps,
  IVW_beta = mrj$ivw_beta,
  IVW_SE = mrj$ivw_SE,
  IVW_P = mrj$ivw_p_value,
  IVW_Bonferroni_significant = mr_sig,
  MR_heterogeneity_P = mrj$ivw_heterogeneity_p_value,
  MR_Egger_intercept_P = mrj$mr_egger_intercept_p_value,
  MR_Egger_intercept_FDR_P = mrj$mr_egger_intercept_fdr_p_value,
  leave_one_out_sign_change = mrj$leave_one_out_sign_change,
  rg_estimable = ldsc_estimable,
  rg = rg$rg,
  rg_SE = rg$rg_se,
  rg_lower_95CI = rg$rg_lower_95CI,
  rg_upper_95CI = rg$rg_upper_95CI,
  rg_P = rg$rg_p,
  rg_Bonferroni_significant = ldsc_sig,
  rg_FDR_significant = fdr,
  outcome_h2 = rg$outcome_univariate_h2,
  outcome_h2_Z = rg$outcome_univariate_h2_z,
  LDSC_readiness = rg$h2_readiness_category,
  LDSC_status = rg$run_status,
  mr_ldsc_pattern = pattern,
  check.names = FALSE
)
write_tsv(comparison, "results/combined/finngen_mr_ldsc_comparison.tsv")

patterns <- data.frame(
  Outcome = comparison$Outcome,
  MR_Bonferroni_significant = comparison$IVW_Bonferroni_significant,
  LDSC_Bonferroni_significant = ifelse(comparison$rg_estimable, comparison$rg_Bonferroni_significant, NA),
  LDSC_estimable = comparison$rg_estimable,
  MR_pleiotropy_warning = as.logical(mrj$directional_pleiotropy_evidence),
  MR_heterogeneity_warning = as.logical(mrj$heterogeneity_evidence),
  LDSC_h2_warning = ifelse(rg$h2_readiness_category == "A — STRONGER_H2_SIGNAL", "FALSE", "TRUE"),
  pattern = pattern,
  check.names = FALSE
)
write_tsv(patterns, "results/combined/finngen_mr_ldsc_result_patterns.tsv")

combined_dissertation <- data.frame(
  Outcome = rg$outcome,
  `MR IVW beta` = fmt_num(mrj$ivw_beta, 4),
  `MR P` = fmt_p(mrj$ivw_p_value),
  `MR interpretation flag` = ifelse(mr_sig, "Primary MR evidence", "No primary MR evidence"),
  `LDSC rg` = fmt_num(rg$rg),
  `LDSC 95% CI` = ifelse(estimable, paste0(fmt_num(rg$rg_lower_95CI), " to ", fmt_num(rg$rg_upper_95CI)), "—"),
  `LDSC P` = fmt_p(rg$rg_p),
  `LDSC interpretation flag` = ifelse(!estimable, "Not estimable",
                                       ifelse(ldsc_sig, "Significant genetic correlation", "No significant genetic correlation")),
  check.names = FALSE
)
write_tsv(combined_dissertation, "results/combined/finngen_mr_ldsc_dissertation_table.tsv")

counts <- table(factor(pattern, levels = c("BOTH_SIGNIFICANT", "LDSC_SIGNIFICANT_MR_NOT_SIGNIFICANT",
                                          "MR_SIGNIFICANT_LDSC_NOT_SIGNIFICANT", "NEITHER_SIGNIFICANT", "LDSC_NOT_ESTIMABLE")))
names_for <- function(category) {
  x <- rg$outcome[pattern == category]
  if (length(x)) paste(x, collapse = "; ") else "None"
}
count_lines <- c(
  "MR–LDSC key counts (recalculated from authoritative source tables)",
  sprintf("outcomes planned = %d", nrow(comparison)),
  sprintf("MR Bonferroni significant = %d", sum(mr_sig)),
  paste("MR Bonferroni significant outcomes =", if (sum(mr_sig)) paste(rg$outcome[mr_sig], collapse = "; ") else "None"),
  sprintf("LDSC Bonferroni significant = %d", sum(ldsc_sig)),
  paste("LDSC Bonferroni significant outcomes =", paste(rg$outcome[ldsc_sig], collapse = "; ")),
  sprintf("both significant = %d", counts[["BOTH_SIGNIFICANT"]]),
  paste("both significant outcomes =", names_for("BOTH_SIGNIFICANT")),
  sprintf("LDSC significant / MR not significant = %d", counts[["LDSC_SIGNIFICANT_MR_NOT_SIGNIFICANT"]]),
  paste("LDSC significant / MR not significant outcomes =", names_for("LDSC_SIGNIFICANT_MR_NOT_SIGNIFICANT")),
  sprintf("neither significant = %d", counts[["NEITHER_SIGNIFICANT"]]),
  paste("neither significant outcomes =", names_for("NEITHER_SIGNIFICANT")),
  sprintf("LDSC not estimable = %d", counts[["LDSC_NOT_ESTIMABLE"]]),
  paste("LDSC not estimable outcomes =", names_for("LDSC_NOT_ESTIMABLE"))
)
writeLines(count_lines, "results/combined/mr_ldsc_key_counts.txt")

# Categorical overview: effect scales are not mixed.
overview <- rbind(
  data.frame(Outcome = rg$outcome, Evidence = "MR evidence\n(IVW Bonferroni)",
             Category = ifelse(mr_sig, "Significant", "Not significant"),
             Label = paste0("P=", fmt_p(mrj$ivw_p_value))),
  data.frame(Outcome = rg$outcome, Evidence = "LDSC genetic correlation\n(rg Bonferroni)",
             Category = ifelse(!estimable, "Not estimable", ifelse(ldsc_sig, "Significant", "Not significant")),
             Label = ifelse(estimable, paste0("rg=", fmt_num(rg$rg)), "Not estimable"))
)
overview$Outcome <- factor(overview$Outcome, levels = rev(display_order))
overview$Evidence <- factor(overview$Evidence, levels = c("MR evidence\n(IVW Bonferroni)", "LDSC genetic correlation\n(rg Bonferroni)"))
p_overview <- ggplot(overview, aes(x = Evidence, y = Outcome, fill = Category)) +
  geom_tile(colour = "white", linewidth = 1) + geom_text(aes(label = Label), size = 3.1) +
  scale_fill_manual(values = c("Significant" = "#a9c6d9", "Not significant" = "#eeeeee", "Not estimable" = "#d9d9d9")) +
  labs(x = NULL, y = NULL, fill = NULL, title = "MR and LDSC evidence overview",
       caption = "MR and LDSC answer different scientific questions; their effect estimates are not directly comparable.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(face = "bold"), legend.position = "bottom", plot.title.position = "plot")
dir.create("results/combined/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/combined/figures/mr_ldsc_comparison_overview.pdf", p_overview, width = 8.5, height = 6.4)
ggsave("results/combined/figures/mr_ldsc_comparison_overview.png", p_overview, width = 8.5, height = 6.4, dpi = 300)

writeLines(c(
  "LDSC estimator/configuration distinction for HDL and LDL",
  "Standalone --h2 used LDSC's automatic two-step estimator with cutoff 30.",
  "Unconstrained --rg used one-step IRWLS for its pairwise trait-specific h2 quantities.",
  "This distinction produced different trait-specific h2 estimates for HDL and LDL.",
  "Standalone h2 restricted to the exact asthma-overlap SNP sets reproduced the original univariate h2 estimates.",
  "SNP loss therefore did not explain the difference.",
  "Parser accuracy, reference LD scores, regression weights, environment, commit, and unconstrained intercept handling were validated.",
  "The authoritative rg estimates were retained unchanged; this is documented as an estimator/configuration distinction, not overstated as an error."
), "results/ldsc/summary/ldsc_estimator_note.txt")

writeLines(c(
  "FinnGen laboratory-trait LDSC methodological note",
  "FinnGen quantitative laboratory GWAS were generated using mixed-model analyses.",
  "Nominal-N scaling can affect absolute LDSC h2 estimates, so univariate h2 is treated primarily as QC here.",
  "Genetic correlation (rg) is the primary cross-trait LDSC quantity of interest.",
  "No sample-size alteration or summary-statistic rescaling was performed."
), "results/ldsc/summary/ldsc_mixed_model_methodological_note.txt")

strong_pos <- rg$outcome[which.max(ifelse(estimable, rg$rg, -Inf))]
strong_neg <- rg$outcome[which.min(ifelse(estimable, rg$rg, Inf))]
summary_lines <- c(
  "Machine-generated MR–LDSC numerical results summary",
  sprintf("No MR IVW association survived Bonferroni correction (%d of 12).", sum(mr_sig)),
  sprintf("Seven genetic correlations survived Bonferroni correction (%d of 12 planned outcomes).", sum(ldsc_sig)),
  sprintf("The strongest positive genetic correlation was asthma–%s (rg=%s).", strong_pos, rg$rg[rg$outcome == strong_pos]),
  sprintf("The strongest negative genetic correlation was asthma–%s (rg=%s).", strong_neg, rg$rg[rg$outcome == strong_neg]),
  "Eosinophilic disease rg was not estimable because its univariate LDSC h2 estimate was non-positive.",
  "These statements are descriptive and do not assign causal meaning."
)
writeLines(summary_lines, "results/combined/mr_ldsc_results_summary.txt")

# Final internal consistency checks against source values and flags.
stopifnot(nrow(comparison) == 12, !anyDuplicated(comparison$Outcome))
stopifnot(identical(comparison$IVW_beta, mrj$ivw_beta), identical(comparison$IVW_P, mrj$ivw_p_value))
stopifnot(identical(comparison$rg, rg$rg), identical(comparison$rg_P, rg$rg_p))
stopifnot(comparison$LDSC_status[comparison$Outcome == "Eosinophilic disease"] == "NOT_ESTIMABLE_LOW_H2")
stopifnot(sum(mr_sig) == 0, sum(ldsc_sig) == 7)
stopifnot(sum(pattern == "LDSC_SIGNIFICANT_MR_NOT_SIGNIFICANT") == 7)

writeLines(c(
  "Final MR–LDSC reporting validation",
  "exactly 12 outcomes = PASS",
  "duplicate outcome rows = 0",
  "MR values match authoritative script-26 summary = PASS",
  "LDSC values and flags match authoritative script-33 summary = PASS",
  "eosinophilic disease remains NOT_ESTIMABLE_LOW_H2 = PASS",
  "MR Bonferroni significant count = 0",
  "LDSC Bonferroni significant count = 7",
  "effect scales are separated in the overview figure = PASS",
  "no MR or LDSC analysis function is called by script 36 = PASS"
), "results/combined/mr_ldsc_validation.txt")

cat("validation=PASS outcomes=12 mr_bonferroni=", sum(mr_sig),
    " ldsc_bonferroni=", sum(ldsc_sig), "\n", sep = "")
