#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Consolidate existing FinnGen MR results. No MR or sensitivity estimator is run here.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(arg)) sub("^--file=", "", arg[1]) else ""
script_dir <- if (nzchar(script_file)) dirname(normalizePath(script_file)) else getwd()
root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
p <- function(...) file.path(root, ...)

files <- list(
  primary_all = p("results", "mr", "finngen", "results", "primary", "finngen_primary_mr_all_methods.tsv"),
  ivw = p("results", "mr", "finngen", "results", "primary", "finngen_primary_ivw_summary.tsv"),
  consistency = p("results", "mr", "finngen", "results", "primary", "finngen_mr_method_consistency.tsv"),
  sensitivity = p("results", "mr", "finngen", "sensitivity", "finngen_mr_sensitivity_summary.tsv"),
  heterogeneity = p("results", "mr", "finngen", "sensitivity", "finngen_mr_heterogeneity.tsv"),
  egger = p("results", "mr", "finngen", "sensitivity", "finngen_mr_egger_intercept.tsv"),
  loo = p("results", "mr", "finngen", "sensitivity", "finngen_leave_one_out_summary.tsv"),
  extraction = p("results", "mr", "finngen", "qc", "finngen_outcome_extraction_grch38_qc.tsv"),
  harmonisation = p("results", "mr", "finngen", "qc", "finngen_harmonisation_qc.tsv")
)
missing_files <- unlist(files)[!file.exists(unlist(files))]
if (length(missing_files)) stop("Missing source files: ", paste(missing_files, collapse = ", "), call. = FALSE)

out_dir <- p("results", "mr", "finngen", "summary")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
read_id <- function(path) fread(path, na.strings = c("NA", ""), colClasses = list(character = "outcome_id"))
pa <- read_id(files$primary_all); ivw <- read_id(files$ivw); con <- read_id(files$consistency)
sens <- read_id(files$sensitivity); het <- read_id(files$heterogeneity); eg <- read_id(files$egger)
loo <- read_id(files$loo); ext <- read_id(files$extraction); harm <- read_id(files$harmonisation)

sources <- list(ivw = ivw, consistency = con, sensitivity = sens, extraction = ext, harmonisation = harm)
for (nm in names(sources)) {
  x <- sources[[nm]]
  if (nrow(x) != 12L || uniqueN(x$outcome_id) != 12L) stop(nm, " must contain exactly 12 unique outcomes", call. = FALSE)
}
if (nrow(pa) != 60L || anyDuplicated(pa[, .(outcome_id, method_id)])) stop("Five-method primary table is incomplete or duplicated", call. = FALSE)
if (nrow(het) != 24L || anyDuplicated(het[, .(outcome_id, method_id)])) stop("Heterogeneity table is incomplete or duplicated", call. = FALSE)
if (nrow(eg) != 12L || uniqueN(eg$outcome_id) != 12L || nrow(loo) != 12L || uniqueN(loo$outcome_id) != 12L) {
  stop("Egger-intercept or leave-one-out summary is incomplete", call. = FALSE)
}

# Validate agreement between redundant saved outputs before consolidation.
ivw_from_all <- pa[method_id == "mr_ivw"]
cmp <- merge(ivw[, .(outcome_id, beta, SE, p_value, nsnp)],
             ivw_from_all[, .(outcome_id, beta_all = beta, SE_all = SE, p_all = p_value, nsnp_all = nsnp)], by = "outcome_id")
if (nrow(cmp) != 12L || any(cmp$beta != cmp$beta_all) || any(cmp$SE != cmp$SE_all) ||
    any(cmp$p_value != cmp$p_all) || any(cmp$nsnp != cmp$nsnp_all)) stop("Primary estimates disagree across script 24 outputs", call. = FALSE)
scmp <- merge(ivw[, .(outcome_id, beta, p_value)], sens[, .(outcome_id, ivw_beta, ivw_p_value)], by = "outcome_id")
if (any(scmp$beta != scmp$ivw_beta) || any(scmp$p_value != scmp$ivw_p_value)) stop("Primary estimates disagree with script 25 summary", call. = FALSE)
ecmp <- merge(eg[, .(outcome_id, egger_intercept, SE, p_value, p_egger_fdr)],
              sens[, .(outcome_id, e2 = egger_intercept, se2 = egger_intercept_SE,
                       p2 = egger_intercept_p_value, fdr2 = p_egger_fdr)], by = "outcome_id")
if (any(ecmp$egger_intercept != ecmp$e2) || any(ecmp$SE != ecmp$se2) ||
    any(ecmp$p_value != ecmp$p2) || any(ecmp$p_egger_fdr != ecmp$fdr2)) stop("Egger-intercept outputs disagree", call. = FALSE)

master <- Reduce(function(x, y) merge(x, y, by = "outcome_id", all = TRUE, sort = FALSE), list(
  ivw[, .(outcome, outcome_id, outcome_type, final_mr_ready_snps = nsnp,
          ivw_beta = beta, ivw_SE = SE, ivw_lower_95CI = beta_lower_95CI,
          ivw_upper_95CI = beta_upper_95CI, ivw_p_value = p_value,
          OR, OR_lower_95CI, OR_upper_95CI,
          ivw_p_bonferroni = p_bonferroni, ivw_bonferroni_significant = bonferroni_significant,
          ivw_p_fdr_bh = p_fdr_bh, ivw_fdr_significant = fdr_significant)],
  ext[, .(outcome_id, finngen_release = release, initial_asthma_instruments = requested_snps,
          finngen_matched_snps = matched_snps, finngen_coverage_percent = coverage_percent)],
  harm[, .(outcome_id, snps_entering_harmonisation, harmonisation_final_snps = final_mr_keep_snps,
           harmonisation_retention_percent = retention_percentage)],
  sens[, .(outcome_id, ivw_heterogeneity_Q = ivw_Q, ivw_heterogeneity_Q_df = ivw_Q_df,
           ivw_heterogeneity_p_value, mr_egger_intercept = egger_intercept,
           mr_egger_intercept_SE = egger_intercept_SE,
           mr_egger_intercept_p_value = egger_intercept_p_value,
           mr_egger_intercept_fdr_p_value = p_egger_fdr,
           leave_one_out_sign_change, most_influential_leave_one_out_snp,
           maximum_absolute_leave_one_out_beta_change)],
  con[, .(outcome_id, all_five_methods_same_direction = all_five_same_direction)]
))
master[, `:=`(
  heterogeneity_evidence = ivw_heterogeneity_p_value < 0.05,
  directional_pleiotropy_evidence = mr_egger_intercept_p_value < 0.05,
  directional_pleiotropy_fdr = mr_egger_intercept_fdr_p_value < 0.05
)]

classify <- function(ivw_p, het_flag, egger_flag, loo_flag) {
  if (!is.na(ivw_p) && ivw_p < 0.05) return("Not classified: nominal IVW association")
  labels <- character()
  if (isTRUE(het_flag) && !isTRUE(egger_flag)) labels <- c(labels, "B")
  if (isTRUE(egger_flag)) labels <- c(labels, "C")
  if (isTRUE(loo_flag)) labels <- c(labels, "D")
  if (!length(labels)) labels <- "A"
  paste(labels, collapse = "; ")
}
master[, interpretation_qc := mapply(classify, ivw_p_value, heterogeneity_evidence,
                                      directional_pleiotropy_evidence, leave_one_out_sign_change)]
master[, outcome_order := match(outcome_id, ivw$outcome_id)]
setorder(master, outcome_order)
if (nrow(master) != 12L || uniqueN(master$outcome_id) != 12L) stop("Final master does not contain 12 unique outcomes", call. = FALSE)
if (any(master$initial_asthma_instruments != 5254L) || any(master$final_mr_ready_snps != master$harmonisation_final_snps)) {
  stop("Instrument-count validation failed", call. = FALSE)
}
essential <- c("ivw_beta", "ivw_SE", "ivw_p_value", "ivw_heterogeneity_Q", "ivw_heterogeneity_p_value",
               "mr_egger_intercept", "mr_egger_intercept_SE", "mr_egger_intercept_p_value")
if (anyNA(master[, ..essential])) stop("Final master contains missing essential results", call. = FALSE)
if (anyNA(master[outcome_type == "binary", .(OR, OR_lower_95CI, OR_upper_95CI)]) ||
    any(!is.na(master[outcome_type == "continuous", OR]))) stop("Binary/continuous OR validation failed", call. = FALSE)
fwrite(master[, !c("outcome_order", "harmonisation_final_snps")],
       file.path(out_dir, "finngen_mr_final_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")

# Human-readable dissertation table: strings are deliberately rounded here only.
fmt <- function(x, digits = 5L) formatC(x, format = "g", digits = digits)
fmt_p <- function(x) formatC(x, format = "g", digits = 3L)
dissertation <- master[, .(
  Outcome = outcome,
  `n SNPs` = final_mr_ready_snps,
  `IVW beta` = fmt(ivw_beta),
  `IVW SE` = fmt(ivw_SE),
  `95% CI` = paste0(fmt(ivw_lower_95CI), " to ", fmt(ivw_upper_95CI)),
  `IVW p-value` = fmt_p(ivw_p_value),
  `OR (95% CI)` = fifelse(outcome_type == "binary",
    sprintf("%.4f (%.4f to %.4f)", OR, OR_lower_95CI, OR_upper_95CI), NA_character_),
  `IVW heterogeneity p-value` = fmt_p(ivw_heterogeneity_p_value),
  `MR-Egger intercept` = fmt(mr_egger_intercept),
  `MR-Egger intercept p-value` = fmt_p(mr_egger_intercept_p_value),
  `Leave-one-out sign change` = fifelse(leave_one_out_sign_change, "Yes", "No"),
  `Interpretation QC` = interpretation_qc
)]
fwrite(dissertation, file.path(out_dir, "finngen_mr_dissertation_table.tsv"), sep = "\t", quote = FALSE, na = "NA")

# Five-method comparison, retaining numerical estimates and explicit directions.
method_order <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode")
method_prefix <- c(mr_ivw = "ivw", mr_egger_regression = "mr_egger", mr_weighted_median = "weighted_median",
                   mr_simple_mode = "simple_mode", mr_weighted_mode = "weighted_mode")
comparison <- unique(pa[, .(outcome, outcome_id, outcome_type)])
for (id in method_order) {
  prefix <- method_prefix[[id]]
  m <- pa[method_id == id, .(outcome_id, beta, p_value)]
  setnames(m, c("beta", "p_value"), paste0(prefix, c("_beta", "_p_value")))
  m[, (paste0(prefix, "_direction")) := fifelse(get(paste0(prefix, "_beta")) > 0, "positive",
    fifelse(get(paste0(prefix, "_beta")) < 0, "negative", "zero"))]
  comparison <- merge(comparison, m, by = "outcome_id", all.x = TRUE, sort = FALSE)
}
comparison[, outcome_order := match(outcome_id, ivw$outcome_id)]
setorder(comparison, outcome_order)
setcolorder(comparison, c("outcome", "outcome_id", "outcome_type", setdiff(names(comparison), c("outcome", "outcome_id", "outcome_type", "outcome_order")), "outcome_order"))
if (anyNA(comparison)) stop("Five-method comparison contains missing results", call. = FALSE)
fwrite(comparison[, !"outcome_order"], file.path(out_dir, "finngen_mr_five_method_comparison.tsv"), sep = "\t", quote = FALSE, na = "NA")

# Plot helpers.
theme_dissertation <- theme_bw(base_size = 11) + theme(
  panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
  axis.title.y = element_blank(), plot.title = element_text(face = "bold"),
  strip.text = element_text(face = "bold"), legend.position = "none"
)
save_pair <- function(plot, stem, width, height) {
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in")
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = 300)
}
plot_data <- master[, .(
  outcome = factor(outcome, levels = rev(outcome)), outcome_type,
  estimate = fifelse(outcome_type == "binary", OR, ivw_beta),
  lower = fifelse(outcome_type == "binary", OR_lower_95CI, ivw_lower_95CI),
  upper = fifelse(outcome_type == "binary", OR_upper_95CI, ivw_upper_95CI),
  scale = fifelse(outcome_type == "binary", "Binary outcomes: odds ratio", "Continuous outcomes: beta")
)]
refs <- data.table(scale = c("Binary outcomes: odds ratio", "Continuous outcomes: beta"), null = c(1, 0))
all_plot <- ggplot(plot_data, aes(estimate, outcome)) +
  geom_vline(data = refs, aes(xintercept = null), linetype = 2, colour = "grey45", inherit.aes = FALSE) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.18, linewidth = 0.45) +
  geom_point(shape = 21, fill = "white", size = 2.5, stroke = 0.7) +
  facet_wrap(~scale, scales = "free", ncol = 1) +
  labs(title = "IVW estimates for FinnGen outcomes", x = "Effect estimate (panel-specific scale)") + theme_dissertation
save_pair(all_plot, "finngen_ivw_all_outcomes", 8.5, 7.5)

binary_plot <- ggplot(plot_data[outcome_type == "binary"], aes(estimate, outcome)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.18) +
  geom_point(shape = 21, fill = "white", size = 2.7) +
  labs(title = "IVW estimates: binary FinnGen outcomes", x = "Odds ratio (95% CI)") + theme_dissertation
save_pair(binary_plot, "finngen_ivw_binary_outcomes", 8, 4.5)
continuous_plot <- ggplot(plot_data[outcome_type == "continuous"], aes(estimate, outcome)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.18) +
  geom_point(shape = 21, fill = "white", size = 2.7) +
  labs(title = "IVW estimates: continuous FinnGen outcomes", x = "Beta (95% CI)") + theme_dissertation
save_pair(continuous_plot, "finngen_ivw_continuous_outcomes", 8, 5.5)

qc_long <- melt(master[, .(outcome,
  `Heterogeneity` = heterogeneity_evidence,
  `Egger intercept` = directional_pleiotropy_evidence,
  `Egger intercept FDR` = directional_pleiotropy_fdr,
  `Leave-one-out sign reversal` = leave_one_out_sign_change)],
  id.vars = "outcome", variable.name = "diagnostic", value.name = "evidence")
qc_long[, `:=`(
  outcome = factor(outcome, levels = rev(master$outcome)),
  diagnostic = factor(diagnostic, levels = c("Heterogeneity", "Egger intercept", "Egger intercept FDR", "Leave-one-out sign reversal")),
  label = fifelse(evidence, "Yes", "No")
)]
qc_plot <- ggplot(qc_long, aes(diagnostic, outcome, fill = evidence)) +
  geom_tile(colour = "white", linewidth = 0.7) + geom_text(aes(label = label), size = 3.2) +
  scale_fill_manual(values = c(`FALSE` = "grey92", `TRUE` = "grey55")) +
  labs(title = "FinnGen MR sensitivity overview", x = NULL, y = NULL,
       caption = "Descriptive diagnostic flags; no results were excluded") +
  theme_minimal(base_size = 10) + theme(panel.grid = element_blank(), legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1), plot.title = element_text(face = "bold"))
save_pair(qc_plot, "finngen_sensitivity_overview", 8.5, 6)

# Key counts and factual summary are recalculated from the consolidated data.
categories <- list(
  "Nominal IVW p < 0.05" = master[ivw_p_value < 0.05, outcome],
  "Bonferroni-significant outcomes" = master[ivw_bonferroni_significant == TRUE, outcome],
  "FDR-significant outcomes" = master[ivw_fdr_significant == TRUE, outcome],
  "Significant IVW heterogeneity" = master[heterogeneity_evidence == TRUE, outcome],
  "Significant Egger intercept" = master[directional_pleiotropy_evidence == TRUE, outcome],
  "FDR-significant Egger intercept" = master[directional_pleiotropy_fdr == TRUE, outcome],
  "Leave-one-out sign reversal" = master[leave_one_out_sign_change == TRUE, outcome],
  "All five methods same beta direction" = master[all_five_methods_same_direction == TRUE, outcome]
)
count_lines <- c("Total outcomes = 12", unlist(lapply(names(categories), function(nm) {
  vals <- categories[[nm]]
  c(paste0(nm, " = ", length(vals), "/12"), paste0("  Outcomes: ", if (length(vals)) paste(vals, collapse = "; ") else "None"))
})))
writeLines(count_lines, file.path(out_dir, "finngen_mr_key_counts.txt"))

known <- list(nominal = 0L, bonferroni = 0L, fdr = 0L, heterogeneity = 10L, egger = 5L,
              egger_fdr = 5L, loo = 3L)
observed <- c(nominal = length(categories[[1]]), bonferroni = length(categories[[2]]),
              fdr = length(categories[[3]]), heterogeneity = length(categories[[4]]),
              egger = length(categories[[5]]), egger_fdr = length(categories[[6]]), loo = length(categories[[7]]))
discrepancies <- names(known)[observed[names(known)] != unlist(known)]
expected_egger <- c("Allergic rhinitis", "Bronchiectasis", "Eosinophilic disease", "AST", "Triglycerides (fasting)")
expected_loo <- c("Eosinophilic disease", "NAFLD", "Triglycerides (fasting)")
if (!setequal(categories[[5]], expected_egger)) discrepancies <- c(discrepancies, "Egger outcome names")
if (!setequal(categories[[7]], expected_loo)) discrepancies <- c(discrepancies, "leave-one-out outcome names")

summary_lines <- c(
  "FinnGen MR numerical-results summary",
  "",
  sprintf("None of the 12 IVW estimates reached nominal p < 0.05."),
  sprintf("No IVW association survived Bonferroni or BH FDR correction."),
  sprintf("Significant IVW heterogeneity was detected for %d of 12 outcomes.", observed["heterogeneity"]),
  sprintf("MR-Egger intercept evidence of possible directional pleiotropy occurred for %d of 12 outcomes; %d remained significant after BH FDR correction.", observed["egger"], observed["egger_fdr"]),
  sprintf("Leave-one-out sign reversals occurred for %d of 12 outcomes.", observed["loo"]),
  "These are descriptive statistical diagnostics; no outcome or SNP was excluded.",
  paste0("Validation against known completed results: ", if (length(discrepancies)) paste("DISCREPANCY -", paste(unique(discrepancies), collapse = "; ")) else "PASS; no discrepancies detected.")
)
writeLines(summary_lines, file.path(out_dir, "finngen_mr_results_summary.txt"))

# Final artifact and label checks.
expected_outputs <- c(
  file.path(out_dir, c("finngen_mr_final_summary.tsv", "finngen_mr_dissertation_table.tsv",
                       "finngen_mr_five_method_comparison.tsv", "finngen_mr_key_counts.txt", "finngen_mr_results_summary.txt")),
  file.path(fig_dir, paste0(rep(c("finngen_ivw_all_outcomes", "finngen_ivw_binary_outcomes",
                                  "finngen_ivw_continuous_outcomes", "finngen_sensitivity_overview"), each = 2),
                               rep(c(".pdf", ".png"), 4)))
)
if (any(!file.exists(expected_outputs)) || any(file.info(expected_outputs)$size == 0)) stop("One or more summary outputs are missing or empty", call. = FALSE)
if (length(levels(plot_data$outcome)) != 12L || !setequal(c(0, 1), refs$null)) stop("Plot label or null-line validation failed", call. = FALSE)

message("Final outcomes: ", nrow(master))
message("Nominal IVW / Bonferroni / FDR significant: ", observed["nominal"], " / ", observed["bonferroni"], " / ", observed["fdr"])
message("Heterogeneity / Egger intercept / leave-one-out sign reversal: ", observed["heterogeneity"], " / ", observed["egger"], " / ", observed["loo"])
message("Validation discrepancies: ", if (length(discrepancies)) paste(unique(discrepancies), collapse = "; ") else "none")
message("Summary directory: ", out_dir)
