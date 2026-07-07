#!/usr/bin/env Rscript

library(data.table)
library(TwoSampleMR)
library(ggplot2)

# ============================================================
# 24_run_mr_analysis.R
#
# Primary MR direction:
#   Asthma -> configured continuous outcome
#
# Supervisor plan:
#   Use already LD-clumped asthma signals as MR instruments and do
#   not clump again using TwoSampleMR.
#
# Purpose:
#   Run the main TwoSampleMR estimates for asthma as the exposure
#   and a configured continuous outcome.
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)
source(file.path(script_dir, "mr_config_helpers.R"))

mr_config <- read_mr_outcome_config(project_dir)
outcome_label <- mr_config$outcome_label
outcome_name <- mr_config$outcome_name
outcome_type <- mr_config$outcome_type

primary_dir <- file.path(project_dir, "results", "mr", paste0("asthma_to_", outcome_label))
harmonised_file <- file.path(primary_dir, "harmonised", paste0("asthma_", outcome_label, "_harmonised.tsv"))
analysis_dir <- file.path(primary_dir, "analysis")
figures_dir <- file.path(primary_dir, "figures")

raw_results_file <- file.path(analysis_dir, paste0("asthma_", outcome_label, "_mr_results.tsv"))
publication_results_file <- file.path(analysis_dir, paste0("asthma_", outcome_label, "_mr_results_publication.tsv"))
scatter_png_file <- file.path(figures_dir, paste0("asthma_", outcome_label, "_mr_scatter_plot.png"))
scatter_pdf_file <- file.path(figures_dir, paste0("asthma_", outcome_label, "_mr_scatter_plot.pdf"))

minimum_usable_instruments <- as.integer(Sys.getenv("MR_MIN_USABLE_INSTRUMENTS", unset = "10"))

if (is.na(minimum_usable_instruments) || minimum_usable_instruments < 1) {
  stop("MR_MIN_USABLE_INSTRUMENTS must be a positive integer.", call. = FALSE)
}

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

required_columns <- c(
  "SNP",
  "beta.exposure",
  "se.exposure",
  "effect_allele.exposure",
  "other_allele.exposure",
  "beta.outcome",
  "se.outcome",
  "effect_allele.outcome",
  "other_allele.outcome",
  "exposure",
  "outcome",
  "id.exposure",
  "id.outcome",
  "mr_keep"
)

stop_if_missing(harmonised_file, paste("Harmonised asthma-", outcome_name, " file", sep = ""))
dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

message("Primary MR direction: Asthma -> ", outcome_name)
print_mr_config(mr_config, primary_dir)
message("Reading harmonised data: ", harmonised_file)
harmonised <- fread(harmonised_file)

missing_columns <- setdiff(required_columns, names(harmonised))

if (length(missing_columns) > 0) {
  stop("Harmonised file is missing required column(s): ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

total_harmonised_rows <- nrow(harmonised)
excluded_mr_keep_false <- nrow(harmonised[mr_keep == FALSE])
analysis_dat <- harmonised[mr_keep == TRUE]
usable_instruments <- uniqueN(analysis_dat$SNP)

message("Total harmonised rows: ", total_harmonised_rows)
message("Rows excluded because mr_keep == FALSE: ", excluded_mr_keep_false)
message("Final usable instruments: ", usable_instruments)

if (usable_instruments < minimum_usable_instruments) {
  stop(
    "Only ", usable_instruments, " usable MR instruments remain. ",
    "This is below MR_MIN_USABLE_INSTRUMENTS = ", minimum_usable_instruments, ".",
    call. = FALSE
  )
}

duplicate_rows <- analysis_dat[duplicated(SNP) | duplicated(SNP, fromLast = TRUE)]

if (nrow(duplicate_rows) > 0) {
  stop("Duplicate SNP rows found after mr_keep filtering: ", paste(sort(unique(duplicate_rows$SNP)), collapse = ", "), call. = FALSE)
}

main_methods <- c(
  "mr_ivw",
  "mr_egger_regression",
  "mr_weighted_median",
  "mr_weighted_mode",
  "mr_simple_mode"
)

message("Running main MR methods: ", paste(main_methods, collapse = ", "))

mr_results <- mr(dat = as.data.frame(analysis_dat), method_list = main_methods)
setDT(mr_results)

if (nrow(mr_results) == 0) {
  stop("TwoSampleMR returned no MR results.", call. = FALSE)
}

mr_results[, beta_lci95 := b - 1.96 * se]
mr_results[, beta_uci95 := b + 1.96 * se]

# Continuous outcomes are reported on their original beta scale. Odds ratios are
# not calculated unless a future config row explicitly supports a binary scale.
publication_results <- mr_results[
  ,
  .(
    method,
    number_of_snps = nsnp,
    beta_estimate = b,
    standard_error = se,
    p_value = pval,
    odds_ratio = NA_real_,
    beta_lci95,
    beta_uci95,
    odds_ratio_lci95 = NA_real_,
    odds_ratio_uci95 = NA_real_,
    odds_ratio_note = paste("Not calculated:", outcome_name, "is a", outcome_type, "outcome")
  )
]

fwrite(mr_results, raw_results_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(publication_results, publication_results_file, sep = "\t", quote = FALSE, na = "NA")

message("Written raw MR results: ", raw_results_file)
message("Written publication MR results: ", publication_results_file)

scatter_plots <- mr_scatter_plot(mr_results = as.data.frame(mr_results), dat = as.data.frame(analysis_dat))

if (length(scatter_plots) > 0) {
  ggsave(scatter_png_file, plot = scatter_plots[[1]], width = 8, height = 6, dpi = 300)
  ggsave(scatter_pdf_file, plot = scatter_plots[[1]], width = 8, height = 6)
  message("Written MR scatter plot PNG: ", scatter_png_file)
  message("Written MR scatter plot PDF: ", scatter_pdf_file)
} else {
  warning("TwoSampleMR did not return a scatter plot.", call. = FALSE)
}

message("")
message("Main MR analysis complete. Sensitivity analyses have not been run.")
message("Usable instruments analysed: ", usable_instruments)
message("Concise MR summary:")
print(publication_results[, .(method, number_of_snps, beta_estimate, standard_error, p_value, beta_lci95, beta_uci95)])
message("")
message("Scientific caution: these are preliminary main MR estimates and do not prove causality.")
