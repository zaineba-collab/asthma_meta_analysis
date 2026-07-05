#!/usr/bin/env Rscript

library(data.table)
library(TwoSampleMR)
library(ggplot2)

# ============================================================
# 25_mr_sensitivity_analysis.R
#
# Primary MR direction:
#   Asthma -> BMI
#
# Supervisor plan:
#   Use already LD-clumped asthma signals as MR instruments and do
#   not clump again using TwoSampleMR.
#
# Purpose:
#   Run sensitivity analyses for the asthma -> BMI TwoSampleMR
#   analysis using the already harmonised dataset.
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

primary_dir <- file.path(project_dir, "results", "mr", "asthma_to_bmi")
harmonised_file <- file.path(primary_dir, "harmonised", "asthma_bmi_harmonised.tsv")
sensitivity_dir <- file.path(primary_dir, "sensitivity")
figures_dir <- file.path(primary_dir, "figures")

heterogeneity_file <- file.path(sensitivity_dir, "asthma_bmi_heterogeneity.tsv")
egger_intercept_file <- file.path(sensitivity_dir, "asthma_bmi_egger_intercept.tsv")
single_snp_file <- file.path(sensitivity_dir, "asthma_bmi_single_snp.tsv")
leave_one_out_file <- file.path(sensitivity_dir, "asthma_bmi_leave_one_out.tsv")
sensitivity_summary_file <- file.path(sensitivity_dir, "asthma_bmi_sensitivity_summary.tsv")

forest_png_file <- file.path(figures_dir, "asthma_bmi_single_snp_forest_plot.png")
forest_pdf_file <- file.path(figures_dir, "asthma_bmi_single_snp_forest_plot.pdf")
loo_png_file <- file.path(figures_dir, "asthma_bmi_leave_one_out_plot.png")
loo_pdf_file <- file.path(figures_dir, "asthma_bmi_leave_one_out_plot.pdf")
funnel_png_file <- file.path(figures_dir, "asthma_bmi_funnel_plot.png")
funnel_pdf_file <- file.path(figures_dir, "asthma_bmi_funnel_plot.pdf")

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

get_metric <- function(dt, method_pattern, column_name) {
  if (!column_name %in% names(dt)) {
    return(NA_real_)
  }

  rows <- grep(method_pattern, dt$method, ignore.case = TRUE)

  if (length(rows) == 0) {
    return(NA_real_)
  }

  as.numeric(dt[[column_name]][rows[1]])
}

save_plot_pair <- function(plot_object, png_file, pdf_file, width, height) {
  if (is.null(plot_object)) {
    return("plot object was NULL")
  }

  tryCatch(
    {
      ggsave(png_file, plot = plot_object, width = width, height = height, dpi = 300)
      ggsave(pdf_file, plot = plot_object, width = width, height = height)
      NA_character_
    },
    error = function(e) conditionMessage(e)
  )
}

stop_if_missing(harmonised_file, "Harmonised asthma-BMI file")
dir.create(sensitivity_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

message("Primary MR direction: Asthma -> BMI")
message("Using already LD-clumped asthma signals. No additional TwoSampleMR clumping performed.")
message("Reading harmonised data: ", harmonised_file)
harmonised <- fread(harmonised_file)

missing_columns <- setdiff(required_columns, names(harmonised))

if (length(missing_columns) > 0) {
  stop("Harmonised file is missing required column(s): ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

analysis_dat <- harmonised[mr_keep == TRUE]
usable_instruments <- uniqueN(analysis_dat$SNP)

message("Usable SNPs with mr_keep == TRUE: ", usable_instruments)

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

analysis_df <- as.data.frame(analysis_dat)

message("Running heterogeneity analysis.")
heterogeneity <- mr_heterogeneity(dat = analysis_df, method_list = c("mr_egger_regression", "mr_ivw"))
setDT(heterogeneity)
fwrite(heterogeneity, heterogeneity_file, sep = "\t", quote = FALSE, na = "NA")

message("Running MR-Egger intercept test for directional pleiotropy.")
egger_intercept <- mr_pleiotropy_test(analysis_df)
setDT(egger_intercept)
fwrite(egger_intercept, egger_intercept_file, sep = "\t", quote = FALSE, na = "NA")

message("Running single-SNP analysis.")
single_snp <- mr_singlesnp(analysis_df)
setDT(single_snp)
fwrite(single_snp, single_snp_file, sep = "\t", quote = FALSE, na = "NA")

message("Running leave-one-out analysis.")
leave_one_out <- mr_leaveoneout(analysis_df)
setDT(leave_one_out)
fwrite(leave_one_out, leave_one_out_file, sep = "\t", quote = FALSE, na = "NA")

ivw_q <- get_metric(heterogeneity, "Inverse variance weighted", "Q")
ivw_q_pval <- get_metric(heterogeneity, "Inverse variance weighted", "Q_pval")
egger_q <- get_metric(heterogeneity, "MR Egger", "Q")
egger_q_pval <- get_metric(heterogeneity, "MR Egger", "Q_pval")

egger_intercept_value <- if ("egger_intercept" %in% names(egger_intercept)) as.numeric(egger_intercept$egger_intercept[1]) else NA_real_
egger_intercept_se <- if ("se" %in% names(egger_intercept)) as.numeric(egger_intercept$se[1]) else NA_real_
egger_intercept_pval <- if ("pval" %in% names(egger_intercept)) as.numeric(egger_intercept$pval[1]) else NA_real_

heterogeneity_significant <- any(c(ivw_q_pval, egger_q_pval) < 0.05, na.rm = TRUE)
directional_pleiotropy_significant <- isTRUE(egger_intercept_pval < 0.05)

sensitivity_summary <- data.table(
  number_of_snps_analysed = usable_instruments,
  ivw_heterogeneity_q = ivw_q,
  ivw_heterogeneity_p_value = ivw_q_pval,
  mr_egger_heterogeneity_q = egger_q,
  mr_egger_heterogeneity_p_value = egger_q_pval,
  mr_egger_intercept = egger_intercept_value,
  mr_egger_intercept_se = egger_intercept_se,
  mr_egger_intercept_p_value = egger_intercept_pval,
  heterogeneity_p_lt_0_05 = heterogeneity_significant,
  directional_pleiotropy_p_lt_0_05 = directional_pleiotropy_significant,
  notes = paste(
    "Sensitivity analyses should be interpreted alongside the main asthma -> BMI MR estimates;",
    "these tests assess evidence of heterogeneity and directional pleiotropy,",
    "but do not prove or disprove causality."
  )
)

fwrite(sensitivity_summary, sensitivity_summary_file, sep = "\t", quote = FALSE, na = "NA")

plot_failures <- character()

message("Creating single-SNP forest plot.")
forest_plots <- mr_forest_plot(as.data.frame(single_snp))
forest_height <- min(max(8, nrow(single_snp) * 0.08), 60)
forest_error <- if (length(forest_plots) > 0) {
  save_plot_pair(forest_plots[[1]], forest_png_file, forest_pdf_file, 9, forest_height)
} else {
  "mr_forest_plot returned no plots"
}

if (!is.na(forest_error)) {
  plot_failures <- c(plot_failures, paste("forest plot:", forest_error))
}

message("Creating leave-one-out plot.")
loo_plots <- mr_leaveoneout_plot(as.data.frame(leave_one_out))
loo_error <- if (length(loo_plots) > 0) {
  save_plot_pair(loo_plots[[1]], loo_png_file, loo_pdf_file, 9, 7)
} else {
  "mr_leaveoneout_plot returned no plots"
}

if (!is.na(loo_error)) {
  plot_failures <- c(plot_failures, paste("leave-one-out plot:", loo_error))
}

message("Creating funnel plot.")
funnel_plots <- mr_funnel_plot(as.data.frame(single_snp))
funnel_error <- if (length(funnel_plots) > 0) {
  save_plot_pair(funnel_plots[[1]], funnel_png_file, funnel_pdf_file, 8, 6)
} else {
  "mr_funnel_plot returned no plots"
}

if (!is.na(funnel_error)) {
  plot_failures <- c(plot_failures, paste("funnel plot:", funnel_error))
}

message("")
message("Sensitivity analysis complete.")
message("SNPs analysed: ", usable_instruments)
message("Written heterogeneity results: ", heterogeneity_file)
message("Written MR-Egger intercept results: ", egger_intercept_file)
message("Written single-SNP results: ", single_snp_file)
message("Written leave-one-out results: ", leave_one_out_file)
message("Written sensitivity summary: ", sensitivity_summary_file)
message("IVW Q = ", round(ivw_q, 4), ", p = ", signif(ivw_q_pval, 4))
message("MR-Egger Q = ", round(egger_q, 4), ", p = ", signif(egger_q_pval, 4))
message(
  "MR-Egger intercept = ", signif(egger_intercept_value, 4),
  ", SE = ", signif(egger_intercept_se, 4),
  ", p = ", signif(egger_intercept_pval, 4)
)

if (heterogeneity_significant) {
  message("There is evidence of heterogeneity at p < 0.05.")
} else {
  message("There is no strong evidence of heterogeneity at p < 0.05.")
}

if (directional_pleiotropy_significant) {
  message("There is evidence of directional pleiotropy at p < 0.05.")
} else {
  message("There is no strong evidence of directional pleiotropy at p < 0.05.")
}

if (length(plot_failures) == 0) {
  message("All requested sensitivity plots were saved.")
} else {
  warning("Some sensitivity plots failed to save: ", paste(plot_failures, collapse = "; "), call. = FALSE)
}

message("Caution: these sensitivity analyses do not prove or disprove causality.")
