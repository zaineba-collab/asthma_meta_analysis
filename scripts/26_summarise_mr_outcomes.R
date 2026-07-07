#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

# ============================================================
# 26_summarise_mr_outcomes.R
#
# Purpose:
#   Combine completed asthma -> outcome MR results into reproducible
#   multi-outcome summary tables and a forest plot.
#
# Supervisor plan:
#   Use already LD-clumped asthma signals. This script only reads existing
#   MR outputs and does not perform clumping, harmonisation, outcome
#   extraction, or MR estimation.
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

config_file <- file.path(project_dir, "config", "mr_outcomes.tsv")
summary_dir <- file.path(project_dir, "results", "mr", "summary")

all_methods_file <- file.path(summary_dir, "asthma_multi_outcome_mr_all_methods.tsv")
ivw_summary_file <- file.path(summary_dir, "asthma_multi_outcome_ivw_summary.tsv")
diagnostics_file <- file.path(summary_dir, "asthma_multi_outcome_diagnostics.tsv")
publication_file <- file.path(summary_dir, "asthma_multi_outcome_publication_table.tsv")
missing_file <- file.path(summary_dir, "asthma_multi_outcome_missing_results.tsv")
forest_png_file <- file.path(summary_dir, "asthma_multi_outcome_ivw_forest_plot.png")
forest_pdf_file <- file.path(summary_dir, "asthma_multi_outcome_ivw_forest_plot.pdf")

required_config_cols <- c("outcome_label", "outcome_name", "outcome_id", "outcome_type")
required_mr_cols <- c(
  "method",
  "number_of_snps",
  "beta_estimate",
  "standard_error",
  "p_value",
  "beta_lci95",
  "beta_uci95"
)
required_sensitivity_cols <- c(
  "number_of_snps_analysed",
  "ivw_heterogeneity_q",
  "ivw_heterogeneity_p_value",
  "mr_egger_heterogeneity_q",
  "mr_egger_heterogeneity_p_value",
  "mr_egger_intercept",
  "mr_egger_intercept_se",
  "mr_egger_intercept_p_value",
  "heterogeneity_p_lt_0_05",
  "directional_pleiotropy_p_lt_0_05"
)

warn_records <- list()

add_warning <- function(label, issue, detail = NA_character_) {
  warn_records[[length(warn_records) + 1L]] <<- data.table(
    outcome_label = label,
    issue = issue,
    detail = detail
  )
  warning(paste(label, issue, detail), call. = FALSE)
}

missing_dt <- function(label, outcome_name, outcome_id, reason, path = NA_character_) {
  data.table(
    outcome_label = label,
    outcome_name = outcome_name,
    outcome_id = outcome_id,
    status = "missing_or_incomplete",
    reason = reason,
    path = path
  )
}

empty_missing_results <- function() {
  data.table(
    outcome_label = character(),
    outcome_name = character(),
    outcome_id = character(),
    status = character(),
    reason = character(),
    path = character()
  )
}

as_scalar <- function(dt, metric_name) {
  if (is.null(dt) || nrow(dt) == 0 || !"metric" %in% names(dt) || !"value" %in% names(dt)) {
    return(NA_character_)
  }

  value <- dt[metric == metric_name, value]

  if (length(value) == 0) {
    return(NA_character_)
  }

  as.character(value[1])
}

as_scalar_pattern <- function(dt, pattern) {
  if (is.null(dt) || nrow(dt) == 0 || !"metric" %in% names(dt) || !"value" %in% names(dt)) {
    return(NA_character_)
  }

  value <- dt[grepl(pattern, metric), value]

  if (length(value) == 0) {
    return(NA_character_)
  }

  as.character(value[1])
}

num_or_na <- function(x) {
  suppressWarnings(as.numeric(x))
}

logical_or_na <- function(x) {
  if (length(x) == 0 || is.na(x)) {
    return(NA)
  }

  if (is.logical(x)) {
    return(x[1])
  }

  tolower(as.character(x[1])) %in% c("true", "t", "1", "yes", "y")
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "fg", flag = "#"))
}

fmt_p <- function(x) {
  vapply(x, function(value) {
    if (is.na(value)) {
      return("NA")
    }

    if (value == 0) {
      return("0 (underflow)")
    }

    if (abs(value) < 0.001) {
      return(formatC(value, digits = 3, format = "e"))
    }

    formatC(value, digits = 3, format = "fg", flag = "#")
  }, character(1))
}

if (!file.exists(config_file)) {
  stop("MR outcome config file not found: ", config_file, call. = FALSE)
}

dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

config <- fread(config_file)
missing_config_cols <- setdiff(required_config_cols, names(config))

if (length(missing_config_cols) > 0) {
  stop("Config missing required column(s): ", paste(missing_config_cols, collapse = ", "), call. = FALSE)
}

duplicate_labels <- config[duplicated(outcome_label) | duplicated(outcome_label, fromLast = TRUE)]

if (nrow(duplicate_labels) > 0) {
  stop("Duplicate outcome labels in config: ", paste(unique(duplicate_labels$outcome_label), collapse = ", "), call. = FALSE)
}

completed <- list()
ivw_rows <- list()
diagnostic_rows <- list()
missing_rows <- list()

for (i in seq_len(nrow(config))) {
  cfg <- config[i]
  label <- cfg$outcome_label
  outcome_name <- cfg$outcome_name
  outcome_id <- cfg$outcome_id
  outcome_type <- cfg$outcome_type

  base_dir <- file.path(project_dir, "results", "mr", paste0("asthma_to_", label))
  mr_file <- file.path(base_dir, "analysis", paste0("asthma_", label, "_mr_results_publication.tsv"))
  sensitivity_file <- file.path(base_dir, "sensitivity", paste0("asthma_", label, "_sensitivity_summary.tsv"))
  harmonisation_file <- file.path(base_dir, "mr_harmonisation_diagnostics.tsv")
  outcome_diag_file <- file.path(base_dir, paste0(label, "_outcome_diagnostics.tsv"))

  if (!file.exists(mr_file)) {
    missing_rows[[length(missing_rows) + 1L]] <- missing_dt(label, outcome_name, outcome_id, "MR publication results file missing", mr_file)
    add_warning(label, "MR publication results file missing", mr_file)
    next
  }

  if (!file.exists(sensitivity_file)) {
    missing_rows[[length(missing_rows) + 1L]] <- missing_dt(label, outcome_name, outcome_id, "Sensitivity summary file missing", sensitivity_file)
    add_warning(label, "Sensitivity summary file missing", sensitivity_file)
    next
  }

  mr <- tryCatch(fread(mr_file), error = function(e) e)
  sensitivity <- tryCatch(fread(sensitivity_file), error = function(e) e)

  if (inherits(mr, "error")) {
    missing_rows[[length(missing_rows) + 1L]] <- missing_dt(label, outcome_name, outcome_id, conditionMessage(mr), mr_file)
    add_warning(label, "Could not read MR publication results file", conditionMessage(mr))
    next
  }

  if (inherits(sensitivity, "error")) {
    missing_rows[[length(missing_rows) + 1L]] <- missing_dt(label, outcome_name, outcome_id, conditionMessage(sensitivity), sensitivity_file)
    add_warning(label, "Could not read sensitivity summary file", conditionMessage(sensitivity))
    next
  }

  missing_mr_cols <- setdiff(required_mr_cols, names(mr))
  missing_sensitivity_cols <- setdiff(required_sensitivity_cols, names(sensitivity))

  if (length(missing_mr_cols) > 0) {
    missing_rows[[length(missing_rows) + 1L]] <- missing_dt(label, outcome_name, outcome_id, paste("MR file missing columns:", paste(missing_mr_cols, collapse = ", ")), mr_file)
    add_warning(label, "MR file missing required columns", paste(missing_mr_cols, collapse = ", "))
    next
  }

  if (length(missing_sensitivity_cols) > 0) {
    missing_rows[[length(missing_rows) + 1L]] <- missing_dt(label, outcome_name, outcome_id, paste("Sensitivity file missing columns:", paste(missing_sensitivity_cols, collapse = ", ")), sensitivity_file)
    add_warning(label, "Sensitivity file missing required columns", paste(missing_sensitivity_cols, collapse = ", "))
    next
  }

  duplicate_methods <- mr[duplicated(method) | duplicated(method, fromLast = TRUE)]

  if (nrow(duplicate_methods) > 0) {
    missing_rows[[length(missing_rows) + 1L]] <- missing_dt(label, outcome_name, outcome_id, paste("Duplicate method rows:", paste(unique(duplicate_methods$method), collapse = ", ")), mr_file)
    add_warning(label, "Duplicate outcome-method rows", paste(unique(duplicate_methods$method), collapse = ", "))
    next
  }

  ivw <- mr[method == "Inverse variance weighted"]

  if (nrow(ivw) != 1) {
    missing_rows[[length(missing_rows) + 1L]] <- missing_dt(label, outcome_name, outcome_id, paste("Expected exactly one IVW row, found", nrow(ivw)), mr_file)
    add_warning(label, "Invalid IVW row count", as.character(nrow(ivw)))
    next
  }

  outcome_diag <- if (file.exists(outcome_diag_file)) fread(outcome_diag_file) else NULL
  harmonisation_diag <- if (file.exists(harmonisation_file)) fread(harmonisation_file) else NULL

  if (is.null(outcome_diag)) {
    add_warning(label, "Outcome diagnostics file missing", outcome_diag_file)
  }

  if (is.null(harmonisation_diag)) {
    add_warning(label, "Harmonisation diagnostics file missing", harmonisation_file)
  }

  diag_outcome_id <- as_scalar(outcome_diag, "outcome_id")
  diag_outcome_name <- as_scalar(outcome_diag, "outcome_name")

  if (!is.na(diag_outcome_id) && !identical(diag_outcome_id, outcome_id)) {
    missing_rows[[length(missing_rows) + 1L]] <- missing_dt(label, outcome_name, outcome_id, paste("Outcome diagnostics ID mismatch:", diag_outcome_id), outcome_diag_file)
    add_warning(label, "Outcome diagnostics ID mismatch", paste("config =", outcome_id, "diagnostics =", diag_outcome_id))
    next
  }

  if (!is.na(diag_outcome_name) && !identical(diag_outcome_name, outcome_name)) {
    add_warning(label, "Outcome diagnostics name differs from config", paste("config =", outcome_name, "diagnostics =", diag_outcome_name))
  }

  completed[[length(completed) + 1L]] <- cbind(
    cfg[, .(outcome_label, outcome_name, outcome_id, outcome_type)],
    mr[
      ,
      .(
        method,
        number_of_snps = num_or_na(number_of_snps),
        beta_estimate = num_or_na(beta_estimate),
        standard_error = num_or_na(standard_error),
        beta_lci95 = num_or_na(beta_lci95),
        beta_uci95 = num_or_na(beta_uci95),
        p_value = num_or_na(p_value)
      )
    ]
  )

  ivw_rows[[length(ivw_rows) + 1L]] <- cbind(
    cfg[, .(outcome_label, outcome_name, outcome_id, outcome_type)],
    ivw[
      ,
      .(
        number_of_snps = num_or_na(number_of_snps),
        beta_estimate = num_or_na(beta_estimate),
        standard_error = num_or_na(standard_error),
        beta_lci95 = num_or_na(beta_lci95),
        beta_uci95 = num_or_na(beta_uci95),
        p_value = num_or_na(p_value)
      )
    ],
    sensitivity[
      1,
      .(
        ivw_heterogeneity_q = num_or_na(ivw_heterogeneity_q),
        ivw_heterogeneity_p_value = num_or_na(ivw_heterogeneity_p_value),
        mr_egger_heterogeneity_q = num_or_na(mr_egger_heterogeneity_q),
        mr_egger_heterogeneity_p_value = num_or_na(mr_egger_heterogeneity_p_value),
        mr_egger_intercept = num_or_na(mr_egger_intercept),
        mr_egger_intercept_se = num_or_na(mr_egger_intercept_se),
        mr_egger_intercept_p_value = num_or_na(mr_egger_intercept_p_value),
        heterogeneity_p_lt_0_05 = logical_or_na(heterogeneity_p_lt_0_05),
        directional_pleiotropy_p_lt_0_05 = logical_or_na(directional_pleiotropy_p_lt_0_05)
      )
    ]
  )

  diagnostic_rows[[length(diagnostic_rows) + 1L]] <- cfg[
    ,
    .(
      outcome_label,
      outcome_name,
      outcome_id,
      outcome_type,
      total_asthma_exposure_snps = num_or_na(as_scalar(outcome_diag, "total_asthma_exposure_snps")),
      outcome_rows_returned = num_or_na(as_scalar_pattern(outcome_diag, "_outcome_rows_returned$")),
      direct_matches = num_or_na(as_scalar(outcome_diag, "direct_matches")),
      proxy_matches = num_or_na(as_scalar(outcome_diag, "proxy_matches")),
      unmatched_snps = num_or_na(as_scalar(outcome_diag, "unmatched_queried_snps")),
      outcome_coverage_percentage = num_or_na(as_scalar(outcome_diag, "percentage_outcome_coverage")),
      harmonised_rows = num_or_na(as_scalar(harmonisation_diag, "snps_retained_after_harmonisation")),
      mr_keep_true = num_or_na(as_scalar(harmonisation_diag, "mr_keep_true_rows")),
      mr_keep_false = num_or_na(as_scalar(harmonisation_diag, "mr_keep_false_rows")),
      final_usable_instruments = num_or_na(as_scalar(harmonisation_diag, "final_usable_mr_instruments")),
      duplicate_exact_outcome_rows = num_or_na(as_scalar(outcome_diag, "duplicate_exact_outcome_rows")),
      duplicate_outcome_snp_rows = num_or_na(as_scalar(outcome_diag, "duplicate_outcome_snp_rows")),
      additional_twosamplemr_clumping_performed = as_scalar(outcome_diag, "additional_twosamplemr_clumping_performed")
    )
  ]
}

all_methods <- rbindlist(completed, fill = TRUE)
ivw_summary <- rbindlist(ivw_rows, fill = TRUE)
diagnostics <- rbindlist(diagnostic_rows, fill = TRUE)
missing_results <- if (length(missing_rows) > 0) {
  rbindlist(missing_rows, fill = TRUE)
} else {
  empty_missing_results()
}
warnings_dt <- rbindlist(warn_records, fill = TRUE)

if (nrow(all_methods) == 0 || nrow(ivw_summary) == 0) {
  fwrite(missing_results, missing_file, sep = "\t", quote = FALSE, na = "NA")
  stop("No completed MR outcomes could be summarised. See: ", missing_file, call. = FALSE)
}

completed_outcomes <- nrow(ivw_summary)
bonferroni_threshold <- 0.05 / completed_outcomes

ivw_summary[
  ,
  `:=`(
    number_of_completed_outcomes = completed_outcomes,
    exploratory_bonferroni_threshold = bonferroni_threshold,
    ivw_p_lt_0_05 = p_value < 0.05,
    ivw_p_lt_bonferroni_threshold = p_value < bonferroni_threshold
  )
]

publication_table <- ivw_summary[
  ,
  .(
    Outcome = outcome_name,
    `N instruments` = number_of_snps,
    `IVW beta` = fmt_num(beta_estimate),
    `IVW SE` = fmt_num(standard_error),
    `95% CI` = paste0(fmt_num(beta_lci95), " to ", fmt_num(beta_uci95)),
    `IVW p-value` = fmt_p(p_value),
    `Heterogeneity p-value` = fmt_p(ivw_heterogeneity_p_value),
    `MR-Egger intercept p-value` = fmt_p(mr_egger_intercept_p_value),
    `Evidence of heterogeneity` = ifelse(heterogeneity_p_lt_0_05 %in% TRUE, "Yes", ifelse(is.na(heterogeneity_p_lt_0_05), "Unknown", "No")),
    `Evidence of directional pleiotropy` = ifelse(directional_pleiotropy_p_lt_0_05 %in% TRUE, "Yes", ifelse(is.na(directional_pleiotropy_p_lt_0_05), "Unknown", "No"))
  )
]

setorder(all_methods, outcome_label, method)
setorder(ivw_summary, p_value)
setorder(diagnostics, outcome_label)
setorder(publication_table, Outcome)

fwrite(all_methods, all_methods_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(ivw_summary, ivw_summary_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(diagnostics, diagnostics_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(publication_table, publication_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(missing_results, missing_file, sep = "\t", quote = FALSE, na = "NA")

forest_data <- copy(ivw_summary)
forest_data[, outcome_name := factor(outcome_name, levels = rev(outcome_name[order(beta_estimate)]))]

forest_plot <- ggplot(forest_data, aes(x = beta_estimate, y = outcome_name)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_errorbar(aes(xmin = beta_lci95, xmax = beta_uci95), orientation = "y", width = 0.18, linewidth = 0.6) +
  geom_point(size = 2.2) +
  labs(
    x = "IVW beta estimate",
    y = NULL,
    title = "Asthma -> outcome MR estimates",
    subtitle = "Inverse variance weighted estimates with 95% confidence intervals"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

plot_height <- max(4.5, 0.45 * nrow(forest_data) + 1.5)
ggsave(forest_png_file, plot = forest_plot, width = 8, height = plot_height, dpi = 300)
ggsave(forest_pdf_file, plot = forest_plot, width = 8, height = plot_height)

included_labels <- ivw_summary$outcome_label
missing_labels <- setdiff(config$outcome_label, included_labels)
ivw_p05 <- ivw_summary[ivw_p_lt_0_05 == TRUE, outcome_label]
ivw_bonf <- ivw_summary[ivw_p_lt_bonferroni_threshold == TRUE, outcome_label]
heterogeneity_labels <- ivw_summary[heterogeneity_p_lt_0_05 == TRUE, outcome_label]
pleiotropy_labels <- ivw_summary[directional_pleiotropy_p_lt_0_05 == TRUE, outcome_label]

message("")
message("Multi-outcome MR summary complete.")
message("Configured outcomes: ", nrow(config))
message("Completed outcomes included: ", completed_outcomes)
message("Outcome labels included: ", paste(included_labels, collapse = ", "))
message("Missing or incomplete outcomes: ", ifelse(length(missing_labels) == 0, "none", paste(missing_labels, collapse = ", ")))
message("Rows in all-methods table: ", nrow(all_methods))
message("Exploratory Bonferroni threshold: ", signif(bonferroni_threshold, 4))
message("Outcomes with IVW p < 0.05: ", ifelse(length(ivw_p05) == 0, "none", paste(ivw_p05, collapse = ", ")))
message("Outcomes passing exploratory Bonferroni threshold: ", ifelse(length(ivw_bonf) == 0, "none", paste(ivw_bonf, collapse = ", ")))
message("Outcomes with evidence of heterogeneity: ", ifelse(length(heterogeneity_labels) == 0, "none", paste(heterogeneity_labels, collapse = ", ")))
message("Outcomes with evidence of directional pleiotropy: ", ifelse(length(pleiotropy_labels) == 0, "none", paste(pleiotropy_labels, collapse = ", ")))
message("Written all-methods table: ", all_methods_file)
message("Written IVW summary: ", ivw_summary_file)
message("Written diagnostics: ", diagnostics_file)
message("Written publication table: ", publication_file)
message("Written missing-results diagnostics: ", missing_file)
message("Written forest plot PNG: ", forest_png_file)
message("Written forest plot PDF: ", forest_pdf_file)

if (nrow(warnings_dt) > 0) {
  message("Warnings were recorded during summarisation:")
  print(warnings_dt)
}

message("Caution: these summaries do not prove causality and should be interpreted alongside heterogeneity and pleiotropy results.")
