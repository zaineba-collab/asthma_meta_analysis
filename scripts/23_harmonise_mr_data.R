#!/usr/bin/env Rscript

library(data.table)
library(TwoSampleMR)

# ============================================================
# 23_harmonise_mr_data.R
#
# Primary MR direction:
#   Asthma -> configured outcome
#
# Supervisor plan:
#   Use already LD-clumped asthma signals as MR instruments and do
#   not clump again using TwoSampleMR.
#
# Purpose:
#   Harmonise asthma exposure instruments with configured outcome associations.
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)
source(file.path(script_dir, "mr_config_helpers.R"))

mr_config <- read_mr_outcome_config(project_dir)
outcome_label <- mr_config$outcome_label
outcome_name <- mr_config$outcome_name

primary_dir <- file.path(project_dir, "results", "mr", paste0("asthma_to_", outcome_label))
harmonised_dir <- file.path(primary_dir, "harmonised")

exposure_file <- file.path(primary_dir, "asthma_exposure_twosamplemr.tsv")
outcome_file <- file.path(primary_dir, paste0(outcome_label, "_outcome_twosamplemr.tsv"))
harmonised_file <- file.path(harmonised_dir, paste0("asthma_", outcome_label, "_harmonised.tsv"))
diagnostics_file <- file.path(primary_dir, "mr_harmonisation_diagnostics.tsv")

minimum_usable_instruments <- as.integer(Sys.getenv("MR_MIN_USABLE_INSTRUMENTS", unset = "10"))

if (is.na(minimum_usable_instruments) || minimum_usable_instruments < 1) {
  stop("MR_MIN_USABLE_INSTRUMENTS must be a positive integer.", call. = FALSE)
}

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

count_true <- function(x) sum(x %in% TRUE, na.rm = TRUE)
count_false <- function(x) sum(x %in% FALSE, na.rm = TRUE)

stop_if_missing(exposure_file, "Asthma exposure file")
stop_if_missing(outcome_file, paste(outcome_name, "outcome file"))
dir.create(harmonised_dir, showWarnings = FALSE, recursive = TRUE)

message("Primary MR direction: Asthma -> ", outcome_name)
print_mr_config(mr_config, primary_dir)
message("Reading asthma exposure data: ", exposure_file)
exposure <- fread(exposure_file)

message("Reading ", outcome_name, " outcome data: ", outcome_file)
outcome <- fread(outcome_file)

required_exposure_cols <- c(
  "SNP",
  "beta.exposure",
  "se.exposure",
  "effect_allele.exposure",
  "other_allele.exposure",
  "exposure",
  "id.exposure"
)

required_outcome_cols <- c(
  "SNP",
  "beta.outcome",
  "se.outcome",
  "effect_allele.outcome",
  "other_allele.outcome",
  "outcome",
  "id.outcome"
)

missing_exposure <- setdiff(required_exposure_cols, names(exposure))
missing_outcome <- setdiff(required_outcome_cols, names(outcome))

if (length(missing_exposure) > 0) {
  stop("Exposure file missing column(s): ", paste(missing_exposure, collapse = ", "), call. = FALSE)
}

if (length(missing_outcome) > 0) {
  stop("Outcome file missing column(s): ", paste(missing_outcome, collapse = ", "), call. = FALSE)
}

duplicate_exposure <- exposure[duplicated(SNP) | duplicated(SNP, fromLast = TRUE)]
duplicate_outcome <- outcome[duplicated(SNP) | duplicated(SNP, fromLast = TRUE)]

if (nrow(duplicate_exposure) > 0) {
  stop("Duplicate asthma exposure SNP rows found: ", paste(sort(unique(duplicate_exposure$SNP)), collapse = ", "), call. = FALSE)
}

if (nrow(duplicate_outcome) > 0) {
  stop("Duplicate ", outcome_name, " outcome SNP rows found: ", paste(sort(unique(duplicate_outcome$SNP)), collapse = ", "), call. = FALSE)
}

exposure_snps_before <- uniqueN(exposure$SNP)
outcome_rows_before <- nrow(outcome)
outcome_snps_before <- uniqueN(outcome$SNP)

message("Exposure SNPs before harmonisation: ", exposure_snps_before)
message(outcome_name, " outcome rows before harmonisation: ", outcome_rows_before)
message("Unique ", outcome_name, " outcome SNPs before harmonisation: ", outcome_snps_before)

harmonised <- harmonise_data(
  exposure_dat = as.data.frame(exposure),
  outcome_dat = as.data.frame(outcome),
  action = 2
)

setDT(harmonised)

snps_retained_after_harmonisation <- uniqueN(harmonised$SNP)
exposure_snps_without_outcome <- exposure_snps_before - outcome_snps_before
outcome_rows_not_harmonised <- outcome_rows_before - nrow(harmonised)

mr_keep_true <- if ("mr_keep" %in% names(harmonised)) count_true(harmonised$mr_keep) else NA_integer_
mr_keep_false <- if ("mr_keep" %in% names(harmonised)) count_false(harmonised$mr_keep) else NA_integer_

palindromic_removed <- if (all(c("palindromic", "mr_keep") %in% names(harmonised))) {
  count_true(harmonised$palindromic %in% TRUE & harmonised$mr_keep %in% FALSE)
} else {
  NA_integer_
}

ambiguous_removed <- if (all(c("ambiguous", "mr_keep") %in% names(harmonised))) {
  count_true(harmonised$ambiguous %in% TRUE & harmonised$mr_keep %in% FALSE)
} else {
  NA_integer_
}

final_usable_instruments <- if ("mr_keep" %in% names(harmonised)) {
  uniqueN(harmonised[mr_keep == TRUE, SNP])
} else {
  snps_retained_after_harmonisation
}

removed_during_harmonisation <- snps_retained_after_harmonisation - final_usable_instruments

diagnostics <- data.table(
  metric = c(
    "primary_mr_direction",
    "supervisor_plan",
    "exposure_snps_before_harmonisation",
    paste0(outcome_label, "_outcome_rows_before_harmonisation"),
    paste0("unique_", outcome_label, "_outcome_snps_before_harmonisation"),
    paste0("exposure_snps_without_", outcome_label, "_outcome_before_harmonisation"),
    "outcome_rows_not_harmonised",
    "snps_retained_after_harmonisation",
    "snps_removed_during_harmonisation",
    "mr_keep_true_rows",
    "mr_keep_false_rows",
    "palindromic_snps_removed",
    "ambiguous_snps_removed",
    "final_usable_mr_instruments",
    "additional_twosamplemr_clumping_performed"
  ),
  value = c(
    paste("Asthma ->", outcome_name),
    "Use already LD-clumped asthma signals; do not clump again using TwoSampleMR",
    exposure_snps_before,
    outcome_rows_before,
    outcome_snps_before,
    exposure_snps_without_outcome,
    outcome_rows_not_harmonised,
    snps_retained_after_harmonisation,
    removed_during_harmonisation,
    mr_keep_true,
    mr_keep_false,
    palindromic_removed,
    ambiguous_removed,
    final_usable_instruments,
    "No"
  )
)

fwrite(harmonised, harmonised_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(diagnostics, diagnostics_file, sep = "\t", quote = FALSE, na = "NA")

message("Written harmonised data: ", harmonised_file)
message("Written harmonisation diagnostics: ", diagnostics_file)
message("SNPs retained after harmonisation: ", snps_retained_after_harmonisation)
message("Rows marked mr_keep = TRUE: ", mr_keep_true)
message("Rows marked mr_keep = FALSE: ", mr_keep_false)
message("Final usable MR instruments: ", final_usable_instruments)

if (final_usable_instruments < minimum_usable_instruments) {
  stop(
    "Only ", final_usable_instruments, " usable MR instruments remain after harmonisation. ",
    "This is below MR_MIN_USABLE_INSTRUMENTS = ", minimum_usable_instruments, ".",
    call. = FALSE
  )
}

message("Harmonisation diagnostic step complete. Final MR interpretation has not been run.")
