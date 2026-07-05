#!/usr/bin/env Rscript

library(data.table)
library(TwoSampleMR)

# ============================================================
# 22_prepare_mr_outcomes.R
#
# Primary MR direction:
#   Asthma -> BMI
#
# Supervisor plan:
#   Use already LD-clumped asthma signals as MR instruments and do
#   not clump again using TwoSampleMR.
#
# Purpose:
#   Extract BMI outcome associations from OpenGWAS for the already
#   PLINK LD-clumped asthma exposure SNPs.
#
# Inputs:
#   results/mr/asthma_to_bmi/asthma_exposure_twosamplemr.tsv
#
# Outputs:
#   results/mr/asthma_to_bmi/bmi_outcome_twosamplemr.tsv
#   results/mr/asthma_to_bmi/bmi_outcome_diagnostics.tsv
#   results/mr/asthma_to_bmi/bmi_outcome_proxy_summary.tsv
#   results/mr/asthma_to_bmi/bmi_outcome_unmatched_snps.txt
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

primary_dir <- file.path(project_dir, "results", "mr", "asthma_to_bmi")
exposure_file <- file.path(primary_dir, "asthma_exposure_twosamplemr.tsv")
outcome_file <- file.path(primary_dir, "bmi_outcome_twosamplemr.tsv")
diagnostics_file <- file.path(primary_dir, "bmi_outcome_diagnostics.tsv")
proxy_summary_file <- file.path(primary_dir, "bmi_outcome_proxy_summary.tsv")
unmatched_file <- file.path(primary_dir, "bmi_outcome_unmatched_snps.txt")
non_rsid_file <- file.path(primary_dir, "bmi_outcome_non_rsid_snps.tsv")
partial_outcome_file <- file.path(primary_dir, "bmi_outcome_twosamplemr_partial.tsv")

bmi_outcome_id <- "ieu-a-835"
bmi_outcome_name <- "BMI"

batch_size <- as.integer(Sys.getenv("OPENGWAS_BATCH_SIZE", unset = "25"))
max_attempts <- as.integer(Sys.getenv("OPENGWAS_MAX_ATTEMPTS", unset = "3"))
use_proxies <- tolower(Sys.getenv("OPENGWAS_PROXIES", unset = "true")) %in%
  c("true", "t", "1", "yes", "y")
force_opengwas <- tolower(Sys.getenv("FORCE_OPENGWAS_EXTRACT", unset = "false")) %in%
  c("true", "t", "1", "yes", "y")

if (is.na(batch_size) || batch_size < 1) {
  stop("OPENGWAS_BATCH_SIZE must be a positive integer.", call. = FALSE)
}

if (is.na(max_attempts) || max_attempts < 1) {
  stop("OPENGWAS_MAX_ATTEMPTS must be a positive integer.", call. = FALSE)
}

batch_dir <- file.path(primary_dir, "opengwas_bmi_batches", paste0("batch_size_", batch_size))

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

split_into_batches <- function(x, n) {
  split(x, ceiling(seq_along(x) / n))
}

batch_file <- function(batch_number) {
  file.path(batch_dir, sprintf("bmi_%s_batch_%03d.rds", bmi_outcome_id, batch_number))
}

combine_saved_batches <- function(files) {
  files <- files[file.exists(files)]

  if (length(files) == 0) {
    return(data.table())
  }

  rbindlist(lapply(files, function(path) {
    x <- readRDS(path)
    setDT(x)
    x
  }), fill = TRUE)
}

query_batch <- function(batch_snps, proxies_enabled) {
  result <- extract_outcome_data(
    snps = batch_snps,
    outcomes = bmi_outcome_id,
    proxies = proxies_enabled,
    rsq = 0.8,
    align_alleles = 1,
    palindromes = 1,
    maf_threshold = 0.3
  )

  if (is.null(result)) {
    return(data.table())
  }

  setDT(result)
  result
}

write_diagnostics <- function(outcome, snps, is_rsid) {
  setDT(outcome)
  outcome[, outcome := bmi_outcome_name]

  matched_snps <- unique(outcome$SNP)
  matched_snps <- matched_snps[!is.na(matched_snps) & matched_snps != ""]
  snps_to_query <- snps[is_rsid]
  unmatched_snps <- setdiff(snps_to_query, matched_snps)

  direct_match_count <- nrow(outcome[is.na(proxy.outcome) | proxy.outcome == FALSE])
  proxy_match_count <- nrow(outcome[proxy.outcome == TRUE])
  duplicate_exact_rows <- sum(duplicated(outcome))
  duplicate_snp_rows <- nrow(outcome) - uniqueN(outcome$SNP)
  outcome_coverage_percent <- length(matched_snps) / length(snps_to_query) * 100

  diagnostics <- data.table(
    metric = c(
      "primary_mr_direction",
      "supervisor_plan",
      "total_asthma_exposure_snps",
      "non_rsid_snps_skipped",
      "rsid_snps_queried",
      "bmi_outcome_rows_returned",
      "unique_bmi_outcome_snps_returned",
      "direct_matches",
      "proxy_matches",
      "unmatched_queried_snps",
      "percentage_outcome_coverage",
      "duplicate_exact_outcome_rows",
      "duplicate_outcome_snp_rows",
      "additional_twosamplemr_clumping_performed"
    ),
    value = c(
      "Asthma -> BMI",
      "Use already LD-clumped asthma signals; do not clump again using TwoSampleMR",
      length(snps),
      sum(!is_rsid),
      length(snps_to_query),
      nrow(outcome),
      length(matched_snps),
      direct_match_count,
      proxy_match_count,
      length(unmatched_snps),
      round(outcome_coverage_percent, 3),
      duplicate_exact_rows,
      duplicate_snp_rows,
      "No"
    )
  )

  proxy_columns <- intersect(
    c(
      "SNP",
      "target_snp.outcome",
      "proxy_snp.outcome",
      "proxy.outcome",
      "target_a1.outcome",
      "target_a2.outcome",
      "proxy_a1.outcome",
      "proxy_a2.outcome",
      "effect_allele.outcome",
      "other_allele.outcome",
      "eaf.outcome",
      "beta.outcome",
      "se.outcome",
      "pval.outcome"
    ),
    names(outcome)
  )

  proxy_summary <- outcome[proxy.outcome == TRUE, ..proxy_columns]

  fwrite(outcome, outcome_file, sep = "\t", quote = FALSE, na = "NA")
  fwrite(diagnostics, diagnostics_file, sep = "\t", quote = FALSE, na = "NA")
  fwrite(proxy_summary, proxy_summary_file, sep = "\t", quote = FALSE, na = "NA")
  fwrite(data.table(SNP = unmatched_snps), unmatched_file, col.names = FALSE)

  message("Written BMI outcome file: ", outcome_file)
  message("Written BMI outcome diagnostics: ", diagnostics_file)
  message("Written BMI proxy summary: ", proxy_summary_file)
  message("Written unmatched SNPs: ", unmatched_file)
  message("Direct BMI matches: ", direct_match_count)
  message("Proxy BMI matches: ", proxy_match_count)
  message("Unmatched asthma SNPs: ", length(unmatched_snps))
  message("Outcome coverage: ", round(outcome_coverage_percent, 3), "%")
}

stop_if_missing(exposure_file, "Asthma exposure file")
dir.create(primary_dir, showWarnings = FALSE, recursive = TRUE)

message("Primary MR direction: Asthma -> BMI")
message("Using already LD-clumped asthma signals. No additional TwoSampleMR clumping performed.")
message("Reading asthma exposure instruments: ", exposure_file)

exposure <- fread(exposure_file)

if (!"SNP" %in% names(exposure)) {
  stop("Exposure file is missing required SNP column.", call. = FALSE)
}

snps <- unique(exposure$SNP)
snps <- snps[!is.na(snps) & snps != ""]

if (length(snps) == 0) {
  stop("No SNPs found in exposure file.", call. = FALSE)
}

is_rsid <- grepl("^rs[0-9]+$", snps)

if (any(!is_rsid)) {
  fwrite(data.table(SNP = snps[!is_rsid]), non_rsid_file, sep = "\t", quote = FALSE, na = "NA")
  warning("Some exposure SNP IDs are not rsIDs and will be skipped. Written: ", non_rsid_file, call. = FALSE)
}

snps_to_query <- snps[is_rsid]

if (length(snps_to_query) == 0) {
  stop("No rsID-style SNPs found in the exposure file.", call. = FALSE)
}

if (file.exists(outcome_file) && !force_opengwas) {
  message("Existing BMI outcome file found; reusing it. Set FORCE_OPENGWAS_EXTRACT=true to re-query OpenGWAS.")
  outcome <- fread(outcome_file)
  write_diagnostics(outcome, snps, is_rsid)
  quit(save = "no", status = 0)
}

dir.create(batch_dir, showWarnings = FALSE, recursive = TRUE)

message("Extracting BMI outcome data from OpenGWAS outcome: ", bmi_outcome_id)
message("SNPs to query: ", length(snps_to_query))
message("Batch size: ", batch_size)
message("Batch cache: ", batch_dir)
message("Proxy lookup: ", ifelse(use_proxies, "enabled", "disabled"))

batches <- split_into_batches(snps_to_query, batch_size)
batch_files <- vapply(seq_along(batches), batch_file, character(1))
missing_batches <- !file.exists(batch_files)

if (any(missing_batches) && identical(Sys.getenv("OPENGWAS_JWT"), "")) {
  stop(
    "OPENGWAS_JWT is not set and some OpenGWAS batches still need to be queried. ",
    "Set OPENGWAS_JWT or reuse the existing BMI outcome file.",
    call. = FALSE
  )
}

for (i in seq_along(batches)) {
  current_file <- batch_files[[i]]

  if (file.exists(current_file)) {
    message("Skipping saved batch ", i, " of ", length(batches), ": ", current_file)
    next
  }

  message("Querying batch ", i, " of ", length(batches), " (", length(batches[[i]]), " SNPs)")
  batch_result <- NULL
  batch_error <- NULL

  for (attempt in seq_len(max_attempts)) {
    batch_result <- tryCatch(
      query_batch(batches[[i]], proxies_enabled = use_proxies),
      error = function(e) {
        batch_error <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(batch_result)) {
      break
    }

    message("Batch ", i, " failed on attempt ", attempt, " of ", max_attempts, ": ", batch_error)
    Sys.sleep(10)
  }

  if (is.null(batch_result) && use_proxies) {
    message("Retrying batch ", i, " without proxies.")

    for (attempt in seq_len(max_attempts)) {
      batch_result <- tryCatch(
        query_batch(batches[[i]], proxies_enabled = FALSE),
        error = function(e) {
          batch_error <<- conditionMessage(e)
          NULL
        }
      )

      if (!is.null(batch_result)) {
        break
      }

      message("Batch ", i, " without proxies failed on attempt ", attempt, " of ", max_attempts, ": ", batch_error)
      Sys.sleep(10)
    }
  }

  if (is.null(batch_result)) {
    partial <- combine_saved_batches(batch_files)

    if (nrow(partial) > 0) {
      fwrite(partial, partial_outcome_file, sep = "\t", quote = FALSE, na = "NA")
      message("Partial progress written: ", partial_outcome_file)
    }

    stop("OpenGWAS batch ", i, " failed after retries. Last error: ", batch_error, call. = FALSE)
  }

  saveRDS(batch_result, current_file)
  message("Saved batch ", i, ": ", current_file, " (", nrow(batch_result), " rows)")
  Sys.sleep(1)
}

outcome <- combine_saved_batches(batch_files)

if (nrow(outcome) == 0) {
  stop("No BMI outcome rows were returned from OpenGWAS.", call. = FALSE)
}

write_diagnostics(outcome, snps, is_rsid)
