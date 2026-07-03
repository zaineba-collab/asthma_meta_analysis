#!/usr/bin/env Rscript

library(data.table)
library(TwoSampleMR)

# ============================================================
# 22_prepare_mr_outcomes.R
#
# Purpose:
#   Prepare BMI outcome association data for Mendelian randomisation
#   using the independent asthma lead SNPs as instruments.
#
# Input:
#   results/mr/asthma_exposure_twosamplemr.tsv
#
# Output:
#   results/mr/bmi_outcome_twosamplemr.tsv
#
# Notes:
#   - This script uses TwoSampleMR::extract_outcome_data() to query
#     the IEU OpenGWAS database.
#   - The default BMI outcome is ieu-a-2, a commonly used BMI GWAS
#     dataset in OpenGWAS/TwoSampleMR examples.
#   - OpenGWAS may require authentication. If so, set OPENGWAS_JWT
#     before running this script.
#   - OpenGWAS usually expects rsID-style SNP identifiers. Your current
#     asthma instruments are chromosome:position:ref:alt IDs, so if the
#     query returns no matches, the next step is to map those IDs to rsIDs.
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

mr_dir <- file.path(project_dir, "results", "mr")
exposure_file <- file.path(mr_dir, "asthma_exposure_twosamplemr.tsv")
outcome_file <- file.path(mr_dir, "bmi_outcome_twosamplemr.tsv")

# OpenGWAS outcome ID for BMI. Change this if you decide to use a different
# BMI/adiposity GWAS, such as body fat percentage or waist-hip ratio.
bmi_outcome_id <- "ieu-a-835"
bmi_outcome_name <- "BMI"

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

stop_if_missing(exposure_file, "Asthma exposure file")
dir.create(mr_dir, showWarnings = FALSE, recursive = TRUE)

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

if (any(!grepl("^rs[0-9]+$", snps))) {
  warning(
    "Some exposure SNP IDs are not rsIDs. OpenGWAS may not find these variants. ",
    "If no outcome rows are returned, map chromosome:position:ref:alt IDs to rsIDs first.",
    call. = FALSE
  )
}

if (identical(Sys.getenv("OPENGWAS_JWT"), "")) {
  stop(
    "OPENGWAS_JWT is not set. OpenGWAS now requires a JWT token for this query. ",
    "Create a token at https://api.opengwas.io/ and set OPENGWAS_JWT before rerunning.",
    call. = FALSE
  )
}

message("Extracting BMI outcome data from OpenGWAS outcome: ", bmi_outcome_id)

outcome <- extract_outcome_data(
  snps = snps,
  outcomes = bmi_outcome_id,
  proxies = TRUE,
  rsq = 0.8,
  align_alleles = 1,
  palindromes = 1,
  maf_threshold = 0.3
)

if (is.null(outcome) || nrow(outcome) == 0) {
  stop(
    "No BMI outcome rows were returned from OpenGWAS. ",
    "This is likely because the exposure SNP IDs are not rsIDs.",
    call. = FALSE
  )
}

setDT(outcome)

# Give the outcome a simple human-readable name for downstream tables.
outcome[, outcome := bmi_outcome_name]

fwrite(outcome, outcome_file, sep = "\t", quote = FALSE, na = "NA")

message("Written: ", outcome_file)
message("Outcome SNPs returned: ", nrow(outcome))
