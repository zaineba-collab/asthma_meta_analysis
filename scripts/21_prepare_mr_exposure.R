#!/usr/bin/env Rscript

library(data.table)

# ============================================================
# 21_prepare_mr_exposure.R
#
# Primary MR direction:
#   Asthma -> BMI
#
# Supervisor plan:
#   Use already LD-clumped asthma signals as MR instruments and do
#   not clump again using TwoSampleMR.
#
# Purpose:
#   Build the asthma exposure dataset for TwoSampleMR directly from
#   the corrected rsID-based PLINK LD-clumped lead SNPs.
#
# Inputs:
#   results/ld_clumping/independent_lead_snps.tsv
#   results/gwama/asthma_fixed_for_plink_clumping.txt
#   results/gwama/asthma_meta.out
#
# Output:
#   results/mr/asthma_to_bmi/asthma_exposure_twosamplemr.tsv
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

lead_file <- file.path(project_dir, "results", "ld_clumping", "independent_lead_snps.tsv")
clump_input_file <- file.path(project_dir, "results", "gwama", "asthma_fixed_for_plink_clumping.txt")
gwama_file <- file.path(project_dir, "results", "gwama", "asthma_meta.out")
output_dir <- file.path(project_dir, "results", "mr", "asthma_to_bmi")
output_file <- file.path(output_dir, "asthma_exposure_twosamplemr.tsv")
diagnostics_file <- file.path(output_dir, "asthma_exposure_diagnostics.tsv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

assert_no_extra_clumping <- function() {
  this_script <- normalizePath(sys.frame(1)$ofile %||% script_file, mustWork = FALSE)

  if (!file.exists(this_script)) {
    return(invisible(TRUE))
  }

  txt <- readLines(this_script, warn = FALSE)
  forbidden <- c(
    paste0("clump", "_data("),
    paste0("extract", "_instruments("),
    paste0("ld", "_clump("),
    paste0("ieugwasr", "::")
  )

  hits <- unlist(lapply(forbidden, function(pattern) grep(pattern, txt, fixed = TRUE, value = TRUE)))

  if (length(hits) > 0) {
    stop(
      "Forbidden clumping/OpenGWAS instrument-selection call found in script 21. ",
      "Supervisor plan requires already PLINK-clumped asthma signals only.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

`%||%` <- function(x, y) if (is.null(x) || is.na(x) || identical(x, "")) y else x

assert_no_extra_clumping()

stop_if_missing(lead_file, "Independent lead SNP table")
stop_if_missing(clump_input_file, "rsID PLINK clumping input")
stop_if_missing(gwama_file, "GWAMA meta-analysis output")

message("Primary MR direction: Asthma -> BMI")
message("Using already LD-clumped asthma signals. No additional TwoSampleMR clumping performed.")
message("Reading corrected PLINK LD-clumped asthma lead SNPs: ", lead_file)

lead <- fread(lead_file)
clump_map <- fread(clump_input_file)

required_lead_cols <- c("SNP", "p_value")
required_map_cols <- c("SNP", "original_SNP")

missing_lead_cols <- setdiff(required_lead_cols, names(lead))
missing_map_cols <- setdiff(required_map_cols, names(clump_map))

if (length(missing_lead_cols) > 0) {
  stop("Lead SNP table missing column(s): ", paste(missing_lead_cols, collapse = ", "), call. = FALSE)
}

if (length(missing_map_cols) > 0) {
  stop("Clumping input missing column(s): ", paste(missing_map_cols, collapse = ", "), call. = FALSE)
}

duplicate_leads <- lead[duplicated(SNP) | duplicated(SNP, fromLast = TRUE)]

if (nrow(duplicate_leads) > 0) {
  stop("Duplicate lead SNPs found: ", paste(sort(unique(duplicate_leads$SNP)), collapse = ", "), call. = FALSE)
}

lead_map <- merge(
  lead[, .(SNP, lead_p_value = p_value)],
  clump_map[, .(SNP, original_snp.exposure = original_SNP)],
  by = "SNP",
  all.x = TRUE,
  sort = FALSE
)

missing_mapping <- lead_map[is.na(original_snp.exposure)]

if (nrow(missing_mapping) > 0) {
  stop("Missing original GWAMA ID mapping for ", nrow(missing_mapping), " lead SNP(s).", call. = FALSE)
}

message("Reading GWAMA asthma exposure columns: ", gwama_file)

gwama <- fread(
  gwama_file,
  select = c(
    "rs_number",
    "reference_allele",
    "other_allele",
    "eaf",
    "beta",
    "se",
    "p-value",
    "n_samples"
  )
)

gwama <- gwama[
  rs_number %in% lead_map$original_snp.exposure,
  .(
    original_snp.exposure = rs_number,
    beta.exposure = beta,
    se.exposure = se,
    effect_allele.exposure = other_allele,
    other_allele.exposure = reference_allele,
    eaf.exposure = eaf,
    pval.exposure = `p-value`,
    samplesize.exposure = n_samples,
    exposure = "Asthma",
    id.exposure = "asthma_gwama"
  )
]

exposure <- merge(
  lead_map,
  gwama,
  by = "original_snp.exposure",
  all.x = TRUE,
  sort = FALSE
)

missing_gwama <- exposure[is.na(beta.exposure) | is.na(se.exposure)]

if (nrow(missing_gwama) > 0) {
  stop("Missing GWAMA effect estimates for ", nrow(missing_gwama), " lead SNP(s).", call. = FALSE)
}

exposure[eaf.exposure == -9, eaf.exposure := NA]
exposure[samplesize.exposure == -9, samplesize.exposure := NA]
exposure[, lead_p_value := NULL]

setcolorder(
  exposure,
  c(
    "SNP",
    "original_snp.exposure",
    "beta.exposure",
    "se.exposure",
    "effect_allele.exposure",
    "other_allele.exposure",
    "eaf.exposure",
    "pval.exposure",
    "samplesize.exposure",
    "exposure",
    "id.exposure"
  )
)

setorder(exposure, pval.exposure)

diagnostics <- data.table(
  metric = c(
    "primary_mr_direction",
    "supervisor_plan",
    "ld_clumped_asthma_snps_read",
    "asthma_exposure_snps_created",
    "additional_twosamplemr_clumping_performed",
    "source_ld_clumped_file"
  ),
  value = c(
    "Asthma -> BMI",
    "Use already LD-clumped asthma signals; do not clump again using TwoSampleMR",
    nrow(lead),
    nrow(exposure),
    "No",
    lead_file
  )
)

fwrite(exposure, output_file, sep = "\t", na = "NA", quote = FALSE)
fwrite(diagnostics, diagnostics_file, sep = "\t", na = "NA", quote = FALSE)

message("Written: ", output_file)
message("Written diagnostics: ", diagnostics_file)
message("Corrected LD-clumped asthma SNPs read: ", nrow(lead))
message("Asthma exposure SNPs created: ", nrow(exposure))
message("Confirmation: no additional TwoSampleMR clumping performed.")
