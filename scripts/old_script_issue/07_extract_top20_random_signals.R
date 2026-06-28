#!/usr/bin/env Rscript

# ============================================================
# 07_extract_top20_random_signals.R
#
# NOTE:
#   This script is retained for a future/random-effects sensitivity analysis.
#   It should only be run if results/gwama/asthma_meta_random.out exists.
#   The current corrected GWAMA run produces results/gwama/asthma_meta.out.
# Purpose:
#   Extract the 20 most statistically significant SNPs from the
#   random-effects GWAMA meta-analysis. These results are used to
#   summarise the strongest association signals while accounting
#   for between-study heterogeneity.
#
# Inputs:
#   ../results/gwama/asthma_meta_random.out
#
# Outputs:
#   ../results/gwama/genome_wide_significant_random.txt
#   ../results/gwama/top20_random_signals.txt
#
# Output columns:
#   rs_number        = single nucleotide polymorphism (SNP) identifier
#   reference_allele = reference allele
#   other_allele     = alternate/non-reference allele
#   beta             = random-effects beta estimate
#   se               = standard error
#   p-value          = association p-value
#   z                = z statistic
#   q_statistic      = Cochran's Q heterogeneity statistic
#   q_p-value        = heterogeneity p-value
#   i2               = I-squared heterogeneity estimate
#   n_studies        = number of contributing studies
#
# Notes:
#   - Genome-wide significant variants are defined as p-value < 5e-8.
#   - The top 20 table is sorted from most significant to least
#     significant using ascending p-value.
#   - Only the columns listed above are retained for a compact
#     random-effects summary table.
#   - This script reproduces the terminal workflow:
#       1. Filter asthma_meta_random.out to genome-wide significant
#          variants.
#       2. Sort the filtered results by p-value.
#       3. Keep the first 20 rows.
#
# Version notes:
#   - Paths are anchored to this script's location, so it can be
#     run from the project root with:
#       Rscript scripts/07_extract_top20_random_signals.R
#   - The script stops with a clear error if any required column is
#     missing from asthma_meta_random.out.
# ============================================================

library(data.table)

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

gwama_results_dir <- file.path(project_dir, "results", "gwama")
input_file <- file.path(gwama_results_dir, "asthma_meta_random.out")
genome_wide_file <- file.path(gwama_results_dir, "genome_wide_significant_random.txt")
output_file <- file.path(gwama_results_dir, "top20_random_signals.txt")

required_columns <- c(
  "rs_number",
  "reference_allele",
  "other_allele",
  "beta",
  "se",
  "p-value",
  "z",
  "q_statistic",
  "q_p-value",
  "i2",
  "n_studies"
)

message("Reading: ", input_file)
results <- fread(input_file)

missing_columns <- setdiff(required_columns, names(results))
if (length(missing_columns) > 0) {
  stop(
    "Missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

genome_wide_significant <- results[`p-value` < 5e-8]

fwrite(
  genome_wide_significant,
  genome_wide_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message(
  "Written: ",
  genome_wide_file,
  " | rows: ",
  nrow(genome_wide_significant)
)

top20 <- genome_wide_significant[
  order(`p-value`),
  ..required_columns
][1:20]

fwrite(top20, output_file, sep = "\t", quote = FALSE, na = "NA")

message("Written: ", output_file, " | rows: ", nrow(top20))
