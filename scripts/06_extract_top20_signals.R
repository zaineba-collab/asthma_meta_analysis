#!/usr/bin/env Rscript

# ============================================================
# 06_extract_top20_fixed_signals.R
#
# Purpose:
#   Extract genome-wide significant SNPs and the 20 most
#   statistically significant SNPs from the current GWAMA
#   meta-analysis output. These results are used for quality
#   control and to summarise the strongest association signals
#   before downstream pruning and locus definition.
#
# Inputs:
#   ../results/gwama/asthma_meta.out
#
# Outputs:
#   ../results/gwama/genome_wide_significant.txt
#   ../results/gwama/top20_signals.txt
#
# Output columns:
#   rs_number        = single nucleotide polymorphism (SNP) identifier
#   reference_allele = reference allele
#   other_allele     = alternate/non-reference allele
#   beta             = fixed-effect beta estimate
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
#   - Top 20 variants are selected from the genome-wide significant
#     variants and sorted by ascending p-value.
#   - Only the columns listed above are retained for a compact
#     fixed-effect summary table.
#
# Version notes:
#   - Paths are anchored to this script's location, so it can be
#     run from the project root with:
#       Rscript scripts/06_extract_top20_signals.R
#   - The script stops with a clear error if any required column is
#     missing from asthma_meta.out.
# ============================================================

library(data.table)

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

gwama_results_dir <- file.path(project_dir, "results", "gwama")
input_file <- file.path(gwama_results_dir, "asthma_meta.out")
genome_wide_file <- file.path(gwama_results_dir, "genome_wide_significant.txt")
output_file <- file.path(gwama_results_dir, "top20_signals.txt")

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
