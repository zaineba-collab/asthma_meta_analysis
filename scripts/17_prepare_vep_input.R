#!/usr/bin/env Rscript

library(data.table)

# ============================================================
# 17_prepare_vep_input.R
#
# Purpose:
#   Prepare independent lead SNPs from PLINK clumping for Ensembl
#   Variant Effect Predictor (VEP) annotation.
#
# Input:
#   results/ld_clumping/independent_lead_snps.tsv
#
# Outputs:
#   results/annotation/independent_lead_snps_vep_input.txt
#   results/annotation/top20_independent_lead_snps_vep_input.txt
#   results/annotation/independent_lead_snps_ids.txt
#
# Notes:
#   - The current clumping workflow uses variant IDs in this format:
#       chromosome:position:reference_allele:alternate_allele
#     for example:
#       11:867462:C:T
#   - VEP accepts a six-column tab-delimited input:
#       chromosome, start, end, allele, strand, identifier
#   - The allele column is written as REF/ALT, e.g. C/T.
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

lead_snps_file <- file.path(
  project_dir,
  "results",
  "ld_clumping",
  "independent_lead_snps.tsv"
)

annotation_dir <- file.path(project_dir, "results", "annotation")
dir.create(annotation_dir, showWarnings = FALSE, recursive = TRUE)

vep_file <- file.path(annotation_dir, "independent_lead_snps_vep_input.txt")
top20_vep_file <- file.path(annotation_dir, "top20_independent_lead_snps_vep_input.txt")
id_file <- file.path(annotation_dir, "independent_lead_snps_ids.txt")

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

stop_if_missing(lead_snps_file, "Independent lead SNP table")

message("Reading independent lead SNPs: ", lead_snps_file)
lead <- fread(lead_snps_file)

required_cols <- c("SNP", "chromosome", "position", "p_value")
missing_cols <- setdiff(required_cols, names(lead))
if (length(missing_cols) > 0) {
  stop(
    "Lead SNP table is missing required column(s): ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

lead[, chromosome := as.character(chromosome)]
lead[, position := as.integer(position)]
lead[, p_value := as.numeric(p_value)]

lead <- lead[
  chromosome %in% as.character(1:22) &
    !is.na(position) &
    !is.na(p_value)
]

if (nrow(lead) == 0) {
  stop("No autosomal lead SNPs with valid position and p-value were found.", call. = FALSE)
}

# Parse the current SNP identifier so VEP receives the allele change directly.
# This avoids relying on a separate PLINK BIM file and makes the step easier to
# rerun on a different machine.
snp_parts <- tstrsplit(lead$SNP, ":", fixed = TRUE)
if (length(snp_parts) != 4) {
  stop(
    "SNP IDs must use chromosome:position:reference_allele:alternate_allele format.",
    call. = FALSE
  )
}

lead[, parsed_chromosome := snp_parts[[1]]]
lead[, parsed_position := as.integer(snp_parts[[2]])]
lead[, reference_allele := toupper(snp_parts[[3]])]
lead[, alternate_allele := toupper(snp_parts[[4]])]

invalid_ids <- lead[
  parsed_chromosome != chromosome |
    parsed_position != position |
    is.na(parsed_position) |
    !grepl("^[ACGT]+$", reference_allele) |
    !grepl("^[ACGT]+$", alternate_allele)
]

if (nrow(invalid_ids) > 0) {
  stop(
    "Some SNP IDs could not be safely parsed for VEP. First invalid ID: ",
    invalid_ids$SNP[1],
    call. = FALSE
  )
}

lead[, allele := paste0(reference_allele, "/", alternate_allele)]
setorder(lead, p_value)

vep_input <- lead[, .(
  chromosome,
  start = position,
  end = position,
  allele,
  strand = "+",
  identifier = SNP
)]

lead_ids <- lead[, .(SNP)]

fwrite(vep_input, vep_file, sep = "\t", quote = FALSE, col.names = FALSE)
fwrite(vep_input[1:min(.N, 20)], top20_vep_file, sep = "\t", quote = FALSE, col.names = FALSE)
fwrite(lead_ids, id_file, sep = "\t", quote = FALSE, col.names = FALSE)

message("Written: ", vep_file)
message("Written: ", top20_vep_file)
message("Written: ", id_file)
message("Lead SNPs prepared for VEP: ", nrow(vep_input))
