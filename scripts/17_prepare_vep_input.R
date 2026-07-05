#!/usr/bin/env Rscript

library(data.table)

# ============================================================
# 17_prepare_vep_input.R
#
# Purpose:
#   Prepare rsID-based independent lead SNPs from PLINK clumping
#   for Ensembl Variant Effect Predictor (VEP) annotation.
#
# Inputs:
#   results/ld_clumping/independent_lead_snps.tsv
#   data/reference_panel/eur/1000G_EUR_biallelic.bim
#
# Outputs:
#   results/annotation/independent_lead_snps_vep_input.txt
#   results/annotation/top20_independent_lead_snps_vep_input.txt
#   results/annotation/independent_lead_snps_ids.txt
#
# Notes:
#   - The cleaned clumping pipeline now uses rsIDs as SNP identifiers.
#   - VEP six-column input still needs an allele change, so this script
#     looks up the alleles for each rsID in the PLINK BIM reference file.
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

lead_snps_file <- file.path(project_dir, "results", "ld_clumping", "independent_lead_snps.tsv")
bim_file <- file.path(project_dir, "data", "reference_panel", "eur", "1000G_EUR_biallelic.bim")

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
stop_if_missing(bim_file, "PLINK reference BIM")

message("Reading independent lead SNPs: ", lead_snps_file)
lead <- fread(lead_snps_file)

required_cols <- c("SNP", "chromosome", "position", "p_value")
missing_cols <- setdiff(required_cols, names(lead))
if (length(missing_cols) > 0) {
  stop("Lead SNP table is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

lead[, chromosome := as.character(chromosome)]
lead[, position := as.integer(position)]
lead[, p_value := as.numeric(p_value)]

lead <- lead[
  grepl("^rs[0-9]+$", SNP) &
    chromosome %in% as.character(1:22) &
    !is.na(position) &
    !is.na(p_value)
]

if (nrow(lead) == 0) {
  stop("No autosomal rsID lead SNPs with valid position and p-value were found.", call. = FALSE)
}

message("Filtering BIM alleles for lead rsIDs...")
rsid_file <- tempfile(pattern = "vep_lead_rsids_", fileext = ".txt")
writeLines(lead$SNP, rsid_file)

bim_filter_cmd <- sprintf(
  "awk 'NR==FNR { ids[$1]; next } ($2 in ids)' %s %s",
  shQuote(rsid_file),
  shQuote(bim_file)
)

bim <- fread(
  cmd = bim_filter_cmd,
  header = FALSE,
  col.names = c("bim_chromosome", "SNP", "genetic_distance", "bim_position", "A1", "A2")
)[
  ,
  .(SNP, bim_chromosome, bim_position, A1, A2)
]

unlink(rsid_file)

bim[, bim_chromosome := as.character(bim_chromosome)]
bim[, bim_position := as.integer(bim_position)]
bim[, A1 := toupper(A1)]
bim[, A2 := toupper(A2)]

lead <- merge(lead, bim, by = "SNP", all.x = TRUE, sort = FALSE)

missing_alleles <- lead[is.na(A1) | is.na(A2)]
if (nrow(missing_alleles) > 0) {
  stop("Missing BIM allele lookup for ", nrow(missing_alleles), " lead SNP(s).", call. = FALSE)
}

position_mismatch <- lead[
  chromosome != bim_chromosome |
    position != bim_position
]
if (nrow(position_mismatch) > 0) {
  warning("BIM position differs from clump position for ", nrow(position_mismatch), " SNP(s). Using clump positions.")
}

lead[, allele := paste0(A1, "/", A2)]
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
