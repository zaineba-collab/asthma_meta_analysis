#!/usr/bin/env Rscript

library(data.table)

# ============================================================
# 16_summarise_plink_clumping.R
#
# Purpose:
#   Summarise PLINK LD clumping results from the asthma
#   meta-analysis run.
#
# Inputs:
#   results/ld_clumping/asthma_fixed_plink_clumped.clumps
#   results/ld_clumping/asthma_fixed_plink_clumped.clumps.missing_id
#   results/ld_clumping/asthma_fixed_plink_clumped.log
#
# Outputs:
#   results/ld_clumping/plink_clumping_summary.txt
#   results/ld_clumping/independent_lead_snps.tsv
#
# Notes:
#   - Each row in the PLINK .clumps file is one independent lead locus.
#   - TOTAL is the number of additional variants assigned to that clump.
#     TOTAL == 0 means the lead SNP had no additional clumped SNPs.
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

clump_dir <- file.path(project_dir, "results", "ld_clumping")

clump_file <- file.path(clump_dir, "asthma_fixed_plink_clumped.clumps")
missing_file <- file.path(clump_dir, "asthma_fixed_plink_clumped.clumps.missing_id")
log_file <- file.path(clump_dir, "asthma_fixed_plink_clumped.log")

summary_file <- file.path(clump_dir, "plink_clumping_summary.txt")
lead_snps_file <- file.path(clump_dir, "independent_lead_snps.tsv")

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

extract_log_integer <- function(log_lines, pattern, replacement) {
  line <- grep(pattern, log_lines, value = TRUE)
  if (length(line) == 0) {
    return(NA_integer_)
  }

  as.integer(sub(replacement, "\\1", line[1]))
}

stop_if_missing(clump_file, "PLINK clumps file")

message("Reading PLINK clumps: ", clump_file)
clumps <- fread(clump_file)

required_cols <- c("#CHROM", "POS", "ID", "P", "TOTAL")
missing_cols <- setdiff(required_cols, names(clumps))
if (length(missing_cols) > 0) {
  stop(
    "Clumps file is missing required column(s): ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

# Ensure numeric summary fields are interpreted consistently across machines.
clumps[, POS := as.integer(POS)]
clumps[, P := as.numeric(P)]
clumps[, TOTAL := as.integer(TOTAL)]

setorder(clumps, P)

# The lead-SNP table is the compact reproducible input for downstream locus
# summaries or annotation: one row per independent LD clump.
lead_snps <- clumps[, .(
  SNP = ID,
  chromosome = `#CHROM`,
  position = POS,
  p_value = P
)]

# Basic counts from the clumps file.
n_clumps <- nrow(clumps)
n_clumps_with_additional_snps <- clumps[TOTAL > 0, .N]
n_clumps_without_additional_snps <- clumps[TOTAL == 0, .N]
n_additional_snps <- sum(clumps$TOTAL, na.rm = TRUE)

# Optional run diagnostics from PLINK sidecar files. These are useful for
# reproducibility, but the summary still works if the sidecar files are absent.
n_missing_variant_ids <- if (file.exists(missing_file)) {
  length(readLines(missing_file, warn = FALSE))
} else {
  NA_integer_
}

n_index_candidates <- NA_integer_
if (file.exists(log_file)) {
  log_lines <- readLines(log_file, warn = FALSE)
  n_index_candidates <- extract_log_integer(
    log_lines,
    "clumps formed from [0-9]+ index candidates",
    ".*formed from ([0-9]+) index candidates.*"
  )
}

summary <- data.table(
  Metric = c(
    "Independent lead loci",
    "Total PLINK LD clumps",
    "Clumps containing additional SNPs",
    "Clumps containing only the lead SNP",
    "Additional SNPs assigned to clumps",
    "PLINK index candidates",
    "Missing variant IDs",
    "Strongest lead SNP",
    "Strongest lead SNP p-value"
  ),
  Value = c(
    n_clumps,
    n_clumps,
    n_clumps_with_additional_snps,
    n_clumps_without_additional_snps,
    n_additional_snps,
    n_index_candidates,
    n_missing_variant_ids,
    if (n_clumps > 0) clumps$ID[1] else NA_character_,
    if (n_clumps > 0) clumps$P[1] else NA_real_
  )
)

fwrite(summary, summary_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(lead_snps, lead_snps_file, sep = "\t", quote = FALSE, na = "NA")

print(summary)
message("Written: ", summary_file)
message("Written: ", lead_snps_file)
