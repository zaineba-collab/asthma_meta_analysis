#!/usr/bin/env Rscript

# ============================================================
# 20_build_eur_reference_panel.R
# 
#
# Purpose:
#   Build the full European 1000 Genomes Phase 3 PLINK LD
#   reference panel for LD clumping.
#
# Inputs:
#   reference_panel/raw/all_phase3.pgen.zst
#   reference_panel/raw/all_phase3.pvar.zst
#   reference_panel/raw/phase3_corrected.psam
#
# Outputs:
#   reference_panel/raw/all_phase3.pgen
#   reference_panel/raw/eur_samples_iid.txt
#   reference_panel/eur/1000G_EUR_biallelic.bed
#   reference_panel/eur/1000G_EUR_biallelic.bim
#   reference_panel/eur/1000G_EUR_biallelic.fam
#   reference_panel/eur/1000G_EUR_biallelic.log
#
# Notes:
#   - This script keeps all 1000 Genomes samples where SuperPop == "EUR".
#   - PLINK2 performs the decompression and reference panel build.
#   - Do not commit reference_panel/raw/ or reference_panel/eur/.
#     This repository currently ignores reference_panel/.
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

plink2 <- "C:/Users/ZainebAHMED(Student)/OneDrive - Birkbeck, University of London/Documents/Tools/plink2/plink2.exe"

raw_dir <- file.path(project_dir, "reference_panel", "raw")
eur_dir <- file.path(project_dir, "reference_panel", "eur")

pgen_zst <- file.path(raw_dir, "all_phase3.pgen.zst")
pvar_zst <- file.path(raw_dir, "all_phase3.pvar.zst")
psam_file <- file.path(raw_dir, "phase3_corrected.psam")
pgen_file <- file.path(raw_dir, "all_phase3.pgen")
eur_keep_file <- file.path(raw_dir, "eur_samples_iid.txt")
out_prefix <- file.path("..", "eur", "1000G_EUR_biallelic")
out_prefix_abs <- file.path(eur_dir, "1000G_EUR_biallelic")

dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(eur_dir, showWarnings = FALSE, recursive = TRUE)

cat("Project directory: ", project_dir, "\n", sep = "")
cat("Raw reference directory: ", raw_dir, "\n", sep = "")
cat("European output directory: ", eur_dir, "\n", sep = "")

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

run_plink2 <- function(args, wd) {
  cat("\nRunning PLINK2:\n  ", plink2, " ", paste(args, collapse = " "), "\n", sep = "")
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(wd)

  status <- system2(plink2, args = args, stdout = TRUE, stderr = TRUE, wait = TRUE)
  cat(paste(status, collapse = "\n"), "\n", sep = "")

  exit_code <- attr(status, "status")
  if (!is.null(exit_code) && exit_code != 0) {
    stop("PLINK2 failed with exit code ", exit_code, call. = FALSE)
  }
}

# Check that PLINK2 and all required source files are available.
stop_if_missing(plink2, "PLINK2 executable")
stop_if_missing(pgen_zst, "Compressed PGEN source")
stop_if_missing(pvar_zst, "Compressed PVAR source")
stop_if_missing(psam_file, "Corrected PSAM sample file")

pgen_has_magic <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 2) {
    return(FALSE)
  }

  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  identical(readBin(con, what = "raw", n = 2), as.raw(c(0x6c, 0x1b)))
}

# Remove a likely corrupted decompressed PGEN. The previous failure was caused
# by invalid first bytes, so checking the PLINK2 PGEN magic number is more
# reliable than checking file size.
if (file.exists(pgen_file)) {
  if (!pgen_has_magic(pgen_file)) {
    cat("\nRemoving invalid all_phase3.pgen; first bytes do not match PLINK2 PGEN magic.\n")
    unlink(pgen_file)
  } else {
    cat("\nExisting all_phase3.pgen has valid PLINK2 PGEN magic; keeping it.\n")
  }
}

# Decompress the PGEN using PLINK2's output argument, not shell redirection.
if (!file.exists(pgen_file)) {
  run_plink2(
    args = c("--zst-decompress", "all_phase3.pgen.zst", "all_phase3.pgen"),
    wd = raw_dir
  )
}

stop_if_missing(pgen_file, "Decompressed PGEN")

# Create the PLINK keep file from all European 1000 Genomes samples.
psam <- read.table(
  psam_file,
  header = TRUE,
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_cols <- c("#IID", "SuperPop")
missing_cols <- setdiff(required_cols, names(psam))
if (length(missing_cols) > 0) {
  stop("PSAM is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

eur_samples <- psam[psam$SuperPop == "EUR", "#IID", drop = TRUE]
write.table(
  eur_samples,
  file = eur_keep_file,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("\nEUR sample count: ", length(eur_samples), "\n", sep = "")
if (length(eur_samples) != 633) {
  warning("Expected 633 EUR samples, but found ", length(eur_samples), call. = FALSE)
}

# Build the all-European biallelic SNP reference panel. The PVAR remains
# compressed, and PLINK2 reads it directly with '--pfile all_phase3 vzs'.
run_plink2(
  args = c(
    "--pfile", "all_phase3", "vzs",
    "--psam", "phase3_corrected.psam",
    "--keep", "eur_samples_iid.txt",
    "--snps-only", "just-acgt",
    "--max-alleles", "2",
    "--make-bed",
    "--out", out_prefix
  ),
  wd = raw_dir
)

# Confirm the expected PLINK binary output files were created.
expected_outputs <- paste0(out_prefix_abs, c(".bed", ".bim", ".fam", ".log"))
created <- file.exists(expected_outputs)
output_check <- data.frame(file = expected_outputs, created = created)

cat("\nReference panel output check:\n")
print(output_check, row.names = FALSE)

if (!all(created[1:3])) {
  stop("One or more required PLINK reference panel files were not created.", call. = FALSE)
}

cat("\nFinished building the all-EUR 1000 Genomes PLINK reference panel.\n")
