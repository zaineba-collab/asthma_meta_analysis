#!/usr/bin/env Rscript

# ============================================================
# Project: Asthma GWAS Meta-analysis
# Script:
# Author: Zaineb Ahmed
#
# Purpose:
#   Inspect the raw asthma GWAS summary statistics before
#   harmonisation.
#
# Inputs:
#   data/raw/*
#
# Outputs:
#   logs/asthma_input_inventory.tsv
#
# What this script records:
#   - Original column names for each raw file.
#   - Likely genome build inferred from file names and columns.
#   - Whether each file appears to contain BETA, SE, OR, CI, and
#     p-value fields.
#
# Version notes:
#   - This is an audit-only script. It does not modify raw data.
#   - Genome build calls are heuristic and should be confirmed
#     against source documentation where possible.
# ============================================================

library(data.table)

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

raw_dir <- file.path(project_dir, "data", "raw")
log_dir <- file.path(project_dir, "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(raw_dir, full.names = TRUE)

audit_one <- function(file) {
  dt <- fread(file, nrows = 5)
  cols <- names(dt)
  cols_lower <- tolower(cols)
  
  data.table(
    file = basename(file),
    columns = paste(cols, collapse = " | "),
    likely_build = ifelse(grepl("grch38|buildgrch38|hg38", basename(file), ignore.case = TRUE), "GRCh38",
                   ifelse(grepl("grch37|buildgrch37|hg19|b37|position_b37", paste(c(basename(file), cols), collapse = " "), ignore.case = TRUE), "GRCh37",
                   "unknown")),
    has_beta = any(cols_lower %in% c("beta", "effect", "b")),
    has_se = any(cols_lower %in% c("se", "standard_error", "se_gc")),
    has_or = any(cols_lower %in% c("or", "odds_ratio")),
    has_ci = any(cols_lower %in% c("or_95l", "or_95u", "ci_lower", "ci_upper", "or_l95", "or_u95")),
    has_p = any(cols_lower %in% c("p", "p_value", "p_gc", "p_nospa"))
  )
}

audit <- rbindlist(lapply(files, audit_one), fill = TRUE)

fwrite(audit, file.path(log_dir, "asthma_input_inventory.tsv"), sep = "\t")

print(audit)
cat("\nWritten to ", file.path(log_dir, "asthma_input_inventory.tsv"), "\n", sep = "")
