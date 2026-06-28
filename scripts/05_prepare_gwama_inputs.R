library(data.table)

# ============================================================
# 05_prepare_gwama_inputs.R
#
# Purpose:
#   Convert final hg38 harmonised asthma GWAS files into
#   GWAMA quantitative-trait input format using BETA and SE.
#
# Inputs:
#   data/final_hg38/*gwama_ready.tsv
#
# Outputs:
#   config/gwama_input/*.gwama.txt
#   config/gwama_input/gwama.in
#
# Notes:
#   - The harmonisation script has already converted OR effects
#     to log(OR), so BETA/SE are now on a consistent scale.
#   This script keeps only the columns required by the current
#   GWAMA run: MARKERNAME, EA, NEA, BETA, and SE.
#   - gwama.in contains the list of generated GWAMA input files.
#
# Version notes:
#   - Requires the final_hg38 directory to be assembled before
#     running this script.
#   - Does not re-filter variants; filtering is expected to have
#     happened during harmonisation and liftover.
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

final_dir <- file.path(project_dir, "data", "final_hg38")
gwama_dir <- file.path(project_dir, "config", "gwama_input")

dir.create(gwama_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(final_dir, pattern = "gwama_ready.tsv$", full.names = TRUE)

cat("Found", length(files), "final hg38 harmonised files.\n")

for (f in files) {
  dt <- fread(f)
  
  out <- dt[, .(
    MARKERNAME,
    EA,
    NEA,
    BETA,
    SE
  )]
  
  outfile <- file.path(
    gwama_dir,
    paste0(tools::file_path_sans_ext(basename(f)), ".gwama.txt")
  )
  
  fwrite(out, outfile, sep = "\t", quote = FALSE, na = "NA")
  
  cat("Written:", basename(outfile), "| rows:", nrow(out), "\n")
}

gwama_files <- list.files(gwama_dir, pattern = "\\.gwama\\.txt$", full.names = FALSE)

writeLines(gwama_files, file.path(gwama_dir, "gwama.in"))

cat("\nCreated gwama.in with", length(gwama_files), "studies.\n")
