library(data.table)

# ============================================================
# 03_prepare_liftover.R
#
# Purpose:
#   Prepare GRCh37 variants for coordinate liftover to GRCh38.
#
# Inputs:
#   data/harmonised/GCST90029018_buildGRCh37.gwama_ready.tsv
#   data/harmonised/Shrine_30552067_moderate-severe_asthma.txt.gwama_ready.tsv
#
# Outputs:
#   ../liftover/*.hg19.bed
#
# Notes:
#   LiftOver tools usually require BED-like coordinates:
#     chromosome, start, end, variant ID
#
#   BED coordinates are 0-based, so:
#     start = POS - 1
#     end   = POS
#
# Version notes:
#   - This script creates BED input only. The actual coordinate
#     conversion is performed outside R with the liftover tool and
#     the hg19ToHg38.over.chain.gz chain file.
#   - MARKERNAME is carried as the BED name field so lifted
#     coordinates can be merged back onto harmonised rows later.
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

harmonised_dir <- file.path(project_dir, "data", "harmonised")
liftover_dir <- file.path(project_dir, "liftover")

dir.create(liftover_dir, showWarnings = FALSE, recursive = TRUE)

grch37_files <- c(
  "GCST90029018_buildGRCh37.gwama_ready.tsv",
  "Shrine_30552067_moderate-severe_asthma.txt.gwama_ready.tsv"
)

for (f in grch37_files) {
  message("Preparing liftover input for: ", f)
  
  dt <- fread(file.path(harmonised_dir, f))
  
  bed <- dt[, .(
    chrom = paste0("chr", CHR),
    start = POS - 1,
    end = POS,
    name = MARKERNAME
  )]
  
  outfile <- file.path(
    liftover_dir,
    paste0(sub(".gwama_ready.tsv", "", f, fixed = TRUE), ".hg19.bed")
  )
  
  fwrite(bed, outfile, sep = "\t", col.names = FALSE, quote = FALSE)
  
  message("Written: ", outfile, " | rows: ", nrow(bed))
}
