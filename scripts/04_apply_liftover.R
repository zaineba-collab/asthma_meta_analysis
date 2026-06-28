library(data.table)

# ============================================================
# 04_apply_liftover.R
#
# Purpose:
#   Replace GRCh37 positions with lifted GRCh38 positions for
#   the two studies that required coordinate conversion.
#
# Inputs:
#   data/harmonised/*.gwama_ready.tsv for the GRCh37 studies
#   liftover/*.hg38.bed files produced by the liftover tool
#
# Outputs:
#   data/harmonised/*GRCh38.lifted.gwama_ready.tsv
#   data/harmonised/*hg38.lifted.gwama_ready.tsv
#
# Notes:
#   - Merges lifted coordinates back by MARKERNAME.
#   - Drops variants lifted to non-standard/random/alternate
#     contigs because GWAMA-ready files use standard chromosome
#     codes.
#   - Encodes chrX, chrY, chrM/chrMT as 23, 24, and 25 if present.
#
# Version notes:
#   - Paths are anchored to this script's location, so it can be
#     run from the project root with:
#       Rscript scripts/04_apply_liftover.R
#   - BUILD is updated to GRCh38 in the lifted outputs.
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

harmonised_dir <- file.path(project_dir, "data", "harmonised")
liftover_dir <- file.path(project_dir, "liftover")

apply_liftover <- function(harmonised_file, lifted_bed, output_file) {
  message("Applying liftover to: ", harmonised_file)
  
  dt <- fread(file.path(harmonised_dir, harmonised_file))
  
  bed <- fread(
    file.path(liftover_dir, lifted_bed),
    col.names = c("chr38", "start38", "end38", "MARKERNAME")
  )
  
  standard_chr <- paste0("chr", c(1:22, "X", "Y", "M", "MT"))
  dropped_nonstandard <- bed[!chr38 %in% standard_chr, .N]
  if (dropped_nonstandard > 0) {
    message("Dropping non-standard lifted contigs: ", dropped_nonstandard)
  }
  bed <- bed[chr38 %in% standard_chr]
  
  bed[, chr_num := sub("^chr", "", chr38)]
  bed[chr_num == "X", chr_num := "23"]
  bed[chr_num == "Y", chr_num := "24"]
  bed[chr_num %in% c("M", "MT"), chr_num := "25"]
  bed[, CHR := as.integer(chr_num)]
  bed[, chr_num := NULL]
  bed[, POS := as.integer(end38)]
  
  dt[, c("CHR", "POS") := NULL]
  
  merged <- merge(dt, bed[, .(MARKERNAME, CHR, POS)],
                  by = "MARKERNAME",
                  all.x = FALSE,
                  all.y = FALSE)
  
  merged[, BUILD := "GRCh38"]
  
  setcolorder(
    merged,
    c("MARKERNAME", "EA", "NEA", "BETA", "SE", "P", "N",
      "CHR", "POS", "BUILD", "SOURCE_FILE")
  )
  
  fwrite(
    merged,
    file.path(harmonised_dir, output_file),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  
  message("Written: ", output_file, " | rows: ", nrow(merged))
}

apply_liftover(
  harmonised_file = "GCST90029018_buildGRCh37.gwama_ready.tsv",
  lifted_bed = "GCST90029018_buildGRCh37.hg38.bed",
  output_file = "GCST90029018_buildGRCh38.lifted.gwama_ready.tsv"
)

apply_liftover(
  harmonised_file = "Shrine_30552067_moderate-severe_asthma.txt.gwama_ready.tsv",
  lifted_bed = "Shrine_30552067_moderate-severe_asthma.txt.hg38.bed",
  output_file = "Shrine_30552067_moderate-severe_asthma.hg38.lifted.gwama_ready.tsv"
)
