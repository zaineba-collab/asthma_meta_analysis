library(data.table)

# ============================================================
# 12_prepare_sentinel_input.R
#
# Purpose:
#   Prepare the current GWAMA meta-analysis output for distance-based
#   sentinel SNP selection after the corrected Shrine file rerun.
#
# Why this is needed:
#   My supervisor-provided sentinel-selection script expects
#   columns named:
#     SNP, chromosome, position, p_value_lrt, beta1, se1
#
#   GWAMA output uses different names and does not include
#   chromosome/position directly, so this script merges GWAMA
#   results with the hg38 coordinates from the final harmonised
#   files in data/final_hg38.
#
# Inputs:
#   results/gwama/asthma_meta.out
#   data/final_hg38/*gwama_ready.tsv
#
# Output:
#   results/pruning/asthma_fixed_for_sentinel_selection.txt
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

gwama_dir <- file.path(project_dir, "results", "gwama")
final_dir <- file.path(project_dir, "data", "final_hg38")
pruning_dir <- file.path(project_dir, "results", "pruning")

dir.create(pruning_dir, showWarnings = FALSE, recursive = TRUE)

gwama_file <- file.path(gwama_dir, "asthma_meta.out")
output_file <- file.path(pruning_dir, "asthma_fixed_for_sentinel_selection.txt")

message("Reading GWAMA results...")

gwas <- fread(
  gwama_file,
  select = c("rs_number", "beta", "se", "p-value")
)

setnames(
  gwas,
  old = c("rs_number", "beta", "se", "p-value"),
  new = c("SNP", "beta1", "se1", "p_value_lrt")
)

message("Collecting hg38 SNP coordinates from final harmonised files...")

coord_files <- list.files(final_dir, pattern = "gwama_ready.tsv$", full.names = TRUE)

coords <- rbindlist(
  lapply(coord_files, function(f) {
    message("Reading coordinates from: ", basename(f))
    fread(f, select = c("MARKERNAME", "CHR", "POS"))
  }),
  fill = TRUE
)

setnames(coords, old = c("MARKERNAME", "CHR", "POS"),
         new = c("SNP", "chromosome", "position"))

coords <- unique(coords)

message("Merging GWAMA results with coordinates...")

dt <- merge(gwas, coords, by = "SNP", all.x = FALSE, all.y = FALSE)

dt[, chromosome := as.integer(chromosome)]
dt[, position := as.integer(position)]
dt[, p_value_lrt := as.numeric(p_value_lrt)]
dt[, beta1 := as.numeric(beta1)]
dt[, se1 := as.numeric(se1)]

dt <- dt[
  chromosome %in% 1:22 &
    !is.na(position) &
    !is.na(p_value_lrt) &
    p_value_lrt > 0 &
    p_value_lrt <= 1
]

setcolorder(dt, c("SNP", "chromosome", "position", "p_value_lrt", "beta1", "se1"))

fwrite(dt, output_file, sep = "\t", quote = FALSE, na = "NA")

message("Written: ", output_file)
message("Rows: ", nrow(dt))
