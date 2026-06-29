library(data.table)

# ============================================================
# 13_sentinel_selection.R
#
# Purpose:
#   Perform distance-based sentinel SNP selection using the
#   GWAMA meta-analysis results.
#
# Method:
#   1. Keep SNPs passing p < 5e-9.
#   2. Sort SNPs by p-value.
#   3. Select the most significant SNP as a sentinel.
#   4. Remove all SNPs within ±1 Mb of that sentinel.
#   5. Repeat until no SNPs remain.
#
# Input:
#   results/pruning/asthma_fixed_for_sentinel_selection.txt
#
# Output:
#   results/pruning/asthma_fixed_sentinels_clean.txt
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

input_file <- file.path(
  project_dir,
  "results",
  "pruning",
  "asthma_fixed_for_sentinel_selection.txt"
)

output_dir <- file.path(project_dir, "results", "pruning")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_file <- file.path(output_dir, "asthma_fixed_sentinels_clean.txt")

WIDTH <- 1e6
PTHRESH <- 5e-9
chromosomes <- 1:22

clean_cols <- c(
  "SNP",
  "chromosome",
  "position",
  "p_value_lrt",
  "beta1",
  "se1"
)

selectSentinels <- function(dt, width_bp = 1e6) {
  dt <- copy(dt)
  setorder(dt, pvalue)

  sent <- dt[0]

  while (nrow(dt) > 0) {
    lead <- dt[1]

    sent <- rbind(sent, lead, use.names = TRUE, fill = TRUE)

    dt <- dt[
      chromosome != lead$chromosome |
        position < (lead$position - width_bp) |
        position > (lead$position + width_bp)
    ]
  }

  sent
}

message("Reading sentinel input: ", input_file)

dt <- fread(
  input_file,
  select = c("SNP", "chromosome", "position", "p_value_lrt", "beta1", "se1")
)

dt[, chromosome := as.integer(chromosome)]
dt[, position := as.integer(position)]
dt[, pvalue := as.numeric(p_value_lrt)]
dt[, beta1 := as.numeric(beta1)]
dt[, se1 := as.numeric(se1)]

dt <- dt[
  chromosome %in% chromosomes &
    !is.na(position) &
    !is.na(pvalue) &
    pvalue > 0 &
    pvalue < PTHRESH
]

message("Variants passing p < ", PTHRESH, ": ", nrow(dt))

if (nrow(dt) == 0) {
  stop("No variants passed p < ", PTHRESH, call. = FALSE)
}

message("Selecting sentinels using ±", WIDTH / 1e6, " Mb window...")

sentinels <- selectSentinels(dt, WIDTH)

sentinels_clean <- sentinels[, ..clean_cols]
setorder(sentinels_clean, chromosome, position)

fwrite(sentinels_clean, output_file, sep = "\t", quote = FALSE, na = "NA")

message("Sentinel file written to: ", output_file)
message("Number of sentinel SNPs: ", nrow(sentinels_clean))