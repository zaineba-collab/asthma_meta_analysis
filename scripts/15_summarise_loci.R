# ============================================================
# Script: 15_summarise_loci.R
#
# Purpose:
#   Summarise the merged asthma loci generated after
#   distance-based sentinel SNP selection and locus merging.
#
# Input:
#   results/loci/asthma_fixed_sentinel_merge_loci.txt
#
# Output:
#   results/loci/loci_summary.txt
#
# Notes:
#   - Total loci represents the number of merged genomic regions.
#   - Locus size is calculated from the merged +/- 1 Mb sentinel windows.
#   - The maximum sentinel SNP count shows the largest number of
#     sentinel variants merged into a single locus.
#   - These results are used to describe the scale of independent
#     asthma-associated genomic regions identified after pruning.
# ============================================================
library(data.table)

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

input_file <- file.path(
  project_dir,
  "results",
  "loci",
  "asthma_fixed_sentinel_merge_loci.txt"
)

summary_file <- file.path(
  project_dir,
  "results",
  "loci",
  "loci_summary.txt"
)

message("Reading loci...")

loci <- fread(input_file)

sentinels_by_chr <- loci[
  ,
  .(sentinel_snps = sum(n_sentinels_in_locus)),
  by = chrom
][order(-sentinel_snps)]

top_chr <- sentinels_by_chr[1]

# Create a concise summary of the merged loci.
# These metrics are useful for the dissertation Results section.
summary <- data.table(
  Metric = c(
    "Total loci",
    "Largest locus (bp)",
    "Median locus size (bp)",
    "Mean locus size (bp)",
    "Maximum sentinel SNPs within one locus",
    "Chromosome with most sentinel SNPs",
    "Sentinel SNPs on top chromosome"
  ),
  Value = c(
    nrow(loci),
    max(loci$locus_size),
    median(loci$locus_size),
    round(mean(loci$locus_size)),
    max(loci$n_sentinels_in_locus),
    top_chr$chrom,
    top_chr$sentinel_snps
  )
)

# Save the locus summary so it can be reported and reused later.
fwrite(summary, summary_file, sep="\t")

print(summary)
message("Summary written to ", summary_file)
