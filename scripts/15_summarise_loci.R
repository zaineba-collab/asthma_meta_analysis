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

summary <- data.table(
  Metric = c(
    "Total loci",
    "Largest locus (bp)",
    "Median locus size (bp)",
    "Mean locus size (bp)",
    "Maximum sentinel SNPs within one locus"
  ),
  Value = c(
    nrow(loci),
    max(loci$locus_size),
    median(loci$locus_size),
    round(mean(loci$locus_size)),
    max(loci$n_sentinels_in_locus)
  )
)

fwrite(summary, summary_file, sep="\t")

print(summary)
message("Summary written to ", summary_file)