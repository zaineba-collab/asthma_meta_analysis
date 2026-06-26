library(data.table)

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

clump_dir <- file.path(project_dir, "results", "ld_clumping")

clump_file <- file.path(clump_dir, "asthma_fixed_plink_clumped.clumps")
missing_file <- file.path(clump_dir, "asthma_fixed_plink_clumped.clumps.missing_id")
log_file <- file.path(clump_dir, "asthma_fixed_plink_clumped.log")
filtered_clump_input <- file.path(
  project_dir,
  "results",
  "gwama",
  "asthma_fixed_for_plink_clumping_p5e-9.txt"
)

summary_file <- file.path(clump_dir, "plink_clumping_summary.txt")
top20_file <- file.path(clump_dir, "top20_plink_clumped_lead_snps.txt")
lead_snps_file <- file.path(clump_dir, "independent_lead_snps.txt")

clumps <- fread(clump_file)

n_clumps <- nrow(clumps)
n_missing <- if (file.exists(missing_file)) length(readLines(missing_file)) else 0
n_significant_candidates <- if (file.exists(filtered_clump_input)) {
  nrow(fread(filtered_clump_input))
} else {
  NA_integer_
}
n_index_candidates <- NA_integer_
if (file.exists(log_file)) {
  log_lines <- readLines(log_file, warn = FALSE)
  index_line <- grep("clumps formed from [0-9]+ index candidates", log_lines, value = TRUE)
  if (length(index_line) > 0) {
    n_index_candidates <- as.integer(sub(".*formed from ([0-9]+) index candidates.*", "\\1", index_line[1]))
  }
}

setorder(clumps, P)

lead_snps <- clumps[, .(
  CHR = `#CHROM`,
  BP = POS,
  SNP = ID,
  P,
  TOTAL
)]

summary <- data.table(
  Metric = c(
    "Significant candidate variants in clumping input",
    "PLINK index candidates after missing IDs removed",
    "Independent lead SNPs",
    "Total PLINK LD clumps",
    "Missing variant IDs",
    "Strongest clumped SNP",
    "Strongest clumped p-value"
  ),
  Value = c(
    n_significant_candidates,
    n_index_candidates,
    n_clumps,
    n_clumps,
    n_missing,
    clumps$ID[1],
    clumps$P[1]
  )
)

top20 <- clumps[1:20, .(
  CHROM = `#CHROM`,
  POS,
  ID,
  P,
  TOTAL
)]

fwrite(summary, summary_file, sep = "\t", quote = FALSE)
fwrite(lead_snps, lead_snps_file, sep = "\t", quote = FALSE)
fwrite(top20, top20_file, sep = "\t", quote = FALSE)

print(summary)
message("Written: ", summary_file)
message("Written: ", lead_snps_file)
message("Written: ", top20_file)
