library(data.table)

# ============================================================
# 09_summarise_heterogeneity.R
#
# Purpose:
#   Summarise heterogeneity across the corrected fixed-effect
#   GWAMA meta-analysis by reporting:
#     - total number of analysed SNPs,
#     - variants with high heterogeneity (I² ≥ 0.75),
#     - variants with significant Cochran's Q test (Q p < 0.05),
#     - distribution of I² values.
#
#   These summaries are used as quality-control metrics to assess
#   the consistency of genetic effects across contributing studies
#   before downstream interpretation.
#
# Inputs:
#   results/gwama/asthma_meta.out
#
# Outputs:
#   results/gwama/high_heterogeneity_i2_75.txt
#   results/gwama/significant_heterogeneity_qp.txt
#   results/gwama/heterogeneity_summary.txt
#   results/gwama/i2_distribution.txt
#
# Notes:
#   I2 >= 0.75 is treated as high heterogeneity.
#   Q p-value < 0.05 indicates statistically significant
#   evidence of heterogeneity between studies.
#   The high-I2 and significant-Q files are regenerated directly
#   from asthma_meta.out, so they match the corrected Shrine run.
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

gwama_results_dir <- file.path(project_dir, "results", "gwama")

input_file <- file.path(gwama_results_dir, "asthma_meta.out")
high_i2_file <- file.path(gwama_results_dir, "high_heterogeneity_i2_75.txt")
significant_q_file <- file.path(gwama_results_dir, "significant_heterogeneity_qp.txt")

summary_file <- file.path(gwama_results_dir, "heterogeneity_summary.txt")
i2_distribution_file <- file.path(gwama_results_dir, "i2_distribution.txt")

i2_threshold <- 0.75
q_p_threshold <- 0.05

message("Reading GWAMA output for heterogeneity columns...")
results <- fread(
  input_file,
  select = c("rs_number", "q_p-value", "q_statistic", "i2")
)

message("Filtering high I2 SNP file...")
high_i2 <- results[i2 >= i2_threshold, .(rs_number)]

message("Filtering significant Q p-value SNP file...")
significant_q <- results[`q_p-value` < q_p_threshold, .(rs_number)]

fwrite(high_i2, high_i2_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(significant_q, significant_q_file, sep = "\t", quote = FALSE, na = "NA")

total_snps <- nrow(results)
n_high_i2 <- nrow(high_i2) 
n_significant_q <- nrow(significant_q) 

summary_table <- data.table(
  metric = c(
    "Total SNPs in GWAMA output",
    "SNPs with I2 >= 0.75",
    "Percentage of SNPs with I2 >= 0.75",
    "SNPs with Q p-value < 0.05",
    "Percentage of SNPs with Q p-value < 0.05"
  ),
  value = c(
    total_snps,
    n_high_i2,
    round((n_high_i2 / total_snps) * 100, 3),
    n_significant_q,
    round((n_significant_q / total_snps) * 100, 3)
  )
)

i2_distribution <- data.table(
  statistic = c(
    "Minimum I2",
    "1st quartile I2",
    "Median I2",
    "Mean I2",
    "3rd quartile I2",
    "Maximum I2"
  ),
  value = c(
    min(results$i2, na.rm = TRUE),
    quantile(results$i2, 0.25, na.rm = TRUE),
    median(results$i2, na.rm = TRUE),
    mean(results$i2, na.rm = TRUE),
    quantile(results$i2, 0.75, na.rm = TRUE),
    max(results$i2, na.rm = TRUE)
  )
)

fwrite(summary_table, summary_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(i2_distribution, i2_distribution_file, sep = "\t", quote = FALSE, na = "NA")

message("Written: ", summary_file)
message("Written: ", i2_distribution_file)

print(summary_table)
print(i2_distribution)
