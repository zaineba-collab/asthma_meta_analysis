library(data.table)

# ============================================================
# 09_summarise_heterogeneity.R
#
# Purpose:
#   Summarise SNP heterogeneity in the fixed-effect GWAMA
#   meta-analysis results.
#
# Inputs:
#   results/gwama/asthma_meta_fixed.out
#   results/gwama/high_heterogeneity_i2_75.txt
#   results/gwama/significant_heterogeneity_qp.txt
#
# Outputs:
#   results/gwama/fixed_effect_heterogeneity_summary.txt
#   results/gwama/fixed_effect_i2_distribution.txt
#
# Notes:
#   I2 >= 0.75 is treated as high heterogeneity.
#   Q p-value < 0.05 indicates statistically significant
#   evidence of heterogeneity between studies.
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

gwama_results_dir <- file.path(project_dir, "results", "gwama")

fixed_file <- file.path(gwama_results_dir, "asthma_meta_fixed.out")
high_i2_file <- file.path(gwama_results_dir, "high_heterogeneity_i2_75.txt")
significant_q_file <- file.path(gwama_results_dir, "significant_heterogeneity_qp.txt")

summary_file <- file.path(gwama_results_dir, "fixed_effect_heterogeneity_summary.txt")
i2_distribution_file <- file.path(gwama_results_dir, "fixed_effect_i2_distribution.txt")

i2_threshold <- 0.75
q_p_threshold <- 0.05

message("Reading full fixed-effect GWAMA output for total SNP count and I2 distribution...")
fixed <- fread(
  fixed_file,
  select = c("rs_number", "q_p-value", "i2")
)

message("Reading high I2 SNP file...")
high_i2 <- fread(high_i2_file, select = "rs_number")

message("Reading significant Q p-value SNP file...")
significant_q <- fread(significant_q_file, select = "rs_number")

total_snps <- nrow(fixed)
n_high_i2 <- nrow(high_i2) - 1
n_significant_q <- nrow(significant_q) - 1

# Safety check: if fread has already treated the first row as header,
# do not subtract one.
if (!"rs_number" %in% high_i2$rs_number[1]) {
  n_high_i2 <- nrow(high_i2)
}

if (!"rs_number" %in% significant_q$rs_number[1]) {
  n_significant_q <- nrow(significant_q)
}

summary_table <- data.table(
  metric = c(
    "Total SNPs in fixed-effect GWAMA output",
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
    min(fixed$i2, na.rm = TRUE),
    quantile(fixed$i2, 0.25, na.rm = TRUE),
    median(fixed$i2, na.rm = TRUE),
    mean(fixed$i2, na.rm = TRUE),
    quantile(fixed$i2, 0.75, na.rm = TRUE),
    max(fixed$i2, na.rm = TRUE)
  )
)

fwrite(summary_table, summary_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(i2_distribution, i2_distribution_file, sep = "\t", quote = FALSE, na = "NA")

message("Written: ", summary_file)
message("Written: ", i2_distribution_file)

print(summary_table)
print(i2_distribution)

