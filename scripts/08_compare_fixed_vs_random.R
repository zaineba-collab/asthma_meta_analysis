library(data.table)

# ============================================================
# 08_compare_fixed_vs_random.R
#
# Purpose:
#   Compare the top GWAMA fixed-effect and random-effects SNPs.
#   This helps assess whether the strongest association signals
#   are robust across both meta-analysis models.
#
# Inputs:
#   results/gwama/top20_fixed_signals.txt
#   results/gwama/top20_random_signals.txt
#
# Output:
#   results/gwama/top20_fixed_random_comparison.txt
# ============================================================

gwama_dir <- "../results/gwama"

fixed <- fread(file.path(gwama_dir, "top20_fixed_signals.txt"))
random <- fread(file.path(gwama_dir, "top20_random_signals.txt"))

comparison <- merge(
  fixed,
  random,
  by = "rs_number",
  all = TRUE,
  suffixes = c("_fixed", "_random")
)

comparison[, in_fixed := !is.na(`p-value_fixed`)]
comparison[, in_random := !is.na(`p-value_random`)]

outfile <- file.path(gwama_dir, "top20_fixed_random_comparison.txt")

fwrite(comparison, outfile, sep = "\t", quote = FALSE, na = "NA")

message("Written: ", outfile)
message("Rows: ", nrow(comparison))