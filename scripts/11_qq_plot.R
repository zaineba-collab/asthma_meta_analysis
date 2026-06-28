library(data.table)
library(ggplot2)

# ============================================================
# 11_qq_plot.R
#
# Purpose:
#   Create a QQ plot from the fixed-effect GWAMA meta-analysis
#   p-values and estimate genomic inflation factor lambda GC.
#
#
# Input file:
#   results/gwama/asthma_meta.out
#
# Output files:
#   results/gwama/qq_plot.png
#   results/gwama/lambda_gc.txt
# ============================================================

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_file <- sub("--file=", "", script_arg)
script_dir <- if (!is.na(script_file)) {
  dirname(normalizePath(script_file))
} else {
  getwd()
}
project_dir <- dirname(script_dir)

gwama_dir <- file.path(project_dir, "results", "gwama")
input_file <- file.path(gwama_dir, "asthma_meta.out")
plot_file <- file.path(gwama_dir, "qq_plot.png")
lambda_file <- file.path(gwama_dir, "lambda_gc.txt")

message("Reading GWAMA p-values...")
dt <- fread(input_file, select = "p-value")
setnames(dt, "p-value", "P")

dt <- dt[!is.na(P) & P > 0 & P <= 1]

message("Calculating expected and observed p-values...")
setorder(dt, P)

dt[, expected := -log10(seq_len(.N) / (.N + 1))]
dt[, observed := -log10(P)]

# Genomic inflation factor lambda GC
# Convert p-values to chi-square statistics with 1 degree of freedom.
chisq <- qchisq(1 - dt$P, df = 1)
lambda_gc <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)

lambda_table <- data.table(
  metric = "Lambda GC",
  value = lambda_gc
)

fwrite(lambda_table, lambda_file, sep = "\t", quote = FALSE)

message("Lambda GC: ", round(lambda_gc, 4))

message("Creating QQ plot...")

p <- ggplot(dt, aes(x = expected, y = observed)) +
  geom_point(size = 0.35, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(
    title = "Asthma GWAS meta-analysis QQ plot",
    subtitle = paste0("GWAMA; lambda GC = ", round(lambda_gc, 3)),
    x = expression(Expected ~ -log[10](P)),
    y = expression(Observed ~ -log[10](P))
  ) +
  theme_minimal(base_size = 12)

ggsave(plot_file, p, width = 6, height = 6, dpi = 300)

message("Written: ", plot_file)
message("Written: ", lambda_file)
