library(data.table)
library(ggplot2)

# ============================================================
# 10_manhattan_plot.R
#
# Purpose:
#   Create a Manhattan plot from the fixed-effect GWAMA
#   meta-analysis results.
#
# Input:
#   results/gwama/asthma_meta_fixed.out
#
# Output:
#   results/gwama/manhattan_fixed_effect.png
#
# Notes:
#   The GWAMA output does not contain chromosome/position, so
#   this script uses SNP positions from the harmonised final hg38
#   files and merges them onto the GWAMA results by rs_number.
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

gwama_dir <- file.path(project_dir, "results", "gwama")
final_dir <- file.path(project_dir, "data", "final_hg38")

gwama_file <- file.path(gwama_dir, "asthma_meta_fixed.out")
plot_file <- file.path(gwama_dir, "manhattan_fixed_effect.png")

message("Reading GWAMA p-values...")
gwas <- fread(
  gwama_file,
  select = c("rs_number", "p-value")
)

setnames(gwas, old = c("rs_number", "p-value"), new = c("SNP", "P"))

message("Collecting SNP coordinates from final hg38 harmonised files...")

coord_files <- list.files(final_dir, pattern = "gwama_ready.tsv$", full.names = TRUE)

coords <- rbindlist(
  lapply(coord_files, function(f) {
    message("Reading coordinates from: ", basename(f))
    fread(f, select = c("MARKERNAME", "CHR", "POS"))
  }),
  fill = TRUE
)

setnames(coords, "MARKERNAME", "SNP")
coords <- unique(coords)

message("Merging GWAMA results with coordinates...")
dt <- merge(gwas, coords, by = "SNP", all.x = FALSE, all.y = FALSE)

dt <- dt[!is.na(P) & P > 0 & !is.na(CHR) & !is.na(POS)]
dt[, CHR := as.integer(CHR)]
dt <- dt[CHR %in% 1:22]

message("Preparing cumulative chromosome positions...")
setorder(dt, CHR, POS)

chr_lengths <- dt[, .(chr_len = max(POS, na.rm = TRUE)), by = CHR]
chr_lengths[, offset := cumsum(shift(chr_len, fill = 0))]
dt <- merge(dt, chr_lengths[, .(CHR, offset)], by = "CHR")
dt[, BPcum := POS + offset]
dt[, logp := -log10(P)]

axis_df <- dt[, .(center = (min(BPcum) + max(BPcum)) / 2), by = CHR]

message("Creating Manhattan plot...")

p <- ggplot(dt, aes(x = BPcum, y = logp, color = as.factor(CHR %% 2))) +
  geom_point(size = 0.4, alpha = 0.6) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
  scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center) +
  scale_color_manual(values = c("grey30", "grey60")) +
  labs(
    title = "Asthma GWAS meta-analysis Manhattan plot",
    x = "Chromosome",
    y = expression(-log[10](P))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(plot_file, p, width = 14, height = 6, dpi = 300)

message("Written: ", plot_file)
message("Variants plotted: ", nrow(dt))
