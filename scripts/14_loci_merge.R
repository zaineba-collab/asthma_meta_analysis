# ============================================================
# 14_loci_merge.R
#
# Purpose:
#   Merge nearby sentinel SNPs into broader genomic loci.
#
# Method:
#   Each sentinel SNP is expanded by +/- 1 Mb.
#   Overlapping windows are merged into loci.
#
# Input:
#   results/pruning/asthma_fixed_sentinels_clean.txt
#
# Output:
#   results/loci/asthma_fixed_sentinel_merge_loci.txt
# ============================================================
sink(stderr())

library(data.table)
library(dplyr)

suppressMessages(library(GenomicRanges))
suppressMessages(library(tidyverse))

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

trait <- "asthma_fixed"

input_file <- file.path(
  project_dir,
  "results",
  "pruning",
  "asthma_fixed_sentinels_clean.txt"
)

output_dir <- file.path(project_dir, "results", "loci")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_file <- file.path(output_dir, paste0(trait, "_sentinel_merge_loci.txt"))

width <- 1e6
gap <- 0

message("Reading sentinel file: ", input_file)

sentinels <- fread(input_file) %>%
  arrange(chromosome, position) %>%
  mutate(
    idx = row_number(),
    start = pmax(position - width, 1),
    end = position + width
  )

message("Number of sentinel SNPs: ", nrow(sentinels))

gr <- makeGRangesFromDataFrame(
  sentinels,
  seqnames.field = "chromosome",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE
)

merged_gr <- GenomicRanges::reduce(gr, min.gapwidth = gap + 1)

loci <- data.table(
  locus = paste0("locus", seq_along(merged_gr)),
  chrom = as.character(seqnames(merged_gr)),
  locus_start = start(merged_gr),
  locus_end = end(merged_gr),
  locus_size = IRanges::width(merged_gr)
)

sentinel_dt <- as.data.table(sentinels)

assignments <- findOverlaps(gr, merged_gr) %>%
  as_tibble() %>%
  mutate(
    SNP = sentinel_dt$SNP[queryHits],
    p_value_lrt = sentinel_dt$p_value_lrt[queryHits],
    beta1 = sentinel_dt$beta1[queryHits],
    se1 = sentinel_dt$se1[queryHits],
    locus = paste0("locus", subjectHits)
  )

lead_snps <- assignments %>%
  group_by(locus) %>%
  arrange(p_value_lrt, .by_group = TRUE) %>%
  summarise(
    first = dplyr::first(SNP),
    lead_p = dplyr::first(p_value_lrt),
    beta1 = dplyr::first(beta1),
    se1 = dplyr::first(se1),
    n_sentinels_in_locus = n(),
    merge.list = paste(unique(SNP), collapse = " "),
    .groups = "drop"
  )

loci <- loci %>%
  left_join(lead_snps, by = "locus") %>%
  arrange(as.integer(chrom), locus_start)

write.table(loci, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

sink()
message("Locus merge table written to: ", output_file)
message("Number of merged loci: ", nrow(loci))
