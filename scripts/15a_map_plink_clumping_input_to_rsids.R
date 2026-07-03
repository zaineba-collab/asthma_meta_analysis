#!/usr/bin/env Rscript

library(data.table)

# ============================================================
# 15a_map_plink_clumping_input_to_rsids.R
#
# Purpose:
#   Convert the GWAMA/PLINK clumping input from coordinate-style
#   variant IDs to rsIDs so the IDs match the 1000 Genomes European
#   PLINK reference panel BIM file.
#
# Input:
#   results/gwama/asthma_fixed_for_plink_clumping.txt
#   data/reference_panel/eur/1000G_EUR_biallelic.bim
#
# Outputs:
#   results/gwama/asthma_fixed_for_plink_clumping_rsids.txt
#   results/gwama/asthma_fixed_for_plink_clumping_unmatched.tsv
#   results/gwama/asthma_fixed_for_plink_clumping_rsid_mapping_summary.txt
#
# Notes:
#   - Input SNP IDs are expected to be CHR:POS:ALLELE1:ALLELE2.
#   - BIM alleles can be in either order, so A1/A2 and A2/A1 are both
#     considered valid allele matches.
#   - Because the downstream PLINK command uses --clump-p1 5e-9 and
#     --clump-p2 5e-9, this script maps only variants with pvalue <= 5e-9.
#     Those are the only variants PLINK will use for this clumping run.
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

gwama_dir <- file.path(project_dir, "results", "gwama")

clump_input_file <- file.path(gwama_dir, "asthma_fixed_for_plink_clumping.txt")
bim_file <- file.path(project_dir, "data", "reference_panel", "eur", "1000G_EUR_biallelic.bim")

rsid_clump_file <- file.path(gwama_dir, "asthma_fixed_for_plink_clumping_rsids.txt")
unmatched_file <- file.path(gwama_dir, "asthma_fixed_for_plink_clumping_unmatched.tsv")
summary_file <- file.path(gwama_dir, "asthma_fixed_for_plink_clumping_rsid_mapping_summary.txt")

clump_p_threshold <- 5e-9

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

stop_if_missing(clump_input_file, "GWAMA clumping input")
stop_if_missing(bim_file, "PLINK reference BIM")

message("Reading clumping input: ", clump_input_file)
clump <- fread(clump_input_file)

required_clump_cols <- c("SNP", "pvalue")
missing_clump_cols <- setdiff(required_clump_cols, names(clump))
if (length(missing_clump_cols) > 0) {
  stop(
    "Clumping input is missing required column(s): ",
    paste(missing_clump_cols, collapse = ", "),
    call. = FALSE
  )
}

clump[, pvalue := as.numeric(pvalue)]
total_variants <- nrow(clump)

clump_candidates <- clump[!is.na(pvalue) & pvalue <= clump_p_threshold]
message("Total variants in clumping input: ", total_variants)
message("Variants with pvalue <= ", clump_p_threshold, ": ", nrow(clump_candidates))

if (nrow(clump_candidates) == 0) {
  stop("No variants pass the clumping p-value threshold.", call. = FALSE)
}

parts <- tstrsplit(clump_candidates$SNP, ":", fixed = TRUE)
if (length(parts) < 4) {
  parts <- c(parts, rep(list(rep(NA_character_, nrow(clump_candidates))), 4 - length(parts)))
}

clump_candidates[, original_SNP := SNP]
clump_candidates[, CHR := parts[[1]]]
clump_candidates[, POS := suppressWarnings(as.integer(parts[[2]]))]
clump_candidates[, input_A1 := toupper(parts[[3]])]
clump_candidates[, input_A2 := toupper(parts[[4]])]

valid_candidates <- clump_candidates[
  !is.na(CHR) &
    !is.na(POS) &
    grepl("^[ACGT]+$", input_A1) &
    grepl("^[ACGT]+$", input_A2)
]

invalid_candidates <- clump_candidates[!valid_candidates, on = "original_SNP"]
if (nrow(invalid_candidates) > 0) {
  invalid_candidates[, reason := "invalid_coordinate_id"]
}

message("Valid coordinate-style variants to map: ", nrow(valid_candidates))

positions_to_keep <- unique(valid_candidates[, .(CHR, POS)])

message("Reading BIM file: ", bim_file)
bim <- fread(
  bim_file,
  header = FALSE,
  select = c(1, 2, 4, 5, 6),
  col.names = c("CHR", "rsID", "POS", "bim_A1", "bim_A2"),
  colClasses = c("character", "character", "integer", "character", "character")
)

bim[, bim_A1 := toupper(bim_A1)]
bim[, bim_A2 := toupper(bim_A2)]

message("Restricting BIM to candidate chromosome/position pairs...")
bim <- bim[positions_to_keep, on = .(CHR, POS), nomatch = 0]
message("BIM rows at candidate positions: ", nrow(bim))

message("Matching by position and allele orientation...")
joined <- merge(
  valid_candidates,
  bim,
  by = c("CHR", "POS"),
  allow.cartesian = TRUE
)

allele_matched <- joined[
  (input_A1 == bim_A1 & input_A2 == bim_A2) |
    (input_A1 == bim_A2 & input_A2 == bim_A1)
]

setorder(allele_matched, original_SNP, rsID)

ambiguous_matches <- allele_matched[
  ,
  .N,
  by = original_SNP
][N > 1, .N]

matched_one <- allele_matched[
  ,
  .SD[1],
  by = original_SNP
]

matched_output <- matched_one[
  ,
  .(
    SNP = rsID,
    original_SNP,
    pvalue
  )
]

unmatched_valid <- valid_candidates[
  !matched_one,
  on = "original_SNP"
][
  ,
  .(
    original_SNP,
    pvalue,
    CHR,
    POS,
    input_A1,
    input_A2,
    reason = "no_position_allele_match_in_bim"
  )
]

unmatched_invalid <- invalid_candidates[
  ,
  .(
    original_SNP,
    pvalue,
    CHR,
    POS,
    input_A1,
    input_A2,
    reason
  )
]

unmatched <- rbindlist(list(unmatched_valid, unmatched_invalid), fill = TRUE)

setorder(matched_output, pvalue)
setorder(unmatched, pvalue)

variants_evaluated <- nrow(clump_candidates)
variants_matched <- nrow(matched_output)
variants_unmatched <- nrow(unmatched)
percent_matched <- 100 * variants_matched / variants_evaluated

summary <- data.table(
  Metric = c(
    "Total variants in clumping input",
    paste0("Variants evaluated for mapping (pvalue <= ", clump_p_threshold, ")"),
    "Variants matched to rsID",
    "Variants not matched",
    "Percentage matched among evaluated variants",
    "Ambiguous allele matches resolved by first rsID"
  ),
  Value = c(
    total_variants,
    variants_evaluated,
    variants_matched,
    variants_unmatched,
    sprintf("%.2f", percent_matched),
    ambiguous_matches
  )
)

fwrite(matched_output, rsid_clump_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(unmatched, unmatched_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(summary, summary_file, sep = "\t", quote = FALSE, na = "NA")

message("Written rsID clumping input: ", rsid_clump_file)
message("Written unmatched variants: ", unmatched_file)
message("Written mapping summary: ", summary_file)
message("Matched variants: ", variants_matched, " / ", variants_evaluated,
        " (", sprintf("%.2f", percent_matched), "%)")

message("\nPLINK command to rerun LD clumping:")
message(
  "plink2 \\\n",
  "  --bfile data/reference_panel/eur/1000G_EUR_biallelic \\\n",
  "  --clump results/gwama/asthma_fixed_for_plink_clumping_rsids.txt \\\n",
  "  --clump-id-field SNP \\\n",
  "  --clump-p-field pvalue \\\n",
  "  --clump-p1 5e-9 \\\n",
  "  --clump-p2 5e-9 \\\n",
  "  --clump-r2 0.1 \\\n",
  "  --clump-kb 99999 \\\n",
  "  --out results/ld_clumping/asthma_fixed_plink_clumped_rsids"
)
