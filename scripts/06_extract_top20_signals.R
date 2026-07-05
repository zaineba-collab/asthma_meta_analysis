#!/usr/bin/env Rscript

# ============================================================
# 06_extract_top20_fixed_signals.R
#
# Purpose:
#   Extract genome-wide significant SNPs and the 20 most
#   statistically significant SNPs from the current GWAMA
#   meta-analysis output. These results are used for quality
#   control and to summarise the strongest association signals
#   before downstream pruning and locus definition.
#
# Inputs:
#   ../results/gwama/asthma_meta.out
#
# Outputs:
#   ../results/gwama/genome_wide_significant.txt
#   ../results/gwama/top20_signals.txt
#   ../results/gwama/asthma_fixed_for_plink_clumping.txt
#   ../results/gwama/asthma_fixed_for_plink_clumping_extract_snps.txt
#   ../results/gwama/asthma_fixed_for_plink_clumping_unmatched.tsv
#   ../results/gwama/asthma_fixed_for_plink_clumping_mapping_summary.txt
#
# Output columns:
#   rs_number        = single nucleotide polymorphism (SNP) identifier
#   reference_allele = reference allele
#   other_allele     = alternate/non-reference allele
#   beta             = fixed-effect beta estimate
#   se               = standard error
#   p-value          = association p-value
#   z                = z statistic
#   q_statistic      = Cochran's Q heterogeneity statistic
#   q_p-value        = heterogeneity p-value
#   i2               = I-squared heterogeneity estimate
#   n_studies        = number of contributing studies
#
# Notes:
#   - Genome-wide significant variants are defined as p-value < 5e-8.
#   - Top 20 variants are selected from the genome-wide significant
#     variants and sorted by ascending p-value.
#   - PLINK LD clumping uses a stricter p-value threshold of 5e-9.
#   - GWAMA variant IDs are chromosome:position:allele1:allele2, while
#     the 1000 Genomes PLINK reference panel uses rsIDs in the BIM file.
#     This script maps the clumping input to rsIDs before writing it.
#   - Only the columns listed above are retained for a compact
#     fixed-effect summary table.
#
# Version notes:
#   - Paths are anchored to this script's location, so it can be
#     run from the project root with:
#       Rscript scripts/06_extract_top20_signals.R
#   - The script stops with a clear error if any required column is
#     missing from asthma_meta.out.
# ============================================================

library(data.table)

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

gwama_results_dir <- file.path(project_dir, "results", "gwama")
input_file <- file.path(gwama_results_dir, "asthma_meta.out")
bim_file <- file.path(project_dir, "data", "reference_panel", "eur", "1000G_EUR_biallelic.bim")
genome_wide_file <- file.path(gwama_results_dir, "genome_wide_significant.txt")
output_file <- file.path(gwama_results_dir, "top20_signals.txt")
plink_clump_file <- file.path(gwama_results_dir, "asthma_fixed_for_plink_clumping.txt")
plink_extract_file <- file.path(gwama_results_dir, "asthma_fixed_for_plink_clumping_extract_snps.txt")
plink_unmatched_file <- file.path(gwama_results_dir, "asthma_fixed_for_plink_clumping_unmatched.tsv")
plink_mapping_summary_file <- file.path(gwama_results_dir, "asthma_fixed_for_plink_clumping_mapping_summary.txt")

genome_wide_threshold <- 5e-8
plink_clump_threshold <- 5e-9

required_columns <- c(
  "rs_number",
  "reference_allele",
  "other_allele",
  "beta",
  "se",
  "p-value",
  "z",
  "q_statistic",
  "q_p-value",
  "i2",
  "n_studies"
)

message("Reading: ", input_file)
results <- fread(input_file)

missing_columns <- setdiff(required_columns, names(results))
if (length(missing_columns) > 0) {
  stop(
    "Missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

genome_wide_significant <- results[`p-value` < genome_wide_threshold]

fwrite(
  genome_wide_significant,
  genome_wide_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message(
  "Written: ",
  genome_wide_file,
  " | rows: ",
  nrow(genome_wide_significant)
)

top20 <- genome_wide_significant[
  order(`p-value`),
  ..required_columns
][1:20]

fwrite(top20, output_file, sep = "\t", quote = FALSE, na = "NA")

message("Written: ", output_file, " | rows: ", nrow(top20))

message("Preparing rsID-based PLINK clumping input...")

if (!file.exists(bim_file)) {
  stop("PLINK reference BIM not found: ", bim_file, call. = FALSE)
}

clump_candidates <- results[
  !is.na(`p-value`) & `p-value` <= plink_clump_threshold,
  .(
    original_SNP = rs_number,
    pvalue = `p-value`
  )
]

message("Variants with p-value <= ", plink_clump_threshold, ": ", nrow(clump_candidates))

rsid_candidates <- clump_candidates[grepl("^rs[0-9]+$", original_SNP)]
coordinate_candidates <- clump_candidates[!grepl("^rs[0-9]+$", original_SNP)]

message("Candidate variants already using rsIDs: ", nrow(rsid_candidates))
message("Candidate variants using coordinate-style IDs: ", nrow(coordinate_candidates))

parts <- tstrsplit(coordinate_candidates$original_SNP, ":", fixed = TRUE)
if (length(parts) < 4) {
  parts <- c(parts, rep(list(rep(NA_character_, nrow(coordinate_candidates))), 4 - length(parts)))
}

coordinate_candidates[, CHR := parts[[1]]]
coordinate_candidates[, POS := suppressWarnings(as.integer(parts[[2]]))]
coordinate_candidates[, input_A1 := toupper(parts[[3]])]
coordinate_candidates[, input_A2 := toupper(parts[[4]])]

valid_coordinate_candidates <- coordinate_candidates[
  !is.na(CHR) &
    !is.na(POS) &
    grepl("^[ACGT]+$", input_A1) &
    grepl("^[ACGT]+$", input_A2)
]

invalid_coordinate_candidates <- coordinate_candidates[
  !valid_coordinate_candidates,
  on = "original_SNP"
]

if (nrow(invalid_coordinate_candidates) > 0) {
  invalid_coordinate_candidates[, reason := "invalid_coordinate_id"]
}

candidate_positions <- unique(valid_coordinate_candidates[, .(CHR, POS)])

message("Filtering reference BIM to candidate rsIDs and positions: ", bim_file)
rsid_key_file <- tempfile(pattern = "plink_clump_rsids_", fileext = ".txt")
position_key_file <- tempfile(pattern = "plink_clump_positions_", fileext = ".txt")
writeLines(rsid_candidates$original_SNP, rsid_key_file)
writeLines(paste(candidate_positions$CHR, candidate_positions$POS, sep = ":"), position_key_file)

bim_filter_cmd <- sprintf(
  "awk 'FILENAME==ARGV[1] { ids[$1]; next } FILENAME==ARGV[2] { pos[$1]; next } (($2 in ids) || (($1 \":\" $4) in pos))' %s %s %s",
  shQuote(rsid_key_file),
  shQuote(position_key_file),
  shQuote(bim_file)
)

bim <- fread(
  cmd = bim_filter_cmd,
  header = FALSE,
  col.names = c("CHR", "rsID", "genetic_distance", "POS", "bim_A1", "bim_A2")
)[
  ,
  .(CHR, rsID, POS, bim_A1, bim_A2)
]

unlink(rsid_key_file)
unlink(position_key_file)

bim[, CHR := as.character(CHR)]
bim[, POS := as.integer(POS)]
bim[, bim_A1 := toupper(bim_A1)]
bim[, bim_A2 := toupper(bim_A2)]

message("BIM rows at candidate IDs/positions: ", nrow(bim))

direct_rsid_matches <- merge(
  rsid_candidates,
  unique(bim[, .(rsID)]),
  by.x = "original_SNP",
  by.y = "rsID",
  all = FALSE,
  sort = FALSE
)[
  ,
  .(
    SNP = original_SNP,
    original_SNP,
    pvalue
  )
]

joined <- merge(
  valid_coordinate_candidates,
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

coordinate_matches <- allele_matched[
  ,
  .SD[1],
  by = original_SNP
]

coordinate_rsid_matches <- coordinate_matches[
  ,
  .(
    SNP = rsID,
    original_SNP,
    pvalue
  )
]

plink_clump_input <- rbindlist(
  list(direct_rsid_matches, coordinate_rsid_matches),
  fill = TRUE
)

duplicate_rsid_rows <- nrow(plink_clump_input) - uniqueN(plink_clump_input$SNP)
setorder(plink_clump_input, SNP, pvalue)
plink_clump_input <- plink_clump_input[
  ,
  .SD[1],
  by = SNP
]

unmatched_rsid <- rsid_candidates[
  !direct_rsid_matches,
  on = "original_SNP"
][
  ,
  .(
    original_SNP,
    pvalue,
    CHR = NA_character_,
    POS = NA_integer_,
    input_A1 = NA_character_,
    input_A2 = NA_character_,
    reason = "rsid_not_found_in_bim"
  )
]

unmatched_valid <- valid_coordinate_candidates[
  !coordinate_matches,
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

unmatched_invalid <- invalid_coordinate_candidates[
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

unmatched <- rbindlist(list(unmatched_rsid, unmatched_valid, unmatched_invalid), fill = TRUE)

setorder(plink_clump_input, pvalue)
setorder(unmatched, pvalue)

mapping_summary <- data.table(
  Metric = c(
    "Total variants in GWAMA output",
    paste0("Variants evaluated for PLINK clumping input (p-value <= ", plink_clump_threshold, ")"),
    "Variants matched to rsID",
    "Variants not matched",
    "Percentage matched among evaluated variants",
    "Duplicate rsID rows collapsed",
    "Ambiguous matches resolved by first rsID"
  ),
  Value = c(
    nrow(results),
    nrow(clump_candidates),
    nrow(plink_clump_input),
    nrow(unmatched),
    sprintf("%.2f", 100 * nrow(plink_clump_input) / nrow(clump_candidates)),
    duplicate_rsid_rows,
    ambiguous_matches
  )
)

fwrite(plink_clump_input, plink_clump_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(plink_clump_input[, .(SNP)], plink_extract_file, sep = "\t", quote = FALSE, col.names = FALSE)
fwrite(unmatched, plink_unmatched_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(mapping_summary, plink_mapping_summary_file, sep = "\t", quote = FALSE, na = "NA")

message("Written: ", plink_clump_file, " | rows: ", nrow(plink_clump_input))
message("Written: ", plink_extract_file, " | rows: ", nrow(plink_clump_input))
message("Written: ", plink_unmatched_file, " | rows: ", nrow(unmatched))
message("Written: ", plink_mapping_summary_file)
message(
  "Matched ",
  nrow(plink_clump_input),
  " of ",
  nrow(clump_candidates),
  " PLINK clumping candidates to rsIDs."
)

message("\nPLINK command:")
message(
  "plink2 \\\n",
  "  --bfile data/reference_panel/eur/1000G_EUR_biallelic \\\n",
  "  --extract results/gwama/asthma_fixed_for_plink_clumping_extract_snps.txt \\\n",
  "  --clump results/gwama/asthma_fixed_for_plink_clumping.txt \\\n",
  "  --clump-id-field SNP \\\n",
  "  --clump-p-field pvalue \\\n",
  "  --clump-p1 5e-9 \\\n",
  "  --clump-p2 5e-9 \\\n",
  "  --clump-r2 0.1 \\\n",
  "  --clump-kb 99999 \\\n",
  "  --out results/ld_clumping/asthma_fixed_plink_clumped"
)
