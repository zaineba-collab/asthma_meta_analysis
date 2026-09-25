#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

# Harmonise the fixed asthma exposure with corrected local FinnGen outcomes.
# This script does not clump, use proxies, query external services, or run MR.

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(arg)) sub("^--file=", "", arg[1]) else ""
script_dir <- if (nzchar(script_file)) dirname(normalizePath(script_file)) else getwd()
root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
p <- function(...) file.path(root, ...)
fail <- function(...) stop(..., call. = FALSE)

exposure_file <- p("results", "mr", "finngen", "exposure", "asthma_exposure_5254_grch38.tsv")
outcome_dir <- p("results", "mr", "finngen", "outcomes_grch38")
config_file <- p("config", "finngen_mr_outcomes.tsv")
harm_dir <- p("results", "mr", "finngen", "harmonised")
ready_dir <- file.path(harm_dir, "mr_ready")
excluded_dir <- file.path(harm_dir, "excluded")
qc_dir <- p("results", "mr", "finngen", "qc")
qc_file <- file.path(qc_dir, "finngen_harmonisation_qc.tsv")
direction_file <- file.path(qc_dir, "harmonisation_direction_audit.tsv")
conflict_file <- file.path(qc_dir, "conflicting_variant_harmonisation_status.tsv")
error_file <- file.path(qc_dir, "finngen_harmonisation_errors.tsv")
dir.create(ready_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(excluded_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(exposure_file)) fail("Corrected exposure missing: ", exposure_file)
if (!file.exists(config_file)) fail("FinnGen MR config missing: ", config_file)
exposure <- fread(exposure_file)
cfg <- fread(config_file, na.strings = c("NA", ""), colClasses = list(character = "outcome_id"))
required_e <- c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure",
                "other_allele.exposure", "pval.exposure", "exposure")
if (length(setdiff(required_e, names(exposure)))) fail("Exposure is missing TwoSampleMR fields")
if (nrow(exposure) != 5254L || uniqueN(exposure$SNP) != 5254L) fail("Exposure instrument set is not 5,254 unique SNPs")
if (anyNA(exposure[, .(beta.exposure, se.exposure, effect_allele.exposure, other_allele.exposure, pval.exposure)])) fail("Exposure has missing required association fields")
if (any(exposure$se.exposure <= 0)) fail("Exposure contains non-positive SE")

# harmonise_data() requires the column to exist. NA explicitly represents the
# genuinely unavailable asthma EAF; no frequency value is inferred or invented.
exposure[, `:=`(eaf.exposure = NA_real_, id.exposure = "asthma_gwama")]
instrument_snps <- exposure$SNP
complement <- c(A = "T", T = "A", C = "G", G = "C")
is_true <- function(x) !is.na(x) & x
same_num <- function(x, y) !is.na(x) & !is.na(y) & abs(x - y) <= pmax(1e-12, abs(y) * 1e-10)

qc_rows <- list(); direction_rows <- list(); conflict_rows <- list(); errors <- list()
completed <- logical(nrow(cfg))
conflict_ids <- c("rs12939565", "rs2055456", "rs146054809", "rs370943675",
                  "rs558418140", "rs560335706", "rs571546412")

for (i in seq_len(nrow(cfg))) {
  meta <- cfg[i]
  outcome_file <- file.path(outcome_dir, paste0(meta$output_stub, "_outcome.tsv"))
  message("Harmonising Asthma -> ", meta$outcome, " with action = 2")
  tryCatch({
    if (!file.exists(outcome_file)) fail("Corrected outcome missing: ", outcome_file)
    outcome <- fread(outcome_file, colClasses = list(character = "outcome_id"))
    required_o <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome",
                    "other_allele.outcome", "eaf.outcome", "pval.outcome", "outcome")
    if (length(setdiff(required_o, names(outcome)))) fail(meta$outcome, " outcome lacks TwoSampleMR fields")
    if (anyDuplicated(outcome$SNP)) fail(meta$outcome, " corrected outcome has duplicate SNPs")
    if (any(!outcome$SNP %in% instrument_snps)) fail(meta$outcome, " contains SNPs outside the exposure set")
    outcome[, id.outcome := as.character(meta$outcome_id)]
    outcome[, `:=`(
      beta.outcome.original = beta.outcome,
      effect_allele.outcome.original = effect_allele.outcome,
      other_allele.outcome.original = other_allele.outcome,
      eaf.outcome.original = eaf.outcome
    )]

    h <- as.data.table(TwoSampleMR::harmonise_data(as.data.frame(exposure), as.data.frame(outcome), action = 2))
    if (nrow(h) != nrow(outcome) || uniqueN(h$SNP) != nrow(h)) fail(meta$outcome, ": unexpected row loss/duplication in harmonise_data")
    required_flags <- c("palindromic", "ambiguous", "remove", "mr_keep")
    if (length(setdiff(required_flags, names(h)))) fail(meta$outcome, ": harmonisation QC flags missing")
    h[, outcome_beta_flipped := same_num(beta.outcome, -beta.outcome.original) & !same_num(beta.outcome, beta.outcome.original)]
    h[, outcome_alleles_reordered := effect_allele.outcome.original == other_allele.outcome &
                                      other_allele.outcome.original == effect_allele.outcome]
    h[, outcome_strand_complemented :=
        complement[effect_allele.outcome.original] == effect_allele.outcome &
        complement[other_allele.outcome.original] == other_allele.outcome]
    h[, exclusion_reason := fcase(
      mr_keep %in% TRUE, NA_character_,
      palindromic %in% TRUE & ambiguous %in% TRUE, "ambiguous_palindromic_missing_exposure_eaf",
      remove %in% TRUE, "allele_incompatible",
      ambiguous %in% TRUE, "ambiguous",
      is.na(beta.exposure) | is.na(beta.outcome) | is.na(se.exposure) | is.na(se.outcome), "missing_beta_or_se",
      default = "mr_keep_false_other"
    )]

    out_h <- file.path(harm_dir, paste0(meta$output_stub, "_harmonised.tsv"))
    out_ready <- file.path(ready_dir, paste0(meta$output_stub, "_mr_ready.tsv"))
    out_excluded <- file.path(excluded_dir, paste0(meta$output_stub, "_excluded.tsv"))
    fwrite(h, out_h, sep = "\t", quote = FALSE, na = "NA")
    ready <- h[mr_keep %in% TRUE]
    excluded <- h[!mr_keep %in% TRUE, .(SNP, effect_allele.exposure, other_allele.exposure,
      effect_allele.outcome, other_allele.outcome, effect_allele.outcome.original,
      other_allele.outcome.original, palindromic, ambiguous, remove, mr_keep, exclusion_reason)]
    fwrite(ready, out_ready, sep = "\t", quote = FALSE, na = "NA")
    fwrite(excluded, out_excluded, sep = "\t", quote = FALSE, na = "NA")

    if (anyDuplicated(ready$SNP) || anyNA(ready[, .(beta.exposure, beta.outcome, se.exposure, se.outcome)]) ||
        any(ready$se.exposure <= 0) || any(ready$se.outcome <= 0) || any(!ready$mr_keep) ||
        any(ready$effect_allele.exposure != ready$effect_allele.outcome)) {
      fail(meta$outcome, ": final MR-ready validation failed")
    }

    qc_rows[[i]] <- data.table(
      outcome = meta$outcome, outcome_id = as.character(meta$outcome_id), outcome_type = meta$outcome_type,
      corrected_finngen_matched_snps = nrow(outcome), snps_entering_harmonisation = nrow(outcome),
      snps_successfully_harmonised = nrow(h), palindromic_snps = sum(is_true(h$palindromic)),
      palindromic_unresolved_missing_exposure_eaf = sum(is_true(h$palindromic) & is_true(h$ambiguous)),
      ambiguous_snps = sum(is_true(h$ambiguous)), allele_incompatible_snps = sum(is_true(h$remove)),
      allele_flips = sum(is_true(h$outcome_beta_flipped)),
      strand_complements = sum(is_true(h$outcome_strand_complemented)),
      removed_snps = sum(!is_true(h$mr_keep)), final_mr_keep_snps = sum(is_true(h$mr_keep)),
      retention_percentage = sum(is_true(h$mr_keep)) / nrow(h) * 100,
      action = 2L, exposure_eaf_available = FALSE
    )

    # Direction audit: ten unchanged, ten sign-flipped, and ten palindromic
    # rows per outcome where available. Original FinnGen fields are retained.
    h[, audit_category := NA_character_]
    no_flip_idx <- which(!h$outcome_beta_flipped & !h$palindromic)
    if (length(no_flip_idx)) h[no_flip_idx[seq_len(min(10L, length(no_flip_idx)))], audit_category := "no_beta_flip"]
    flip_idx <- which(h$outcome_beta_flipped & !h$palindromic)
    if (length(flip_idx)) h[flip_idx[seq_len(min(10L, length(flip_idx)))], audit_category := "outcome_beta_flipped"]
    pal_idx <- which(h$palindromic)
    if (length(pal_idx)) h[pal_idx[seq_len(min(10L, length(pal_idx)))], audit_category := "palindromic"]
    direction_rows[[i]] <- h[!is.na(audit_category), .(
      outcome, outcome_id, audit_category, SNP,
      exposure_effect_allele = effect_allele.exposure,
      exposure_other_allele = other_allele.exposure,
      original_finngen_effect_allele = effect_allele.outcome.original,
      original_finngen_other_allele = other_allele.outcome.original,
      original_finngen_beta = beta.outcome.original,
      harmonised_effect_allele = effect_allele.outcome,
      harmonised_other_allele = other_allele.outcome,
      harmonised_beta = beta.outcome,
      beta_flipped = outcome_beta_flipped,
      strand_complemented = outcome_strand_complemented,
      palindromic, ambiguous, remove, retained = mr_keep
    )]

    present_conflicts <- h[SNP %in% conflict_ids]
    conflict_rows[[i]] <- rbindlist(list(
      present_conflicts[, .(outcome = meta$outcome, outcome_id = as.character(meta$outcome_id), SNP,
        enters_harmonisation = TRUE, retained = mr_keep %in% TRUE, removed = !mr_keep %in% TRUE,
        absent_from_corrected_outcome = FALSE, palindromic, ambiguous, remove, exclusion_reason)],
      data.table(outcome = meta$outcome, outcome_id = as.character(meta$outcome_id),
        SNP = setdiff(conflict_ids, outcome$SNP), enters_harmonisation = FALSE, retained = FALSE,
        removed = FALSE, absent_from_corrected_outcome = TRUE, palindromic = NA,
        ambiguous = NA, remove = NA, exclusion_reason = "absent_from_corrected_outcome")
    ), fill = TRUE)
    completed[i] <- TRUE
  }, error = function(e) {
    errors[[length(errors) + 1L]] <<- data.table(outcome = meta$outcome,
      outcome_id = as.character(meta$outcome_id), error = conditionMessage(e))
    warning(meta$outcome, " failed: ", conditionMessage(e), call. = FALSE)
  })
}

fwrite(rbindlist(qc_rows, fill = TRUE), qc_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(rbindlist(direction_rows, fill = TRUE), direction_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(rbindlist(conflict_rows, fill = TRUE), conflict_file, sep = "\t", quote = FALSE, na = "NA")
if (length(errors)) {
  fwrite(rbindlist(errors), error_file, sep = "\t", quote = FALSE, na = "NA")
} else {
  fwrite(data.table(outcome = character(), outcome_id = character(), error = character()), error_file, sep = "\t")
}

message("Successfully harmonised outcomes: ", sum(completed), "/12")
message("Harmonisation QC: ", qc_file)
if (!all(completed)) quit(status = 1L)
