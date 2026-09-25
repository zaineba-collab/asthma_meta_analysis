#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

# Extract the fixed 5,254 asthma instruments from local FinnGen summary
# statistics. This script does not clump, harmonise, use proxies, query an API,
# or run MR. FinnGen documentation defines the files as GRCh38 and beta/af_alt
# as referring to ALT (the effect allele).

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else ""
script_dir <- if (nzchar(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

exposure_file <- file.path(project_dir, "results/mr/finngen/exposure/asthma_exposure_5254.tsv")
config_file <- file.path(project_dir, "config/finngen_mr_outcomes.tsv")
outcome_dir <- file.path(project_dir, "results/mr/finngen/outcomes")
unmatched_dir <- file.path(outcome_dir, "unmatched")
qc_dir <- file.path(project_dir, "results/mr/finngen/qc")
audit_file <- file.path(qc_dir, "finngen_file_structure_audit.tsv")
qc_file <- file.path(qc_dir, "finngen_outcome_extraction_qc.tsv")
abnormal_file <- file.path(qc_dir, "finngen_outcome_abnormal_rows.tsv")

dir.create(unmatched_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

fail <- function(...) stop(..., call. = FALSE)
if (!file.exists(exposure_file)) fail("Exposure file missing: ", exposure_file)
if (!file.exists(config_file)) fail("FinnGen outcome config missing: ", config_file)
if (Sys.which("gzip") == "") fail("gzip executable is required but unavailable")

cfg <- fread(config_file, na.strings = c("NA", ""), colClasses = list(character = "outcome_id"))
required_cfg <- c("output_stub", "outcome", "outcome_id", "release", "outcome_type",
                  "file", "cases", "controls", "sample_size")
if (nrow(cfg) != 12L) fail("Config must contain exactly 12 outcomes; found ", nrow(cfg))
if (length(setdiff(required_cfg, names(cfg)))) fail("Config is missing required columns")
if (anyDuplicated(cfg$output_stub) || anyDuplicated(cfg$outcome_id)) fail("Config contains duplicate output stubs or outcome IDs")

cfg[, absolute_file := file.path(project_dir, file)]
missing_gz <- cfg$file[!file.exists(cfg$absolute_file)]
missing_tbi_idx <- !file.exists(paste0(cfg$absolute_file, ".tbi"))
missing_tbi <- if (any(missing_tbi_idx)) paste0(cfg$file[missing_tbi_idx], ".tbi") else character()
if (length(missing_gz) || length(missing_tbi)) {
  fail("Missing required FinnGen input(s): ", paste(c(missing_gz, missing_tbi), collapse = ", "))
}

message("Checking gzip integrity for all 12 FinnGen files...")
bad_gzip <- character()
for (f in cfg$absolute_file) {
  status <- system2("gzip", c("-t", shQuote(f)), stdout = FALSE, stderr = FALSE)
  if (!identical(status, 0L)) bad_gzip <- c(bad_gzip, f)
}
if (length(bad_gzip)) fail("Gzip integrity check failed: ", paste(bad_gzip, collapse = ", "))

exposure <- fread(exposure_file)
required_exposure <- c("SNP", "chromosome", "position", "effect_allele.exposure", "other_allele.exposure")
if (length(setdiff(required_exposure, names(exposure)))) fail("Exposure file lacks required SNP/position/allele columns")
if (nrow(exposure) != 5254L || uniqueN(exposure$SNP) != 5254L) {
  fail("Exposure must contain exactly 5,254 rows and 5,254 unique SNPs")
}
if (anyNA(exposure[, ..required_exposure])) fail("Exposure has missing SNP, position, or allele fields")
exposure[, chromosome := as.character(chromosome)]
exposure[, position := as.integer(position)]
exposure[, `:=`(effect_allele.exposure = toupper(effect_allele.exposure),
                other_allele.exposure = toupper(other_allele.exposure))]

read_header <- function(f) {
  con <- gzfile(f, open = "rt")
  on.exit(close(con))
  sub("^#", "", readLines(con, n = 1L, warn = FALSE))
}

header_rows <- vector("list", nrow(cfg))
headers <- vector("list", nrow(cfg))
for (i in seq_len(nrow(cfg))) {
  h <- strsplit(read_header(cfg$absolute_file[i]), "\t", fixed = TRUE)[[1]]
  headers[[i]] <- h
  needed <- c("chrom", "pos", "ref", "alt", "beta", "sebeta", "pval", "af_alt")
  missing_cols <- setdiff(needed, h)
  if (length(missing_cols)) fail(cfg$outcome[i], " is missing columns: ", paste(missing_cols, collapse = ", "))
  header_rows[[i]] <- data.table(
    outcome = cfg$outcome[i], outcome_id = cfg$outcome_id[i], release = cfg$release[i],
    file = cfg$file[i], header = paste(h, collapse = " | "),
    chromosome_column = "#chrom", position_column = "pos", reference_allele_column = "ref",
    alternate_allele_column = "alt", rsid_column = if ("rsids" %in% h) "rsids" else NA_character_,
    beta_column = "beta", se_column = "sebeta", pvalue_column = "pval",
    allele_frequency_column = "af_alt", genome_build = "GRCh38",
    beta_effect_allele = "ALT", has_rsids = "rsids" %in% h,
    extraction_match = if ("rsids" %in% h) "exact rsID token + exact exposure allele pair; FinnGen row supplies GRCh38 coordinate" else "FinnGen-R13 rsID-confirmed GRCh38 coordinate + exact REF/ALT (R12 has no rsids column)"
  )
}
fwrite(rbindlist(header_rows), audit_file, sep = "\t", quote = FALSE, na = "NA")

# Neither tabix nor an installed HTS-capable R library is available in the
# current environment. The PLINK clumping table positions come from its BIM
# and are not GRCh38 for many instruments, so R13 is streamed by exact rsID.
# An rsID-confirmed R13 variant map then supplies GRCh38 coordinates/alleles
# for the R12 files, which omit rsids entirely.
rsid_file <- tempfile("asthma_rsids_", fileext = ".txt")
writeLines(exposure$SNP, rsid_file)
on.exit(unlink(rsid_file), add = TRUE)

shell_quote <- function(x) shQuote(x)
stream_rsids <- function(gwas_file, header) {
  awk <- paste0(
    "awk 'BEGIN{FS=OFS=\"\\t\"} NR==FNR{want[$1]=1;next} FNR==1{next} ",
    "{n=split($5,a,/[,;|[:space:]]+/); for(i=1;i<=n;i++) if(a[i] in want){print;break}}' ",
    shell_quote(rsid_file), " <(gzip -cd ", shell_quote(gwas_file), ")"
  )
  cmd <- paste("bash -c", shell_quote(awk))
  fread(cmd = cmd, header = FALSE, col.names = header, showProgress = FALSE)
}

stream_positions <- function(gwas_file, header, region_file) {
  awk <- paste0(
    "awk 'BEGIN{FS=OFS=\"\\t\"} NR==FNR{k[$1 SUBSEP $2]=1;next} ",
    "FNR==1{next} (($1 SUBSEP $2) in k){print}' ",
    shell_quote(region_file), " <(gzip -cd ", shell_quote(gwas_file), ")"
  )
  cmd <- paste("bash -c", shell_quote(awk))
  fread(cmd = cmd, header = FALSE, col.names = header, showProgress = FALSE)
}

rsid_tokens <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  unique(strsplit(x, "[,;|[:space:]]+", perl = TRUE)[[1]])
}

# BMI is used only as a local R13 variant dictionary: exact rsID -> GRCh38
# chromosome, position, REF and ALT. Association values still come from each
# outcome's own file.
map_i <- which(cfg$outcome_id == "BMI_IRN")
map_raw <- stream_rsids(cfg$absolute_file[map_i], headers[[map_i]])
map_raw[, `:=`(chrom = as.character(chrom), pos = as.integer(pos), ref = toupper(ref), alt = toupper(alt), row_id = .I)]
map_hits <- map_raw[, .(SNP = unlist(lapply(rsids, rsid_tokens))), by = row_id]
map_hits <- map_hits[nzchar(SNP) & SNP %in% exposure$SNP]
r13_variant_candidates <- merge(map_hits, map_raw[, .(row_id, chrom, pos, ref, alt)], by = "row_id", allow.cartesian = TRUE)
r13_variant_candidates <- merge(r13_variant_candidates,
                                exposure[, .(SNP, exp_ea = effect_allele.exposure, exp_oa = other_allele.exposure)],
                                by = "SNP", allow.cartesian = TRUE)
map_conflicts <- r13_variant_candidates[!((ref == exp_ea & alt == exp_oa) | (ref == exp_oa & alt == exp_ea))]
r13_variant_map <- r13_variant_candidates[(ref == exp_ea & alt == exp_oa) | (ref == exp_oa & alt == exp_ea)]
r13_variant_map <- unique(r13_variant_map[, .(SNP, chrom, pos, ref, alt)])
map_dups <- r13_variant_map[, .N, by = SNP][N > 1L]
if (nrow(map_dups)) fail("R13 rsID dictionary has multiple variants for: ", paste(map_dups$SNP, collapse = ", "))
r12_region_file <- tempfile("asthma_r13_grch38_regions_", fileext = ".tsv")
fwrite(unique(r13_variant_map[, .(chrom, pos)]), r12_region_file, sep = "\t", col.names = FALSE)
on.exit(unlink(r12_region_file), add = TRUE)

qc_rows <- list()
abnormal_rows <- if (nrow(map_conflicts)) list(map_conflicts[, .(
  outcome = "R13 variant dictionary", outcome_id = "BMI_IRN", SNP,
  issue = "requested_rsid_found_but_allele_pair_differs",
  chromosome = chrom, position = pos, ref, alt,
  exposure_effect_allele = exp_ea, exposure_other_allele = exp_oa
)]) else list()
completed <- logical(nrow(cfg))

for (i in seq_len(nrow(cfg))) {
  meta <- cfg[i]
  message("Extracting ", meta$outcome, " from ", meta$file)
  tryCatch({
    raw <- if ("rsids" %in% headers[[i]]) {
      stream_rsids(meta$absolute_file, headers[[i]])
    } else {
      stream_positions(meta$absolute_file, headers[[i]], r12_region_file)
    }
    raw[, `:=`(chrom = as.character(chrom), pos = as.integer(pos),
               ref = toupper(ref), alt = toupper(alt))]
    raw[, row_id := .I]

    if ("rsids" %in% names(raw)) {
      hits <- raw[, .(SNP = unlist(lapply(rsids, rsid_tokens))), by = row_id]
      hits <- hits[nzchar(SNP) & SNP %in% exposure$SNP]
      candidates <- merge(hits, raw, by = "row_id", allow.cartesian = TRUE)
      candidates <- merge(candidates,
                          exposure[, .(SNP, exp_ea = effect_allele.exposure, exp_oa = other_allele.exposure)],
                          by = "SNP", allow.cartesian = TRUE)
      candidates <- candidates[(ref == exp_ea & alt == exp_oa) | (ref == exp_oa & alt == exp_ea)]
    } else {
      candidates <- merge(
        raw,
        r13_variant_map,
        by = c("chrom", "pos", "ref", "alt"), allow.cartesian = TRUE
      )
    }

    match_counts <- candidates[, .N, by = SNP]
    ambiguous <- match_counts[N > 1L, SNP]
    duplicate_matches <- length(ambiguous)
    if (length(ambiguous)) {
      abnormal_rows[[length(abnormal_rows) + 1L]] <- candidates[SNP %in% ambiguous,
        .(outcome = meta$outcome, outcome_id = meta$outcome_id, SNP,
          issue = "multiple_valid_finngen_rows", chromosome = chrom, position = pos,
          ref, alt, beta, sebeta, pval)]
    }
    matched <- candidates[!SNP %in% ambiguous]

    result <- matched[, .(
      SNP, chromosome.outcome = chrom, position.outcome = pos,
      beta.outcome = as.numeric(beta), se.outcome = as.numeric(sebeta),
      effect_allele.outcome = alt, other_allele.outcome = ref,
      eaf.outcome = as.numeric(af_alt), pval.outcome = as.numeric(pval),
      outcome = meta$outcome, outcome_id = meta$outcome_id,
      outcome_type = meta$outcome_type, finngen_release = meta$release
    )]
    if (meta$outcome_type == "binary") {
      result[, `:=`(ncase.outcome = as.integer(meta$cases),
                    ncontrol.outcome = as.integer(meta$controls),
                    samplesize.outcome = as.integer(meta$sample_size))]
    } else {
      result[, samplesize.outcome := as.integer(meta$sample_size)]
    }

    issues <- rbindlist(list(
      result[is.na(beta.outcome), .(outcome, outcome_id, SNP, issue = "missing_beta", chromosome = chromosome.outcome, position = position.outcome, ref = other_allele.outcome, alt = effect_allele.outcome, beta = beta.outcome, sebeta = se.outcome, pval = pval.outcome)],
      result[is.na(se.outcome) | se.outcome <= 0, .(outcome, outcome_id, SNP, issue = "missing_or_nonpositive_se", chromosome = chromosome.outcome, position = position.outcome, ref = other_allele.outcome, alt = effect_allele.outcome, beta = beta.outcome, sebeta = se.outcome, pval = pval.outcome)],
      result[is.na(pval.outcome) | pval.outcome < 0 | pval.outcome > 1, .(outcome, outcome_id, SNP, issue = "invalid_pvalue", chromosome = chromosome.outcome, position = position.outcome, ref = other_allele.outcome, alt = effect_allele.outcome, beta = beta.outcome, sebeta = se.outcome, pval = pval.outcome)],
      result[!grepl("^[ACGT]+$", effect_allele.outcome) | !grepl("^[ACGT]+$", other_allele.outcome), .(outcome, outcome_id, SNP, issue = "invalid_allele", chromosome = chromosome.outcome, position = position.outcome, ref = other_allele.outcome, alt = effect_allele.outcome, beta = beta.outcome, sebeta = se.outcome, pval = pval.outcome)]
    ), fill = TRUE)
    if (nrow(issues)) abnormal_rows[[length(abnormal_rows) + 1L]] <- issues
    if (anyDuplicated(result$SNP)) fail(meta$outcome, ": unexpected duplicated SNP rows after validation")

    result[, input_order__ := match(SNP, exposure$SNP)]
    setorder(result, input_order__)
    result[, input_order__ := NULL]
    unmatched <- exposure[!result, on = "SNP", .(SNP, chromosome, position)]
    fwrite(result, file.path(outcome_dir, paste0(meta$output_stub, "_outcome.tsv")), sep = "\t", quote = FALSE, na = "NA")
    fwrite(unmatched, file.path(unmatched_dir, paste0(meta$output_stub, "_unmatched_snps.tsv")), sep = "\t", quote = FALSE, na = "NA")

    qc_rows[[i]] <- data.table(
      outcome = meta$outcome, outcome_id = meta$outcome_id, release = meta$release,
      outcome_type = meta$outcome_type, requested_snps = 5254L, matched_snps = nrow(result),
      unmatched_snps = 5254L - nrow(result), coverage_percent = nrow(result) / 5254 * 100,
      duplicate_matches = duplicate_matches, missing_beta = sum(is.na(result$beta.outcome)),
      missing_se = sum(is.na(result$se.outcome)),
      missing_effect_allele = sum(is.na(result$effect_allele.outcome) | result$effect_allele.outcome == ""),
      missing_other_allele = sum(is.na(result$other_allele.outcome) | result$other_allele.outcome == ""),
      missing_pvalue = sum(is.na(result$pval.outcome)), missing_eaf = sum(is.na(result$eaf.outcome)),
      cases = as.integer(meta$cases), controls = as.integer(meta$controls), sample_size = as.integer(meta$sample_size),
      extraction_method = "single-pass gzip streaming fallback"
    )
    completed[i] <- TRUE
  }, error = function(e) {
    warning(meta$outcome, " failed: ", conditionMessage(e), call. = FALSE)
    qc_rows[[i]] <<- data.table(
      outcome = meta$outcome, outcome_id = meta$outcome_id, release = meta$release,
      outcome_type = meta$outcome_type, requested_snps = 5254L, matched_snps = NA_integer_,
      unmatched_snps = NA_integer_, coverage_percent = NA_real_, duplicate_matches = NA_integer_,
      missing_beta = NA_integer_, missing_se = NA_integer_, missing_effect_allele = NA_integer_,
      missing_other_allele = NA_integer_, missing_pvalue = NA_integer_, missing_eaf = NA_integer_,
      cases = as.integer(meta$cases), controls = as.integer(meta$controls), sample_size = as.integer(meta$sample_size),
      extraction_method = paste0("FAILED: ", conditionMessage(e))
    )
  })
}

fwrite(rbindlist(qc_rows, fill = TRUE), qc_file, sep = "\t", quote = FALSE, na = "NA")
if (length(abnormal_rows)) {
  fwrite(rbindlist(abnormal_rows, fill = TRUE), abnormal_file, sep = "\t", quote = FALSE, na = "NA")
} else {
  fwrite(data.table(outcome = character(), outcome_id = character(), SNP = character(), issue = character()),
         abnormal_file, sep = "\t", quote = FALSE, na = "NA")
}

message("Completed outcomes: ", sum(completed), "/", length(completed))
message("Structure audit: ", audit_file)
message("Extraction QC: ", qc_file)
if (!all(completed)) quit(status = 1L)
