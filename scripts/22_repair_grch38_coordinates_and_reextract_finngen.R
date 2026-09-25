#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

# Repair exposure coordinates without changing the fixed instrument set or
# association statistics, then re-extract local FinnGen outcomes. No clumping,
# proxies, harmonisation, API queries, or MR are performed.

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(arg)) sub("^--file=", "", arg[1]) else ""
script_dir <- if (nzchar(script_file)) dirname(normalizePath(script_file)) else getwd()
root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
p <- function(...) file.path(root, ...)
fail <- function(...) stop(..., call. = FALSE)

exposure_file <- p("results", "mr", "finngen", "exposure", "asthma_exposure_5254.tsv")
corrected_exposure_file <- p("results", "mr", "finngen", "exposure", "asthma_exposure_5254_grch38.tsv")
lead_file <- p("results", "ld_clumping", "independent_lead_snps.tsv")
clump_map_file <- p("results", "gwama", "asthma_fixed_for_plink_clumping.txt")
config_file <- p("config", "finngen_mr_outcomes.tsv")
old_qc_file <- p("results", "mr", "finngen", "qc", "finngen_outcome_extraction_qc.tsv")
out_dir <- p("results", "mr", "finngen", "outcomes_grch38")
unmatched_dir <- file.path(out_dir, "unmatched")
qc_dir <- p("results", "mr", "finngen", "qc")
dir.create(unmatched_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

audit_file <- file.path(qc_dir, "exposure_coordinate_build_audit.tsv")
unresolved_file <- file.path(qc_dir, "exposure_grch38_unresolved_snps.tsv")
coordinate_qc_file <- file.path(qc_dir, "exposure_coordinate_repair_qc.tsv")
new_qc_file <- file.path(qc_dir, "finngen_outcome_extraction_grch38_qc.tsv")
comparison_file <- file.path(qc_dir, "extraction_coverage_comparison.tsv")
conflict_file <- file.path(qc_dir, "conflicting_variant_resolution.tsv")
abnormal_file <- file.path(qc_dir, "finngen_outcome_extraction_grch38_abnormal.tsv")

required <- c(exposure_file, lead_file, clump_map_file, config_file, old_qc_file)
if (any(!file.exists(required))) fail("Missing required input(s): ", paste(required[!file.exists(required)], collapse = ", "))
if (Sys.which("gzip") == "") fail("gzip is required")

exposure <- fread(exposure_file)
lead <- fread(lead_file)
clump_map <- fread(clump_map_file)
cfg <- fread(config_file, na.strings = c("NA", ""), colClasses = list(character = "outcome_id"))
if (nrow(exposure) != 5254L || uniqueN(exposure$SNP) != 5254L) fail("Original exposure is not 5,254 unique SNPs")
if (!identical(exposure$SNP, lead$SNP)) fail("Exposure and authoritative PLINK SNP order/set differ")
exposure[, `:=`(effect_allele.exposure = toupper(effect_allele.exposure), other_allele.exposure = toupper(other_allele.exposure))]
instrument_snps <- exposure$SNP

rsid_file <- tempfile("asthma_rsids_", fileext = ".txt")
writeLines(exposure$SNP, rsid_file)
on.exit(unlink(rsid_file), add = TRUE)
sq <- function(x) shQuote(x)
tokens <- function(x) {
  if (is.na(x) || !nzchar(x) || x == ".") return(character())
  unique(strsplit(x, "[,;|[:space:]]+", perl = TRUE)[[1]])
}
allele_ok <- function(ref, alt, ea, oa) (ref == ea & alt == oa) | (ref == oa & alt == ea)

stream_final_rsid <- function(f) {
  cmd0 <- paste0("awk 'BEGIN{FS=OFS=\"\\t\"} NR==FNR{w[$1]=1;next} FNR==1{next} ($1 in w){print $1,$8,$9,$2,$3}' ",
                 sq(rsid_file), " ", sq(f))
  x <- fread(cmd = cmd0, header = FALSE, col.names = c("SNP", "chrom", "pos", "a1", "a2"), showProgress = FALSE)
  x[, `:=`(chrom = as.character(chrom), pos = as.integer(pos), a1 = toupper(a1), a2 = toupper(a2))]
  x
}

read_gz_header <- function(f) {
  con <- gzfile(f, "rt"); on.exit(close(con))
  strsplit(sub("^#", "", readLines(con, n = 1L, warn = FALSE)), "\t", fixed = TRUE)[[1]]
}
stream_finngen_rsids <- function(f, header) {
  awk <- paste0("awk 'BEGIN{FS=OFS=\"\\t\"} NR==FNR{w[$1]=1;next} FNR==1{next} ",
                "{n=split($5,a,/[,;|[:space:]]+/);for(i=1;i<=n;i++)if(a[i] in w){print;break}}' ",
                sq(rsid_file), " <(gzip -cd ", sq(f), ")")
  fread(cmd = paste("bash -c", sq(awk)), header = FALSE, col.names = header, showProgress = FALSE)
}
stream_finngen_targets <- function(f, header, region_file, has_rsids) {
  if (has_rsids) {
    awk <- paste0("awk 'BEGIN{FS=OFS=\"\\t\"} NR==FNR{if(NF==1)w[$1]=1;else k[$1 SUBSEP $2]=1;next} FNR==1{next} ",
                  "{hit=(($1 SUBSEP $2) in k);n=split($5,a,/[,;|[:space:]]+/);for(i=1;i<=n;i++)if(a[i] in w)hit=1;if(hit)print}' ",
                  sq(region_file), " <(gzip -cd ", sq(f), ")")
  } else {
    awk <- paste0("awk 'BEGIN{FS=OFS=\"\\t\"} NR==FNR{k[$1 SUBSEP $2]=1;next} FNR==1{next} (($1 SUBSEP $2) in k){print}' ",
                  sq(region_file), " <(gzip -cd ", sq(f), ")")
  }
  fread(cmd = paste("bash -c", sq(awk)), header = FALSE, col.names = header, showProgress = FALSE)
}

# Project-first mapping sources. These two files were explicitly lifted to
# GRCh38 by the current asthma pipeline; BUILD=unknown files are not used.
trusted_files <- c(
  p("data", "final_hg38", "GCST90029018_buildGRCh38.lifted.gwama_ready.tsv"),
  p("data", "final_hg38", "Shrine_30552067_moderate-severe_asthma.hg38.lifted.gwama_ready.tsv")
)
if (any(!file.exists(trusted_files))) fail("Explicitly lifted mapping source missing")
lifted <- rbindlist(lapply(trusted_files, stream_final_rsid), fill = TRUE)
lifted <- merge(lifted, exposure[, .(SNP, ea = effect_allele.exposure, oa = other_allele.exposure)], by = "SNP")
lifted <- unique(lifted[allele_ok(a1, a2, ea, oa), .(SNP, chrom, pos, a1, a2, coordinate_source = "current_asthma_explicit_liftover")])

# Coordinate-style GWAMA IDs were generated from post-harmonisation variant
# identifiers. They recover variants whose GWAMA key was not an rsID.
cm <- merge(exposure[, .(SNP)], clump_map[, .(SNP, original_SNP)], by = "SNP", all.x = TRUE, sort = FALSE)
coord <- cm[grepl(":", original_SNP)]
parts <- tstrsplit(coord$original_SNP, ":", fixed = TRUE)
coord[, `:=`(chrom = as.character(parts[[1]]), pos = as.integer(parts[[2]]),
             a1 = toupper(parts[[3]]), a2 = toupper(parts[[4]]), coordinate_source = "current_gwama_coordinate_id")]
coord <- coord[, .(SNP, chrom, pos, a1, a2, coordinate_source)]

# Prefer explicitly lifted evidence when present; otherwise use a unique GWAMA
# coordinate ID. FinnGen exact rsID+allele is only a gap-filling validator.
pipeline <- rbindlist(list(lifted, coord), fill = TRUE)
chosen_lifted <- lifted[, .SD[1], by = SNP]
chosen_coord <- coord[!chosen_lifted, on = "SNP"][, .SD[1], by = SNP]
chosen <- rbindlist(list(chosen_lifted, chosen_coord), fill = TRUE)

bmi_i <- which(cfg$outcome_id == "BMI_IRN")
bmi_header <- read_gz_header(p(cfg$file[bmi_i]))
bmi_raw <- stream_finngen_rsids(p(cfg$file[bmi_i]), bmi_header)
bmi_raw[, `:=`(chrom = as.character(chrom), pos = as.integer(pos), ref = toupper(ref), alt = toupper(alt), row_id = .I)]
bmi_hits <- bmi_raw[, .(SNP = unlist(lapply(rsids, tokens))), by = row_id][SNP %in% instrument_snps]
bmi_map <- merge(bmi_hits, bmi_raw[, .(row_id, chrom, pos, ref, alt)], by = "row_id", allow.cartesian = TRUE)
bmi_map <- merge(bmi_map, exposure[, .(SNP, ea = effect_allele.exposure, oa = other_allele.exposure)], by = "SNP")
bmi_map <- unique(bmi_map[allele_ok(ref, alt, ea, oa), .(SNP, chrom, pos, a1 = ref, a2 = alt, coordinate_source = "FinnGen_R13_exact_rsid_alleles")])
bmi_map <- bmi_map[, .SD[1], by = SNP]
# Validate every overlapping project mapping against FinnGen. Exact rsID plus
# compatible alleles is definitive for FinnGen's documented GRCh38 build and
# overrides discordant coordinate-style GWAMA IDs (which can originate from
# inputs whose BUILD remained unknown). Agreeing project mappings retain their
# project provenance; FinnGen also fills otherwise unresolved gaps.
overlap <- merge(chosen[, .(SNP, project_chrom = chrom, project_pos = pos)],
                 bmi_map[, .(SNP, fg_chrom = chrom, fg_pos = pos)], by = "SNP")
discordant_ids <- overlap[project_chrom != fg_chrom | project_pos != fg_pos, SNP]
chosen <- chosen[!SNP %in% discordant_ids]
supplement <- bmi_map[!chosen, on = "SNP"]
supplement[SNP %in% discordant_ids, coordinate_source := "FinnGen_R13_exact_rsid_alleles_overrode_discordant_unknown_build_GWAMA_ID"]
chosen <- rbindlist(list(chosen, supplement), fill = TRUE)

chosen_counts <- chosen[, .N, by = SNP]
if (any(chosen_counts$N > 1L)) fail("Coordinate selection produced duplicate SNP mappings")
corrected <- merge(exposure, chosen, by = "SNP", all.x = TRUE, sort = FALSE)
corrected[, input_order__ := match(SNP, instrument_snps)]
setorder(corrected, input_order__)
corrected[, `:=`(
  chromosome.plink_reference = chromosome,
  position.plink_reference = position,
  chromosome.grch38 = chrom,
  position.grch38 = pos
)]
corrected[, c("chrom", "pos", "a1", "a2", "input_order__") := NULL]
setcolorder(corrected, c("SNP", "chromosome.plink_reference", "position.plink_reference",
                        "chromosome.grch38", "position.grch38", "coordinate_source",
                        setdiff(names(corrected), c("SNP", "chromosome.plink_reference", "position.plink_reference",
                                                    "chromosome.grch38", "position.grch38", "coordinate_source",
                                                    "chromosome", "position"))))
corrected[, c("chromosome", "position") := NULL]
fwrite(corrected, corrected_exposure_file, sep = "\t", quote = FALSE, na = "NA")
unresolved <- corrected[is.na(position.grch38), .(SNP, chromosome.plink_reference, position.plink_reference,
                                                  effect_allele.exposure, other_allele.exposure)]
fwrite(unresolved, unresolved_file, sep = "\t", quote = FALSE, na = "NA")

# Coordinate build audit: top 20 plus deterministic random 20, retaining
# FinnGen position where exact rsID+alleles were available.
set.seed(220826)
audit_ids <- unique(c(exposure$SNP[1:20], sample(exposure$SNP, 20)))
audit <- corrected[SNP %in% audit_ids, .(SNP, chromosome.plink_reference, position.plink_reference,
                                        chromosome.grch38, position.grch38, coordinate_source,
                                        audit_group = fifelse(match(SNP, instrument_snps) <= 20, "top_association", "deterministic_random"))]
audit <- merge(audit, bmi_map[, .(SNP, finngen_r13_chromosome = chrom, finngen_r13_position = pos)], by = "SNP", all.x = TRUE)
audit[, `:=`(old_equals_new = chromosome.plink_reference == chromosome.grch38 & position.plink_reference == position.grch38,
             new_equals_finngen = chromosome.grch38 == finngen_r13_chromosome & position.grch38 == finngen_r13_position)]
setorder(audit, audit_group, SNP)
fwrite(audit, audit_file, sep = "\t", quote = FALSE, na = "NA")

resolved <- corrected[!is.na(position.grch38)]
coord_qc <- data.table(
  metric = c("coordinate_provenance", "inferred_PLINK_build", "total_instruments", "GRCh38_coordinates_recovered",
             "coordinates_unresolved", "old_position_equals_new", "old_position_differs_new", "chromosome_changed",
             "duplicate_coordinate_allele_combinations", "rs9273374_GRCh38"),
  value = c("PLINK .clumps #CHROM/POS copied from 1000G_EUR_biallelic.bim; repaired from current explicit liftover/GWAMA IDs, then FinnGen exact-rsID gaps",
            "GRCh37/hg19", nrow(corrected), nrow(resolved), nrow(unresolved),
            resolved[chromosome.plink_reference == chromosome.grch38 & position.plink_reference == position.grch38, .N],
            resolved[chromosome.plink_reference != chromosome.grch38 | position.plink_reference != position.grch38, .N],
            resolved[chromosome.plink_reference != chromosome.grch38, .N],
            resolved[, .N, by = .(chromosome.grch38, position.grch38, effect_allele.exposure, other_allele.exposure)][N > 1, .N],
            paste0(corrected[SNP == "rs9273374", chromosome.grch38], ":", corrected[SNP == "rs9273374", position.grch38]))
)
fwrite(coord_qc, coordinate_qc_file, sep = "\t", quote = FALSE, na = "NA")

# One regions file contains corrected coordinates and rsIDs. R13 awk uses both;
# R12 ignores one-field rsID rows and queries all corrected positions.
target_file <- tempfile("asthma_grch38_targets_", fileext = ".tsv")
fwrite(resolved[, .(chromosome.grch38, position.grch38)], target_file, sep = "\t", col.names = FALSE)
write(exposure$SNP, file = target_file, append = TRUE)
on.exit(unlink(target_file), add = TRUE)

qc_rows <- list(); abnormal <- list(); completed <- logical(nrow(cfg))
for (i in seq_len(nrow(cfg))) {
  meta <- cfg[i]
  message("Corrected extraction: ", meta$outcome)
  tryCatch({
    f <- p(meta$file); h <- read_gz_header(f); has_rsids <- "rsids" %in% h
    raw <- stream_finngen_targets(f, h, target_file, has_rsids)
    raw[, `:=`(chrom = as.character(chrom), pos = as.integer(pos), ref = toupper(ref), alt = toupper(alt), row_id = .I)]
    map <- resolved[, .(SNP, chrom = chromosome.grch38, pos = position.grch38,
                        ea = effect_allele.exposure, oa = other_allele.exposure)]
    rsid_matches <- data.table()
    if (has_rsids) {
      hit <- raw[, .(SNP = unlist(lapply(rsids, tokens))), by = row_id][SNP %in% instrument_snps]
      rsid_matches <- merge(hit, raw, by = "row_id", allow.cartesian = TRUE)
      rsid_matches <- merge(rsid_matches, map[, .(SNP, ea, oa)], by = "SNP")
      rsid_matches <- rsid_matches[allele_ok(ref, alt, ea, oa)]
      rsid_matches[, match_method := "rsid"]
    }
    already <- unique(rsid_matches$SNP)
    coord_matches <- merge(raw, map[!SNP %in% already], by = c("chrom", "pos"), allow.cartesian = TRUE)
    coord_matches <- coord_matches[allele_ok(ref, alt, ea, oa)]
    coord_matches[, match_method := "grch38_position_alleles"]
    candidates <- rbindlist(list(rsid_matches, coord_matches), fill = TRUE)
    counts <- candidates[, .N, by = SNP]
    ambiguous <- counts[N > 1L, SNP]
    if (length(ambiguous)) abnormal[[length(abnormal) + 1L]] <- candidates[SNP %in% ambiguous,
      .(outcome = meta$outcome, outcome_id = meta$outcome_id, SNP, issue = "multiple_valid_rows", chrom, pos, ref, alt, match_method)]
    matched <- candidates[!SNP %in% ambiguous]
    result <- matched[, .(SNP, chromosome.outcome = chrom, position.outcome = pos,
      beta.outcome = as.numeric(beta), se.outcome = as.numeric(sebeta),
      effect_allele.outcome = alt, other_allele.outcome = ref,
      eaf.outcome = as.numeric(af_alt), pval.outcome = as.numeric(pval),
      outcome = meta$outcome, outcome_id = meta$outcome_id, outcome_type = meta$outcome_type,
      finngen_release = meta$release, match_method)]
    if (meta$outcome_type == "binary") result[, `:=`(ncase.outcome = as.integer(meta$cases), ncontrol.outcome = as.integer(meta$controls), samplesize.outcome = as.integer(meta$sample_size))]
    else result[, samplesize.outcome := as.integer(meta$sample_size)]
    result[, order__ := match(SNP, instrument_snps)]; setorder(result, order__); result[, order__ := NULL]
    bad <- result[is.na(beta.outcome) | is.na(se.outcome) | se.outcome <= 0 | is.na(pval.outcome) |
                    pval.outcome < 0 | pval.outcome > 1 | !grepl("^[ACGT]+$", effect_allele.outcome) |
                    !grepl("^[ACGT]+$", other_allele.outcome)]
    if (nrow(bad)) abnormal[[length(abnormal) + 1L]] <- bad[, .(outcome, outcome_id, SNP, issue = "invalid_association_row")]
    if (nrow(bad) || anyDuplicated(result$SNP) || any(!result$SNP %in% instrument_snps)) fail(meta$outcome, ": corrected validation failed")
    um <- exposure[!result, on = "SNP", .(SNP)]
    fwrite(result, file.path(out_dir, paste0(meta$output_stub, "_outcome.tsv")), sep = "\t", quote = FALSE, na = "NA")
    fwrite(um, file.path(unmatched_dir, paste0(meta$output_stub, "_unmatched_snps.tsv")), sep = "\t", quote = FALSE, na = "NA")
    qc_rows[[i]] <- data.table(outcome = meta$outcome, outcome_id = meta$outcome_id, release = meta$release,
      outcome_type = meta$outcome_type, requested_snps = 5254L, matched_snps = nrow(result), unmatched_snps = 5254L - nrow(result),
      coverage_percent = nrow(result) / 5254 * 100, duplicate_matches = length(ambiguous),
      rsid_matches = result[match_method == "rsid", .N], coordinate_allele_matches = result[match_method == "grch38_position_alleles", .N],
      missing_beta = sum(is.na(result$beta.outcome)), missing_se = sum(is.na(result$se.outcome)),
      missing_effect_allele = sum(is.na(result$effect_allele.outcome)), missing_other_allele = sum(is.na(result$other_allele.outcome)),
      missing_pvalue = sum(is.na(result$pval.outcome)), missing_eaf = sum(is.na(result$eaf.outcome)),
      cases = as.integer(meta$cases), controls = as.integer(meta$controls), sample_size = as.integer(meta$sample_size))
    completed[i] <- TRUE
  }, error = function(e) {
    warning(meta$outcome, " failed: ", conditionMessage(e), call. = FALSE)
    qc_rows[[i]] <<- data.table(outcome = meta$outcome, outcome_id = meta$outcome_id, release = meta$release,
      outcome_type = meta$outcome_type, requested_snps = 5254L, error = conditionMessage(e))
  })
}
new_qc <- rbindlist(qc_rows, fill = TRUE)
fwrite(new_qc, new_qc_file, sep = "\t", quote = FALSE, na = "NA")
if (length(abnormal)) {
  fwrite(rbindlist(abnormal, fill = TRUE), abnormal_file, sep = "\t", quote = FALSE, na = "NA")
} else {
  fwrite(data.table(outcome = character(), outcome_id = character(), SNP = character(), issue = character()), abnormal_file, sep = "\t")
}

old_qc <- fread(old_qc_file, colClasses = list(character = "outcome_id"))
comparison <- merge(old_qc[, .(outcome_id, outcome, original_matched = matched_snps, original_coverage_percent = coverage_percent)],
                    new_qc[, .(outcome_id, corrected_matched = matched_snps, corrected_coverage_percent = coverage_percent,
                               rsid_matches, coordinate_allele_matches)], by = "outcome_id", all = TRUE)
comparison[, `:=`(gain_in_snps = corrected_matched - original_matched,
                  coverage_percentage_point_change = corrected_coverage_percent - original_coverage_percent)]
setcolorder(comparison, c("outcome", "outcome_id", "original_matched", "corrected_matched", "gain_in_snps",
                         "original_coverage_percent", "corrected_coverage_percent", "coverage_percentage_point_change",
                         "rsid_matches", "coordinate_allele_matches"))
fwrite(comparison, comparison_file, sep = "\t", quote = FALSE, na = "NA")

conflict_ids <- c("rs146054809", "rs370943675", "rs558418140", "rs560335706", "rs571546412", "rs12939565", "rs2055456")
conf <- corrected[SNP %in% conflict_ids, .(SNP, chromosome.grch38, position.grch38,
  exposure_effect_allele = effect_allele.exposure, exposure_other_allele = other_allele.exposure, coordinate_source)]
raw_conf <- bmi_raw[, .(row_id, finngen_chromosome = chrom, finngen_position = pos, ref, alt, rsids)]
raw_conf <- raw_conf[grepl(paste(conflict_ids, collapse = "|"), rsids)]
conf_rows <- raw_conf[, .(SNP = unlist(lapply(rsids, tokens))), by = .(row_id, finngen_chromosome, finngen_position, ref, alt)][SNP %in% conflict_ids]
conf <- merge(conf, conf_rows, by = "SNP", all.x = TRUE, allow.cartesian = TRUE)
conf[, allele_compatible := allele_ok(ref, alt, exposure_effect_allele, exposure_other_allele)]
conf_summary <- conf[, .(
  chromosome.grch38 = chromosome.grch38[1], position.grch38 = position.grch38[1],
  exposure_effect_allele = exposure_effect_allele[1], exposure_other_allele = exposure_other_allele[1],
  coordinate_source = coordinate_source[1],
  finngen_rows_with_rsid = sum(!is.na(finngen_position)),
  compatible_variant_rows = sum(allele_compatible %in% TRUE),
  incompatible_variant_rows = sum(allele_compatible %in% FALSE),
  compatible_variants = paste(unique(paste(finngen_chromosome[allele_compatible %in% TRUE],
                                           finngen_position[allele_compatible %in% TRUE],
                                           ref[allele_compatible %in% TRUE], alt[allele_compatible %in% TRUE], sep = ":")), collapse = ";"),
  incompatible_variants = paste(unique(paste(finngen_chromosome[allele_compatible %in% FALSE],
                                             finngen_position[allele_compatible %in% FALSE],
                                             ref[allele_compatible %in% FALSE], alt[allele_compatible %in% FALSE], sep = ":")), collapse = ";")
), by = SNP]
conf_summary[, resolution := fcase(
  compatible_variant_rows == 1L, "one_unambiguous_valid_variant",
  compatible_variant_rows > 1L, "multiple_valid_variants",
  finngen_rows_with_rsid > 0L, "incompatible_alleles",
  default = "no_valid_outcome_variant"
)]
fwrite(conf_summary, conflict_file, sep = "\t", quote = FALSE, na = "NA")

message("GRCh38 coordinates recovered: ", nrow(resolved), "/5254; unresolved: ", nrow(unresolved))
message("Completed corrected outcomes: ", sum(completed), "/12")
if (!all(completed)) quit(status = 1L)
