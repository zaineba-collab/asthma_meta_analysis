#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

# Corrective GRCh38 VEP preparation. Scripts 17--19 used GRCh37 PLINK
# coordinates with a GRCh38 cache and are superseded for final reporting by
# scripts 38--40. This script never changes the authoritative SNP set or any
# association statistic.

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
p <- function(...) file.path(root, ...)
fail <- function(...) stop(..., call. = FALSE)

lead_file <- p("results", "ld_clumping", "independent_lead_snps.tsv")
exposure_file <- p("results", "mr", "finngen", "exposure", "asthma_exposure_5254_grch38.tsv")
repair_qc_file <- p("results", "mr", "finngen", "qc", "exposure_coordinate_repair_qc.tsv")
base <- p("results", "annotation_grch38")
dirs <- file.path(base, c("input", "raw", "summary", "qc", "figures"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
mapping_file <- file.path(base, "qc", "vep_grch38_input_mapping.tsv")
probe_input <- file.path(base, "input", "lead_snps_grch38_allele_orientation_probe.tsv")
probe_output <- file.path(base, "raw", "lead_snps_grch38_allele_orientation_probe_vep.txt")
final_input <- file.path(base, "input", "lead_snps_grch38_vep_input.tsv")

lead <- fread(lead_file)
exp <- fread(exposure_file, na.strings = c("NA", ""))
qc <- fread(repair_qc_file)
if (nrow(lead) != 5254L || uniqueN(lead$SNP) != 5254L) fail("Authoritative lead set is not 5,254 unique SNPs")
if (nrow(exp) != 5254L || uniqueN(exp$SNP) != 5254L || !identical(lead$SNP, exp$SNP)) fail("Repaired exposure differs from authoritative lead set/order")
if (!identical(as.numeric(lead$p_value), as.numeric(exp$pval.exposure))) fail("Association P-values differ between authoritative lead and repaired exposure")
expected <- setNames(as.character(qc$value), qc$metric)
if (expected[["GRCh38_coordinates_recovered"]] != "5238" || expected[["coordinates_unresolved"]] != "16" ||
    expected[["old_position_differs_new"]] != "3047" || expected[["old_position_equals_new"]] != "2191" ||
    expected[["chromosome_changed"]] != "0") fail("Coordinate repair QC does not match validated expectations")
key <- exp[SNP == "rs9273374"]
if (nrow(key) != 1L || key$chromosome.plink_reference != 6L || key$position.plink_reference != 32626614L ||
    key$chromosome.grch38 != 6L || key$position.grch38 != 32658837L) fail("rs9273374 coordinate validation failed")

comp <- function(x) chartr("ACGT", "TGCA", x)
valid <- function(x) !is.na(x) & grepl("^[ACGT]+$", x)
exp[, `:=`(effect_allele.exposure = toupper(effect_allele.exposure), other_allele.exposure = toupper(other_allele.exposure))]

# Prefer genomic REF/ALT already observed in corrected local FinnGen outcome
# extracts. These tables define outcome effect allele as FinnGen ALT and other
# allele as FinnGen REF (script 22), and are used only for allele orientation.
fg_files <- list.files(p("results", "mr", "finngen", "outcomes_grch38"), pattern = "_outcome\\.tsv$", full.names = TRUE)
fg <- if (length(fg_files)) rbindlist(lapply(fg_files, function(f) {
  x <- fread(f, select = c("SNP", "chromosome.outcome", "position.outcome", "effect_allele.outcome", "other_allele.outcome"))
  x[, .(SNP, corrected_chr_grch38 = as.integer(chromosome.outcome), corrected_pos_grch38 = as.integer(position.outcome),
        REF = toupper(other_allele.outcome), ALT = toupper(effect_allele.outcome), source = "FinnGen_GRCh38_REF_ALT")]
}), fill = TRUE) else data.table()
if (nrow(fg)) fg <- unique(fg[valid(REF) & valid(ALT)])
fg_unique <- if (nrow(fg)) fg[, if (uniqueN(paste(corrected_chr_grch38, corrected_pos_grch38, REF, ALT)) == 1L) .SD[1], by = SNP] else fg

map <- exp[, .(
  SNP,
  original_plink_chr = as.integer(chromosome.plink_reference),
  original_plink_pos_grch37 = as.integer(position.plink_reference),
  corrected_chr_grch38 = as.integer(chromosome.grch38),
  corrected_pos_grch38 = as.integer(position.grch38),
  association_p = as.numeric(pval.exposure),
  ea = effect_allele.exposure,
  oa = other_allele.exposure,
  coordinate_source
)]
map[, coordinate_status := fifelse(!is.na(corrected_chr_grch38) & corrected_chr_grch38 %in% 1:22 &
                                     !is.na(corrected_pos_grch38) & corrected_pos_grch38 > 0,
                                   "RESOLVED_GRCH38", "UNRESOLVED_GRCH38_COORDINATE")]
map[, `:=`(REF = NA_character_, ALT = NA_character_, allele_status = "UNRESOLVED_ALLELES", mapping_source = coordinate_source,
           notes = fifelse(coordinate_status == "UNRESOLVED_GRCH38_COORDINATE", "No validated GRCh38 coordinate; excluded from coordinate-based VEP", "Awaiting genomic REF/ALT resolution"))]

# Direct reference-base retrieval from the indexed Ensembl release-116 GRCh38
# primary assembly. The helper keeps one faidx handle open for all 5,238 loci.
fasta <- p("data", "reference", "ensembl_116_GRCh38_primary_assembly.fa")
ref_helper <- p("scripts", "38_extract_grch38_reference_bases.pl")
if (!file.exists(fasta) || !file.exists(paste0(fasta, ".fai"))) fail("Indexed Ensembl-116 GRCh38 FASTA is unavailable")
ref_cmd <- sprintf("perl %s %s %s", shQuote(ref_helper), shQuote(fasta), shQuote(exposure_file))
refs <- fread(cmd = ref_cmd)
if (nrow(refs) != 5238L || uniqueN(refs$SNP) != 5238L)
  fail("GRCh38 FASTA reference-base extraction did not return exactly 5,238 coordinate-resolved SNPs")
tmp <- merge(map[, .(SNP, ea, oa)], refs, by = "SNP")
tmp[, direct_match := REF_FASTA %in% c(ea, oa), by = SNP]
tmp[, complement_match := REF_FASTA %in% c(comp(ea), comp(oa)), by = SNP]
tmp[, ALT_FASTA := fcase(direct_match & REF_FASTA == ea, oa,
                         direct_match & REF_FASTA == oa, ea,
                         !direct_match & complement_match & REF_FASTA == comp(ea), comp(oa),
                         !direct_match & complement_match & REF_FASTA == comp(oa), comp(ea),
                         default = NA_character_)]
fasta_resolved <- tmp[valid(ALT_FASTA) & ALT_FASTA != REF_FASTA, .(SNP, REF = REF_FASTA, ALT = ALT_FASTA)]
if (nrow(fasta_resolved) < 1L) fail("No SNP alleles were compatible with the GRCh38 reference")
map[fasta_resolved, on = "SNP", `:=`(REF = i.REF, ALT = i.ALT, allele_status = "RESOLVED_ENSEMBL116_GRCH38_FASTA",
                                      mapping_source = paste0(coordinate_source, ";Ensembl_116_GRCh38_primary_assembly_FASTA"),
                                      notes = "Genomic REF from indexed Ensembl release-116 GRCh38 FASTA; ALT selected from association allele pair after strand check")]
map[SNP %in% refs[REF_FASTA == "N", SNP], notes := "GRCh38 primary-assembly reference base unavailable (N); genomic REF/ALT unresolved"]
map[SNP %in% tmp[REF_FASTA %in% c("A", "C", "G", "T") & is.na(ALT_FASTA), SNP],
    notes := "Association allele pair and its strand complement are incompatible with the GRCh38 FASTA reference base; genomic REF/ALT unresolved"]

if (nrow(fg_unique)) {
  z <- merge(map[, .(SNP, corrected_chr_grch38, corrected_pos_grch38)], fg_unique,
             by = c("SNP", "corrected_chr_grch38", "corrected_pos_grch38"))
  if (nrow(z)) {
    map[z, on = "SNP", `:=`(REF = i.REF, ALT = i.ALT, allele_status = "RESOLVED_FINNGEN_GRCH38",
                             mapping_source = paste0(coordinate_source, ";FinnGen_GRCh38_REF_ALT"), notes = "Genomic REF/ALT from local FinnGen corrected outcome extract")]
  }
}

# For every coordinate-resolved SNV, probe A/C/G/T against the matching
# Ensembl-116 GRCh38 FASTA via VEP --lookup_ref. The true reference candidate
# is the single allele for which VEP returns no variant consequence. REF is
# then matched to either the association allele pair or its strand complement;
# the remaining allele is ALT. This never assumes that EA is REF.
# No cache-only fallback is attempted for FASTA-incompatible loci: absence of
# an exact dbSNP co-location cannot establish genomic REF. Retain them as
# allele-unresolved rather than infer an orientation.
need <- map[FALSE]
if (file.exists(probe_input)) unlink(probe_input)
probe <- rbindlist(lapply(seq_len(nrow(need)), function(i) {
  r <- need[i]
  cand <- data.table(ALT = c("A", "C", "G", "T"))
  cand[, .(chromosome = r$corrected_chr_grch38, start = r$corrected_pos_grch38, end = r$corrected_pos_grch38,
           allele = paste0("N/", ALT), strand = "+", identifier = paste(r$SNP, ALT, sep = "|"))]
}), fill = TRUE)
if (nrow(probe)) fwrite(probe, probe_input, sep = "\t", col.names = FALSE, quote = FALSE)

if (nrow(need) > 0L && file.exists(probe_output) && file.info(probe_output)$size > 0) {
  v <- fread(probe_output, skip = "#Uploaded_variation", fill = TRUE)
  names(v) <- sub("^#", "", names(v))
  if (!all(c("Uploaded_variation", "Existing_variation") %in% names(v))) fail("Orientation probe VEP output lacks required columns")
  ids <- tstrsplit(v$Uploaded_variation, "|", fixed = TRUE)
  v[, `:=`(SNP = ids[[1]], probe_ALT = ids[[2]])]
  returned <- unique(v[, .(SNP, probe_ALT)])
  grid <- CJ(SNP = need$SNP, candidate = c("A", "C", "G", "T"), unique = TRUE)
  grid[, returned := paste(SNP, candidate) %in% paste(returned$SNP, returned$probe_ALT)]
  refs <- grid[returned == FALSE, if (.N == 1L) .(REF_FASTA = candidate) else NULL, by = SNP]
  pairs <- merge(need[, .(SNP, ea, oa)], refs, by = "SNP")
  pairs[, `:=`(direct_match = REF_FASTA %in% c(ea, oa), complement_match = REF_FASTA %in% c(comp(ea), comp(oa))), by = SNP]
  pairs[, ALT_FASTA := fcase(direct_match & REF_FASTA == ea, oa,
                             direct_match & REF_FASTA == oa, ea,
                             !direct_match & complement_match & REF_FASTA == comp(ea), comp(oa),
                             !direct_match & complement_match & REF_FASTA == comp(oa), comp(ea),
                             default = NA_character_)]
  resolved <- pairs[valid(ALT_FASTA) & ALT_FASTA != REF_FASTA, .(SNP, REF = REF_FASTA, ALT = ALT_FASTA)]
  # Existing local FinnGen mappings must agree with the reference FASTA.
  check <- merge(map[allele_status == "RESOLVED_FINNGEN_GRCH38", .(SNP, old_REF = REF, old_ALT = ALT)], resolved, by = "SNP")
  if (nrow(check) && any(check$old_REF != check$REF | check$old_ALT != check$ALT)) fail("FinnGen REF/ALT conflicts with Ensembl-116 GRCh38 FASTA")
  if (nrow(resolved)) map[resolved, on = "SNP", `:=`(REF = i.REF, ALT = i.ALT, allele_status = "RESOLVED_ENSEMBL116_GRCH38_FASTA",
                                                       mapping_source = paste0(coordinate_source, ";Ensembl_116_GRCh38_primary_assembly_FASTA"),
                                                       notes = "Genomic REF from Ensembl release-116 GRCh38 FASTA; ALT selected from association allele pair after strand check")]
}

map[, vep_ready := coordinate_status == "RESOLVED_GRCH38" & valid(REF) & valid(ALT) & REF != ALT]
map[coordinate_status == "RESOLVED_GRCH38" & !vep_ready & !valid(ea), notes := "Invalid or unavailable association alleles; VEP allele unresolved"]
map[coordinate_status == "RESOLVED_GRCH38" & !vep_ready & notes == "Awaiting genomic REF/ALT resolution",
    notes := "No unique validated genomic REF/ALT mapping; excluded from final VEP"]
out_map <- map[, .(SNP, original_plink_chr, original_plink_pos_grch37, corrected_chr_grch38, corrected_pos_grch38,
                   REF, ALT, association_p, coordinate_status, allele_status, vep_ready, mapping_source, notes)]
fwrite(out_map, mapping_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(out_map[coordinate_status == "RESOLVED_GRCH38" & vep_ready == FALSE],
       file.path(base, "qc", "vep_grch38_allele_unresolved.tsv"), sep = "\t", quote = FALSE, na = "NA")

ready <- out_map[vep_ready == TRUE, .(chromosome = corrected_chr_grch38, start = corrected_pos_grch38,
                                      end = corrected_pos_grch38, allele = paste0(REF, "/", ALT), strand = "+", identifier = SNP)]
if (anyDuplicated(ready$identifier) || anyDuplicated(ready[, .(chromosome, start, allele)])) fail("Duplicate final VEP records detected")
if (ready[identifier == "rs9273374", start] != 32658837L) fail("Final rs9273374 VEP position is not chr6:32658837")
if (nrow(ready)) fwrite(ready, final_input, sep = "\t", col.names = FALSE, quote = FALSE)

cat(sprintf("lead_total=%d\ncoordinate_resolved=%d\nallele_resolved=%d\nvep_ready=%d\nunresolved_coordinate=%d\n",
            nrow(out_map), sum(out_map$coordinate_status == "RESOLVED_GRCH38"), sum(!is.na(out_map$REF) & !is.na(out_map$ALT)),
            sum(out_map$vep_ready), sum(out_map$coordinate_status == "UNRESOLVED_GRCH38_COORDINATE")))
