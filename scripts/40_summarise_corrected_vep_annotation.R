#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

# Summarise only the corrected GRCh38 VEP run and update only annotation-related
# dissertation evidence. Scripts 17--19 remain preserved but are superseded.

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
p <- function(...) file.path(root, ...)
fail <- function(...) stop(..., call. = FALSE)
base <- p("results", "annotation_grch38")
mapping_file <- file.path(base, "qc", "vep_grch38_input_mapping.tsv")
raw_file <- file.path(base, "raw", "lead_snps_grch38_vep_output.txt")
summary_dir <- file.path(base, "summary")
qc_dir <- file.path(base, "qc")
dis <- p("results", "dissertation", "results_chapter")

map <- fread(mapping_file, na.strings = c("NA", ""))
v <- fread(raw_file, skip = "#Uploaded_variation", fill = TRUE)
names(v) <- sub("^#", "", names(v))
required <- c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Consequence", "Existing_variation", "IMPACT", "SYMBOL")
if (!all(required %in% names(v))) fail("Corrected VEP output missing columns: ", paste(setdiff(required, names(v)), collapse = ", "))
setnames(v, "Uploaded_variation", "SNP")
v[, SNP := as.character(SNP)]
annot_ids <- unique(v$SNP)
if (nrow(map) != 5254L || uniqueN(map$SNP) != 5254L) fail("Mapping audit is not one row per 5,254 lead SNPs")
if (!all(annot_ids %in% map$SNP)) fail("VEP returned non-authoritative IDs")

severity <- c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant",
              "stop_lost", "start_lost", "transcript_amplification", "feature_elongation", "feature_truncation",
              "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant", "splice_donor_5th_base_variant",
              "splice_region_variant", "splice_donor_region_variant", "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant",
              "start_retained_variant", "stop_retained_variant", "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant",
              "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant",
              "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant", "regulatory_region_ablation",
              "regulatory_region_amplification", "regulatory_region_variant", "intergenic_variant", "sequence_variant")
rank_consequence <- function(x) {
  toks <- strsplit(as.character(x), "[,&]", perl = TRUE)
  vapply(toks, function(z) { m <- match(z, severity); if (all(is.na(m))) 999L else min(m, na.rm = TRUE) }, integer(1))
}
v[, consequence_rank__ := rank_consequence(Consequence)]
v[, impact_rank__ := match(IMPACT, c("HIGH", "MODERATE", "LOW", "MODIFIER", "-"))]
v[is.na(impact_rank__), impact_rank__ := 99L]
v[, gene_label__ := fifelse(!is.na(SYMBOL) & SYMBOL != "-", SYMBOL, fifelse(!is.na(Gene) & Gene != "-", Gene, NA_character_))]
v[, gene_missing__ := is.na(gene_label__)]
setorder(v, SNP, consequence_rank__, impact_rank__, gene_missing__, Feature)
best <- v[, .SD[1], by = SNP]
best[, Location_chr := as.integer(sub(":.*", "", Location))]
best[, Location_pos := as.integer(sub("-.*", "", sub(".*:", "", Location)))]
best <- merge(map, best[, .(SNP, vep_chr = Location_chr, vep_pos = Location_pos, vep_allele = Allele,
                            Gene = gene_label__, Most_severe_consequence = Consequence, Impact = IMPACT,
                            Feature, Existing_variation)], by = "SNP", all.x = TRUE, sort = FALSE)
best[, annotation_status := fcase(
  coordinate_status == "UNRESOLVED_GRCH38_COORDINATE", "COORDINATE_UNRESOLVED",
  vep_ready != TRUE, "ALLELE_UNRESOLVED",
  !SNP %in% annot_ids, "REJECTED_OR_NOT_RETURNED_BY_VEP",
  is.na(Gene), "ANNOTATED_NO_GENE_ASSIGNMENT",
  default = "SUCCESSFULLY_ANNOTATED"
)]
best[, coordinate_match := !is.na(vep_chr) & vep_chr == corrected_chr_grch38 & vep_pos == corrected_pos_grch38]
if (any(best[SNP %in% annot_ids]$coordinate_match != TRUE)) fail("At least one VEP output coordinate differs from corrected GRCh38 mapping")
if (best[SNP == "rs9273374", corrected_pos_grch38] != 32658837L) fail("rs9273374 corrected coordinate lost")

completion_by_snp <- best[, .(SNP, coordinate_status, allele_status, vep_ready, annotation_status,
                               corrected_chr_grch38, corrected_pos_grch38, REF, ALT, coordinate_match)]
fwrite(completion_by_snp, file.path(qc_dir, "vep_annotation_completion_by_snp.tsv"), sep = "\t", na = "NA", quote = FALSE)
n_coord <- sum(map$coordinate_status == "RESOLVED_GRCH38")
n_ready <- sum(map$vep_ready == TRUE)
n_ann <- uniqueN(v$SNP)
n_failed <- n_ready - n_ann
metrics <- data.table(metric = c("total_lead_SNPs", "GRCh38_coordinate_resolved", "VEP_input_ready", "successfully_annotated",
                                 "failed_or_not_returned", "coordinate_unresolved", "allele_unresolved_among_coordinate_resolved",
                                 "successfully_annotated_percent_of_all", "successfully_annotated_percent_of_VEP_ready"),
                      value = c(5254, n_coord, n_ready, n_ann, n_failed, 5254 - n_coord, n_coord - n_ready,
                                100 * n_ann / 5254, if (n_ready) 100 * n_ann / n_ready else NA_real_))
fwrite(metrics, file.path(qc_dir, "vep_annotation_completion_qc.tsv"), sep = "\t", na = "NA", quote = FALSE)

collapse_consequence <- function(x) {
  fcase(grepl("intron", x), "intronic", grepl("intergenic", x), "intergenic",
        grepl("upstream|downstream", x), "upstream/downstream", grepl("regulatory|TFBS|TF_binding", x), "regulatory",
        grepl("UTR", x), "UTR", grepl("synonymous", x), "synonymous", grepl("missense|protein_altering", x), "missense",
        grepl("splice", x), "splice-related", default = "other")
}
raw_summary <- best[SNP %in% annot_ids, .N, by = .(consequence = Most_severe_consequence)][order(-N, consequence)]
raw_summary[, `:=`(summary_type = "raw_most_severe_VEP_term", percentage = 100 * N / sum(N),
                   collapse_rule = "One most-severe VEP consequence per annotated lead SNP; VEP terminology retained")]
collapsed <- best[SNP %in% annot_ids, .(consequence = collapse_consequence(Most_severe_consequence))][, .N, by = consequence][order(-N, consequence)]
collapsed[, `:=`(summary_type = "dissertation_collapsed_group", percentage = 100 * N / sum(N),
                 collapse_rule = "intron->intronic; intergenic->intergenic; upstream/downstream combined; regulatory/TFBS combined; UTR combined; synonymous; missense/protein_altering; splice terms combined; remainder=other")]
cons <- rbindlist(list(raw_summary, collapsed), use.names = TRUE)[, .(summary_type, consequence, count = N, percentage, collapse_rule)]
fwrite(cons, file.path(summary_dir, "vep_consequence_summary.tsv"), sep = "\t", na = "NA", quote = FALSE)

gene_rows <- unique(v[gene_label__ %in% gene_label__[!is.na(gene_label__)], .(SNP, gene = gene_label__, Consequence)])
gene_rows <- merge(gene_rows, map[, .(SNP, association_p)], by = "SNP")
gene_summary <- gene_rows[, {
  z <- .SD[order(association_p)][1]
  rr <- rank_consequence(Consequence); severe <- Consequence[which.min(rr)]
  .(number_of_lead_SNP_annotations = uniqueN(SNP), strongest_associated_SNP = z$SNP,
    strongest_association_P = z$association_p, most_severe_consequence_observed = severe)
}, by = gene][order(-number_of_lead_SNP_annotations, strongest_association_P, gene)]
fwrite(gene_summary, file.path(summary_dir, "vep_gene_summary.tsv"), sep = "\t", na = "NA", quote = FALSE)

top20_ids <- map[order(association_p)]$SNP[1:20]
top20 <- best[match(top20_ids, SNP), .(SNP, Chromosome = corrected_chr_grch38, `GRCh38 position` = corrected_pos_grch38,
                                      `Association P` = association_p, REF, ALT, Gene,
                                      `Most severe consequence` = Most_severe_consequence, Impact, Feature,
                                      `Annotation status` = annotation_status)]
fwrite(top20, file.path(summary_dir, "top_20_annotated_lead_snps.tsv"), sep = "\t", na = "NA", quote = FALSE)

set.seed(220826)
eligible_random <- setdiff(best[SNP %in% annot_ids]$SNP, top20_ids)
random_ids <- unique(best[SNP %in% eligible_random, .SD[sample(.N, 1)], by = corrected_chr_grch38]$SNP)
if (length(random_ids) < 20L) random_ids <- unique(c(random_ids, sample(setdiff(eligible_random, random_ids), 20L - length(random_ids))))
key_ids <- unique(c("rs9273374", top20_ids, random_ids))
key_qc <- best[match(key_ids, SNP), .(SNP, audit_group = fcase(SNP == "rs9273374", "special_validation",
                                                              SNP %in% top20_ids, "top_20_association", default = "deterministic_random_across_chromosomes"),
                                        GRCh38_chromosome = corrected_chr_grch38, GRCh38_position = corrected_pos_grch38,
                                        REF, ALT, VEP_gene = Gene, most_severe_consequence = Most_severe_consequence,
                                        impact = Impact, VEP_location_chromosome = vep_chr, VEP_location_position = vep_pos,
                                        coordinate_match, annotation_status)]
fwrite(key_qc, file.path(qc_dir, "vep_key_variant_validation.tsv"), sep = "\t", na = "NA", quote = FALSE)

# Provenance is intentionally explicit and the superseded output is retained.
writeLines(c(
  "VEP coordinate correction provenance", "",
  "The final PLINK clumping reference positions were inherited from the 1000 Genomes European panel and are GRCh37/hg19.",
  "Script 17 used those positions directly and script 18 ran VEP with --assembly GRCh38; therefore the original results/annotation output was coordinate-incompatible and invalid.",
  "Validated GRCh38 coordinates were subsequently recovered in results/mr/finngen/exposure/asthma_exposure_5254_grch38.tsv.",
  sprintf("The corrected VEP run submitted %d variants with validated GRCh38 coordinates and genomic REF/ALT; %d remained coordinate unresolved and %d additional coordinate-resolved variants remained allele unresolved.", n_ready, 5254 - n_coord, n_coord - n_ready),
  "VEP was rerun offline using Ensembl VEP/cache 116 and GRCh38. Old annotation outputs were excluded from dissertation reporting and remain preserved only for provenance."
), file.path(qc_dir, "vep_coordinate_correction_provenance.txt"))

# Figure 3.4: collapsed one-SNP/one-most-severe-consequence distribution.
figdat <- collapsed[order(N)]
pdf_path <- file.path(dis, "figures", "Figure_3_4_VEP_consequence_distribution.pdf")
png_path <- file.path(dis, "figures", "Figure_3_4_VEP_consequence_distribution.png")
plot_fig <- function() {
  par(mar = c(5, 12, 4, 2) + 0.1)
  barplot(figdat$N, names.arg = figdat$consequence, horiz = TRUE, las = 1, col = "#4C78A8", border = NA,
          xlab = "Independent lead SNPs", main = "GRCh38-corrected VEP consequence distribution")
}
pdf(pdf_path, width = 9, height = 6.5); plot_fig(); dev.off()
png(png_path, width = 1800, height = 1300, res = 200); plot_fig(); dev.off()

# Preserve the previous placeholder separately, then activate corrected Table 3.4.
active_table <- file.path(dis, "tables", "Table_3_4_top_annotated_lead_variants.tsv")
placeholder <- file.path(dis, "tables", "Table_3_4_top_annotated_lead_variants_SUPERSEDED_COORDINATE_MISMATCH.tsv")
if (file.exists(active_table) && !file.exists(placeholder) && any(grepl("NOT_REPORTABLE_COORDINATE_BUILD_MISMATCH", readLines(active_table))))
  file.copy(active_table, placeholder, overwrite = FALSE)
fwrite(top20, active_table, sep = "\t", na = "NA", quote = FALSE)

# Annotation-only evidence-pack updates.
inventory_file <- file.path(dis, "results_inventory.tsv")
inventory <- fread(inventory_file)
i <- which(inventory$section == "3.6")
inventory[i, `:=`(analysis_stage = "Corrected GRCh38 VEP functional annotation",
                   operation = "Annotated coordinate- and allele-resolved independent lead variants with local VEP",
                   design = "Validated GRCh38 coordinates; Ensembl VEP/cache 116; one most-severe consequence per SNP for summaries",
                   source_files = "results/annotation_grch38/raw/lead_snps_grch38_vep_output.txt; results/annotation_grch38/qc/vep_annotation_completion_qc.tsv",
                   key_result = sprintf("%d of 5254 lead SNPs VEP-ready; %d successfully annotated", n_ready, n_ann),
                   statistical_framework = "Descriptive counts and percentages only",
                   recommended_main_table = "Table 3.4", recommended_main_figure = "Figure 3.4",
                   immediate_conclusion = "Corrected GRCh38 annotation characterises the functional context of resolved lead SNPs.",
                   logical_link = "The annotated lead variants characterised independent signals; the same lead-SNP set underpinned MR instrument preparation.",
                   main_text_or_appendix = "Main text")]
fwrite(inventory, inventory_file, sep = "\t", na = "NA")

num_file <- file.path(dis, "results_numbers_master.tsv")
nums <- fread(num_file)
nums <- nums[section != "3.6"]
newnums <- data.table(section = "3.6", metric = c("Lead SNPs", "Validated GRCh38 coordinates", "VEP-ready SNPs", "Successfully annotated SNPs", "Coordinate unresolved SNPs", "Annotation completion"),
                      value = c(5254, n_coord, n_ready, n_ann, 5254 - n_coord, 100 * n_ann / 5254),
                      unit = c(rep("SNPs", 5), "percent"),
                      source_file = c("results/ld_clumping/independent_lead_snps.tsv", "results/annotation_grch38/qc/vep_grch38_input_mapping.tsv",
                                      "results/annotation_grch38/qc/vep_grch38_input_mapping.tsv", "results/annotation_grch38/qc/vep_annotation_completion_qc.tsv",
                                      "results/annotation_grch38/qc/vep_grch38_input_mapping.tsv", "results/annotation_grch38/qc/vep_annotation_completion_qc.tsv"),
                      source_column = c("rows", "coordinate_status", "vep_ready", "successfully_annotated", "coordinate_status", "successfully_annotated_percent_of_all"),
                      row_or_contrast = c("all", "RESOLVED_GRCH38", "TRUE", "metric", "UNRESOLVED_GRCH38_COORDINATE", "metric"),
                      formatted_value = c("5,254", format(n_coord, big.mark = ","), format(n_ready, big.mark = ","), format(n_ann, big.mark = ","),
                                          as.character(5254 - n_coord), sprintf("%.2f%%", 100 * n_ann / 5254)), notes = "Corrected GRCh38 VEP annotation")
nums <- rbindlist(list(nums, newnums), fill = TRUE)
setorder(nums, section)
fwrite(nums, num_file, sep = "\t", na = "NA")

fig_index_file <- file.path(dis, "figure_index.tsv"); fi <- fread(fig_index_file)
fi[Number == "Figure 3.4", `:=`(Title = "GRCh38-corrected VEP consequence distribution of independent lead variants",
                                 Recommendation = "Essential main figure", File = "figures/Figure_3_4_VEP_consequence_distribution.pdf", Status = "READY")]
fwrite(fi, fig_index_file, sep = "\t", na = "NA")
tab_index_file <- file.path(dis, "table_index.tsv"); ti <- fread(tab_index_file)
ti[Number == "Table 3.4", `:=`(Title = "Top functionally annotated independent lead variants",
                                File = "tables/Table_3_4_top_annotated_lead_variants.tsv", Status = "READY_CORRECTED_GRCH38")]
fwrite(ti, tab_index_file, sep = "\t", na = "NA")

replace_line <- function(path, prefix, replacement) {
  x <- readLines(path, warn = FALSE); hit <- startsWith(x, prefix); if (!any(hit)) fail("Caption entry not found: ", prefix)
  x[hit] <- replacement; writeLines(x, path)
}
replace_line(file.path(dis, "captions", "results_figure_captions.txt"), "Figure 3.4.",
             sprintf("Figure 3.4. GRCh38-corrected Ensembl VEP consequence distribution for %s successfully annotated independent asthma lead SNPs. Categories represent one most-severe consequence per SNP and are shown as descriptive counts.", format(n_ann, big.mark = ",")))
replace_line(file.path(dis, "captions", "results_table_captions.txt"), "Table 3.4.",
             "Table 3.4. Top 20 independent asthma lead SNPs ranked by association P-value and annotated using validated GRCh38 coordinates with Ensembl VEP/cache release 116. Unavailable gene or consequence assignments are retained as NA.")

blue_file <- file.path(dis, "results_writing_blueprint.md")
blue <- readLines(blue_file, warn = FALSE)
s <- grep("^## 3\\.6 ", blue); e <- grep("^## 3\\.7", blue)
if (length(s) != 1L || !length(e) || e[1] <= s) fail("Could not locate Section 3.6 in Results blueprint")
topcats <- paste(sprintf("%s: %d (%.1f%%)", head(raw_summary$consequence, 5), head(raw_summary$N, 5), head(raw_summary$percentage, 5)), collapse = "; ")
section <- c("## 3.6 Functional annotation of independent lead variants", "",
  "**Operation carried out**", "- Ensembl VEP/cache release 116 was run offline on the independent lead SNPs with validated GRCh38 coordinates and resolved genomic REF/ALT alleles.", "",
  "**Design features that need mentioning**", "- The authoritative set remained 5,254 SNPs. GRCh37 PLINK positions were not submitted; unresolved coordinates or alleles were retained in QC and excluded from coordinate-based VEP.", "",
  "**Results/numbers that MUST appear**",
  sprintf("- %s of 5,254 SNPs had validated GRCh38 coordinates; %s were VEP-ready; %s were successfully annotated (%.2f%% of all lead SNPs).", format(n_coord, big.mark = ","), format(n_ready, big.mark = ","), format(n_ann, big.mark = ","), 100 * n_ann / 5254),
  sprintf("- The five most common raw most-severe consequence categories were %s.", topcats),
  "- Report the top variants and gene assignments exactly as supported in Table 3.4; do not infer biological function.", "",
  "**Statistical analysis**", "- Descriptive counts and percentages; one most-severe VEP consequence per successfully annotated SNP.", "",
  "**Recommended table/figure**", "- Table 3.4 and Figure 3.4.", "",
  "**Immediate conclusion**", "- Corrected annotation describes the genomic consequence context of the coordinate- and allele-resolved independent association signals without functional interpretation.", "",
  "**Logical link to next analysis**", "- The annotated lead variants were subsequently used to characterise the functional context of the independent association signals, while the lead SNP set also formed the basis of downstream MR instrument preparation.", "")
blue <- c(blue[seq_len(s - 1L)], section, blue[e[1]:length(blue)])
writeLines(blue, blue_file)

cat(sprintf("total=%d\ncoordinate_resolved=%d\nvep_ready=%d\nannotated=%d\nfailed=%d\nunresolved=%d\ncompletion_percent=%.4f\n",
            nrow(map), n_coord, n_ready, n_ann, n_failed, 5254 - n_coord, 100 * n_ann / 5254))
