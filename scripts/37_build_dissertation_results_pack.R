#!/usr/bin/env Rscript

# Build a dissertation Results evidence pack from existing outputs only.
# No GWAMA, clumping, annotation, MR, LDSC, harmonisation, or statistical test is run.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", arg))
root <- dirname(dirname(script_path))
setwd(root)

out <- "results/dissertation/results_chapter"
tables_dir <- file.path(out, "tables")
figures_dir <- file.path(out, "figures")
captions_dir <- file.path(out, "captions")
appendix_dir <- file.path(out, "appendix_candidates")
invisible(lapply(c(out, tables_dir, figures_dir, captions_dir, appendix_dir), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

read_tsv <- function(path) fread(path, sep = "\t", header = TRUE, na.strings = "NA")
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
copy_as <- function(source, target) {
  stopifnot(file.exists(source))
  ok <- file.copy(source, target, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
  stopifnot(ok)
}
fmt_p <- function(x) vapply(x, function(v) {
  if (is.na(v)) "NA" else if (v < 0.001) format(v, digits = 3, scientific = TRUE) else format(round(v, 3), trim = TRUE)
}, character(1))
fmt_n <- function(x, digits = 4) ifelse(is.na(x), "NA", format(signif(x, digits), trim = TRUE, scientific = FALSE))
bool <- function(x) tolower(as.character(x)) %in% c("true", "yes")

# ------------------------------ authoritative source reads
gw_top <- read_tsv("results/gwama/top20_signals.txt")
lambda <- read_tsv("results/gwama/lambda_gc.txt")
hetero <- read_tsv("results/gwama/heterogeneity_summary.txt")
i2_dist <- read_tsv("results/gwama/i2_distribution.txt")
map_summary <- read_tsv("results/gwama/asthma_fixed_for_plink_clumping_mapping_summary.txt")
sentinels <- read_tsv("results/pruning/asthma_fixed_sentinels_clean.txt")
loci <- read_tsv("results/loci/asthma_fixed_sentinel_merge_loci.txt")
loci_summary <- read_tsv("results/loci/loci_summary.txt")
clump_summary <- read_tsv("results/ld_clumping/plink_clumping_summary.txt")
leads <- read_tsv("results/ld_clumping/independent_lead_snps.tsv")
lead_annotation <- read_tsv("results/annotation/lead_snp_annotation_table.tsv")
vep_clean <- read_tsv("results/annotation/independent_lead_snps_vep_clean_summary.txt")
coord_repair <- read_tsv("results/mr/finngen/qc/exposure_coordinate_repair_qc.tsv")
extract_qc <- read_tsv("results/mr/finngen/qc/finngen_outcome_extraction_grch38_qc.tsv")
harm_qc <- read_tsv("results/mr/finngen/qc/finngen_harmonisation_qc.tsv")
mr <- read_tsv("results/mr/finngen/summary/finngen_mr_final_summary.tsv")
sens <- read_tsv("results/mr/finngen/sensitivity/finngen_mr_sensitivity_summary.tsv")
h2 <- read_tsv("results/ldsc/qc/ldsc_h2_summary.tsv")
h2_ready <- read_tsv("results/ldsc/qc/ldsc_h2_readiness.tsv")
rg <- read_tsv("results/ldsc/qc/ldsc_rg_summary.tsv")
combined <- read_tsv("results/combined/finngen_mr_ldsc_comparison.tsv")
combined_diss <- read_tsv("results/combined/finngen_mr_ldsc_dissertation_table.tsv")

metric <- function(dt, label) dt[Metric == label, Value][1]
metric_lower <- function(dt, label) dt[tolower(metric) == tolower(label), value][1]
coord_metric <- function(label) coord_repair[metric == label, value][1]

# ------------------------------ core validations and derived reporting counts
n_studies <- 12L
n_variants <- as.integer(hetero[metric == "Total SNPs in GWAMA output", value])
n_conventional <- nrow(read_tsv("results/gwama/genome_wide_significant.txt"))
n_project <- as.integer(map_summary[Metric == "Variants evaluated for PLINK clumping input (p-value <= 5e-09)", Value])
n_sentinels <- nrow(sentinels)
n_loci <- nrow(loci)
n_leads <- nrow(leads)
top <- gw_top[which.min(`p-value`)]
lambda_gc <- as.numeric(lambda$value[1])
n_i2 <- as.integer(hetero[metric == "SNPs with I2 >= 0.75", value])
pct_i2 <- as.numeric(hetero[metric == "Percentage of SNPs with I2 >= 0.75", value])
n_q <- as.integer(hetero[metric == "SNPs with Q p-value < 0.05", value])
pct_q <- as.numeric(hetero[metric == "Percentage of SNPs with Q p-value < 0.05", value])
n_rsid <- as.integer(map_summary[Metric == "Variants matched to rsID", Value])
n_unmatched <- as.integer(map_summary[Metric == "Variants not matched", Value])
map_pct <- as.numeric(map_summary[Metric == "Percentage matched among evaluated variants", Value])
n_dup <- as.integer(map_summary[Metric == "Duplicate rsID rows collapsed", Value])

stopifnot(n_studies == 12L, n_variants == 34316972L, n_project == 39933L,
          n_sentinels == 1563L, n_loci == 405L, n_leads == 5254L,
          top$rs_number == "rs9273374",
          abs(top$`p-value` / 5.22e-253 - 1) < 1e-12)
stopifnot(nrow(mr) == 12L, nrow(rg) == 12L, nrow(combined) == 12L)
stopifnot(sum(bool(mr$ivw_bonferroni_significant)) == 0L)
stopifnot(sum(bool(rg$bonferroni12_significant), na.rm = TRUE) == 7L)
stopifnot(rg[outcome == "Eosinophilic disease", run_status] == "NOT_ESTIMABLE_LOW_H2")

# VEP audit: the active output was prepared with PLINK GRCh37 positions but run as GRCh38.
vep_input_n <- length(readLines("results/annotation/independent_lead_snps_vep_input.txt"))
vep_unique_n <- uniqueN(vep_clean$SNP)
vep_gene_n <- uniqueN(vep_clean[!is.na(SYMBOL) & SYMBOL != "-", SNP])
plink_build <- coord_metric("inferred_PLINK_build")
vep_build_mismatch <- identical(plink_build, "GRCh37/hg19") &&
  any(grepl("--assembly GRCh38", readLines("scripts/18_run_vep_annotation.sh", warn = FALSE), fixed = TRUE))
stopifnot(vep_build_mismatch, vep_input_n == 5252L, vep_unique_n == 5252L)

# ------------------------------ main tables
table_3_1 <- data.table(
  Metric = c("Number of input GWAS", "Variants analysed", "Variants with P < 5e-8",
             "Project-significant variants (P <= 5e-9)", "Top SNP", "Top chromosome",
             "Top position (corrected GRCh38 coordinate)", "Top beta", "Top SE",
             "Top association P", "lambda_GC", "Variants with I2 >= 0.75",
             "Percentage with I2 >= 0.75", "Variants with Cochran Q P < 0.05",
             "Percentage with Cochran Q P < 0.05", "Project signal-selection threshold"),
  Result = c(n_studies, n_variants, n_conventional, n_project, top$rs_number, 6, 32658837,
             top$beta, top$se, top$`p-value`, lambda_gc, n_i2, pct_i2, n_q, pct_q, "P <= 5e-9")
)
write_tsv(table_3_1, file.path(tables_dir, "Table_3_1_GWAMA_summary.tsv"))

table_3_2 <- data.table(
  Metric = c("Project-significant candidate variants", "Distance-based sentinel SNPs", "Merged loci",
             "Median locus size (bp)", "Mean locus size (bp)", "Largest locus (bp)",
             "Maximum sentinel SNPs within one locus"),
  Result = c(n_project, n_sentinels, n_loci,
             metric(loci_summary, "Median locus size (bp)"), metric(loci_summary, "Mean locus size (bp)"),
             metric(loci_summary, "Largest locus (bp)"), metric(loci_summary, "Maximum sentinel SNPs within one locus"))
)
write_tsv(table_3_2, file.path(tables_dir, "Table_3_2_signal_refinement_summary.tsv"))

top_loci <- loci[order(lead_p)][1:10, .(Locus = locus, Chromosome = chrom, Start = locus_start,
                                       End = locus_end, Lead_sentinel_SNP = first, P = lead_p,
                                       Number_of_sentinels = n_sentinels_in_locus)]
write_tsv(top_loci, file.path(tables_dir, "Table_3_2b_top_loci.tsv"))

table_3_3 <- data.table(
  Metric = c("Input variants", "rsID matched", "rsID unmatched", "Mapping percentage",
             "Duplicate rsID rows collapsed", "Clumping r2", "Clumping window (kb)",
             "P1", "P2", "Final independent lead SNPs"),
  Result = c(n_project, n_rsid, n_unmatched, map_pct, n_dup, 0.1, 99999, 5e-9, 5e-9, n_leads)
)
write_tsv(table_3_3, file.path(tables_dir, "Table_3_3_LD_clumping_summary.tsv"))

# Annotation values are deliberately withheld because the active VEP coordinates are build-mismatched.
top20_assoc <- leads[order(p_value)][1:20]
table_3_4 <- top20_assoc[, .(SNP, Chromosome = chromosome, Position = position, P = p_value)]
table_3_4[, `:=`(Gene = NA_character_, Most_severe_consequence = NA_character_, Impact = NA_character_,
                 Annotation_status = "NOT_REPORTABLE_COORDINATE_BUILD_MISMATCH")]
write_tsv(table_3_4, file.path(tables_dir, "Table_3_4_top_annotated_lead_variants.tsv"))
write_tsv(data.table(
  Check = c("Independent lead SNPs", "rsID SNPs submitted to VEP", "Unique SNPs returned by VEP",
            "PLINK coordinate build established later", "VEP assembly setting", "Dissertation reporting decision"),
  Result = c(n_leads, vep_input_n, vep_unique_n, plink_build, "GRCh38",
             "WITHHOLD functional annotations pending coordinate-corrected VEP")
), file.path(tables_dir, "Table_3_4a_VEP_annotation_audit.tsv"))

harm_idx <- match(extract_qc$outcome_id, harm_qc$outcome_id)
stopifnot(!anyNA(harm_idx))
table_3_5 <- data.table(
  Outcome = extract_qc$outcome,
  FinnGen_matched_SNPs = extract_qc$matched_snps,
  Entering_harmonisation = harm_qc$snps_entering_harmonisation[harm_idx],
  Palindromic_ambiguous_removed = harm_qc$palindromic_unresolved_missing_exposure_eaf[harm_idx],
  Final_MR_ready_SNPs = harm_qc$final_mr_keep_snps[harm_idx],
  Retention_percent = harm_qc$retention_percentage[harm_idx]
)
write_tsv(table_3_5, file.path(tables_dir, "Table_3_5_MR_instrument_harmonisation_summary.tsv"))

table_3_6 <- mr[, .(
  Outcome = outcome, n_SNPs = final_mr_ready_snps, IVW_beta = ivw_beta, SE = ivw_SE,
  Lower_95CI = ivw_lower_95CI, Upper_95CI = ivw_upper_95CI, P = ivw_p_value,
  OR = OR, OR_lower_95CI = OR_lower_95CI, OR_upper_95CI = OR_upper_95CI,
  Bonferroni_P = ivw_p_bonferroni, FDR_P = ivw_p_fdr_bh
)]
write_tsv(table_3_6, file.path(tables_dir, "Table_3_6_primary_MR_results.tsv"))

table_3_7 <- sens[, .(
  Outcome = outcome, IVW_heterogeneity_P = ivw_heterogeneity_p_value,
  Heterogeneity_evidence = heterogeneity_evidence, Egger_intercept = egger_intercept,
  Egger_intercept_P = egger_intercept_p_value, Egger_FDR_P = p_egger_fdr,
  Directional_pleiotropy_flag = directional_pleiotropy_evidence,
  LOO_sign_reversal = leave_one_out_sign_change,
  Most_influential_SNP = most_influential_leave_one_out_snp
)]
write_tsv(table_3_7, file.path(tables_dir, "Table_3_7_MR_sensitivity_summary.tsv"))

table_3_8 <- rg[, .(
  Outcome = outcome, rg, SE = rg_se, Lower_95CI = rg_lower_95CI, Upper_95CI = rg_upper_95CI,
  P = rg_p, Bonferroni_P = p_bonferroni_12, FDR_P = p_fdr_bh_11,
  Cross_trait_intercept = cross_trait_intercept, SNPs, Status = run_status
)]
write_tsv(table_3_8, file.path(tables_dir, "Table_3_8_LDSC_genetic_correlations.tsv"))
write_tsv(combined_diss, file.path(tables_dir, "Table_3_9_MR_LDSC_comparison.tsv"))

# ------------------------------ figures
# A publication-oriented Manhattan plot from a deterministic systematic sample plus all P <= 5e-9 rows.
# The source itself contains duplicate coordinate mappings; the first occurrence per SNP is retained consistently.
awk_cmd <- paste(
  "awk -F '\\t' 'BEGIN{OFS=\"\\t\"} NR==1 || $4<=5e-9 || NR%500==0 {print $1,$2,$3,$4}'",
  shQuote("results/pruning/asthma_fixed_for_sentinel_selection.txt")
)
man <- fread(cmd = awk_cmd, header = TRUE, col.names = c("SNP", "CHR", "POS", "P"))
man <- man[CHR %in% 1:22 & !is.na(POS) & !is.na(P) & P > 0]
setorder(man, SNP)
man <- man[, .SD[1], by = SNP]
setorder(man, CHR, POS)
chr_span <- man[, .(chr_len = as.numeric(max(POS))), by = CHR][order(CHR)]
chr_span[, offset := cumsum(shift(chr_len, fill = 0))]
man <- merge(man, chr_span[, .(CHR, offset)], by = "CHR")
man[, `:=`(BPcum = as.numeric(POS) + offset, logP = -log10(P))]
axis <- man[, .(centre = min(BPcum) + (max(BPcum) - min(BPcum))/2), by = CHR][order(CHR)]
top_plot <- man[SNP == "rs9273374"][which.max(logP)]
p_man <- ggplot(man, aes(BPcum, logP, colour = factor(CHR %% 2))) +
  geom_point(size = 0.45, alpha = 0.65) +
  geom_hline(yintercept = -log10(5e-9), linetype = 2, linewidth = 0.45) +
  geom_text(data = top_plot, aes(label = SNP), colour = "black", vjust = -0.7, size = 3.2) +
  scale_x_continuous(breaks = axis$centre, labels = axis$CHR, expand = expansion(mult = c(0.01, 0.01))) +
  scale_colour_manual(values = c("grey25", "grey65")) +
  labs(title = "Asthma GWAS meta-analysis", x = "Chromosome", y = expression(-log[10](P)),
       caption = "Dashed line: project signal-selection threshold, P = 5e-9. All variants at or below the threshold and a systematic background sample are displayed.") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.text.x = element_text(size = 8), plot.title.position = "plot")
ggsave(file.path(figures_dir, "Figure_3_1_Manhattan_plot.pdf"), p_man, width = 11.5, height = 5.8)
ggsave(file.path(figures_dir, "Figure_3_1_Manhattan_plot.png"), p_man, width = 11.5, height = 5.8, dpi = 300)

copy_as("results/gwama/qq_plot.png", file.path(figures_dir, "Figure_3_2_QQ_plot.png"))

flow <- rbind(
  data.table(Pathway = "Distance-based locus definition", Stage = c("Project-significant\nvariants", "Sentinel SNPs", "Merged loci"),
             Count = c(n_project, n_sentinels, n_loci), Order = 1:3),
  data.table(Pathway = "European-reference LD clumping", Stage = c("Project-significant\nvariants", "rsID-mapped\ncandidates", "Independent lead\nSNPs"),
             Count = c(n_project, n_rsid, n_leads), Order = 1:3)
)
flow[, Stage := factor(Stage, levels = unique(Stage))]
p_flow <- ggplot(flow, aes(Order, Count)) +
  geom_line(linewidth = 0.65, colour = "grey40") + geom_point(size = 3, colour = "#355f7a") +
  geom_text(aes(label = format(Count, big.mark = ",", scientific = FALSE)), vjust = -0.8, size = 3.5) +
  facet_wrap(~Pathway, scales = "free_x") +
  scale_x_continuous(breaks = 1:3, labels = c("Initial", "Refined", "Final")) +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0.02, 0.15))) +
  labs(title = "Refinement of asthma association signals", x = NULL, y = "Number of variants/features",
       caption = "Distance-based sentinel/locus definition and LD clumping were parallel refinement paths from the project-significant set.") +
  theme_classic(base_size = 11) + theme(strip.background = element_rect(fill = "grey92"), plot.title.position = "plot")
ggsave(file.path(figures_dir, "Figure_3_3_signal_refinement_flow.pdf"), p_flow, width = 9, height = 5.2)
ggsave(file.path(figures_dir, "Figure_3_3_signal_refinement_flow.png"), p_flow, width = 9, height = 5.2, dpi = 300)

# Figure 3.4 is intentionally not generated from build-mismatched VEP annotations.
writeLines(c(
  "Figure 3.4 VEP consequence distribution was not generated.",
  "Reason: the active VEP input used PLINK positions later validated as GRCh37/hg19,",
  "whereas VEP was run against GRCh38. A corrected annotation run is required before",
  "functional consequences can be presented as dissertation results. No annotation was rerun in this reporting task."
), file.path(figures_dir, "Figure_3_4_VEP_consequence_distribution_NOT_GENERATED.txt"))

for (ext in c("pdf", "png")) {
  copy_as(file.path("results/mr/finngen/summary/figures", paste0("finngen_ivw_all_outcomes.", ext)),
          file.path(figures_dir, paste0("Figure_3_5_MR_IVW_forest.", ext)))
  copy_as(file.path("results/ldsc/summary/figures", paste0("ldsc_asthma_rg_forest.", ext)),
          file.path(figures_dir, paste0("Figure_3_6_LDSC_genetic_correlation_forest.", ext)))
  copy_as(file.path("results/combined/figures", paste0("mr_ldsc_comparison_overview.", ext)),
          file.path(figures_dir, paste0("Figure_3_7_MR_LDSC_comparison.", ext)))
  copy_as(file.path("results/ldsc/summary/figures", paste0("ldsc_h2_qc.", ext)),
          file.path(appendix_dir, paste0("Appendix_figure_LDSC_h2_QC.", ext)))
}

# ------------------------------ indexes
figure_index <- data.table(
  Number = paste("Figure", paste0("3.", 1:7)),
  Title = c("Manhattan plot of asthma GWAS meta-analysis", "QQ plot of asthma GWAS meta-analysis",
            "Refinement of association signals to sentinel/locus and independent-lead sets",
            "Functional consequence distribution of independent lead variants",
            "Mendelian randomisation estimates for asthma and FinnGen outcomes",
            "LDSC genetic correlations between asthma and FinnGen outcomes",
            "Comparison of MR and LDSC evidence across outcomes"),
  Recommendation = c(rep("Essential main figure", 3), "Blocked—do not use until corrected VEP annotation exists",
                     "Essential main figure", "Essential main figure", "Main or supplementary"),
  File = c("figures/Figure_3_1_Manhattan_plot.pdf", "figures/Figure_3_2_QQ_plot.png",
           "figures/Figure_3_3_signal_refinement_flow.pdf",
           "figures/Figure_3_4_VEP_consequence_distribution_NOT_GENERATED.txt",
           "figures/Figure_3_5_MR_IVW_forest.pdf", "figures/Figure_3_6_LDSC_genetic_correlation_forest.pdf",
           "figures/Figure_3_7_MR_LDSC_comparison.pdf"),
  Status = c("READY", "READY", "READY", "NOT_READY_COORDINATE_BUILD_MISMATCH", "READY", "READY", "READY")
)
write_tsv(figure_index, file.path(out, "figure_index.tsv"))

table_index <- data.table(
  Number = paste("Table", paste0("3.", 1:9)),
  Title = c("Summary of asthma GWAS meta-analysis results", "Summary of sentinel SNP and locus identification",
            "rsID mapping and LD-clumping results", "Top functionally annotated independent lead variants",
            "FinnGen outcome extraction and MR harmonisation", "Primary IVW Mendelian randomisation results",
            "MR sensitivity analyses", "LDSC genetic-correlation results", "Comparison of MR and LDSC findings"),
  File = c("tables/Table_3_1_GWAMA_summary.tsv", "tables/Table_3_2_signal_refinement_summary.tsv",
           "tables/Table_3_3_LD_clumping_summary.tsv", "tables/Table_3_4_top_annotated_lead_variants.tsv",
           "tables/Table_3_5_MR_instrument_harmonisation_summary.tsv", "tables/Table_3_6_primary_MR_results.tsv",
           "tables/Table_3_7_MR_sensitivity_summary.tsv", "tables/Table_3_8_LDSC_genetic_correlations.tsv",
           "tables/Table_3_9_MR_LDSC_comparison.tsv"),
  Status = c(rep("READY", 3), "WITHHELD_ANNOTATIONS_COORDINATE_BUILD_MISMATCH", rep("READY", 5))
)
write_tsv(table_index, file.path(out, "table_index.tsv"))

# ------------------------------ master numerical audit sheet
numbers <- data.table(section = character(), metric = character(), value = character(), units = character(),
                      source_file = character(), source_column = character(), source_row_or_filter = character(),
                      rounded_display_value = character(), notes = character())
add_num <- function(section, metric, value, units, source, column, filter, display, notes = "") {
  numbers <<- rbind(numbers, data.table(section, metric, value = as.character(value), units, source_file = source,
                                        source_column = column, source_row_or_filter = filter,
                                        rounded_display_value = as.character(display), notes))
}
add_num("3.2", "Input GWAS", n_studies, "studies", "results/gwama/asthma_meta.log.out", "Study count", "Study count line", n_studies)
add_num("3.2", "Variants analysed", n_variants, "variants", "results/gwama/heterogeneity_summary.txt", "value", "Total SNPs in GWAMA output", "34,316,972")
add_num("3.2", "Variants with P < 5e-8", n_conventional, "variants", "results/gwama/genome_wide_significant.txt", "rows", "excluding header", format(n_conventional, big.mark = ","))
add_num("3.2", "Project-significant variants", n_project, "variants", "results/gwama/asthma_fixed_for_plink_clumping_mapping_summary.txt", "Value", "P <= 5e-9 row", "39,933")
for (item in list(c("Top SNP", top$rs_number, "SNP", "rs_number"), c("Top beta", top$beta, "beta", "beta"),
                  c("Top SE", top$se, "SE", "se"), c("Top P", top$`p-value`, "P", "p-value")))
  add_num("3.2", item[1], item[2], item[3], "results/gwama/top20_signals.txt", item[4], "minimum P row", if (item[1] == "Top P") "5.22 × 10^-253" else item[2])
add_num("3.2", "Top SNP chromosome", 6, "chromosome", "results/mr/finngen/qc/exposure_coordinate_repair_qc.tsv", "value", "rs9273374_GRCh38", 6,
        "Corrected GRCh38 coordinate; clumping used the older GRCh37 reference coordinate")
add_num("3.2", "Top SNP GRCh38 position", 32658837, "bp", "results/mr/finngen/qc/exposure_coordinate_repair_qc.tsv", "value", "rs9273374_GRCh38", "32,658,837",
        "PLINK/reference-panel coordinate was 32,626,614 (GRCh37)")
add_num("3.2", "lambda_GC", lambda_gc, "ratio", "results/gwama/lambda_gc.txt", "value", "Lambda GC", "1.3171")
add_num("3.3", "Variants with I2 >= 0.75", n_i2, "variants", "results/gwama/heterogeneity_summary.txt", "value", "I2 >= 0.75", "229,395")
add_num("3.3", "Percentage with I2 >= 0.75", pct_i2, "percent", "results/gwama/heterogeneity_summary.txt", "value", "I2 >= 0.75 percentage", "0.668%")
add_num("3.3", "Cochran Q P < 0.05", n_q, "variants", "results/gwama/heterogeneity_summary.txt", "value", "Q p-value < 0.05", "966,173")
add_num("3.3", "Cochran Q P < 0.05 percentage", pct_q, "percent", "results/gwama/heterogeneity_summary.txt", "value", "Q p-value < 0.05 percentage", "2.815%")
for (x in list(c("Sentinel SNPs", n_sentinels), c("Merged loci", n_loci),
               c("Median locus width", metric(loci_summary, "Median locus size (bp)")),
               c("Mean locus width", metric(loci_summary, "Mean locus size (bp)")),
               c("Largest locus", metric(loci_summary, "Largest locus (bp)")),
               c("Maximum sentinels per locus", metric(loci_summary, "Maximum sentinel SNPs within one locus"))))
  add_num("3.4", x[1], x[2], ifelse(grepl("width|locus$", x[1]), "bp", "count"),
          ifelse(x[1] == "Sentinel SNPs", "results/pruning/asthma_fixed_sentinels_clean.txt", "results/loci/loci_summary.txt"),
          ifelse(x[1] == "Sentinel SNPs", "rows", "Value"), x[1], format(as.numeric(x[2]), big.mark = ",", scientific = FALSE))
for (x in list(c("rsID mapping input", n_project), c("rsID matched", n_rsid), c("rsID unmatched", n_unmatched),
               c("Mapping rate", map_pct), c("Duplicate mappings collapsed", n_dup), c("Independent lead SNPs", n_leads)))
  add_num("3.5", x[1], x[2], ifelse(x[1] == "Mapping rate", "percent", "variants"),
          ifelse(x[1] == "Independent lead SNPs", "results/ld_clumping/independent_lead_snps.tsv", "results/gwama/asthma_fixed_for_plink_clumping_mapping_summary.txt"),
          ifelse(x[1] == "Independent lead SNPs", "rows", "Value"), x[1], x[2])
for (x in list(c("Clumping r2", 0.1), c("Clumping window", 99999), c("P1", 5e-9), c("P2", 5e-9)))
  add_num("3.5", x[1], x[2], ifelse(x[1] == "Clumping window", "kb", "parameter"),
          "results/ld_clumping/asthma_fixed_plink_clumped.log", "PLINK option", x[1], x[2])
add_num("3.6", "Lead SNPs", n_leads, "SNPs", "results/ld_clumping/independent_lead_snps.tsv", "rows", "all", "5,254")
add_num("3.6", "VEP input SNPs", vep_input_n, "SNPs", "results/annotation/independent_lead_snps_vep_input.txt", "rows", "all", "5,252", "Not reportable functionally due build mismatch")
add_num("3.6", "VEP returned SNPs", vep_unique_n, "SNPs", "results/annotation/independent_lead_snps_vep_clean_summary.txt", "SNP", "unique", "5,252", "Not reportable functionally due build mismatch")

for (key in c("total_instruments", "GRCh38_coordinates_recovered", "coordinates_unresolved", "old_position_differs_new")) {
  value <- coord_metric(key)
  add_num("3.7.1", key, value, "SNPs", "results/mr/finngen/qc/exposure_coordinate_repair_qc.tsv", "value", key,
          format(as.numeric(value), big.mark = ",", scientific = FALSE))
}
add_num("3.7.1", "GRCh38 coordinate recovery", as.numeric(coord_metric("GRCh38_coordinates_recovered"))/as.numeric(coord_metric("total_instruments"))*100,
        "percent", "results/mr/finngen/qc/exposure_coordinate_repair_qc.tsv", "value", "recovered / total", "99.70%")
for (i in seq_len(nrow(table_3_5))) {
  for (col in setdiff(names(table_3_5), "Outcome"))
    add_num("3.7.1", paste(table_3_5$Outcome[i], col), table_3_5[[col]][i], ifelse(grepl("percent", col), "percent", "SNPs"),
            ifelse(col == "FinnGen_matched_SNPs", "results/mr/finngen/qc/finngen_outcome_extraction_grch38_qc.tsv", "results/mr/finngen/qc/finngen_harmonisation_qc.tsv"),
            col, table_3_5$Outcome[i], ifelse(grepl("percent", col), sprintf("%.2f%%", table_3_5[[col]][i]), table_3_5[[col]][i]))
}
for (i in seq_len(nrow(mr))) {
  for (col in c("final_mr_ready_snps", "ivw_beta", "ivw_SE", "ivw_lower_95CI", "ivw_upper_95CI", "ivw_p_value", "ivw_p_bonferroni", "ivw_p_fdr_bh"))
    add_num("3.7.2", paste(mr$outcome[i], col), mr[[col]][i], col, "results/mr/finngen/summary/finngen_mr_final_summary.tsv", col,
            mr$outcome[i], ifelse(grepl("p", col, ignore.case = TRUE), fmt_p(mr[[col]][i]), fmt_n(mr[[col]][i])))
}
add_num("3.7.2", "Nominal IVW significant", sum(mr$ivw_p_value < 0.05), "outcomes", "results/mr/finngen/summary/finngen_mr_final_summary.tsv", "ivw_p_value", "P < 0.05", 0)
add_num("3.7.2", "Bonferroni IVW significant", sum(bool(mr$ivw_bonferroni_significant)), "outcomes", "results/mr/finngen/summary/finngen_mr_final_summary.tsv", "ivw_bonferroni_significant", "TRUE", 0)
for (i in seq_len(nrow(sens))) {
  for (col in c("ivw_heterogeneity_p_value", "egger_intercept", "egger_intercept_p_value", "p_egger_fdr", "leave_one_out_sign_change"))
    add_num("3.7.3", paste(sens$outcome[i], col), sens[[col]][i], col, "results/mr/finngen/sensitivity/finngen_mr_sensitivity_summary.tsv", col,
            sens$outcome[i], ifelse(grepl("p_value|fdr", col), fmt_p(sens[[col]][i]), as.character(sens[[col]][i])))
}
add_num("3.7.3", "Heterogeneity evidence", sum(bool(sens$heterogeneity_evidence)), "outcomes", "results/mr/finngen/sensitivity/finngen_mr_sensitivity_summary.tsv", "heterogeneity_evidence", "TRUE", 10)
add_num("3.7.3", "Egger intercept P < 0.05", sum(sens$egger_intercept_p_value < 0.05), "outcomes", "results/mr/finngen/sensitivity/finngen_mr_sensitivity_summary.tsv", "egger_intercept_p_value", "P < 0.05", 5)
add_num("3.7.3", "Egger intercept FDR significant", sum(sens$p_egger_fdr < 0.05), "outcomes", "results/mr/finngen/sensitivity/finngen_mr_sensitivity_summary.tsv", "p_egger_fdr", "P < 0.05", 5)
add_num("3.7.3", "LOO sign reversals", sum(bool(sens$leave_one_out_sign_change)), "outcomes", "results/mr/finngen/sensitivity/finngen_mr_sensitivity_summary.tsv", "leave_one_out_sign_change", "TRUE", 3)
for (i in seq_len(nrow(h2))) {
  for (col in c("h2_obs", "h2_se", "h2_z", "h2_p"))
    add_num("3.8.1", paste(h2$trait[i], col), h2[[col]][i], col, "results/ldsc/qc/ldsc_h2_summary.tsv", col, h2$trait[i],
            ifelse(col == "h2_p", fmt_p(h2[[col]][i]), fmt_n(h2[[col]][i])))
}
add_num("3.8.1", "Positive h2 traits", sum(h2$h2_obs > 0), "traits", "results/ldsc/qc/ldsc_h2_summary.tsv", "h2_obs", "> 0", 12)
for (i in seq_len(nrow(rg))) {
  for (col in c("rg", "rg_se", "rg_lower_95CI", "rg_upper_95CI", "rg_p", "p_bonferroni_12", "p_fdr_bh_11", "cross_trait_intercept", "SNPs"))
    add_num("3.8.2", paste(rg$outcome[i], col), rg[[col]][i], col, "results/ldsc/qc/ldsc_rg_summary.tsv", col, rg$outcome[i],
            ifelse(grepl("_p|p_", col), fmt_p(rg[[col]][i]), fmt_n(rg[[col]][i])))
}
add_num("3.8.2", "Planned rg", nrow(rg), "outcomes", "results/ldsc/qc/ldsc_rg_summary.tsv", "rows", "all", 12)
add_num("3.8.2", "Estimable rg", sum(!is.na(rg$rg)), "outcomes", "results/ldsc/qc/ldsc_rg_summary.tsv", "rg", "non-missing", 11)
add_num("3.8.2", "Bonferroni-significant rg", sum(bool(rg$bonferroni12_significant), na.rm = TRUE), "outcomes", "results/ldsc/qc/ldsc_rg_summary.tsv", "bonferroni12_significant", "yes", 7)
for (category in unique(combined$mr_ldsc_pattern))
  add_num("3.9", category, sum(combined$mr_ldsc_pattern == category), "outcomes", "results/combined/finngen_mr_ldsc_comparison.tsv", "mr_ldsc_pattern", category,
          sum(combined$mr_ldsc_pattern == category))
write_tsv(numbers, file.path(out, "results_numbers_master.tsv"))

# ------------------------------ result inventory
inventory <- data.table(
  section = c("3.1", "3.2", "3.3", "3.4", "3.5", "3.6", "3.7.1", "3.7.2", "3.7.3", "3.8.1", "3.8.2", "3.9"),
  analysis = c("Workflow overview", "Fixed-effect GWAMA", "Between-study heterogeneity", "Sentinel and locus definition",
               "rsID mapping and LD clumping", "VEP functional annotation audit", "MR instrument preparation",
               "Primary FinnGen MR", "MR sensitivity analyses", "Univariate LDSC QC", "LDSC genetic correlation", "MR–LDSC comparison"),
  operation_performed = c("Progressed 12 harmonised asthma GWAS through association and downstream analyses",
                          "Combined study-level effects", "Summarised I2 and Cochran Q", "Selected distance-based sentinels and merged overlapping windows",
                          "Mapped project-significant variants to rsIDs and clumped against European LD", "Audited complete VEP output provenance",
                          "Recovered GRCh38 positions, extracted outcomes, and harmonised alleles", "Estimated 12 IVW associations",
                          "Evaluated heterogeneity, Egger intercepts, and leave-one-out stability", "Estimated observed-scale h2 for 13 traits",
                          "Estimated asthma rg for eligible outcomes", "Joined primary MR and LDSC significance patterns"),
  design_feature = c("Sequential workflow; no pooled raw genotype data", "12 GWAS; fixed-effect; P <= 5e-9 project threshold",
                     "I2 >= 0.75 and Q P < 0.05", "P < 5e-9 and ±1 Mb sentinel windows", "1000G European reference; r2=0.1; 99,999 kb",
                     "Active VEP used GRCh37 PLINK coordinates with GRCh38 cache—invalid for reporting", "5,254 initial instruments; unavailable exposure EAF",
                     "IVW primary; 12-outcome Bonferroni family", "All retained SNPs; sensitivity findings did not alter primary set",
                     "GCST90302886 asthma, not 12-study GWAMA; observed scale", "12 planned; 11 estimable; 0.05/12 primary threshold",
                     "MR and rg kept on distinct inferential/effect scales"),
  authoritative_source_file = c("scripts/01_audit_inputs.R; results/gwama/asthma_meta.log.out", "results/gwama/asthma_meta.out; results/gwama/top20_signals.txt",
                                "results/gwama/heterogeneity_summary.txt", "results/pruning/asthma_fixed_sentinels_clean.txt; results/loci/loci_summary.txt",
                                "results/gwama/asthma_fixed_for_plink_clumping_mapping_summary.txt; results/ld_clumping/independent_lead_snps.tsv",
                                "scripts/17_prepare_vep_input.R; scripts/18_run_vep_annotation.sh; results/mr/finngen/qc/exposure_coordinate_repair_qc.tsv",
                                "results/mr/finngen/qc/exposure_coordinate_repair_qc.tsv; results/mr/finngen/qc/finngen_harmonisation_qc.tsv",
                                "results/mr/finngen/summary/finngen_mr_final_summary.tsv", "results/mr/finngen/sensitivity/finngen_mr_sensitivity_summary.tsv",
                                "results/ldsc/qc/ldsc_h2_summary.tsv", "results/ldsc/qc/ldsc_rg_summary.tsv", "results/combined/finngen_mr_ldsc_comparison.tsv"),
  main_result = c("12 GWAS progressed to GWAMA and downstream refinement", sprintf("%s variants; top SNP %s", format(n_variants, big.mark=","), top$rs_number),
                  sprintf("%s variants had I2 >= 0.75", format(n_i2, big.mark=",")), sprintf("%s sentinels; %s loci", n_sentinels, n_loci),
                  sprintf("%s mapped candidates; %s independent lead SNPs", n_rsid, n_leads), "Functional results withheld because coordinate builds mismatch",
                  "5,238/5,254 GRCh38 coordinates recovered; 623–652 MR-ready SNPs", "No nominal or Bonferroni-significant IVW results",
                  "10 heterogeneity flags; 5 FDR Egger flags; 3 LOO reversals", "12/13 positive h2; eosinophilic disease non-positive",
                  "7 Bonferroni-significant rg results", "0 MR versus 7 LDSC Bonferroni-significant outcomes"),
  statistical_result = c("Not applicable", sprintf("P=%s; lambda_GC=%.4f", top$`p-value`, lambda_gc), sprintf("%.3f%%; Q P<0.05 for %.3f%%", pct_i2, pct_q),
                         "P < 5e-9 selection", "P1=P2=5e-9", "No valid consequence statistics reported", "Retention approximately 80.43–80.60%",
                         "Minimum IVW P remained >=0.05", "Egger intercept FDR P < 0.05 for five outcomes", "Asthma h2=0.0423 (SE=0.0040)",
                         "Bonferroni threshold 0.0041667", "Seven LDSC_SIGNIFICANT_MR_NOT_SIGNIFICANT"),
  recommended_main_table = c("None", "Table 3.1", "Table 3.1", "Table 3.2", "Table 3.3", "Table 3.4a audit only",
                             "Table 3.5", "Table 3.6", "Table 3.7", "Appendix h2 table", "Table 3.8", "Table 3.9"),
  recommended_main_figure = c("Figure 3.3 can serve as compact workflow/refinement visual", "Figures 3.1 and 3.2", "None", "Figure 3.3", "Figure 3.3",
                              "Figure 3.4 blocked", "None", "Figure 3.5", "Supplementary diagnostics", "Appendix h2 QC", "Figure 3.6", "Figure 3.7 optional"),
  immediate_conclusion = c("The workflow yielded defined association, MR, and LDSC result sets.", "The strongest association was on chromosome 6.",
                           "A minority of variants met the high-I2 threshold.", "Distance-based refinement produced 405 merged loci.",
                           "LD clumping yielded 5,254 independent lead SNPs.", "Current annotations cannot support dissertation claims.",
                           "Coordinate repair enabled outcome matching; ambiguous palindromes were removed.", "No primary IVW estimate reached nominal significance.",
                           "Sensitivity analyses identified heterogeneity and directional pleiotropy flags.", "Eosinophilic disease did not support stable rg estimation.",
                           "Seven outcomes met the primary rg threshold.", "Seven rg findings occurred without primary MR significance."),
  logical_next_step = c("Association meta-analysis was summarised next.", "Association signals were refined and heterogeneity summarised.",
                        "Significant signals were grouped into sentinel-defined loci.", "Signals were mapped and LD-clumped.", "Independent lead SNPs entered annotation and MR instrument preparation.",
                        "Corrected VEP annotation is required outside this no-new-analysis task.", "Harmonised instruments entered primary MR.", "Sensitivity diagnostics assessed robustness.",
                        "Genome-wide genetic correlation was then evaluated.", "Eligible traits proceeded to rg.", "MR and LDSC patterns were compared.", "Results proceed to separate Discussion."),
  main_text_or_appendix = c("Main text", "Main text", "Main text", "Main text", "Main text", "Main-text audit limitation; invalid annotations excluded",
                            "Main text", "Main text", "Main text summary; full outputs appendix", "Brief main text; full table appendix", "Main text", "Main or supplementary")
)
write_tsv(inventory, file.path(out, "results_inventory.tsv"))

# ------------------------------ captions
fig_caps <- c(
  "Figure 3.1. Manhattan plot of the fixed-effect asthma GWAS meta-analysis. Points show −log10(P) by chromosome. The dashed horizontal line marks the project signal-selection threshold of P = 5 × 10⁻⁹. All variants meeting that threshold and a deterministic systematic sample of the remaining variants are displayed; rs9273374 is labelled as the strongest association.",
  sprintf("Figure 3.2. Quantile–quantile plot of asthma GWAS meta-analysis P-values. Observed and expected −log10(P) values are shown with the equality line; genomic inflation was lambda_GC = %.4f.", lambda_gc),
  sprintf("Figure 3.3. Parallel refinement of %s project-significant variants. Distance-based selection yielded %s sentinel SNPs and %s merged loci; European-reference LD clumping of %s rsID-mapped candidates yielded %s independent lead SNPs.", n_project, n_sentinels, n_loci, n_rsid, n_leads),
  "Figure 3.4. Not available. The active VEP output was generated using positions subsequently validated as GRCh37 with a GRCh38 annotation cache; consequence distributions are therefore withheld pending corrected annotation.",
  "Figure 3.5. Inverse-variance weighted Mendelian randomisation estimates for asthma and 12 FinnGen outcomes. Binary outcomes are displayed as odds ratios with 95% confidence intervals and continuous outcomes as beta coefficients with 95% confidence intervals on a separate panel-specific scale. No IVW result met the 12-outcome Bonferroni threshold.",
  "Figure 3.6. LDSC genetic correlations between asthma and FinnGen outcomes. Points denote genetic correlation (rg), horizontal bars denote 95% confidence intervals, and the vertical line marks rg = 0. Filled points met the primary Bonferroni threshold of 0.05/12. Eosinophilic disease was not estimated because its observed-scale h2 point estimate was non-positive.",
  "Figure 3.7. Categorical comparison of primary MR and LDSC evidence across 12 outcomes. MR classification uses IVW Bonferroni significance and LDSC classification uses rg Bonferroni significance. MR and LDSC answer different scientific questions and their effect estimates are not directly comparable."
)
writeLines(fig_caps, file.path(captions_dir, "results_figure_captions.txt"))

table_caps <- c(
  "Table 3.1. Summary of the 12-study fixed-effect asthma GWAS meta-analysis, genomic inflation, and between-study heterogeneity. I2 denotes the proportion of variation attributed to between-study heterogeneity.",
  "Table 3.2. Summary of project-significant variant refinement by distance-based sentinel selection and overlapping-window locus merging. The project signal-selection threshold was P < 5 × 10⁻⁹.",
  "Table 3.3. rsID mapping and European-reference LD clumping of project-significant asthma variants. LD, linkage disequilibrium; r2, squared allelic correlation. The final count denotes independent lead SNPs, not genomic loci.",
  "Table 3.4. Association-ranked independent lead SNPs reserved for functional annotation. Gene, consequence, and impact fields are withheld because the active VEP output has a coordinate-build mismatch; the table must not be interpreted as an annotation result.",
  "Table 3.5. FinnGen outcome extraction and harmonisation of 5,254 initial asthma instruments. Palindromic variants lacking exposure effect-allele frequency were removed as ambiguous; retention is relative to SNPs entering harmonisation.",
  "Table 3.6. Primary inverse-variance weighted Mendelian randomisation estimates for 12 FinnGen outcomes. Beta, effect estimate; SE, standard error; CI, confidence interval; OR, odds ratio for binary outcomes. Bonferroni and Benjamini–Hochberg FDR-adjusted P-values are shown.",
  "Table 3.7. MR sensitivity results across 12 outcomes. Heterogeneity was assessed for the IVW model; the MR-Egger intercept was used as a directional horizontal pleiotropy diagnostic; leave-one-out sign reversal records whether the pooled estimate changed sign.",
  "Table 3.8. LDSC genetic correlations between GCST90302886 asthma and 12 planned FinnGen outcomes. CI, confidence interval; rg, genetic correlation. Bonferroni adjustment retains the 12 planned hypotheses and FDR adjustment uses the 11 estimable results. Eosinophilic disease was not estimable because h2 was non-positive.",
  "Table 3.9. Descriptive comparison of primary IVW MR evidence and LDSC genetic-correlation evidence. Significance flags use the prespecified Bonferroni thresholds; the two effect measures are not directly comparable."
)
writeLines(table_caps, file.path(captions_dir, "results_table_captions.txt"))

# ------------------------------ writing blueprint
bp <- c(
  "# Results chapter writing blueprint", "", "Use the exact values below; keep interpretation factual and reserve biological explanation for Discussion.", "",
  "## 3.1 Overview of the analytical workflow", "", "**Operation carried out**", "- Twelve asthma GWAS were harmonised and meta-analysed, followed by signal refinement, FinnGen MR, and LDSC.",
  "", "**Design features that need mentioning**", "- Fixed-effect GWAMA; parallel distance-based locus definition and European-reference LD clumping.",
  "", "**Results/numbers that MUST appear**", "- 12 input GWAS; 34,316,972 meta-analysed variants; 5,254 LD-clumped independent lead SNPs; 12 FinnGen outcomes.",
  "", "**Statistical analysis**", "- State only the project P ≤ 5 × 10⁻⁹ selection threshold here.", "", "**Recommended table/figure**", "- Short prose; Figure 3.3 provides the refinement overview.",
  "", "**Immediate conclusion**", "- The workflow produced validated association, MR, and LDSC result sets.", "", "**Logical link to next analysis**", "- Begin with the fixed-effect asthma meta-analysis.", "",
  "## 3.2 Asthma GWAS meta-analysis", "", "**Operation carried out**", "- Fixed-effect association estimates were combined across 12 asthma GWAS.",
  "", "**Design features that need mentioning**", "- Distinguish conventional P < 5 × 10⁻⁸ from project signal selection at P ≤ 5 × 10⁻⁹.",
  "", "**Results/numbers that MUST appear**", sprintf("- %s variants; %s at P < 5 × 10⁻⁸; %s at P ≤ 5 × 10⁻⁹.", format(n_variants,big.mark=","), format(n_conventional,big.mark=","), format(n_project,big.mark=",")),
  sprintf("- Top association: %s, chromosome 6, corrected GRCh38 position 32,658,837, beta = %.6f, SE = %.6f, P = 5.22 × 10⁻²⁵³. The PLINK clumping coordinate was 32,626,614 in GRCh37.", top$rs_number, top$beta, top$se),
  sprintf("- lambda_GC = %.4f.", lambda_gc), "", "**Statistical analysis**", "- Fixed-effect beta, SE, z and P; QQ genomic inflation statistic.",
  "", "**Recommended table/figure**", "- Table 3.1; Figures 3.1 and 3.2.", "", "**Immediate conclusion**", "- The strongest association was the chromosome 6 rs9273374 signal.",
  "", "**Logical link to next analysis**", "- Quantify between-study heterogeneity and refine significant signals.", "",
  "## 3.3 Between-study heterogeneity and association signals", "", "**Operation carried out**", "- I2 and Cochran Q summaries were extracted from GWAMA output.",
  "", "**Design features that need mentioning**", "- High heterogeneity was defined in the existing pipeline as I2 ≥ 0.75.",
  "", "**Results/numbers that MUST appear**", sprintf("- %s variants (%.3f%%) had I2 ≥ 0.75.", format(n_i2,big.mark=","), pct_i2), sprintf("- %s variants (%.3f%%) had Cochran Q P < 0.05.", format(n_q,big.mark=","), pct_q),
  "", "**Statistical analysis**", "- Report I2 threshold and Q P-value criterion without adding new tests.", "", "**Recommended table/figure**", "- Table 3.1; no additional heterogeneity figure.",
  "", "**Immediate conclusion**", "- High I2 was observed for 0.668% of analysed variants.", "", "**Logical link to next analysis**", "- Project-significant variants were reduced to sentinel-defined loci.", "",
  "## 3.4 Identification of sentinel variants and associated loci", "", "**Operation carried out**", "- Distance-based sentinel selection and overlapping ±1 Mb window merging were applied.",
  "", "**Design features that need mentioning**", "- P < 5 × 10⁻⁹; ±1 Mb sentinel windows.", "", "**Results/numbers that MUST appear**",
  sprintf("- %s sentinel SNPs and %s merged loci.", format(n_sentinels,big.mark=","), n_loci),
  sprintf("- Median width %s bp; mean width %s bp; largest %s bp; maximum 16 sentinels in one locus.",
          format(as.numeric(metric(loci_summary,"Median locus size (bp)")),big.mark=","), format(as.numeric(metric(loci_summary,"Mean locus size (bp)")),big.mark=","), format(as.numeric(metric(loci_summary,"Largest locus (bp)")),big.mark=",")),
  "", "**Statistical analysis**", "- Association P ranking; descriptive locus-size summaries.", "", "**Recommended table/figure**", "- Table 3.2 and optional Table 3.2b; Figure 3.3.",
  "", "**Immediate conclusion**", "- Distance-based refinement grouped 1,563 sentinels into 405 loci.", "", "**Logical link to next analysis**", "- A separate rsID-mapping and LD-clumping path defined independent lead SNPs.", "",
  "## 3.5 LD-based refinement of independent association signals", "", "**Operation carried out**", "- Project-significant variants were mapped to rsIDs and clumped using a European reference panel.",
  "", "**Design features that need mentioning**", "- P1=P2=5 × 10⁻⁹, r2=0.1, window=99,999 kb.", "", "**Results/numbers that MUST appear**",
  sprintf("- %s evaluated; %s rsID matched; %s unmatched; %.2f%% matched; %s duplicate rows collapsed.", format(n_project,big.mark=","), format(n_rsid,big.mark=","), format(n_unmatched,big.mark=","), map_pct, n_dup),
  sprintf("- %s independent lead SNPs.", format(n_leads,big.mark=",")), "", "**Statistical analysis**", "- PLINK LD clumping parameters and result counts.",
  "", "**Recommended table/figure**", "- Table 3.3; Figure 3.3.", "", "**Immediate conclusion**", "- LD clumping yielded 5,254 independent lead SNPs; do not call these 5,254 loci.",
  "", "**Logical link to next analysis**", "- Lead SNPs were intended for functional annotation and MR instrument preparation.", "",
  "## 3.6 Functional annotation of independent lead variants", "", "**Operation carried out**", "- The active VEP provenance was audited; no annotation was rerun.",
  "", "**Design features that need mentioning**", "- PLINK positions were later validated as GRCh37/hg19, but the active VEP run specified GRCh38.",
  "", "**Results/numbers that MUST appear**", "- 5,254 lead SNPs existed; 5,252 rsID SNPs entered and returned from VEP; functional consequences and genes are withheld because coordinates were build-mismatched.",
  "", "**Statistical analysis**", "- None; this is a provenance/QC finding.", "", "**Recommended table/figure**", "- Table 3.4a audit. Do not use Figure 3.4 or gene/consequence claims until corrected annotation exists.",
  "", "**Immediate conclusion**", "- The current VEP annotations are not dissertation-ready.", "", "**Logical link to next analysis**", "- Use the independently validated coordinate-repaired instruments for FinnGen MR.", "",
  "## 3.7.1 Instrument preparation and harmonisation", "", "**Operation carried out**", "- GRCh38 coordinates were repaired, FinnGen associations extracted, and alleles harmonised.",
  "", "**Design features that need mentioning**", "- Exposure EAF was unavailable; unresolved palindromic variants were removed.", "", "**Results/numbers that MUST appear**",
  "- 5,238/5,254 coordinates recovered (99.70%); 16 unresolved; 3,047 positions changed.", "- R13: 810 matched SNPs; R12: 773–798; final MR-ready: 623–652; retention 80.43–80.60%.",
  "", "**Statistical analysis**", "- Descriptive extraction and harmonisation QC.", "", "**Recommended table/figure**", "- Table 3.5.", "", "**Immediate conclusion**", "- Coordinate repair enabled consistent local FinnGen matching and yielded outcome-specific MR sets.",
  "", "**Logical link to next analysis**", "- The retained instruments entered primary IVW MR.", "",
  "## 3.7.2 Primary MR results", "", "**Operation carried out**", "- IVW estimates were evaluated for 12 FinnGen outcomes.", "", "**Design features that need mentioning**", "- Binary outcomes use OR; continuous outcomes use beta; 12-outcome Bonferroni family.",
  "", "**Results/numbers that MUST appear**", "- 0/12 nominal IVW P < 0.05; 0/12 Bonferroni significant. Report all outcome estimates from Table 3.6.",
  "", "**Statistical analysis**", "- IVW beta, SE, 95% CI, P, Bonferroni P and FDR P.", "", "**Recommended table/figure**", "- Table 3.6; Figure 3.5.",
  "", "**Immediate conclusion**", "- No primary IVW estimate reached nominal significance; do not state that this proves absence of an effect.", "", "**Logical link to next analysis**", "- Sensitivity diagnostics were assessed next.", "",
  "## 3.7.3 Sensitivity analyses", "", "**Operation carried out**", "- IVW heterogeneity, MR-Egger intercept and leave-one-out diagnostics were summarised.",
  "", "**Design features that need mentioning**", "- Diagnostics were retained without SNP exclusion.", "", "**Results/numbers that MUST appear**", "- Heterogeneity: 10/12.",
  paste0("- Egger P < 0.05 and FDR significant: 5/12 (", paste(sens[p_egger_fdr < 0.05, outcome], collapse = "; "), ")."),
  paste0("- Leave-one-out sign reversal: 3/12 (", paste(sens[bool(leave_one_out_sign_change), outcome], collapse = "; "), ")."),
  "", "**Statistical analysis**", "- Q-test P, Egger intercept P and BH-FDR P, and sign reversal.", "", "**Recommended table/figure**", "- Table 3.7; detailed plots in appendix.",
  "", "**Immediate conclusion**", "- Sensitivity analyses identified widespread heterogeneity and directional horizontal pleiotropy flags in several outcomes.",
  "", "**Logical link to next analysis**", "- Genome-wide genetic correlation was assessed to quantify shared genetic architecture without assigning causality.", "",
  "## 3.8.1 Univariate LDSC quality control", "", "**Operation carried out**", "- Observed-scale SNP h2 was estimated for the LDSC asthma dataset and 12 outcomes.",
  "", "**Design features that need mentioning**", "- Asthma was GCST90302886: 48,623 cases, 290,722 controls, N=339,345; it was not the 12-study GWAMA.",
  "", "**Results/numbers that MUST appear**", "- Asthma h2=0.0423, SE=0.0040, Z=10.575; 12/13 traits had positive h2.", "- Eosinophilic disease h2=−0.0004, SE=0.0013, Z=−0.308 and was not eligible for stable rg.",
  "", "**Statistical analysis**", "- Univariate LDSC observed-scale h2 and SE; no liability-scale interpretation.", "", "**Recommended table/figure**", "- Brief prose; h2 QC table/figure in appendix.",
  "", "**Immediate conclusion**", "- Eleven outcome correlations were estimable; eosinophilic disease was retained as not estimable.", "", "**Logical link to next analysis**", "- Run-status-eligible traits were summarised using the completed rg results.", "",
  "## 3.8.2 Genetic correlation with asthma", "", "**Operation carried out**", "- Existing asthma–outcome rg estimates were consolidated.", "", "**Design features that need mentioning**", "- 12 planned hypotheses; 11 estimable; primary Bonferroni threshold=0.05/12.", "- The Bonferroni threshold was retained at 0.05/12, reflecting the pre-specified set of 12 planned outcomes, rather than recalculated as 0.05/11 after eosinophilic disease was excluded because of non-positive heritability, to avoid post hoc adjustment of the significance criterion.",
  "", "**Results/numbers that MUST appear**", paste0("- Seven Bonferroni- and FDR-significant outcomes: ", paste(rg[bool(bonferroni12_significant), outcome], collapse = "; "), "."),
  "- Strongest positive: allergic rhinitis rg=0.7114, SE=0.0547, 95% CI 0.6042–0.8186, P=1.3027×10⁻³⁸.", "- Strongest negative: HDL cholesterol rg=−0.1329, SE=0.0324, 95% CI −0.1964 to −0.0694, P=4.1846×10⁻⁵.",
  "", "**Statistical analysis**", "- rg, SE, 95% CI, P, 12-test Bonferroni P and 11-estimate BH-FDR P.", "", "**Recommended table/figure**", "- Table 3.8; Figure 3.6.",
  "", "**Immediate conclusion**", "- Seven outcomes showed significant genetic correlation with asthma.", "", "**Logical link to next analysis**", "- Compare the primary MR and LDSC evidence categories.", "",
  "## 3.9 Comparison of MR and genetic-correlation findings", "", "**Operation carried out**", "- Primary Bonferroni significance patterns were joined across MR and LDSC.",
  "", "**Design features that need mentioning**", "- MR and rg answer distinct questions and effect estimates are not directly compared.", "", "**Results/numbers that MUST appear**", "- MR significant 0/12; LDSC significant 7/12; both significant 0; LDSC significant/MR not significant 7; neither 4; LDSC not estimable 1.",
  paste0("- LDSC significant/MR not significant: ", paste(combined[mr_ldsc_pattern == "LDSC_SIGNIFICANT_MR_NOT_SIGNIFICANT", Outcome], collapse = "; "), "."),
  paste0("- Neither: ", paste(combined[mr_ldsc_pattern == "NEITHER_SIGNIFICANT", Outcome], collapse = "; "), "; not estimable: Eosinophilic disease."),
  "", "**Statistical analysis**", "- Descriptive cross-classification using prespecified Bonferroni flags.", "", "**Recommended table/figure**", "- Table 3.9; Figure 3.7 as main or supplementary.",
  "", "**Immediate conclusion**", "- Seven outcomes showed significant genetic correlation despite no significant primary MR association.", "", "**Logical link to next analysis**", "- End Results and move to a separate Discussion chapter."
)
writeLines(bp, file.path(out, "results_writing_blueprint.md"))

# ------------------------------ appendix plan
appendix_plan <- c(
  "# Appendix and supplementary-material plan", "",
  "Do not duplicate raw GWAS or FinnGen files. Link or cite the compact outputs below.", "",
  "## Recommended appendix tables", "",
  "- Full 405-locus table: `results/loci/asthma_fixed_sentinel_merge_loci.txt`.",
  "- Complete 5,254 independent-lead-SNP table: `results/ld_clumping/independent_lead_snps.tsv`.",
  "- Full MR five-method comparison: `results/mr/finngen/summary/finngen_mr_five_method_comparison.tsv`.",
  "- Full MR sensitivity, single-SNP, and leave-one-out tables under `results/mr/finngen/sensitivity/`.",
  "- Complete LDSC h2 QC table: `results/ldsc/qc/ldsc_h2_summary.tsv` and readiness table.",
  "- LDSC technical/environment audits only if examiner reproducibility detail is required.", "",
  "## Recommended appendix figures", "",
  "- LDSC h2 QC plot copied here as `Appendix_figure_LDSC_h2_QC.*`.",
  "- Outcome-specific MR scatter, funnel, leave-one-out, and single-SNP forest plots under `results/mr/finngen/plots/`.",
  "- MR sensitivity overview from `results/mr/finngen/summary/figures/`.", "",
  "## Annotation hold", "",
  "- Do not include the active full VEP output, gene table, consequence table, or annotation plot as scientific results.",
  "- The active VEP input used PLINK coordinates later validated as GRCh37/hg19, while VEP used GRCh38.",
  "- A corrected VEP run is required in a separately authorised analysis stage before annotation can enter the dissertation.", "",
  "## Technical records normally excluded from the main chapter", "",
  "- Coordinate-repair debugging details, P=0 underflow audit, HDL/LDL estimator audit, installation logs, and checksums.",
  "- Retain these in the project archive; include only concise methodological caveats when needed."
)
writeLines(appendix_plan, file.path(appendix_dir, "appendix_plan.md"))

# ------------------------------ validation report
discussion_terms <- c("confirms the central role", "proves common", "causes", "protective effect", "no causal effect", "shared pathway proven")
generated_text <- c(readLines(file.path(out, "results_writing_blueprint.md")), fig_caps, table_caps)
stopifnot(!any(vapply(discussion_terms, function(term) any(grepl(tolower(term), tolower(generated_text), fixed = TRUE)), logical(1))))
stopifnot(nrow(table_3_5) == 12, nrow(table_3_6) == 12, nrow(table_3_7) == 12,
          nrow(table_3_8) == 12, nrow(combined_diss) == 12)
stopifnot(all(table_3_5$Final_MR_ready_SNPs == mr$final_mr_ready_snps[match(table_3_5$Outcome, mr$outcome)]))
writeLines(c(
  "Dissertation Results evidence-pack validation",
  "Authoritative analytical files were read only; script 37 calls no analysis function.",
  "Archived results/old_position_issue was not read or used.",
  "Input GWAS count = 12; GWAMA variant count = 34,316,972.",
  "Independent lead SNP count = 5,254; terminology retained as independent lead SNPs, not loci.",
  "FinnGen MR outcomes = 12; MR Bonferroni significant = 0.",
  "LDSC planned outcomes = 12; estimable = 11; Bonferroni significant = 7.",
  "Eosinophilic disease rg = NOT_ESTIMABLE_LOW_H2.",
  "LDSC asthma source = GCST90302886, not the 12-study GWAMA.",
  "Binary and continuous MR scales remain in separate panels in Figure 3.5.",
  "DISCUSSION-style prohibited phrases audit = PASS.",
  "MATERIAL DISCREPANCY: active VEP positions were GRCh37/hg19 but VEP assembly was GRCh38; functional annotations withheld.",
  "Figure 3.4 status = NOT GENERATED pending corrected annotation."
), file.path(out, "validation_report.txt"))

cat("results_pack_validation=PASS\n")
cat("vep_annotation_status=WITHHELD_COORDINATE_BUILD_MISMATCH\n")
cat("tables=9 indexed; figures=6 ready plus 1 blocked\n")
