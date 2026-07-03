library(data.table)

# ============================================================
# 19a_publication_annotation_table.R
#
# Purpose:
#   Create a dissertation/publication-ready annotation table for
#   the independent lead SNPs identified by PLINK LD clumping.
#
# Inputs:
#   results/ld_clumping/independent_lead_snps.tsv
#   results/annotation/independent_lead_snps_vep_clean_summary.txt
#
# Output:
#   results/annotation/lead_snp_annotation_table.tsv
#
# Notes:
#   - The lead SNP table has one row per independent clumped locus.
#   - VEP can return multiple annotation rows per SNP because one variant
#     can overlap multiple transcripts, genes, or consequence categories.
#   - This script collapses those VEP annotations so the final table has
#     one row per lead SNP, suitable for reporting and interpretation.
# ============================================================


script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

lead_file <- file.path(project_dir, "results", "ld_clumping", "independent_lead_snps.tsv")

# Cleaned VEP annotation table created by scripts/19_summarise_vep_annotation.R.
vep_clean_file <- file.path(
  project_dir,
  "results",
  "annotation",
  "independent_lead_snps_vep_clean_summary.txt"
)

# Final table for dissertation/results reporting.
output_file <- file.path(
  project_dir,
  "results",
  "annotation",
  "lead_snp_annotation_table.tsv"
)

lead <- fread(lead_file)
vep_clean <- fread(vep_clean_file)

# Collapse VEP annotations to one row per SNP. Multiple gene symbols,
# consequences, or impact terms are combined into comma-separated lists.
vep_one_row <- vep_clean[
    ,
    .(
        Gene_symbol = paste(unique(na.omit(SYMBOL[!is.na(SYMBOL) & SYMBOL != "-"])), collapse = ","),
        Consequence = paste(unique(na.omit(Consequence)), collapse = ","),
        Impact = paste(unique(na.omit(IMPACT)), collapse = ",")
    ),
    by = SNP
]

# Merge the lead SNP table with the collapsed VEP annotations
publication_annotated <- merge(
    lead,
    vep_one_row,
    by = "SNP",
    all.x = TRUE,
    sort = FALSE
)

# Rename the columns for clarity
setnames(
    publication_annotated,
    old = c("SNP", "chromosome", "position", "p_value"),
    new = c("Lead_SNP", "Chromosome", "Position", "P_value")
)

# Arrange columns in the order expected for the final reporting table.
setcolorder(
    publication_annotated,
    c("Lead_SNP", "Chromosome", "Position", "P_value", "Gene_symbol", "Consequence", "Impact")
)

# Sort by association strength so the most significant lead SNPs appear first.
setorder(publication_annotated, P_value)

# Write the annotated lead SNP table to a TSV file
fwrite(publication_annotated, output_file, sep = "\t", na = "NA", quote = FALSE)

message("Written: ", output_file)
message("Rows: ", nrow(publication_annotated))
