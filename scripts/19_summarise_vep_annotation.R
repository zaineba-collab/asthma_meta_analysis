library(data.table)

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

annotation_dir <- file.path(project_dir, "results", "annotation")

vep_file <- file.path(annotation_dir, "independent_lead_snps_vep_output.txt")
clean_file <- file.path(annotation_dir, "independent_lead_snps_vep_clean_summary.txt")
gene_file <- file.path(annotation_dir, "vep_gene_summary.txt")
consequence_file <- file.path(annotation_dir, "vep_consequence_summary.txt")

message("Reading VEP output...")

vep <- fread(vep_file, skip = "#Uploaded_variation")

# Clean column names
names(vep) <- gsub("^#", "", names(vep))

wanted_cols <- c(
  "Uploaded_variation",
  "Location",
  "Allele",
  "Gene",
  "Feature",
  "Feature_type",
  "Consequence",
  "IMPACT",
  "SYMBOL",
  "BIOTYPE",
  "VARIANT_CLASS",
  "GENE_PHENO"
)

available_cols <- intersect(wanted_cols, names(vep))

clean <- vep[, ..available_cols]

# Rename Uploaded_variation to SNP
if ("Uploaded_variation" %in% names(clean)) {
  setnames(clean, "Uploaded_variation", "SNP")
}

# Keep one row per SNP-gene-consequence combination
clean <- unique(clean)

# Gene summary
gene_summary <- clean[
  !is.na(SYMBOL) & SYMBOL != "-",
  .N,
  by = SYMBOL
][order(-N)]

# Consequence summary
consequence_summary <- clean[
  !is.na(Consequence) & Consequence != "-",
  .N,
  by = Consequence
][order(-N)]

fwrite(clean, clean_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(gene_summary, gene_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(consequence_summary, consequence_file, sep = "\t", quote = FALSE, na = "NA")

message("Written: ", clean_file)
message("Written: ", gene_file)
message("Written: ", consequence_file)

message("Rows in clean annotation table: ", nrow(clean))
message("Unique annotated genes: ", nrow(gene_summary))