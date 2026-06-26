library(data.table)

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

input_file <- file.path(project_dir, "results", "ld_clumping", "independent_lead_snps.txt")
bim_file <- file.path(project_dir, "reference_panel", "eur", "1000G_EUR_asthma_snps.bim")
annotation_dir <- file.path(project_dir, "results", "annotation")
dir.create(annotation_dir, showWarnings = FALSE, recursive = TRUE)

vep_file <- file.path(annotation_dir, "independent_lead_snps_vep_input.txt")
top20_vep_file <- file.path(annotation_dir, "top20_independent_lead_snps_vep_input.txt")
rsid_file <- file.path(annotation_dir, "independent_lead_snps_rsids.txt")
vep_install_notes <- file.path(annotation_dir, "vep_install_notes.txt")

lead <- fread(input_file)

required_cols <- c("CHR", "BP", "SNP", "P", "TOTAL")
missing_cols <- setdiff(required_cols, names(lead))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

lead[, CHR := as.integer(CHR)]
lead[, BP := as.integer(BP)]
lead[, P := as.numeric(P)]

lead <- lead[!is.na(CHR) & CHR %in% 1:22 & !is.na(BP)]

message("Reading PLINK BIM alleles...")

bim <- fread(
  bim_file,
  header = FALSE,
  select = c(1, 2, 4, 5, 6),
  col.names = c("CHR", "SNP", "BP", "A1", "A2")
)

bim <- bim[CHR %in% as.character(1:22)]
bim[, CHR := as.integer(CHR)]
bim[, BP := as.integer(BP)]

lead <- merge(
  lead,
  bim,
  by = c("CHR", "BP", "SNP"),
  all.x = TRUE,
  sort = FALSE
)

missing_alleles <- lead[is.na(A1) | is.na(A2)]

if (nrow(missing_alleles) > 0) {
  warning(
    nrow(missing_alleles),
    " lead SNPs were not found in the PLINK BIM file and will use allele N."
  )
}

lead[, allele := fifelse(is.na(A1) | is.na(A2), "N", paste0(A1, "/", A2))]

setorder(lead, P)

vep_input <- lead[, .(
  chromosome = CHR,
  start = BP,
  end = BP,
  allele,
  strand = "+",
  identifier = SNP
)]

rsids <- lead[, .(SNP)]

install_notes <- c(
  "VEP installer note",
  "",
  "The installer stopped while compiling HTSlib because lzma.h was not found.",
  "On this Mac, Homebrew xz is present at /opt/homebrew/opt/xz, so rerun the installer with these paths exported:",
  "",
  "cd tools/ensembl-vep",
  "export PATH=\"/opt/homebrew/opt/xz/bin:$PATH\"",
  "export CPPFLAGS=\"-I/opt/homebrew/opt/xz/include\"",
  "export LDFLAGS=\"-L/opt/homebrew/opt/xz/lib\"",
  "export PKG_CONFIG_PATH=\"/opt/homebrew/opt/xz/lib/pkgconfig:$PKG_CONFIG_PATH\"",
  "export DYLD_LIBRARY_PATH=\"$PWD/htslib:/opt/homebrew/opt/xz/lib:$DYLD_LIBRARY_PATH\"",
  "perl INSTALL.pl",
  "",
  "The DBD::mysql warning is not fatal if you run VEP offline with a local cache."
)

fwrite(vep_input, vep_file, sep = "\t", quote = FALSE, col.names = FALSE)
fwrite(vep_input[1:min(.N, 20)], top20_vep_file, sep = "\t", quote = FALSE, col.names = FALSE)
fwrite(rsids, rsid_file, sep = "\t", quote = FALSE, col.names = FALSE)
writeLines(install_notes, vep_install_notes)

message("Written: ", vep_file)
message("Written: ", top20_vep_file)
message("Written: ", rsid_file)
message("Written: ", vep_install_notes)
message("Lead SNPs prepared for VEP: ", nrow(vep_input))
