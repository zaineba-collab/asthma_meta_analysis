# ============================================================
# 02_harmonise_inputs.R
#
# Purpose:
#   Prepare raw asthma GWAS summary statistics for downstream
#   coordinate harmonisation and GWAMA input creation.
#
#   This script will:
#   1. Read the raw downloaded GWAS files.
#   2. Identify their column formats.
#   3. Standardise key columns.
#   4. Convert OR effects to log(OR) where needed.
#   5. Derive SE from 95% CI for files that do not provide SE.
#   6. Record genome build status.
#   7. Write harmonised *.gwama_ready.tsv files.
#
# Inputs:
#   data/raw/*
#
# Outputs:
#   logs/02_input_format_inventory.tsv
#   data/harmonised/*.gwama_ready.tsv
#
# Notes:
#   - Asthma is a binary trait. All effects are represented
#     internally on the log-odds scale as BETA with matching SE.
#   - GRCh37 files are labelled here but not lifted in this script.
#     Liftover is handled later by 03_prepare_liftover.R and
#     04_apply_liftover.R.
#   - OR-files from GCST9007* use coordinate/allele-based marker
#     names because their native variant IDs are not used here.
#
# Version notes:
#   - Filters to autosomes 1-22, valid p-values, positive SE, and
#     rows with complete marker/allele/effect fields.
#   - Writes one harmonised file per raw study.
# ============================================================

library(data.table)

# ============================================================
# STEP 1: Define project folders
# ============================================================

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

raw_dir <- file.path(project_dir, "data", "raw")
harmonised_dir <- file.path(project_dir, "data", "harmonised")
log_dir <- file.path(project_dir, "logs")

dir.create(harmonised_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# STEP 2: List raw input files
# ============================================================

files <- list.files(raw_dir, full.names = TRUE)

cat("Found", length(files), "raw input files:\n")
print(basename(files))

# ============================================================
# STEP 3: Create a function to read each file safely
# ============================================================

read_gwas_file <- function(file) {
  cat("\nReading:", basename(file), "\n")
  fread(file)
}

# ============================================================
# STEP 4: Create a small test read
# Purpose:
#   Before harmonising all data, check that every file can be
#   opened successfully by R/data.table.
# ============================================================

for (file in files) {
  test <- fread(file, nrows = 2)
  cat("OK:", basename(file), "| columns:", ncol(test), "\n")
}

cat("\nStep 4 complete: all files checked.\n")

# ============================================================
# STEP 5: Classify each file format
# Purpose:
#   Before harmonisation, record whether each file contains
#   BETA/SE data, OR/CI data, genome build information, and
#   whether it may need liftover to hg38.
# ============================================================

classify_file <- function(file) {
  dt <- fread(file, nrows = 5)
  cols <- names(dt)
  cols_lower <- tolower(cols)
  
  likely_build <- ifelse(
    grepl("grch38|buildgrch38|hg38", basename(file), ignore.case = TRUE),
    "GRCh38",
    ifelse(
      grepl("grch37|buildgrch37|hg19|b37|position_b37",
            paste(c(basename(file), cols), collapse = " "),
            ignore.case = TRUE),
      "GRCh37",
      "unknown"
    )
  )
  
  data.table(
    file = basename(file),
    likely_build = likely_build,
    has_beta = any(cols_lower %in% c("beta", "b")),
    has_se = any(cols_lower %in% c("se", "standard_error", "se_gc")),
    has_or = any(cols_lower %in% c("or", "odds_ratio")),
    has_ci = any(cols_lower %in% c("or_95l", "or_95u", "or_l95", "or_u95", "ci_lower", "ci_upper")),
    has_p = any(cols_lower %in% c("p", "p_value", "p_gc", "p_nospa")),
    columns = paste(cols, collapse = " | ")
  )
}

inventory <- rbindlist(lapply(files, classify_file), fill = TRUE)

inventory_file <- file.path(log_dir, "02_input_format_inventory.tsv")
fwrite(inventory, inventory_file, sep = "\t")

cat("\nFormat inventory written to:\n", inventory_file, "\n")
print(inventory)

# ============================================================
# STEP 6: Harmonise files to a common analysis-ready format
#
# Purpose:
#   Create one standard format across all asthma GWAS files.
#
#   Standard output columns:
#     MARKERNAME = SNP / variant identifier, or chr:pos:EA:NEA
#                  for OR-files where that is more consistent
#     EA         = effect allele
#     NEA        = non-effect / other allele
#     BETA       = effect size on log-odds scale
#     SE         = standard error of BETA
#     P          = p-value
#     N          = sample size, if available
#     CHR        = chromosome
#     POS        = base-pair position
#     BUILD      = genome build label
#     SOURCE_FILE = original file name
#
# Notes:
#   - For files with BETA/SE, BETA is used directly.
#   - For files with OR/SE, BETA = log(OR).
#   - Files with known GRCh37 are flagged for later liftover to hg38.
#   - Rows on non-autosomal chromosomes are removed at this stage.
# ============================================================

write_standard <- function(dt, source_file, build) {
  dt[, SOURCE_FILE := source_file]
  dt[, BUILD := build]
  
  dt <- dt[, .(
    MARKERNAME, EA, NEA, BETA, SE, P, N,
    CHR, POS, BUILD, SOURCE_FILE
  )]
  
  dt <- dt[
    !is.na(MARKERNAME) &
      !is.na(EA) &
      !is.na(NEA) &
      !is.na(BETA) &
      !is.na(SE) &
      !is.na(P)
  ]
  
  dt[, CHR := as.integer(CHR)]
  dt[, POS := as.integer(POS)]
  dt[, BETA := as.numeric(BETA)]
  dt[, SE := as.numeric(SE)]
  dt[, P := as.numeric(P)]
  
  dt <- dt[CHR %in% 1:22]
  dt <- dt[P > 0 & P <= 1]
  dt <- dt[SE > 0]
  
  out_file <- file.path(
    harmonised_dir,
    paste0(tools::file_path_sans_ext(source_file), ".gwama_ready.tsv")
  )
  
  fwrite(dt, out_file, sep = "\t", quote = FALSE, na = "NA")
  
  cat("Written:", out_file, "| rows:", nrow(dt), "\n")
}

# ============================================================
# STEP 7: Harmonise BETA/SE files
# ============================================================

# 7.1 495_PheCode
dt <- fread(file.path(raw_dir, "495_PheCode.v1.0.fastGWA.tsv"))
std <- dt[, .(
  MARKERNAME = variant_id,
  EA = effect_allele,
  NEA = other_allele,
  BETA = beta,
  SE = standard_error,
  P = p_value,
  N = N,
  CHR = chromosome,
  POS = base_pair_location
)]
write_standard(std, "495_PheCode.v1.0.fastGWA.tsv", "unknown")

# 7.2 asthma_DISCOVERY
dt <- fread(file.path(raw_dir, "asthma_DISCOVERY.tsv"))
std <- dt[, .(
  MARKERNAME = VARIANT_id,
  EA = effect_allele,
  NEA = other_allele,
  BETA = beta,
  SE = standard_error,
  P = p_value,
  N = N,
  CHR = chromosome,
  POS = base_pair_location
)]
write_standard(std, "asthma_DISCOVERY.tsv", "unknown")

# 7.3 GCST90029018, known GRCh37
dt <- fread(file.path(raw_dir, "GCST90029018_buildGRCh37.tsv"))
std <- dt[, .(
  MARKERNAME = variant_id,
  EA = effect_allele,
  NEA = other_allele,
  BETA = beta,
  SE = standard_error,
  P = p_value,
  N = n,
  CHR = chromosome,
  POS = base_pair_location
)]
write_standard(std, "GCST90029018_buildGRCh37.tsv", "GRCh37")

# 7.4 ukb.allasthma.final.assoc
dt <- fread(file.path(raw_dir, "ukb.allasthma.final.assoc"))
std <- dt[, .(
  MARKERNAME = SNP,
  EA = A1,
  NEA = A2,
  BETA = BETA,
  SE = SE,
  P = P,
  N = NA_integer_,
  CHR = CHR,
  POS = BP
)]
write_standard(std, "ukb.allasthma.final.assoc", "unknown")

# ============================================================
# STEP 8: Harmonise OR/SE files
#
# Purpose:
#   Several GWAS studies report odds ratios rather than beta
#   coefficients.
#
#   GWAMA can use ORs for binary traits, but we convert them
#   to the natural log scale here to maintain a consistent
#   internal representation across all studies.
#
#     BETA = log(OR)
#
#   The SE column is retained as the standard error for the
#   log-odds effect where provided.
# ============================================================

or_files <- c(
  "GCST90077671_buildGRCh38.tsv",
  "GCST90077672_buildGRCh38.tsv",
  "GCST90079456_buildGRCh38.tsv",
  "GCST90080141_buildGRCh38.tsv",
  "GCST90081448_buildGRCh38.tsv"
)

for (f in or_files) {
  dt <- fread(file.path(raw_dir, f))
  
  std <- dt[, .(
    MARKERNAME = paste(chromosome, base_pair_location, effect_allele, other_allele, sep = ":"),
    EA = effect_allele,
    NEA = other_allele,
    BETA = log(odds_ratio),
    SE = standard_error,
    P = p_value,
    N = NA_integer_,
    CHR = chromosome,
    POS = base_pair_location
  )]
  
  write_standard(std, f, "GRCh38")
}

# ============================================================
# STEP 9: Harmonise UKBB.asthma-2.assoc
#
# Scientific rationale:
#   This file reports asthma association results using odds
#   ratios because asthma is a binary phenotype.
#
#   To keep a single internal analysis-ready format, odds ratios
#   are converted onto the log-odds scale:
#
#       BETA = log(OR)
#
#   The file already contains SE and P columns, so these are
#   retained directly.
# ============================================================

dt <- fread(file.path(raw_dir, "UKBB.asthma-2.assoc"))

std <- dt[, .(
  MARKERNAME = SNP,
  EA = A1,
  NEA = A2,
  BETA = log(OR),
  SE = SE,
  P = P,
  N = NA_integer_,
  CHR = CHR,
  POS = POS
)]

write_standard(std, "UKBB.asthma-2.assoc", "unknown")

# ============================================================
# STEP 10: Harmonise HanY_prePMID_asthma_UKBB.txt.gz
#
# Scientific rationale:
#   This file reports odds ratios and 95% confidence intervals,
#   but does not provide a direct SE column.
#
#   We convert OR to log-odds:
#
#       BETA = log(OR)
#
#   and derive SE from the 95% confidence interval:
#
#       SE = (log(OR_95U) - log(OR_95L)) / (2 * 1.96)
#
#   This is the standard approximation when a log-odds confidence
#   interval is available.
# ============================================================

dt <- fread(file.path(raw_dir, "HanY_prePMID_asthma_UKBB.txt.gz"))

std <- dt[, .(
  MARKERNAME = SNP,
  EA = EA,
  NEA = NEA,
  BETA = log(OR),
  SE = (log(OR_95U) - log(OR_95L)) / (2 * 1.96),
  P = P,
  N = N,
  CHR = CHR,
  POS = BP
)]

write_standard(std, "HanY_prePMID_asthma_UKBB.txt.gz", "unknown")

# ============================================================
# STEP 11: Harmonise Shrine moderate-severe asthma file
#
# Scientific rationale:
#   This file contains both beta/SE and OR/CI columns.
#
#   We use beta and SE directly because they are already on the
#   log-odds scale for the coded/effect allele. This avoids
#   recalculating beta from OR when beta is already provided.
#
#   Genome build:
#   The position column is named Position_b37, so this file is
#   labelled GRCh37 and will require liftover to hg38 later.
# ============================================================

dt <- fread(file.path(raw_dir, "Shrine_30552067_moderate-severe_asthma.txt.gz"))

setnames(dt, old = "#SNP", new = "SNP")

std <- dt[, .(
  MARKERNAME = SNP,
  EA = Coded,
  NEA = Non_coded,
  BETA = beta,
  SE = SE_GC,
  P = P_GC,
  N = NA_integer_,
  CHR = Chromosome,
  POS = Position_b37
)]

write_standard(std, "Shrine_30552067_moderate-severe_asthma.txt.gz", "GRCh37")
