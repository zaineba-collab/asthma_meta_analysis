library(data.table)

script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- dirname(script_dir)

lead_file <- file.path(project_dir, "results", "ld_clumping", "independent_lead_snps.tsv")
vep_clean_file <- file.path(
    project_dir,
    "results",
    "annotation",
    "independent_lead_snps_vep_clean_summary.txt"
)
gwama_file <- file.path(project_dir, "results", "gwama", "asthma_meta.out")
output_file <- file.path(project_dir, "results", "mr", "asthma_exposure_twosamplemr.tsv")

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

lead <- fread(lead_file)
vep_clean <- fread(vep_clean_file)
gwama <- fread(gwama_file)

# Map the current chromosome:position:ref:alt lead SNP IDs to dbSNP rsIDs
# reported by VEP. TwoSampleMR/OpenGWAS works much more reliably with rsIDs.
rsid_map <- unique(vep_clean[
    Existing_variation != "-" & !is.na(Existing_variation),
    .(original_snp.exposure = SNP, SNP = Existing_variation)
])

if (nrow(rsid_map) == 0) {
    stop("No rsIDs found in VEP Existing_variation column. Rerun VEP with --check_existing.", call. = FALSE)
}

if (anyDuplicated(rsid_map$original_snp.exposure)) {
    stop("Some lead SNPs map to multiple rsIDs; inspect VEP Existing_variation before continuing.", call. = FALSE)
}

exposure <- gwama[
    rs_number %in% lead$SNP,
    .(
        original_snp.exposure = rs_number,
        beta.exposure = beta,
        se.exposure = se,
        effect_allele.exposure = other_allele,
        other_allele.exposure = reference_allele,
        eaf.exposure = eaf,
        pval.exposure = `p-value`,
        samplesize.exposure = n_samples,
        exposure = "Asthma"
    )
]

exposure <- merge(
    exposure,
    rsid_map,
    by = "original_snp.exposure",
    all.x = TRUE,
    sort = FALSE
)

missing_rsids <- exposure[is.na(SNP)]
if (nrow(missing_rsids) > 0) {
    stop("Missing rsID mapping for ", nrow(missing_rsids), " exposure SNP(s).", call. = FALSE)
}

exposure[eaf.exposure == -9, eaf.exposure := NA]
exposure[samplesize.exposure == -9, samplesize.exposure := NA]

setcolorder(
    exposure,
    c(
        "SNP",
        "original_snp.exposure",
        "beta.exposure",
        "se.exposure",
        "effect_allele.exposure",
        "other_allele.exposure",
        "eaf.exposure",
        "pval.exposure",
        "samplesize.exposure",
        "exposure"
    )
)

setorder(exposure, pval.exposure)

fwrite(exposure, output_file, sep = "\t", na = "NA", quote = FALSE)

message("Written: ", output_file)
message("Exposure SNPs: ", nrow(exposure))
