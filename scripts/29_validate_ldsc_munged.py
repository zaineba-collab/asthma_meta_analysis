#!/usr/bin/env python3
"""Validate munged LDSC files, parse logs, and calculate asthma SNP overlap."""

import csv
import gzip
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
QC = ROOT / "results/ldsc/qc"
LOGS = ROOT / "results/ldsc/logs/munging"
MUNGED = ROOT / "data/ldsc/prepared/munged"

TRAITS = [
    ("asthma", "Asthma"), ("allergic_rhinitis", "Allergic rhinitis"),
    ("bronchiectasis", "Bronchiectasis"), ("eosinophilic_disease", "Eosinophilic disease"),
    ("bmi_irn", "BMI IRN"), ("nafld", "NAFLD"), ("crp", "CRP"),
    ("alt", "ALT"), ("ast", "AST"), ("ggt", "GGT"),
    ("hdl", "HDL cholesterol"), ("ldl", "LDL cholesterol"),
    ("triglycerides_fasting", "Fasting triglycerides"),
]


def load_premunge():
    with (QC / "ldsc_premunging_qc.tsv").open() as handle:
        return {row["trait"]: row for row in csv.DictReader(handle, delimiter="\t")}


def metric(text, pattern, default="NA", flags=0):
    match = re.search(pattern, text, flags)
    return match.group(1) if match else default


def validate(path):
    seen = set()
    issues = []
    count = 0
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header != ["SNP", "A1", "A2", "Z", "N"]:
            issues.append("unexpected_header=" + "|".join(header))
        for lineno, line in enumerate(handle, 2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                issues.append(f"row_{lineno}_field_count")
                continue
            snp, a1, a2, z, n = fields
            count += 1
            if not re.fullmatch(r"rs[0-9]+", snp):
                issues.append(f"row_{lineno}_noncanonical_snp")
            if snp in seen:
                issues.append(f"row_{lineno}_duplicate_snp")
            seen.add(snp)
            if a1 not in "ACGT" or a2 not in "ACGT":
                issues.append(f"row_{lineno}_invalid_allele")
            try:
                if not math.isfinite(float(z)):
                    issues.append(f"row_{lineno}_invalid_z")
                if not math.isfinite(float(n)) or float(n) <= 0:
                    issues.append(f"row_{lineno}_invalid_n")
            except ValueError:
                issues.append(f"row_{lineno}_nonnumeric_z_or_n")
            if len(issues) > 20:
                break
    if count < 500000:
        issues.append("insubstantial_hm3_retention")
    return count, seen, issues


def main():
    prem = load_premunge()
    rows = []
    snp_sets = {}
    for slug, trait in TRAITS:
        pre = prem[trait]
        path = MUNGED / f"{slug}.sumstats.gz"
        log_path = LOGS / f"{slug}.log"
        if pre["status"] != "ready" or not path.exists():
            rows.append({
                "trait": trait, "input_HM3_rows": pre["final_premunge_count"],
                "output_SNPs": "NA", "percent_retained": "NA",
                "variants_removed_missing": "NA", "variants_removed_info": "NA",
                "variants_removed_maf": "NA", "variants_removed_p": "NA",
                "variants_removed_alleles": "NA", "variants_removed_duplicates": "NA",
                "strand_ambiguous_removed": "NA", "mean_chisq": "NA",
                "lambda_gc": "NA", "max_chisq": "NA", "N": pre["N"],
                "munge_exit_status": "not_run", "validation_status": "BLOCKED_PREMUNGING",
                "warning": f"{pre['invalid_P']} P<=0 rows; review results/ldsc/qc/ldsc_invalid_pvalues.tsv",
            })
            continue
        text = log_path.read_text()
        output_count, snps, issues = validate(path)
        snp_sets[slug] = snps
        input_count = int(pre["final_premunge_count"])
        warnings = re.findall(r"^.*WARNING.*$", text, re.MULTILINE)
        if issues:
            warnings.extend(issues)
        rows.append({
            "trait": trait, "input_HM3_rows": input_count, "output_SNPs": output_count,
            "percent_retained": f"{100 * output_count / input_count:.6f}",
            "variants_removed_missing": metric(text, r"Removed ([0-9]+) SNPs with missing values"),
            "variants_removed_info": metric(text, r"Removed ([0-9]+) SNPs with INFO"),
            "variants_removed_maf": metric(text, r"Removed ([0-9]+) SNPs with MAF"),
            "variants_removed_p": metric(text, r"Removed ([0-9]+) SNPs with out-of-bounds p-values"),
            "variants_removed_alleles": metric(text, r"Removed ([0-9]+) SNPs whose alleles did not match"),
            "variants_removed_duplicates": metric(text, r"Removed ([0-9]+) SNPs with duplicated rs numbers"),
            "strand_ambiguous_removed": metric(text, r"Removed ([0-9]+) variants that were not SNPs or were strand-ambiguous"),
            "mean_chisq": metric(text, r"Mean chi\^2 = ([0-9.eE+-]+)"),
            "lambda_gc": metric(text, r"Lambda GC = ([0-9.eE+-]+)"),
            "max_chisq": metric(text, r"Max chi\^2 = ([0-9.eE+-]+)"),
            "N": pre["N"], "munge_exit_status": "0",
            "validation_status": "PASS" if not issues else "FAIL",
            "warning": "; ".join(warnings) if warnings else "",
        })

    fields = ["trait", "input_HM3_rows", "output_SNPs", "percent_retained",
              "variants_removed_missing", "variants_removed_info", "variants_removed_maf",
              "variants_removed_p", "variants_removed_alleles", "variants_removed_duplicates",
              "strand_ambiguous_removed", "mean_chisq", "lambda_gc", "max_chisq", "N",
              "munge_exit_status", "validation_status", "warning"]
    with (QC / "ldsc_munged_summary_qc.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)

    munged_by_trait = {row["trait"]: row for row in rows}
    report_fields = ["Trait", "Raw_variants", "HM3_premunging_variants", "Final_munged_SNPs",
                     "HM3_coverage_percent", "Munging_retention_percent", "Mean_chi_square",
                     "Lambda_GC", "Max_chi_square", "N", "Status"]
    with (QC / "ldsc_13_trait_preparation_munging_status.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=report_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for _, trait in TRAITS:
            pre = prem[trait]
            munged = munged_by_trait[trait]
            writer.writerow({
                "Trait": trait, "Raw_variants": pre["raw_variant_count"],
                "HM3_premunging_variants": pre["final_premunge_count"],
                "Final_munged_SNPs": munged["output_SNPs"],
                "HM3_coverage_percent": f"{float(pre['percent_HM3_coverage']):.6f}",
                "Munging_retention_percent": munged["percent_retained"],
                "Mean_chi_square": munged["mean_chisq"], "Lambda_GC": munged["lambda_gc"],
                "Max_chi_square": munged["max_chisq"], "N": pre["N"],
                "Status": munged["validation_status"],
            })

    asthma = snp_sets.get("asthma", set())
    with (QC / "ldsc_asthma_outcome_snp_overlap.tsv").open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["outcome", "asthma_SNPs", "outcome_SNPs", "shared_SNPs", "percent_asthma_shared", "percent_outcome_shared"])
        for slug, trait in TRAITS[1:]:
            outcome = snp_sets.get(slug)
            if outcome is None:
                writer.writerow([trait, len(asthma), "NA", "NA", "NA", "NA"])
            else:
                shared = len(asthma.intersection(outcome))
                writer.writerow([trait, len(asthma), len(outcome), shared,
                                 f"{100 * shared / len(asthma):.6f}",
                                 f"{100 * shared / len(outcome):.6f}"])

    recovery_summary_path = QC / "finngen_pzero_recovery_summary.tsv"
    if recovery_summary_path.exists():
        with recovery_summary_path.open() as handle:
            recovery = {row["trait"]: row for row in csv.DictReader(handle, delimiter="\t")}
        with (QC / "ldsc_underflow_exclusion_impact.tsv").open("w", newline="") as out:
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            writer.writerow(["trait", "final_munged_SNPs", "underflow_variants_excluded",
                             "excluded_fraction_percent", "whether_excluded_fraction_gt_0.01_percent"])
            for trait, rec in recovery.items():
                final_count = int(munged_by_trait[trait]["output_SNPs"])
                excluded = int(rec["excluded_unrepresentable"])
                fraction = 100 * excluded / (final_count + excluded)
                writer.writerow([trait, final_count, excluded, f"{fraction:.9f}",
                                 "yes" if fraction > 0.01 else "no"])

    audit_path = QC / "finngen_pzero_underflow_audit.tsv"
    with (QC / "finngen_pzero_munged_validation.tsv").open("w", newline="") as out:
        fields = ["trait", "SNP", "BETA", "P_LDSC", "munged_Z", "finite_Z",
                  "sign_matches_BETA", "status"]
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        if audit_path.exists():
            with audit_path.open() as handle:
                recovered = [row for row in csv.DictReader(handle, delimiter="\t")
                             if row["reconstruction_status"] == "RECOVERABLE_NORMAL"]
            munged_values = {}
            for slug, trait in TRAITS:
                if any(row["trait"] == trait for row in recovered):
                    with gzip.open(MUNGED / f"{slug}.sumstats.gz", "rt") as handle:
                        for row in csv.DictReader(handle, delimiter="\t"):
                            munged_values[(trait, row["SNP"])] = float(row["Z"])
            for row in recovered:
                z = munged_values.get((row["trait"], row["SNP"]))
                beta = float(row["beta"])
                finite = z is not None and math.isfinite(z)
                sign_match = finite and ((z > 0) == (beta > 0))
                writer.writerow({"trait": row["trait"], "SNP": row["SNP"],
                                 "BETA": row["beta"], "P_LDSC": row["P_reconstructed"],
                                 "munged_Z": z if z is not None else "", "finite_Z": finite,
                                 "sign_matches_BETA": sign_match,
                                 "status": "PASS" if finite and sign_match else "FAIL"})


if __name__ == "__main__":
    main()
