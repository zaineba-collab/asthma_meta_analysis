#!/usr/bin/env python3
"""Parse and validate the 13 univariate LDSC h2 diagnostic runs."""

from __future__ import annotations

import csv
import math
import re
import statistics
from pathlib import Path

from scipy.stats import norm


ROOT = Path(__file__).resolve().parents[1]
H2_DIR = ROOT / "results/ldsc/h2"
QC_DIR = ROOT / "results/ldsc/qc"

TRAITS = [
    ("asthma", "Asthma"),
    ("allergic_rhinitis", "Allergic rhinitis"),
    ("bronchiectasis", "Bronchiectasis"),
    ("eosinophilic_disease", "Eosinophilic disease"),
    ("bmi_irn", "BMI IRN"),
    ("nafld", "NAFLD"),
    ("crp", "CRP"),
    ("alt", "ALT"),
    ("ast", "AST"),
    ("ggt", "GGT"),
    ("hdl", "HDL cholesterol"),
    ("ldl", "LDL cholesterol"),
    ("triglycerides_fasting", "Fasting triglycerides"),
]


def load_tsv(path, key):
    with path.open() as handle:
        return {row[key]: row for row in csv.DictReader(handle, delimiter="\t")}


def number(text, pattern):
    match = re.search(pattern, text, re.MULTILINE)
    return float(match.group(1)) if match else None


def integer(text, pattern):
    match = re.search(pattern, text, re.MULTILINE)
    return int(match.group(1)) if match else None


def pair(text, label):
    pattern = rf"^{re.escape(label)}:\s+([-+0-9.eE]+)\s+\(([-+0-9.eE]+)\)"
    match = re.search(pattern, text, re.MULTILINE)
    return (float(match.group(1)), float(match.group(2))) if match else (None, None)


def finite(value):
    return value is not None and math.isfinite(value)


def yn(value):
    return "yes" if value else "no"


def main():
    premunging = load_tsv(QC_DIR / "ldsc_premunging_qc.tsv", "trait")
    munging = load_tsv(QC_DIR / "ldsc_munged_summary_qc.tsv", "trait")
    run_status = load_tsv(H2_DIR / "run_status.tsv", "trait")
    summary_rows = []
    readiness_rows = []
    errors = []

    parsed = []
    for slug, trait in TRAITS:
        log_path = H2_DIR / f"{slug}.log"
        status = run_status.get(slug, {})
        exit_status = int(status.get("exit_status", 999))
        if not log_path.exists() or exit_status != 0:
            errors.append([trait, "ldsc_h2", "missing log or nonzero exit", exit_status, "inspect run logs; do not run rg"])
            parsed.append((slug, trait, None))
            continue
        text = log_path.read_text()
        munged_snps = integer(text, r"Read summary statistics for ([0-9]+) SNPs\.")
        ref_merged = integer(text, r"After merging with reference panel LD, ([0-9]+) SNPs remain\.")
        regression_snps = integer(text, r"After merging with regression SNP LD, ([0-9]+) SNPs remain\.")
        h2_obs, h2_se = pair(text, "Total Observed scale h2")
        intercept, intercept_se = pair(text, "Intercept")
        ratio, ratio_se = pair(text, "Ratio")
        lambda_gc = number(text, r"^Lambda GC:\s+([-+0-9.eE]+)")
        mean_chisq = number(text, r"^Mean Chi\^2:\s+([-+0-9.eE]+)")
        warnings = re.findall(r"^.*WARNING.*$", text, re.MULTILINE)
        prior_warning = munging[trait].get("warning", "").strip()
        if prior_warning:
            warnings.append("Premunging LDSC warning: " + prior_warning)
        warning = "; ".join(dict.fromkeys(warnings))

        required_finite = [h2_obs, h2_se, lambda_gc, mean_chisq, intercept, intercept_se]
        valid = (all(finite(x) for x in required_finite) and h2_se > 0 and intercept_se > 0
                 and munged_snps is not None and ref_merged is not None and regression_snps is not None
                 and regression_snps > 0)
        if not valid:
            errors.append([trait, "parse_validation", "missing/nonfinite estimate, nonpositive SE, or missing SNP count", exit_status, "inspect LDSC log; do not run rg"])
            parsed.append((slug, trait, None))
            continue
        h2_z = h2_obs / h2_se
        h2_p = float(2 * norm.sf(abs(h2_z)))
        if not (finite(h2_z) and finite(h2_p) and 0 <= h2_p <= 1):
            errors.append([trait, "derived_statistics", "invalid h2 Z or two-sided P", exit_status, "inspect calculation"])
            parsed.append((slug, trait, None))
            continue
        parsed.append((slug, trait, {
            "munged_snps": munged_snps, "ref_merged": ref_merged,
            "regression_snps": regression_snps, "h2_obs": h2_obs, "h2_se": h2_se,
            "h2_z": h2_z, "h2_p": h2_p, "lambda_gc": lambda_gc,
            "mean_chisq": mean_chisq, "intercept": intercept,
            "intercept_se": intercept_se, "ratio": ratio, "ratio_se": ratio_se,
            "warning": warning, "exit_status": exit_status,
        }))

    successful_counts = [item[2]["regression_snps"] for item in parsed if item[2] is not None]
    median_count = statistics.median(successful_counts) if successful_counts else None

    for slug, trait, result in parsed:
        pre = premunging[trait]
        if result is None:
            summary_rows.append({
                "trait": trait, "trait_type": pre["trait_type"], "release": pre["release"],
                "N": pre["N"], "cases": pre["cases"], "controls": pre["controls"],
                "munged_SNPs": "NA", "regression_SNPs": "NA", "h2_obs": "NA",
                "h2_se": "NA", "h2_z": "NA", "h2_p": "NA", "lambda_gc": "NA",
                "mean_chisq": "NA", "intercept": "NA", "intercept_se": "NA",
                "ratio": "NA", "ratio_se": "NA", "warning": "technical failure",
                "run_status": "FAIL",
            })
            readiness_rows.append({
                "trait": trait, "positive_h2": "NA", "h2_nominally_different_from_zero": "NA",
                "weak_h2_signal": "NA", "very_low_mean_chisq": "NA",
                "intercept_above_1": "NA", "intercept_below_1": "NA",
                "intercept_deviation": "NA", "ldsc_warning_present": "yes",
                "regression_SNP_count_unexpectedly_low": "NA", "rg_numerically_possible": "NO",
                "rg_readiness_category": "E — TECHNICAL_FAILURE",
            })
            continue

        ratio_value = "NA" if result["ratio"] is None else f"{result['ratio']:.17g}"
        ratio_se_value = "NA" if result["ratio_se"] is None else f"{result['ratio_se']:.17g}"
        summary_rows.append({
            "trait": trait, "trait_type": pre["trait_type"], "release": pre["release"],
            "N": pre["N"], "cases": pre["cases"], "controls": pre["controls"],
            "munged_SNPs": result["munged_snps"], "regression_SNPs": result["regression_snps"],
            "h2_obs": f"{result['h2_obs']:.17g}", "h2_se": f"{result['h2_se']:.17g}",
            "h2_z": f"{result['h2_z']:.17g}", "h2_p": f"{result['h2_p']:.17g}",
            "lambda_gc": f"{result['lambda_gc']:.17g}",
            "mean_chisq": f"{result['mean_chisq']:.17g}",
            "intercept": f"{result['intercept']:.17g}",
            "intercept_se": f"{result['intercept_se']:.17g}",
            "ratio": ratio_value, "ratio_se": ratio_se_value,
            "warning": result["warning"], "run_status": "PASS",
        })

        abs_z = abs(result["h2_z"])
        positive = result["h2_obs"] > 0
        major_warning = bool(result["warning"])
        if not positive:
            category = "D — H2_OUT_OF_BOUNDS"
        elif abs_z >= 4 and not major_warning:
            category = "A — STRONGER_H2_SIGNAL"
        elif abs_z >= 2:
            category = "B — MODERATE_H2_SIGNAL"
        else:
            category = "C — WEAK_H2_SIGNAL"
        low_count = median_count is not None and result["regression_snps"] < 0.95 * median_count
        readiness_rows.append({
            "trait": trait, "positive_h2": yn(positive),
            "h2_nominally_different_from_zero": yn(result["h2_p"] < 0.05),
            "weak_h2_signal": yn(abs_z < 2),
            "very_low_mean_chisq": yn(result["mean_chisq"] <= 1.02),
            "intercept_above_1": yn(result["intercept"] > 1),
            "intercept_below_1": yn(result["intercept"] < 1),
            "intercept_deviation": f"{result['intercept'] - 1:.17g}",
            "ldsc_warning_present": yn(major_warning),
            "regression_SNP_count_unexpectedly_low": yn(low_count),
            "rg_numerically_possible": "YES" if positive else "NO",
            "rg_readiness_category": category,
        })

    summary_fields = ["trait", "trait_type", "release", "N", "cases", "controls",
                      "munged_SNPs", "regression_SNPs", "h2_obs", "h2_se", "h2_z", "h2_p",
                      "lambda_gc", "mean_chisq", "intercept", "intercept_se", "ratio", "ratio_se",
                      "warning", "run_status"]
    with (QC_DIR / "ldsc_h2_summary.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=summary_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(summary_rows)

    readiness_fields = ["trait", "positive_h2", "h2_nominally_different_from_zero",
                        "weak_h2_signal", "very_low_mean_chisq", "intercept_above_1",
                        "intercept_below_1", "intercept_deviation", "ldsc_warning_present",
                        "regression_SNP_count_unexpectedly_low", "rg_numerically_possible",
                        "rg_readiness_category"]
    with (QC_DIR / "ldsc_h2_readiness.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=readiness_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(readiness_rows)

    with (QC_DIR / "ldsc_h2_errors.tsv").open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["trait", "stage", "error", "exit_status", "action"])
        writer.writerows(errors)

    if errors:
        raise SystemExit(f"Validation found {len(errors)} error(s); inspect ldsc_h2_errors.tsv")


if __name__ == "__main__":
    main()
