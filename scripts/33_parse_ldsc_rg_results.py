#!/usr/bin/env python3
"""Parse, validate, and multiple-test-correct asthma LDSC rg results."""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RG_DIR = ROOT / "results/ldsc/rg"
QC_DIR = ROOT / "results/ldsc/qc"
PRIMARY_BONF_THRESHOLD = 0.05 / 12
ESTIMABLE_BONF_THRESHOLD = 0.05 / 11

OUTCOMES = [
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

LAB_TRAITS = {"CRP", "ALT", "AST", "GGT", "HDL cholesterol", "LDL cholesterol", "Fasting triglycerides"}


def load_tsv(path, key):
    with path.open() as handle:
        return {row[key]: row for row in csv.DictReader(handle, delimiter="\t")}


def pair(section, label):
    match = re.search(rf"^{re.escape(label)}:\s+([-+0-9.eE]+)\s+\(([-+0-9.eE]+)\)", section, re.MULTILINE)
    return (float(match.group(1)), float(match.group(2))) if match else (None, None)


def scalar(section, label):
    match = re.search(rf"^{re.escape(label)}:\s+([-+0-9.eE]+)", section, re.MULTILINE)
    return float(match.group(1)) if match else None


def section(text, start, end=None):
    begin = text.find(start)
    if begin < 0:
        return ""
    begin += len(start)
    finish = text.find(end, begin) if end else len(text)
    return text[begin: finish if finish >= 0 else len(text)]


def finite(value):
    return value is not None and math.isfinite(value)


def bh_adjust(pvalues):
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    adjusted = [None] * n
    running = 1.0
    for rank_index in range(n - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, pvalues[original_index] * n / rank)
        adjusted[original_index] = min(running, 1.0)
    return adjusted


def main():
    h2_summary = load_tsv(QC_DIR / "ldsc_h2_summary.tsv", "trait")
    h2_readiness = load_tsv(QC_DIR / "ldsc_h2_readiness.tsv", "trait")
    overlap = load_tsv(QC_DIR / "ldsc_asthma_outcome_snp_overlap.tsv", "outcome")
    run_status = load_tsv(RG_DIR / "run_status.tsv", "outcome")
    parsed = {}
    errors = []

    for slug, outcome in OUTCOMES:
        if outcome == "Eosinophilic disease":
            continue
        log_path = RG_DIR / f"{slug}.log"
        exit_status = int(run_status.get(slug, {}).get("exit_status", 999))
        if exit_status != 0 or not log_path.exists():
            errors.append([outcome, "ldsc_rg", exit_status, "missing log or nonzero exit", "inspect logs; retain failed status"])
            continue
        text = log_path.read_text()
        h1_section = section(text, "Heritability of phenotype 1", "Heritability of phenotype 2/2")
        h2_section = section(text, "Heritability of phenotype 2/2", "Genetic Covariance")
        cov_section = section(text, "Genetic Covariance", "Genetic Correlation")
        rg_section = section(text, "Genetic Correlation", "Summary of Genetic Correlation Results")
        asthma_h2, _ = pair(h1_section, "Total Observed scale h2")
        outcome_h2, _ = pair(h2_section, "Total Observed scale h2")
        gencov, gencov_se = pair(cov_section, "Total Observed scale gencov")
        cross_intercept, cross_intercept_se = pair(cov_section, "Intercept")
        rg, rg_se = pair(rg_section, "Genetic Correlation")
        rg_z = scalar(rg_section, "Z-score")
        rg_p = scalar(rg_section, "P")
        snp_match = re.search(r"^([0-9]+) SNPs with valid alleles\.$", text, re.MULTILINE)
        snps = int(snp_match.group(1)) if snp_match else None
        warnings = re.findall(r"^.*WARNING.*$", text, re.MULTILINE)

        required = [asthma_h2, outcome_h2, gencov, gencov_se, cross_intercept,
                    cross_intercept_se, rg, rg_se, rg_z, rg_p]
        valid = all(finite(x) for x in required) and rg_se > 0 and gencov_se > 0 and cross_intercept_se > 0 and snps and 0 <= rg_p <= 1
        if not valid:
            errors.append([outcome, "parse_validation", exit_status, "missing/nonfinite rg fields or invalid SE/P/SNP count", "inspect LDSC log"])
            continue
        expected_snps = int(overlap[outcome]["shared_SNPs"])
        if snps != expected_snps:
            warnings.append(f"RG SNP count {snps} differs from audited overlap {expected_snps}")

        univariate_outcome_h2 = float(h2_summary[outcome]["h2_obs"])
        univariate_asthma_h2 = float(h2_summary["Asthma"]["h2_obs"])
        outcome_abs_diff = abs(outcome_h2 - univariate_outcome_h2)
        outcome_rel_diff = outcome_abs_diff / abs(univariate_outcome_h2) if univariate_outcome_h2 != 0 else math.inf
        asthma_rel_diff = abs(asthma_h2 - univariate_asthma_h2) / abs(univariate_asthma_h2)
        pairwise_discrepancy = outcome_abs_diff > 0.01 and outcome_rel_diff > 0.25
        if pairwise_discrepancy:
            warnings.append(f"PAIRWISE_H2_REVIEW: outcome h2 {outcome_h2} vs univariate {univariate_outcome_h2}")
        if asthma_rel_diff > 0.25:
            warnings.append(f"PAIRWISE_H2_REVIEW: asthma h2 {asthma_h2} vs univariate {univariate_asthma_h2}")
        if outcome == "Bronchiectasis":
            warnings.append("MODERATE_H2_CAUTION")
        if outcome == "ALT":
            warnings.append("MODERATE_H2_CAUTION")
        if outcome in LAB_TRAITS:
            warnings.append("MIXED_MODEL_EFFECTIVE_N_CAVEAT")
        lower, upper = rg - 1.96 * rg_se, rg + 1.96 * rg_se
        outside = abs(rg) > 1
        if lower < -1 or upper > 1:
            warnings.append("RG_CI_EXTENDS_OUTSIDE_THEORETICAL_BOUNDS")

        parsed[outcome] = {
            "outcome": outcome, "rg": rg, "rg_se": rg_se, "rg_z": rg_z, "rg_p": rg_p,
            "genetic_covariance": gencov, "genetic_covariance_se": gencov_se,
            "cross_trait_intercept": cross_intercept,
            "cross_trait_intercept_se": cross_intercept_se,
            "cross_trait_intercept_z": cross_intercept / cross_intercept_se,
            "SNPs": snps, "asthma_h2_pairwise": asthma_h2,
            "outcome_h2_pairwise": outcome_h2,
            "pairwise_h2_material_discrepancy": pairwise_discrepancy,
            "rg_lower_95CI": lower, "rg_upper_95CI": upper,
            "rg_outside_bounds": outside, "warning": "; ".join(dict.fromkeys(warnings)),
            "run_status": "PASS",
        }

    valid_outcomes = [outcome for _, outcome in OUTCOMES if outcome in parsed]
    adjusted = bh_adjust([parsed[outcome]["rg_p"] for outcome in valid_outcomes])
    for outcome, fdr in zip(valid_outcomes, adjusted):
        parsed[outcome]["p_fdr_bh_11"] = fdr
        parsed[outcome]["p_bonferroni_12"] = min(parsed[outcome]["rg_p"] * 12, 1.0)

    raw_fields = ["outcome", "rg", "rg_se", "rg_z", "rg_p", "genetic_covariance",
                  "genetic_covariance_se", "cross_trait_intercept", "cross_trait_intercept_se",
                  "cross_trait_intercept_z", "SNPs", "asthma_h2_pairwise", "outcome_h2_pairwise",
                  "pairwise_h2_material_discrepancy", "warning", "run_status"]
    with (QC_DIR / "ldsc_rg_results_raw.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=raw_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for outcome in valid_outcomes:
            row = parsed[outcome]
            writer.writerow({field: row[field] for field in raw_fields})

    all12_fields = ["outcome", "rg", "rg_se", "rg_z", "rg_p", "univariate_h2", "note", "run_status"]
    with (QC_DIR / "ldsc_rg_all_12_outcomes.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=all12_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for _, outcome in OUTCOMES:
            if outcome == "Eosinophilic disease":
                writer.writerow({"outcome": outcome, "rg": "NA", "rg_se": "NA", "rg_z": "NA",
                                 "rg_p": "NA", "univariate_h2": h2_summary[outcome]["h2_obs"],
                                 "note": "Standard rg not run because univariate LDSC h2 was non-positive",
                                 "run_status": "NOT_ESTIMABLE_LOW_H2"})
            elif outcome in parsed:
                row = parsed[outcome]
                writer.writerow({"outcome": outcome, "rg": row["rg"], "rg_se": row["rg_se"],
                                 "rg_z": row["rg_z"], "rg_p": row["rg_p"],
                                 "univariate_h2": h2_summary[outcome]["h2_obs"], "note": "",
                                 "run_status": "PASS"})
            else:
                writer.writerow({"outcome": outcome, "rg": "NA", "rg_se": "NA", "rg_z": "NA",
                                 "rg_p": "NA", "univariate_h2": h2_summary[outcome]["h2_obs"],
                                 "note": "Technical failure", "run_status": "TECHNICAL_FAILURE"})

    summary_fields = ["outcome", "outcome_type", "rg_estimable", "rg", "rg_se",
                      "rg_lower_95CI", "rg_upper_95CI", "rg_z", "rg_p", "p_bonferroni_12",
                      "p_fdr_bh_11", "bonferroni12_significant", "fdr11_significant",
                      "cross_trait_intercept", "cross_trait_intercept_se", "intercept_z", "SNPs",
                      "outcome_univariate_h2", "outcome_univariate_h2_z", "h2_readiness_category",
                      "rg_outside_bounds", "rg_warning", "run_status"]
    summary_rows = []
    for _, outcome in OUTCOMES:
        h2 = h2_summary[outcome]
        ready = h2_readiness[outcome]
        if outcome == "Eosinophilic disease":
            summary_rows.append({
                "outcome": outcome, "outcome_type": h2["trait_type"], "rg_estimable": "no",
                "rg": "NA", "rg_se": "NA", "rg_lower_95CI": "NA", "rg_upper_95CI": "NA",
                "rg_z": "NA", "rg_p": "NA", "p_bonferroni_12": "NA", "p_fdr_bh_11": "NA",
                "bonferroni12_significant": "NA", "fdr11_significant": "NA",
                "cross_trait_intercept": "NA", "cross_trait_intercept_se": "NA", "intercept_z": "NA",
                "SNPs": "NA", "outcome_univariate_h2": h2["h2_obs"],
                "outcome_univariate_h2_z": h2["h2_z"],
                "h2_readiness_category": ready["rg_readiness_category"], "rg_outside_bounds": "NA",
                "rg_warning": "Non-positive univariate h2; standard rg deliberately not run",
                "run_status": "NOT_ESTIMABLE_LOW_H2",
            })
            continue
        if outcome not in parsed:
            continue
        row = parsed[outcome]
        summary_rows.append({
            "outcome": outcome, "outcome_type": h2["trait_type"], "rg_estimable": "yes",
            "rg": row["rg"], "rg_se": row["rg_se"], "rg_lower_95CI": row["rg_lower_95CI"],
            "rg_upper_95CI": row["rg_upper_95CI"], "rg_z": row["rg_z"], "rg_p": row["rg_p"],
            "p_bonferroni_12": row["p_bonferroni_12"], "p_fdr_bh_11": row["p_fdr_bh_11"],
            "bonferroni12_significant": "yes" if row["rg_p"] < PRIMARY_BONF_THRESHOLD else "no",
            "fdr11_significant": "yes" if row["p_fdr_bh_11"] < 0.05 else "no",
            "cross_trait_intercept": row["cross_trait_intercept"],
            "cross_trait_intercept_se": row["cross_trait_intercept_se"],
            "intercept_z": row["cross_trait_intercept_z"], "SNPs": row["SNPs"],
            "outcome_univariate_h2": h2["h2_obs"], "outcome_univariate_h2_z": h2["h2_z"],
            "h2_readiness_category": ready["rg_readiness_category"],
            "rg_outside_bounds": "yes" if row["rg_outside_bounds"] else "no",
            "rg_warning": row["warning"], "run_status": row["run_status"],
        })
    with (QC_DIR / "ldsc_rg_summary.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=summary_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(summary_rows)

    with (QC_DIR / "ldsc_rg_multiple_testing.txt").open("w") as out:
        out.write(f"planned_hypotheses\t12\nprimary_bonferroni_threshold\t{PRIMARY_BONF_THRESHOLD:.17g}\n")
        out.write(f"estimated_hypotheses\t11\nestimable_only_bonferroni_threshold\t{ESTIMABLE_BONF_THRESHOLD:.17g}\n")
        out.write("primary_adjustment\tp_bonferroni_12=min(rg_p*12,1)\n")
        out.write("fdr_adjustment\tBenjamini-Hochberg across 11 estimated rg p-values only\n")

    with (QC_DIR / "ldsc_rg_errors.tsv").open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["outcome", "stage", "exit_status", "error", "action"])
        writer.writerows(errors)

    if len(summary_rows) != 12:
        errors.append(["ALL", "final_validation", 1, f"Expected 12 summary rows, found {len(summary_rows)}", "inspect parser"])
        raise SystemExit("Final rg summary does not contain exactly 12 rows")
    if errors:
        raise SystemExit(f"RG validation found {len(errors)} error(s)")


if __name__ == "__main__":
    main()
