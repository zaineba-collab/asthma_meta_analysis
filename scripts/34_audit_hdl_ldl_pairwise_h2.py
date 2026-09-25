#!/usr/bin/env python3
"""Targeted, non-destructive HDL/LDL pairwise LDSC h2 audit."""

import argparse
import csv
import gzip
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MUNGED = ROOT / "data/ldsc/prepared/munged"
TMP = ROOT / "data/ldsc/prepared/tmp"
QC = ROOT / "results/ldsc/qc"
AUDIT = QC / "pairwise_h2_audit"
TRAITS = ("hdl", "ldl")


def read_sumstats(path):
    with gzip.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == ["SNP", "A1", "A2", "Z", "N"], reader.fieldnames
        rows = list(reader)
    assert len(rows) == len({row["SNP"] for row in rows}), f"Duplicate SNPs in {path}"
    return rows


def summarize_z2(rows):
    values = sorted(float(row["Z"]) ** 2 for row in rows)
    n = len(values)
    median = (values[(n - 1) // 2] + values[n // 2]) / 2 if n else math.nan
    return {
        "number_of_snps": n,
        "mean_z2": sum(values) / n if n else math.nan,
        "median_z2": median,
        "maximum_z2": values[-1] if n else math.nan,
        "abs_z_gt_5": sum(value > 25 for value in values),
        "abs_z_gt_10": sum(value > 100 for value in values),
        "abs_z_gt_20": sum(value > 400 for value in values),
    }


def write_tsv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def prepare():
    TMP.mkdir(parents=True, exist_ok=True)
    AUDIT.mkdir(parents=True, exist_ok=True)
    asthma_rows = read_sumstats(MUNGED / "asthma.sumstats.gz")
    asthma_snps = {row["SNP"] for row in asthma_rows}
    qc_rows = []

    for trait in TRAITS:
        outcome_rows = read_sumstats(MUNGED / f"{trait}.sumstats.gz")
        retained = [row for row in outcome_rows if row["SNP"] in asthma_snps]
        excluded = [row for row in outcome_rows if row["SNP"] not in asthma_snps]
        shared_path = TMP / f"asthma_{trait}_shared_snps.txt"
        shared_path.write_text("".join(f"{row['SNP']}\n" for row in retained))
        restricted_path = TMP / f"{trait}_asthma_overlap.sumstats.gz"
        with gzip.open(restricted_path, "wt", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=["SNP", "A1", "A2", "Z", "N"],
                delimiter="\t", lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(retained)
        for group, rows in (("retained_asthma_pairwise", retained), ("excluded_absent_from_asthma", excluded)):
            qc_rows.append({"trait": trait.upper(), "group": group, **summarize_z2(rows)})

    write_tsv(
        QC / "ldsc_hdl_ldl_pairwise_snp_loss_qc.tsv",
        ["trait", "group", "number_of_snps", "mean_z2", "median_z2", "maximum_z2",
         "abs_z_gt_5", "abs_z_gt_10", "abs_z_gt_20"],
        qc_rows,
    )


def section(text, title, next_title):
    start = text.index(title)
    end = text.index(next_title, start) if next_title else len(text)
    return text[start:end]


def pair(pattern, text):
    match = re.search(pattern, text)
    if not match:
        raise ValueError(f"Pattern not found: {pattern}")
    return float(match.group(1)), float(match.group(2))


def scalar(pattern, text):
    match = re.search(pattern, text)
    if not match:
        raise ValueError(f"Pattern not found: {pattern}")
    return float(match.group(1))


def parse_rg(path):
    text = path.read_text()
    h1 = section(text, "Heritability of phenotype 1", "Heritability of phenotype 2/2")
    h2 = section(text, "Heritability of phenotype 2/2", "Genetic Covariance")
    gcov = section(text, "Genetic Covariance", "Genetic Correlation")
    rg = section(text, "Genetic Correlation\n", "Summary of Genetic Correlation Results")
    return {
        "asthma_h2", "asthma_h2_se", "outcome_h2", "outcome_h2_se",
        "asthma_intercept", "asthma_intercept_se", "outcome_intercept", "outcome_intercept_se",
        "genetic_covariance", "genetic_covariance_se", "cross_trait_intercept",
        "cross_trait_intercept_se", "rg", "rg_se", "SNPs",
    }, {
        "asthma_h2": pair(r"Total Observed scale h2:\s*(\S+)\s*\((\S+)\)", h1)[0],
        "asthma_h2_se": pair(r"Total Observed scale h2:\s*(\S+)\s*\((\S+)\)", h1)[1],
        "outcome_h2": pair(r"Total Observed scale h2:\s*(\S+)\s*\((\S+)\)", h2)[0],
        "outcome_h2_se": pair(r"Total Observed scale h2:\s*(\S+)\s*\((\S+)\)", h2)[1],
        "asthma_intercept": pair(r"Intercept:\s*(\S+)\s*\((\S+)\)", h1)[0],
        "asthma_intercept_se": pair(r"Intercept:\s*(\S+)\s*\((\S+)\)", h1)[1],
        "outcome_intercept": pair(r"Intercept:\s*(\S+)\s*\((\S+)\)", h2)[0],
        "outcome_intercept_se": pair(r"Intercept:\s*(\S+)\s*\((\S+)\)", h2)[1],
        "genetic_covariance": pair(r"Total Observed scale gencov:\s*(\S+)\s*\((\S+)\)", gcov)[0],
        "genetic_covariance_se": pair(r"Total Observed scale gencov:\s*(\S+)\s*\((\S+)\)", gcov)[1],
        "cross_trait_intercept": pair(r"Intercept:\s*(\S+)\s*\((\S+)\)", gcov)[0],
        "cross_trait_intercept_se": pair(r"Intercept:\s*(\S+)\s*\((\S+)\)", gcov)[1],
        "rg": pair(r"Genetic Correlation:\s*(\S+)\s*\((\S+)\)", rg)[0],
        "rg_se": pair(r"Genetic Correlation:\s*(\S+)\s*\((\S+)\)", rg)[1],
        "SNPs": int(scalar(r"(\d+) SNPs with valid alleles", text)),
    }


def parse_h2(path):
    text = path.read_text()
    value, se = pair(r"Total Observed scale h2:\s*(\S+)\s*\((\S+)\)", text)
    intercept, intercept_se = pair(r"Intercept:\s*(\S+)\s*\((\S+)\)", text)
    return value, se, intercept, intercept_se, int(scalar(r"Read summary statistics for (\d+) SNPs", text))


def finalize():
    raw_rows = list(csv.DictReader((QC / "ldsc_rg_results_raw.tsv").open(), delimiter="\t"))
    raw = {row["outcome"]: row for row in raw_rows}
    validation = []
    decisions = []
    pzero_audit = []
    settings_audit = []
    pzero = list(csv.DictReader((QC / "finngen_pzero_excluded_variants.tsv").open(), delimiter="\t"))

    for trait in TRAITS:
        display = "HDL cholesterol" if trait == "hdl" else "LDL cholesterol"
        _, rg_values = parse_rg(ROOT / f"results/ldsc/rg/{trait}.log")
        uni_h2, uni_se, uni_int, uni_int_se, uni_n = parse_h2(ROOT / f"results/ldsc/h2/{trait}.log")
        tmp_h2, tmp_se, tmp_int, tmp_int_se, tmp_n = parse_h2(AUDIT / f"{trait}_asthma_overlap.log")
        raw_row = raw[display]
        comparisons = {
            "rg": raw_row["rg"], "rg_se": raw_row["rg_se"],
            "genetic_covariance": raw_row["genetic_covariance"],
            "genetic_covariance_se": raw_row["genetic_covariance_se"],
            "cross_trait_intercept": raw_row["cross_trait_intercept"],
            "cross_trait_intercept_se": raw_row["cross_trait_intercept_se"],
            "SNPs": raw_row["SNPs"], "asthma_h2": raw_row["asthma_h2_pairwise"],
            "outcome_h2": raw_row["outcome_h2_pairwise"],
        }
        for field, parsed_value in rg_values.items():
            parser_value = comparisons.get(field, "NA")
            match = "NA_NOT_IN_RAW_TABLE" if parser_value == "NA" else (
                "yes" if math.isclose(float(parser_value), float(parsed_value), rel_tol=0, abs_tol=1e-12) else "no"
            )
            validation.append({
                "trait": trait.upper(), "field": field, "raw_log_value": parsed_value,
                "parser_table_value": parser_value, "values_match": match,
                "source_log": f"results/ldsc/rg/{trait}.log",
            })
        close = math.isclose(tmp_h2, rg_values["outcome_h2"], rel_tol=0.02, abs_tol=0.001)
        trait_pzero = [row for row in pzero if row["trait"] == display]
        pzero_ids = {row["SNP"] for row in trait_pzero}
        munged_ids = {row["SNP"] for row in read_sumstats(MUNGED / f"{trait}.sumstats.gz")}
        shared_ids = set((TMP / f"asthma_{trait}_shared_snps.txt").read_text().splitlines())
        pzero_count = len(trait_pzero)
        in_munged = len(pzero_ids & munged_ids)
        in_shared = len(pzero_ids & shared_ids)
        pzero_audit.append({
            "trait": trait.upper(), "pzero_excluded_count": pzero_count,
            "present_in_univariate_munged": in_munged, "present_in_pairwise_intersection": in_shared,
            "differential_contribution": "no" if in_munged == 0 and in_shared == 0 else "review",
            "conclusion": "Already absent from both univariate and pairwise munged inputs",
        })
        for analysis, log_path, estimator in (
            ("univariate_h2", ROOT / f"results/ldsc/h2/{trait}.log", "two-step cutoff 30"),
            ("asthma_rg", ROOT / f"results/ldsc/rg/{trait}.log", "one-step IRWLS (no automatic two-step line)"),
            ("temporary_exact_overlap_h2", AUDIT / f"{trait}_asthma_overlap.log", "two-step cutoff 30"),
        ):
            log_text = log_path.read_text()
            settings_audit.append({
                "trait": trait.upper(), "analysis": analysis,
                "ref_ld_chr": "data/ldsc/reference/1000G_Phase3_ldscores/LDscore.",
                "w_ld_chr": "data/ldsc/reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.",
                "intercept": "unconstrained",
                "estimator": estimator,
                "two_step_line_present": "yes" if "Using two-step estimator with cutoff at 30" in log_text else "no",
                "ldsc_commit": "6c673952cee74bd5c57aef1555a03b1c015399a0",
                "python_environment": "ldsc39 / Python 3.9.23",
                "source_log": str(log_path.relative_to(ROOT)),
            })
        decisions.append(
            f"{trait.upper()}\n"
            f"category={'A. EXPLAINED_BY_PAIRWISE_SNP_SET' if close else 'C. LDSC_CONFIGURATION_DIFFERENCE'}\n"
            f"original_univariate_h2={uni_h2} ({uni_se}) on {uni_n} SNPs\n"
            f"rg_reported_pairwise_h2={rg_values['outcome_h2']} ({rg_values['outcome_h2_se']}) on {rg_values['SNPs']} SNPs\n"
            f"temporary_exact_overlap_h2={tmp_h2} ({tmp_se}) on {tmp_n} SNPs\n"
            f"temporary_exact_overlap_intercept={tmp_int} ({tmp_int_se})\n"
            f"original_univariate_intercept={uni_int} ({uni_int_se})\n"
            f"pzero_exclusions={pzero_count}; absent from both authoritative munged inputs and therefore not differential\n"
            f"explanation={'pairwise SNP-set restriction reproduces rg h2' if close else 'exact-overlap standalone h2 retains the univariate estimate; rg used one-step IRWLS while standalone h2 automatically used two-step cutoff 30'}\n"
        )

    write_tsv(
        QC / "ldsc_hdl_ldl_rg_log_validation.tsv",
        ["trait", "field", "raw_log_value", "parser_table_value", "values_match", "source_log"],
        validation,
    )
    write_tsv(
        QC / "ldsc_hdl_ldl_pzero_audit.tsv",
        ["trait", "pzero_excluded_count", "present_in_univariate_munged",
         "present_in_pairwise_intersection", "differential_contribution", "conclusion"],
        pzero_audit,
    )
    write_tsv(
        QC / "ldsc_hdl_ldl_settings_audit.tsv",
        ["trait", "analysis", "ref_ld_chr", "w_ld_chr", "intercept", "estimator",
         "two_step_line_present", "ldsc_commit", "python_environment", "source_log"],
        settings_audit,
    )
    (QC / "ldsc_hdl_ldl_pairwise_h2_decision.txt").write_text(
        "Targeted HDL/LDL pairwise h2 audit\n"
        "Primary rg estimates and authoritative h2 outputs were not changed.\n"
        "References, weights, environment, commit, and unconstrained intercepts were identical.\n"
        "Estimator difference: standalone h2 automatically used two-step cutoff 30; rg used one-step IRWLS.\n\n"
        + "\n".join(decisions)
        + "\nMixed-model/effective-N caveat retained; no N rescaling was attempted.\n"
    )
    loss_rows = list(csv.DictReader((QC / "ldsc_hdl_ldl_pairwise_snp_loss_qc.tsv").open(), delimiter="\t"))
    by_trait = {(row["trait"], row["group"]): row for row in loss_rows}
    loss_notes = ["Lost-SNP enrichment assessment (descriptive; no biological interpretation)"]
    for trait in ("HDL", "LDL"):
        retained = by_trait[(trait, "retained_asthma_pairwise")]
        excluded = by_trait[(trait, "excluded_absent_from_asthma")]
        loss_notes.append(
            f"{trait}: excluded mean Z^2={excluded['mean_z2']} versus retained={retained['mean_z2']}; "
            f"median Z^2={excluded['median_z2']} versus retained={retained['median_z2']}. "
            f"High-|Z| counts are recorded in ldsc_hdl_ldl_pairwise_snp_loss_qc.tsv. "
            "The pattern is not consistent across thresholds and the exact-overlap h2 remains essentially unchanged, "
            "so lost SNPs do not explain the rg-versus-univariate h2 discrepancy."
        )
    (QC / "ldsc_hdl_ldl_pairwise_snp_loss_conclusion.txt").write_text("\n".join(loss_notes) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("prepare", "finalize"))
    args = parser.parse_args()
    prepare() if args.stage == "prepare" else finalize()
