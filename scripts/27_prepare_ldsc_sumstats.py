#!/usr/bin/env python3
"""Prepare 13 plaintext HapMap3 summary-statistic files for CBIIT LDSC.

Run from the repository root. Large FinnGen gzip sources are streamed. The
only decompressed working files are HM3-restricted candidates under
data/ldsc/prepared/tmp; full GWAS files are never decompressed to disk.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import heapq
import math
import random
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
HM3_PATH = ROOT / "data/ldsc/reference/w_hm3.snplist"
HM3_MERGE_PATH = ROOT / "data/ldsc/prepared/w_hm3_merge_alleles.tsv"
PREP = ROOT / "data/ldsc/prepared"
RAW_HM3 = PREP / "raw_hm3"
TMP = PREP / "tmp"
QC = ROOT / "results/ldsc/qc"
MAP_PATH = PREP / "finngen_r13_grch38_hm3_variant_map.tsv.gz"
RS_RE = re.compile(r"^rs[0-9]+$")
RS_SPLIT_RE = re.compile(r"[,;|\s]+")
VALID_ALLELES = {"A", "C", "G", "T"}

R13 = [
    ("allergic_rhinitis", "Allergic rhinitis", "data/finngen/r13_disease/finngen_R13_ALLERG_RHINITIS.gz", "case_control", 489568, 15958, 473610),
    ("bronchiectasis", "Bronchiectasis", "data/finngen/r13_disease/finngen_R13_J10_BRONCHIECTASIS.gz", "case_control", 409082, 3096, 405986),
    ("eosinophilic_disease", "Eosinophilic disease", "data/finngen/r13_disease/finngen_R13_EOSINOPHIL_DISEASE.gz", "case_control", 264712, 571, 264141),
    ("bmi_irn", "BMI IRN", "data/finngen/r13_disease/finngen_R13_BMI_IRN.gz", "continuous", 362216, None, None),
    ("nafld", "NAFLD", "data/finngen/r13_disease/finngen_R13_NAFLD.gz", "case_control", 500186, 3943, 496243),
]

R12 = [
    ("crp", "CRP", "data/finngen/lab_values/finngen_R12_CRP_3020460.gz", 319548),
    ("alt", "ALT", "data/finngen/lab_values/finngen_R12_ALT_3006923.gz", 413996),
    ("ast", "AST", "data/finngen/lab_values/finngen_R12_AST_3013721.gz", 124818),
    ("ggt", "GGT", "data/finngen/lab_values/finngen_R12_GGT_3026910.gz", 257813),
    ("hdl", "HDL cholesterol", "data/finngen/lab_values/finngen_R12_HDL_3023602.gz", 361195),
    ("ldl", "LDL cholesterol", "data/finngen/lab_values/finngen_R12_LDL_3001308.gz", 368976),
    ("triglycerides_fasting", "Fasting triglycerides", "data/finngen/lab_values/finngen_R12_TRIGLYCERIDES_3048773.gz", 347199),
]

QC_COLUMNS = [
    "trait", "source_file", "release", "trait_type", "N", "cases", "controls",
    "raw_variant_count", "canonical_rsid_count", "HM3_identifier_matches",
    "HM3_exact_allele_matches", "palindromic_SNPs", "duplicate_SNPs",
    "missing_beta", "missing_SE", "invalid_SE", "missing_P", "invalid_P",
    "missing_A1", "missing_A2", "missing_FRQ", "invalid_FRQ",
    "final_premunge_count", "percent_HM3_coverage", "status", "notes",
    "coordinate_matches", "exact_REF_ALT_matches", "reversed_REF_ALT_candidates",
    "unmapped_rows", "mapping_conflicts",
]


def fail(message: str) -> None:
    raise RuntimeError(message)


def fnum(value: str):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def valid_chr(chrom: str) -> bool:
    return chrom in {str(i) for i in range(1, 23)}


def palindromic(a1: str, a2: str) -> bool:
    return {a1, a2} in ({"A", "T"}, {"C", "G"})


def load_hm3():
    alleles = {}
    pair_seen = set()
    bad = duplicate = duplicate_pair = 0
    with HM3_PATH.open() as handle:
        for lineno, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split()
            if len(fields) != 3:
                fail(f"Malformed HM3 row {lineno}")
            snp, a1, a2 = fields
            a1, a2 = a1.upper(), a2.upper()
            if not RS_RE.fullmatch(snp) or a1 not in VALID_ALLELES or a2 not in VALID_ALLELES:
                bad += 1
            if snp in alleles:
                duplicate += 1
            if (snp, a1, a2) in pair_seen:
                duplicate_pair += 1
            alleles[snp] = (a1, a2)
            pair_seen.add((snp, a1, a2))
    if bad or duplicate or duplicate_pair:
        fail(f"Invalid HM3 reference: bad={bad}, duplicate={duplicate}, duplicate_pair={duplicate_pair}")
    with HM3_MERGE_PATH.open("w") as out:
        out.write("SNP\tA1\tA2\n")
        for snp, (a1, a2) in alleles.items():
            out.write(f"{snp}\t{a1}\t{a2}\n")
    return alleles


def base_qc(trait, source, release, trait_type, n, cases=None, controls=None):
    row = {key: 0 for key in QC_COLUMNS}
    row.update({
        "trait": trait, "source_file": source, "release": release,
        "trait_type": trait_type, "N": n, "cases": cases if cases is not None else "NA",
        "controls": controls if controls is not None else "NA", "status": "pending",
        "notes": "", "percent_HM3_coverage": 0.0,
    })
    return row


def assess_values(qc, snp, a1, a2, beta, se, p, frq, invalid_rows):
    b = fnum(beta)
    s = fnum(se)
    pv = fnum(p)
    fq = fnum(frq) if frq not in (None, "") else None
    if b is None:
        qc["missing_beta"] += 1
    if se in (None, ""):
        qc["missing_SE"] += 1
    elif s is None or s <= 0:
        qc["invalid_SE"] += 1
    if p in (None, ""):
        qc["missing_P"] += 1
    elif pv is None or pv <= 0 or pv > 1:
        qc["invalid_P"] += 1
        invalid_rows.append([qc["trait"], snp, beta, se, p, "missing/non-numeric or P<=0/P>1"])
    if not a1:
        qc["missing_A1"] += 1
    if not a2:
        qc["missing_A2"] += 1
    if frq in (None, ""):
        qc["missing_FRQ"] += 1
    elif fq is None or fq < 0 or fq > 1:
        qc["invalid_FRQ"] += 1
    return b, s, pv, fq


def prepare_asthma(hm3, invalid_rows):
    source = "data/raw/asthma_DISCOVERY.tsv"
    qc = base_qc("Asthma", source, "GCST90302886", "case_control", 339345, 48623, 290722)
    out_path = RAW_HM3 / "asthma_hm3.tsv"
    seen = set()
    with (ROOT / source).open(newline="") as inp, out_path.open("w", newline="") as out:
        reader = csv.DictReader(inp, delimiter="\t")
        required = {"effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value", "VARIANT_id", "N"}
        if not required.issubset(reader.fieldnames or []):
            fail(f"Asthma header missing: {sorted(required - set(reader.fieldnames or []))}")
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "A1", "A2", "BETA", "SE", "P", "FRQ", "N"])
        for row in reader:
            qc["raw_variant_count"] += 1
            snp = row["VARIANT_id"]
            if not RS_RE.fullmatch(snp):
                continue
            qc["canonical_rsid_count"] += 1
            if snp not in hm3:
                continue
            qc["HM3_identifier_matches"] += 1
            a1, a2 = row["effect_allele"].upper(), row["other_allele"].upper()
            b, s, pv, fq = assess_values(qc, snp, a1, a2, row["beta"], row["standard_error"], row["p_value"], row["effect_allele_frequency"], invalid_rows)
            if {a1, a2} != set(hm3[snp]) or a1 not in VALID_ALLELES or a2 not in VALID_ALLELES:
                continue
            qc["HM3_exact_allele_matches"] += 1
            if palindromic(a1, a2):
                qc["palindromic_SNPs"] += 1
            if snp in seen:
                qc["duplicate_SNPs"] += 1
                continue
            seen.add(snp)
            n = fnum(row["N"])
            if b is None or s is None or s <= 0 or pv is None or pv <= 0 or pv > 1 or fq is None or fq < 0 or fq > 1 or n is None or n <= 0:
                continue
            writer.writerow([snp, a1, a2, row["beta"], row["standard_error"], row["p_value"], row["effect_allele_frequency"], int(n)])
            qc["final_premunge_count"] += 1
    qc["percent_HM3_coverage"] = 100 * qc["final_premunge_count"] / len(hm3)
    qc["status"] = "ready" if qc["duplicate_SNPs"] == 0 and qc["invalid_P"] == 0 else "blocked"
    qc["notes"] = "BETA is aligned to effect_allele (A1); variant-specific N preserved (expected 339345). Counts for missing/invalid fields are among direct HM3 identifier matches."
    return qc


def scan_r13(hm3, invalid_rows):
    observations = {}  # (SNP, chrom, pos, ref, alt) -> source bitmask
    ambiguous = []
    qcs = {}
    candidate_header = ["SNP", "chrom", "pos", "ref", "alt", "BETA", "SE", "P", "FRQ", "N"]
    for source_index, (slug, trait, source, trait_type, n, cases, controls) in enumerate(R13):
        qc = base_qc(trait, source, "R13", trait_type, n, cases, controls)
        qcs[slug] = qc
        tmp_path = TMP / f"{slug}_r13_hm3_candidates.tsv"
        with gzip.open(ROOT / source, "rt", newline="") as inp, tmp_path.open("w", newline="") as out:
            required = {"#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"}
            header = inp.readline().rstrip("\r\n").split("\t")
            if not required.issubset(header):
                fail(f"{source} header missing: {sorted(required - set(header))}")
            ix = {name: header.index(name) for name in required}
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            writer.writerow(candidate_header)
            for line in inp:
                fields = line.rstrip("\r\n").split("\t")
                qc["raw_variant_count"] += 1
                rs_field = fields[ix["rsids"]].strip()
                if RS_RE.fullmatch(rs_field):
                    tokens = (rs_field,)
                else:
                    tokens = tuple(x for x in RS_SPLIT_RE.split(rs_field) if RS_RE.fullmatch(x))
                if tokens:
                    qc["canonical_rsid_count"] += 1
                matches = sorted({x for x in tokens if x in hm3})
                if not matches:
                    continue
                qc["HM3_identifier_matches"] += 1
                chrom, pos = fields[ix["#chrom"]], fields[ix["pos"]]
                ref, alt = fields[ix["ref"]].upper(), fields[ix["alt"]].upper()
                if len(matches) != 1:
                    ambiguous.append(["multiple_hm3_rsids_in_one_r13_row", slug, chrom, pos, ref, alt, ",".join(matches)])
                    continue
                snp = matches[0]
                if not valid_chr(chrom) or not pos.isdigit() or ref not in VALID_ALLELES or alt not in VALID_ALLELES:
                    ambiguous.append(["invalid_variant_key", slug, chrom, pos, ref, alt, snp])
                    continue
                if {ref, alt} != set(hm3[snp]):
                    ambiguous.append(["hm3_allele_incompatible", slug, chrom, pos, ref, alt, snp])
                    continue
                qc["HM3_exact_allele_matches"] += 1
                beta, se, pval, frq = fields[ix["beta"]], fields[ix["sebeta"]], fields[ix["pval"]], fields[ix["af_alt"]]
                b, s, pv, fq = assess_values(qc, snp, alt, ref, beta, se, pval, frq, invalid_rows)
                key = (snp, chrom, int(pos), ref, alt)
                observations[key] = observations.get(key, 0) | (1 << source_index)
                writer.writerow([snp, chrom, pos, ref, alt, beta, se, pval, frq, n])
    return observations, ambiguous, qcs


def resolve_map(hm3, observations, ambiguous):
    by_snp = defaultdict(set)
    by_variant = defaultdict(set)
    for snp, chrom, pos, ref, alt in observations:
        variant = (chrom, pos, ref, alt)
        by_snp[snp].add(variant)
        by_variant[variant].add(snp)
    accepted = {}
    conflicts = list(ambiguous)
    for key, mask in observations.items():
        snp, chrom, pos, ref, alt = key
        variant = (chrom, pos, ref, alt)
        if len(by_snp[snp]) == 1 and len(by_variant[variant]) == 1:
            accepted[variant] = (snp, bin(mask).count("1"))
        else:
            reason = "snp_maps_multiple_variants" if len(by_snp[snp]) > 1 else "variant_maps_multiple_snps"
            conflicts.append([reason, "union", chrom, pos, ref, alt, snp])
    with gzip.open(MAP_PATH, "wt", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "chrom", "pos", "ref", "alt", "n_r13_sources"])
        for variant, (snp, count) in sorted(accepted.items(), key=lambda x: (int(x[0][0]), x[0][1], x[0][2], x[0][3])):
            writer.writerow([snp, *variant, count])
    conflict_path = QC / "finngen_grch38_hm3_mapping_conflicts.tsv"
    with conflict_path.open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["conflict_type", "source", "chrom", "pos", "ref", "alt", "SNP_or_SNPs"])
        writer.writerows(conflicts)
    counts = defaultdict(int)
    for _, count in accepted.values():
        counts[count] += 1
    mapped_snps = {snp for snp, _ in accepted.values()}
    true_conflicts = sum(1 for row in conflicts if row[0] != "hm3_allele_incompatible")
    allele_incompatible = sum(1 for row in conflicts if row[0] == "hm3_allele_incompatible")
    with (QC / "finngen_grch38_hm3_map_qc.tsv").open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["hapmap3_reference_snps", len(hm3)])
        writer.writerow(["mapped_unambiguously", len(mapped_snps)])
        writer.writerow(["unmapped", len(hm3) - len(mapped_snps)])
        writer.writerow(["percent_mapped", f"{100 * len(mapped_snps) / len(hm3):.6f}"])
        writer.writerow(["mappings_found_in_all_5_r13_files", counts[5]])
        writer.writerow(["mappings_found_in_fewer_than_5", sum(v for k, v in counts.items() if k < 5)])
        writer.writerow(["conflicting_mapping_rows", true_conflicts])
        writer.writerow(["rejected_hm3_allele_incompatible_rows_across_sources", allele_incompatible])
        writer.writerow(["unique_SNPs", len(mapped_snps)])
        writer.writerow(["unique_variants", len(accepted)])
        writer.writerow(["allele_compatible_with_hm3", len(accepted)])
    return accepted, true_conflicts


def finalize_r13(hm3, accepted, conflict_count, qcs):
    for slug, trait, source, trait_type, n, cases, controls in R13:
        qc = qcs[slug]
        seen = set()
        with (TMP / f"{slug}_r13_hm3_candidates.tsv").open(newline="") as inp, (RAW_HM3 / f"{slug}_hm3.tsv").open("w", newline="") as out:
            reader = csv.DictReader(inp, delimiter="\t")
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            writer.writerow(["SNP", "A1", "A2", "BETA", "SE", "P", "FRQ", "N"])
            for row in reader:
                key = (row["chrom"], int(row["pos"]), row["ref"], row["alt"])
                resolved = accepted.get(key)
                if resolved is None or resolved[0] != row["SNP"]:
                    qc["mapping_conflicts"] += 1
                    continue
                snp = row["SNP"]
                if snp in seen:
                    qc["duplicate_SNPs"] += 1
                    continue
                seen.add(snp)
                b, s, pv, fq = map(fnum, [row["BETA"], row["SE"], row["P"], row["FRQ"]])
                if b is None or s is None or s <= 0 or pv is None or pv <= 0 or pv > 1 or fq is None or fq < 0 or fq > 1:
                    continue
                writer.writerow([snp, row["alt"], row["ref"], row["BETA"], row["SE"], row["P"], row["FRQ"], row["N"]])
                qc["final_premunge_count"] += 1
                if palindromic(row["alt"], row["ref"]):
                    qc["palindromic_SNPs"] += 1
        qc["percent_HM3_coverage"] = 100 * qc["final_premunge_count"] / len(hm3)
        qc["status"] = "ready" if qc["duplicate_SNPs"] == 0 and qc["invalid_P"] == 0 and qc["mapping_conflicts"] == 0 else "blocked"
        qc["notes"] = "ALT is beta-aligned A1; REF is A2; constant phenotype-level N used. Field QC counts are among HM3 identifier/allele candidates."


def prepare_r12(hm3, accepted, invalid_rows):
    qcs = []
    for slug, trait, source, n in R12:
        qc = base_qc(trait, source, "R12", "continuous", n)
        seen = set()
        with gzip.open(ROOT / source, "rt", newline="") as inp, (RAW_HM3 / f"{slug}_hm3.tsv").open("w", newline="") as out:
            required = {"#chrom", "pos", "ref", "alt", "pval", "beta", "sebeta", "af_alt"}
            header = inp.readline().rstrip("\r\n").split("\t")
            if not required.issubset(header):
                fail(f"{source} header missing: {sorted(required - set(header))}")
            ix = {name: header.index(name) for name in required}
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            writer.writerow(["SNP", "A1", "A2", "BETA", "SE", "P", "FRQ", "N"])
            for line in inp:
                fields = line.rstrip("\r\n").split("\t")
                qc["raw_variant_count"] += 1
                chrom, pos = fields[ix["#chrom"]], fields[ix["pos"]]
                ref, alt = fields[ix["ref"]].upper(), fields[ix["alt"]].upper()
                if not valid_chr(chrom) or not pos.isdigit():
                    qc["unmapped_rows"] += 1
                    continue
                key = (chrom, int(pos), ref, alt)
                resolved = accepted.get(key)
                if resolved is None:
                    reversed_key = (chrom, int(pos), alt, ref)
                    if reversed_key in accepted:
                        qc["reversed_REF_ALT_candidates"] += 1
                    qc["unmapped_rows"] += 1
                    continue
                qc["coordinate_matches"] += 1
                qc["exact_REF_ALT_matches"] += 1
                snp = resolved[0]
                qc["canonical_rsid_count"] += 1
                qc["HM3_identifier_matches"] += 1
                qc["HM3_exact_allele_matches"] += 1
                beta, se, pval, frq = fields[ix["beta"]], fields[ix["sebeta"]], fields[ix["pval"]], fields[ix["af_alt"]]
                b, s, pv, fq = assess_values(qc, snp, alt, ref, beta, se, pval, frq, invalid_rows)
                if snp in seen:
                    qc["duplicate_SNPs"] += 1
                    continue
                seen.add(snp)
                if b is None or s is None or s <= 0 or pv is None or pv <= 0 or pv > 1 or fq is None or fq < 0 or fq > 1:
                    continue
                writer.writerow([snp, alt, ref, beta, se, pval, frq, n])
                qc["final_premunge_count"] += 1
                if palindromic(alt, ref):
                    qc["palindromic_SNPs"] += 1
        qc["percent_HM3_coverage"] = 100 * qc["final_premunge_count"] / len(hm3)
        qc["status"] = "ready" if qc["duplicate_SNPs"] == 0 and qc["invalid_P"] == 0 and qc["mapping_conflicts"] == 0 else "blocked"
        qc["notes"] = "Exact GRCh38 chrom/pos/REF/ALT map only; reversed candidates reported but not transformed. ALT is beta-aligned A1; REF is A2."
        qcs.append(qc)
    return qcs


def write_qc(qcs, invalid_rows):
    with (QC / "ldsc_premunging_qc.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=QC_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(qcs)
    with (QC / "ldsc_invalid_pvalues.tsv").open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["trait", "SNP", "beta", "SE", "P", "issue"])
        writer.writerows(invalid_rows)
    with (QC / "ldsc_munging_errors.tsv").open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["trait", "stage", "error", "action"])


def load_finngen_map():
    mapping = {}
    with gzip.open(MAP_PATH, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mapping[(row["chrom"], int(row["pos"]), row["ref"], row["alt"])] = row["SNP"]
    return mapping


def resolve_finngen_p_underflow():
    """Recover FinnGen P=0 from official mlogp for five blocked R12 traits.

    Only normal float64 p-values are accepted. Values beyond the calculated
    normal-float boundary are excluded; no flooring, capping, subnormal value,
    or arbitrary replacement is used. Ordinary P>0 rows retain source pval.
    """
    mapping = load_finngen_map()
    max_normal_mlogp = -math.log10(np.finfo(np.float64).tiny)
    audit_rows = []
    excluded_rows = []
    summaries = []
    validation_rows = []
    trait_updates = {}
    rng = random.Random(20260823)

    blocked_slugs = {"alt", "ggt", "hdl", "ldl", "triglycerides_fasting"}
    for slug, trait, source, n in R12:
        if slug not in blocked_slugs:
            continue
        ordinary_sample = []
        ordinary_seen = 0
        smallest_heap = []  # top 100 mlogp values among nonzero P
        boundary_heap = []  # 100 mlogp values closest to normal-float boundary
        pzero = []
        final_count = 0
        out_path = RAW_HM3 / f"{slug}_hm3.tsv"
        with gzip.open(ROOT / source, "rt", newline="") as inp, out_path.open("w", newline="") as out:
            header = inp.readline().rstrip("\r\n").split("\t")
            required = {"#chrom", "pos", "ref", "alt", "pval", "mlogp", "beta", "sebeta", "af_alt"}
            if not required.issubset(header):
                fail(f"{source} header missing for underflow recovery: {sorted(required - set(header))}")
            ix = {name: header.index(name) for name in required}
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            writer.writerow(["SNP", "A1", "A2", "BETA", "SE", "P", "FRQ", "N"])
            for line in inp:
                fields = line.rstrip("\r\n").split("\t")
                chrom, pos_text = fields[ix["#chrom"]], fields[ix["pos"]]
                if not valid_chr(chrom) or not pos_text.isdigit():
                    continue
                ref, alt = fields[ix["ref"]].upper(), fields[ix["alt"]].upper()
                snp = mapping.get((chrom, int(pos_text), ref, alt))
                if snp is None:
                    continue
                beta_text, se_text = fields[ix["beta"]], fields[ix["sebeta"]]
                p_text, mlogp_text, frq_text = fields[ix["pval"]], fields[ix["mlogp"]], fields[ix["af_alt"]]
                beta, se, pval, mlogp, frq = map(fnum, [beta_text, se_text, p_text, mlogp_text, frq_text])
                if beta is None or se is None or se <= 0 or pval is None or pval < 0 or pval > 1 or frq is None or frq < 0 or frq > 1:
                    continue
                if pval > 0:
                    writer.writerow([snp, alt, ref, beta_text, se_text, p_text, frq_text, n])
                    final_count += 1
                    if mlogp is not None and mlogp >= 0:
                        record = (snp, pval, mlogp)
                        ordinary_seen += 1
                        if len(ordinary_sample) < 100:
                            ordinary_sample.append(record)
                        else:
                            j = rng.randrange(ordinary_seen)
                            if j < 100:
                                ordinary_sample[j] = record
                        heap_record = (mlogp, snp, pval)
                        if len(smallest_heap) < 100:
                            heapq.heappush(smallest_heap, heap_record)
                        elif mlogp > smallest_heap[0][0]:
                            heapq.heapreplace(smallest_heap, heap_record)
                        distance = abs(mlogp - max_normal_mlogp)
                        boundary_record = (-distance, mlogp, snp, pval)
                        if len(boundary_heap) < 100:
                            heapq.heappush(boundary_heap, boundary_record)
                        elif distance < -boundary_heap[0][0]:
                            heapq.heapreplace(boundary_heap, boundary_record)
                    continue

                zwald = beta / se
                status = "INVALID_MLOGP"
                reconstructed = None
                notes = []
                if mlogp is not None and math.isfinite(mlogp) and mlogp > 0:
                    if mlogp <= max_normal_mlogp:
                        candidate = float(np.float64(10.0) ** np.float64(-mlogp))
                        if math.isfinite(candidate) and 0 < candidate <= 1 and candidate >= np.finfo(np.float64).tiny:
                            reconstructed = candidate
                            status = "RECOVERABLE_NORMAL"
                            if abs(-math.log10(candidate) - mlogp) > 1e-12:
                                notes.append("roundtrip_difference_gt_1e-12")
                        else:
                            status = "INVALID_MLOGP"
                    else:
                        status = "BEYOND_NORMAL_FLOAT_RANGE"
                if not math.isfinite(zwald) or se <= 0:
                    notes.append("invalid_wald_diagnostic")
                pzero_row = {
                    "trait": trait, "chrom": chrom, "pos": pos_text, "ref": ref, "alt": alt,
                    "SNP": snp, "beta": beta_text, "sebeta": se_text,
                    "beta_over_se": f"{zwald:.17g}", "abs_beta_over_se": f"{abs(zwald):.17g}",
                    "pval_original": p_text, "mlogp": mlogp_text,
                    "P_reconstructed": f"{reconstructed:.17g}" if reconstructed is not None else "",
                    "reconstruction_status": status, "HM3_status": "exact_GRCh38_REF_ALT_HM3_match",
                    "notes": ";".join(notes),
                }
                audit_rows.append(pzero_row)
                pzero.append((mlogp, abs(zwald), status))
                if status == "RECOVERABLE_NORMAL":
                    writer.writerow([snp, alt, ref, beta_text, se_text, f"{reconstructed:.17g}", frq_text, n])
                    final_count += 1
                else:
                    excluded = dict(pzero_row)
                    excluded["exclusion_reason"] = "excluded_numeric_underflow_unrepresentable" if status == "BEYOND_NORMAL_FLOAT_RANGE" else "excluded_invalid_mlogp"
                    excluded_rows.append(excluded)

        for group, sample in (
            ("ordinary_random_100", ordinary_sample),
            ("smallest_nonzero_100", [(snp, p, m) for m, snp, p in sorted(smallest_heap, reverse=True)]),
            ("near_underflow_boundary", [(snp, p, m) for _, m, snp, p in sorted(boundary_heap, reverse=True)]),
        ):
            differences = [abs(m + math.log10(p)) for _, p, m in sample]
            normal_differences = [abs(m + math.log10(p)) for _, p, m in sample if m <= max_normal_mlogp and p >= np.finfo(np.float64).tiny]
            beyond_count = len(differences) - len(normal_differences)
            strict_pass = bool(normal_differences) and all(d <= 0.0001 for d in normal_differences)
            if group == "ordinary_random_100":
                group_status = "PASS" if strict_pass and beyond_count == 0 else "FAIL"
            elif normal_differences:
                group_status = "PASS" if strict_pass else "FAIL"
            else:
                group_status = "PASS_SUBNORMAL_ROUNDING_DIAGNOSTIC"
            validation_rows.append({
                "trait": trait, "sample_group": group, "n": len(differences),
                "median_absolute_difference": f"{statistics.median(differences):.12g}" if differences else "NA",
                "maximum_absolute_difference": f"{max(differences):.12g}" if differences else "NA",
                "normal_float_domain_n": len(normal_differences),
                "beyond_normal_float_domain_n": beyond_count,
                "tolerance_normal_domain": "0.0001",
                "normal_domain_agreeing_within_tolerance": sum(d <= 0.0001 for d in normal_differences),
                "status": group_status,
            })
        mvalues = [x[0] for x in pzero if x[0] is not None]
        zvalues = [x[1] for x in pzero]
        status_counts = defaultdict(int)
        for _, _, status in pzero:
            status_counts[status] += 1
        summaries.append({
            "trait": trait, "original_P0": len(pzero), "valid_mlogp": len(mvalues),
            "recoverable_normal": status_counts["RECOVERABLE_NORMAL"],
            "beyond_normal_float_range": status_counts["BEYOND_NORMAL_FLOAT_RANGE"],
            "invalid_mlogp": status_counts["INVALID_MLOGP"],
            "recovered_and_retained": status_counts["RECOVERABLE_NORMAL"],
            "excluded_unrepresentable": status_counts["BEYOND_NORMAL_FLOAT_RANGE"] + status_counts["INVALID_MLOGP"],
            "min_mlogp": min(mvalues), "median_mlogp": statistics.median(mvalues), "max_mlogp": max(mvalues),
            "min_abs_Zwald": min(zvalues), "median_abs_Zwald": statistics.median(zvalues), "max_abs_Zwald": max(zvalues),
        })
        trait_updates[trait] = (final_count, len(pzero), status_counts)

    if any(row["status"] == "FAIL" for row in validation_rows):
        fail("FinnGen mlogp validation failed; prepared files must not be used")

    audit_fields = ["trait", "chrom", "pos", "ref", "alt", "SNP", "beta", "sebeta",
                    "beta_over_se", "abs_beta_over_se", "pval_original", "mlogp",
                    "P_reconstructed", "reconstruction_status", "HM3_status", "notes"]
    with (QC / "finngen_pzero_underflow_audit.tsv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=audit_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(audit_rows)
    with (QC / "finngen_mlogp_validation.tsv").open("w", newline="") as out:
        fields = list(validation_rows[0])
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(validation_rows)
    with (QC / "finngen_pzero_recovery_summary.tsv").open("w", newline="") as out:
        fields = list(summaries[0])
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(summaries)
    with (QC / "finngen_pzero_excluded_variants.tsv").open("w", newline="") as out:
        fields = audit_fields + ["exclusion_reason"]
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(excluded_rows)
    (QC / "finngen_pzero_float64_limit.txt").write_text(
        f"numpy_version\t{np.__version__}\n"
        f"float64_tiny_normal\t{np.finfo(np.float64).tiny:.17g}\n"
        f"MAX_NORMAL_MLOGP\t{max_normal_mlogp:.17g}\n"
        "subnormal_values_used\tno\npvalue_floor_or_cap_used\tno\n"
    )

    prem_path = QC / "ldsc_premunging_qc.tsv"
    with prem_path.open() as handle:
        prem_rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in prem_rows:
        if row["trait"] in trait_updates:
            final_count, original_pzero, counts = trait_updates[row["trait"]]
            row["final_premunge_count"] = str(final_count)
            row["percent_HM3_coverage"] = str(100 * final_count / 1187349)
            row["invalid_P"] = "0"
            row["status"] = "ready"
            row["notes"] += (f" Underflow recovery reviewed: {counts['RECOVERABLE_NORMAL']} P=0 rows reconstructed from official mlogp; "
                             f"{counts['BEYOND_NORMAL_FLOAT_RANGE'] + counts['INVALID_MLOGP']} unrepresentable/invalid rows explicitly excluded; no flooring.")
    with prem_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=QC_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(prem_rows)
    print(f"MAX_NORMAL_MLOGP={max_normal_mlogp:.17g}")
    print(f"audited_P0={len(audit_rows)}")
    print(f"recoverable={sum(x['recoverable_normal'] for x in summaries)}")
    print(f"excluded={len(excluded_rows)}")


def write_metadata():
    (QC / "asthma_ldsc_column_mapping.txt").write_text(
        "Asthma LDSC column mapping\n==========================\n\n"
        "Source: data/raw/asthma_DISCOVERY.tsv (GCST90302886)\n"
        "SNP: VARIANT_id (canonical rsIDs retained; coordinate-style IDs excluded)\n"
        "Effect allele / A1: effect_allele\nOther allele / A2: other_allele\n"
        "Beta: beta\nSE: standard_error\nP: p_value\nEffect-allele frequency: effect_allele_frequency\nN: N\n\n"
        "Orientation confirmation: the source uses the GWAS Catalog harmonised summary-statistics schema, in which beta and effect_allele_frequency are aligned to effect_allele. The current upstream script scripts/02_harmonise_inputs.R independently maps effect_allele to EA and beta to BETA without inversion (lines 218-230). Thus A1=effect_allele is the allele for which BETA was estimated.\n"
    )
    (QC / "ldsc_asthma_exposure_provenance.txt").write_text(
        "LDSC asthma exposure provenance\n===============================\n\n"
        "Dataset: GCST90302886 / data/raw/asthma_DISCOVERY.tsv\n"
        "Population: British / European\nCases: 48,623\nControls: 290,722\nTotal N: 339,345\n"
        "Data type: array/imputed genome-wide GWAS\nDirect reconstructed no-MHC HapMap3 coverage in prior audit: approximately 99.7%\n\n"
        "The project-wide 12-input GWAMA is not used for LDSC because it has extensive UK Biobank participant overlap, one multi-ancestry input, one ancestry-unresolved input, five sparse WES datasets, and no defensible combined sample size. LDSC therefore represents this specific high-quality British/European asthma GWAS.\n"
    )
    (QC / "ldsc_reference_paths.txt").write_text(
        "European LD scores (chromosomes 1-22): data/ldsc/reference/1000G_Phase3_ldscores/LDscore.[1-22]\n"
        "European HM3 regression weights (chromosomes 1-22): data/ldsc/reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.[1-22]\n"
        "Reconstructed no-MHC HapMap3 allele list: data/ldsc/reference/w_hm3.snplist\n"
        "Header-bearing LDSC --merge-alleles view: data/ldsc/prepared/w_hm3_merge_alleles.tsv\n"
        "HapMap3 allele-list rows: 1,187,349; SNP duplicates: 0; SNP/allele duplicates: 0; invalid alleles: 0.\n"
    )


def main():
    for path in (RAW_HM3, PREP / "munged", TMP, ROOT / "results/ldsc/logs/munging", QC):
        path.mkdir(parents=True, exist_ok=True)
    hm3 = load_hm3()
    invalid_rows = []
    write_metadata()
    qcs = [prepare_asthma(hm3, invalid_rows)]
    observations, ambiguous, r13_qcs = scan_r13(hm3, invalid_rows)
    accepted, conflict_count = resolve_map(hm3, observations, ambiguous)
    finalize_r13(hm3, accepted, conflict_count, r13_qcs)
    qcs.extend(r13_qcs[slug] for slug, *_ in R13)
    qcs.extend(prepare_r12(hm3, accepted, invalid_rows))
    write_qc(qcs, invalid_rows)
    print(f"HapMap3 SNPs: {len(hm3)}")
    print(f"Accepted FinnGen GRCh38 mappings: {len(accepted)}")
    print(f"Mapping conflict rows: {conflict_count}")
    print(f"Invalid p-value rows: {len(invalid_rows)}")
    for qc in qcs:
        print(f"{qc['trait']}\t{qc['raw_variant_count']}\t{qc['final_premunge_count']}\t{qc['status']}")


if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument("--recover-pzero-only", action="store_true",
                            help="Recreate only five blocked R12 HM3 files using reviewed FinnGen mlogp underflow recovery")
        args = parser.parse_args()
        if args.recover_pzero_only:
            resolve_finngen_p_underflow()
        else:
            main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
