#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 18_run_vep_annotation.sh
#
# Purpose:
#   Run local Ensembl Variant Effect Predictor (VEP) annotation on
#   the independent lead SNP VEP input files prepared by:
#     scripts/17_prepare_vep_input.R
#
# Inputs:
#   results/annotation/top20_independent_lead_snps_vep_input.txt
#   results/annotation/independent_lead_snps_vep_input.txt
#
# Outputs:
#   results/annotation/top20_independent_lead_snps_vep_output.txt
#   results/annotation/independent_lead_snps_vep_output.txt
#
# Notes:
#   - This uses a local command-line VEP install, not the Ensembl web UI.
#   - VEP runs offline against the local cache in ~/.vep.
#   - --check_existing asks VEP to report known co-located variant IDs
#     such as dbSNP rsIDs in the Existing_variation column.
#   - Run scripts/19_summarise_vep_annotation.R after this to create
#     the cleaned gene and consequence summaries.
# ============================================================

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(dirname "$script_dir")"

vep="$project_dir/tools/ensembl-vep/vep"
cache_dir="$HOME/.vep"
annotation_dir="$project_dir/results/annotation"

top20_input="$annotation_dir/top20_independent_lead_snps_vep_input.txt"
all_input="$annotation_dir/independent_lead_snps_vep_input.txt"
top20_output="$annotation_dir/top20_independent_lead_snps_vep_output.txt"
all_output="$annotation_dir/independent_lead_snps_vep_output.txt"

if [[ ! -x "$vep" ]]; then
  echo "ERROR: VEP executable not found or not executable: $vep" >&2
  echo "Install VEP under tools/ensembl-vep, or update the vep path in this script." >&2
  exit 1
fi

if [[ ! -d "$cache_dir" ]]; then
  echo "ERROR: VEP cache directory not found: $cache_dir" >&2
  echo "Install/download the GRCh38 VEP cache before running offline annotation." >&2
  exit 1
fi

if [[ ! -s "$top20_input" ]]; then
  echo "ERROR: Missing top-20 VEP input: $top20_input" >&2
  echo "Run scripts/17_prepare_vep_input.R first." >&2
  exit 1
fi

if [[ ! -s "$all_input" ]]; then
  echo "ERROR: Missing full lead-SNP VEP input: $all_input" >&2
  echo "Run scripts/17_prepare_vep_input.R first." >&2
  exit 1
fi

common_args=(
  --offline
  --cache
  --dir_cache "$cache_dir"
  --species homo_sapiens
  --assembly GRCh38
  --cache_version 116
  --tab
  --symbol
  --gene_phenotype
  --variant_class
  --check_existing
  --force_overwrite
  --no_stats
)

echo "Running VEP on top-20 independent lead SNPs..."
"$vep" "${common_args[@]}" \
  --input_file "$top20_input" \
  --output_file "$top20_output"

echo "Running VEP on all independent lead SNPs..."
"$vep" "${common_args[@]}" \
  --input_file "$all_input" \
  --output_file "$all_output"

echo "Written: $top20_output"
echo "Written: $all_output"
