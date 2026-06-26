#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(dirname "$script_dir")"

vep="$project_dir/tools/ensembl-vep/vep"
cache_dir="$HOME/.vep"
annotation_dir="$project_dir/results/annotation"

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
  --force_overwrite
  --no_stats
)

"$vep" "${common_args[@]}" \
  --input_file "$annotation_dir/top20_independent_lead_snps_vep_input.txt" \
  --output_file "$annotation_dir/top20_independent_lead_snps_vep_output.txt"

"$vep" "${common_args[@]}" \
  --input_file "$annotation_dir/independent_lead_snps_vep_input.txt" \
  --output_file "$annotation_dir/independent_lead_snps_vep_output.txt"

echo "Written: $annotation_dir/top20_independent_lead_snps_vep_output.txt"
echo "Written: $annotation_dir/independent_lead_snps_vep_output.txt"
