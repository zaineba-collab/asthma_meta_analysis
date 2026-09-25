#!/usr/bin/env bash
set -euo pipefail

# Corrected GRCh38 VEP workflow. Scripts 17--19 are retained only as provenance
# and are superseded for final annotation reporting by scripts 38--40.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(dirname "$script_dir")"
vep="$project_dir/tools/ensembl-vep/vep"
cache_dir="${VEP_CACHE_DIR:-$HOME/.vep}"
base="$project_dir/results/annotation_grch38"
probe_input="$base/input/lead_snps_grch38_allele_orientation_probe.tsv"
probe_output="$base/raw/lead_snps_grch38_allele_orientation_probe_vep.txt"
probe_log="$base/raw/lead_snps_grch38_allele_orientation_probe_vep.log"
final_input="$base/input/lead_snps_grch38_vep_input.tsv"
final_output="$base/raw/lead_snps_grch38_vep_output.txt"
final_log="$base/raw/lead_snps_grch38_vep.log"
metadata="$base/qc/vep_run_metadata.txt"
fasta="$project_dir/data/reference/ensembl_116_GRCh38_primary_assembly.fa"

[[ -x "$vep" ]] || { echo "ERROR: VEP executable unavailable: $vep" >&2; exit 1; }
[[ -d "$cache_dir/homo_sapiens/116_GRCh38" ]] || { echo "ERROR: VEP 116 GRCh38 cache unavailable" >&2; exit 1; }
[[ -s "$fasta" ]] || { echo "ERROR: matching Ensembl-116 GRCh38 FASTA unavailable: $fasta" >&2; exit 1; }
mkdir -p "$base"/{input,raw,summary,qc,figures}

rm -f "$probe_output"
Rscript "$script_dir/38_prepare_corrected_grch38_vep_input.R"

common=(--offline --cache --dir_cache "$cache_dir" --species homo_sapiens --assembly GRCh38 --cache_version 116
        --tab --symbol --gene_phenotype --variant_class --check_existing --fork 4 --force_overwrite --no_stats)

if [[ -s "$probe_input" ]]; then
  "$vep" "${common[@]}" --fasta "$fasta" --lookup_ref --input_file "$probe_input" --output_file "$probe_output" >"$probe_log" 2>&1
fi
Rscript "$script_dir/38_prepare_corrected_grch38_vep_input.R"
[[ -s "$final_input" ]] || { echo "ERROR: no final VEP-ready records" >&2; exit 1; }
"$vep" "${common[@]}" --input_file "$final_input" --output_file "$final_output" >"$final_log" 2>&1

{
  "$vep" --help 2>&1 | head -n 1 || true
  echo "cache_version=116"
  echo "assembly=GRCh38 (cache reports GRCh38.p14)"
  echo "reference_fasta=data/reference/ensembl_116_GRCh38_primary_assembly.fa (downloaded from official Ensembl release-116 archive; compressed SHA256 d8c3af0094a7bba6125763bad779ec18a81483c739c6ed122094bdf86c187b92)"
  echo "options=--offline --cache --species homo_sapiens --assembly GRCh38 --cache_version 116 --tab --symbol --gene_phenotype --variant_class --check_existing --fork 4 --no_stats"
  printf 'command=%q ' "$vep" "${common[@]}" --input_file "$final_input" --output_file "$final_output"
  echo
} > "$metadata"

Rscript "$script_dir/40_summarise_corrected_vep_annotation.R"
echo "Corrected VEP output: $final_output"
