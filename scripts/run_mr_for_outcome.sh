#!/usr/bin/env bash
set -euo pipefail

# Run the full supervisor-aligned asthma -> outcome MR workflow.
#
# Supervisor plan:
#   Use already PLINK LD-clumped asthma signals as instruments and do not
#   perform additional clumping in TwoSampleMR.
#
# Usage:
#   bash scripts/run_mr_for_outcome.sh hdl

if [[ $# -lt 1 || -z "${1:-}" ]]; then
  echo "Usage: bash scripts/run_mr_for_outcome.sh <outcome_label>" >&2
  exit 1
fi

outcome_label="$1"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "${script_dir}/.." && pwd)"
config_file="${project_dir}/config/mr_outcomes.tsv"

if [[ ! -f "${config_file}" ]]; then
  echo "MR outcome config file not found: ${config_file}" >&2
  exit 1
fi

config_row="$(awk -F '\t' -v label="${outcome_label}" 'NR > 1 && $1 == label { print; found = 1 } END { if (!found) exit 1 }' "${config_file}")" || {
  echo "Outcome label '${outcome_label}' not found in ${config_file}" >&2
  exit 1
}

outcome_name="$(printf '%s\n' "${config_row}" | awk -F '\t' '{ print $2 }')"
outcome_id="$(printf '%s\n' "${config_row}" | awk -F '\t' '{ print $3 }')"
outcome_type="$(printf '%s\n' "${config_row}" | awk -F '\t' '{ print $4 }')"
output_dir="${project_dir}/results/mr/asthma_to_${outcome_label}"

if [[ -z "${outcome_id}" || "${outcome_id}" == "TO_FILL" ]]; then
  echo "Outcome ID for '${outcome_label}' is TO_FILL in ${config_file}; skipping." >&2
  exit 1
fi

echo "Selected outcome label: ${outcome_label}"
echo "Outcome name: ${outcome_name}"
echo "Outcome ID: ${outcome_id}"
echo "Outcome type: ${outcome_type}"
echo "Output directory: ${output_dir}"
echo "Using already LD-clumped asthma signals. No additional TwoSampleMR clumping performed."

cd "${project_dir}"
export MR_OUTCOME_LABEL="${outcome_label}"

Rscript scripts/21_prepare_mr_exposure.R
Rscript scripts/22_prepare_mr_outcomes.R
Rscript scripts/23_harmonise_mr_data.R
Rscript scripts/24_run_mr_analysis.R
Rscript scripts/25_mr_sensitivity_analysis.R
