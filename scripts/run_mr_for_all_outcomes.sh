#!/usr/bin/env bash
set -euo pipefail

# Run the full asthma -> outcome MR workflow for every configured outcome
# with a completed OpenGWAS ID.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "${script_dir}/.." && pwd)"
config_file="${project_dir}/config/mr_outcomes.tsv"

if [[ ! -f "${config_file}" ]]; then
  echo "MR outcome config file not found: ${config_file}" >&2
  exit 1
fi

awk -F '\t' 'NR > 1 && $3 != "TO_FILL" && $3 != "" { print $1 }' "${config_file}" |
while IFS= read -r outcome_label; do
  if [[ -z "${outcome_label}" ]]; then
    continue
  fi

  echo ""
  echo "=== Running MR workflow for ${outcome_label} ==="
  bash "${script_dir}/run_mr_for_outcome.sh" "${outcome_label}"
done
