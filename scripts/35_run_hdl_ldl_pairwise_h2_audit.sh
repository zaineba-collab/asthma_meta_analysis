#!/usr/bin/env bash
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONDA_SETUP="${HOME}/miniforge3/etc/profile.d/conda.sh"
LDSC="${PROJECT_ROOT}/tools/ldsc/ldsc.py"
TMP="${PROJECT_ROOT}/data/ldsc/prepared/tmp"
OUT="${PROJECT_ROOT}/results/ldsc/qc/pairwise_h2_audit"
REF="${PROJECT_ROOT}/data/ldsc/reference/1000G_Phase3_ldscores/LDscore."
WEIGHTS="${PROJECT_ROOT}/data/ldsc/reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."

# shellcheck disable=SC1090
source "${CONDA_SETUP}"
conda activate ldsc39
mkdir -p "${OUT}"

python --version > "${OUT}/environment.txt"
git -C "${PROJECT_ROOT}/tools/ldsc" branch --show-current >> "${OUT}/environment.txt"
git -C "${PROJECT_ROOT}/tools/ldsc" rev-parse HEAD >> "${OUT}/environment.txt"

failures=0
printf 'trait\texit_status\n' > "${OUT}/run_status.tsv"
for trait in hdl ldl; do
  python "${LDSC}" \
    --h2 "${TMP}/${trait}_asthma_overlap.sumstats.gz" \
    --ref-ld-chr "${REF}" \
    --w-ld-chr "${WEIGHTS}" \
    --out "${OUT}/${trait}_asthma_overlap" > "${OUT}/${trait}_asthma_overlap.console.log" 2>&1
  code=$?
  printf '%s\t%s\n' "${trait}" "${code}" >> "${OUT}/run_status.tsv"
  if [[ ${code} -ne 0 ]]; then failures=$((failures + 1)); fi
done
exit "${failures}"
