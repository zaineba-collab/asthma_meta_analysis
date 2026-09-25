#!/usr/bin/env bash
set -u

# Run unconstrained-intercept, observed-scale univariate LDSC for 13 traits.
# No liability conversion and no genetic-correlation analysis are performed.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONDA_SETUP="${HOME}/miniforge3/etc/profile.d/conda.sh"
LDSC="${PROJECT_ROOT}/tools/ldsc/ldsc.py"
MUNGED_DIR="${PROJECT_ROOT}/data/ldsc/prepared/munged"
REF_PREFIX="${PROJECT_ROOT}/data/ldsc/reference/1000G_Phase3_ldscores/LDscore."
WEIGHT_PREFIX="${PROJECT_ROOT}/data/ldsc/reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
OUTPUT_DIR="${PROJECT_ROOT}/results/ldsc/h2"
RUN_STATUS="${OUTPUT_DIR}/run_status.tsv"

if [[ ! -f "${CONDA_SETUP}" ]]; then
  echo "Missing Conda setup: ${CONDA_SETUP}" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "${CONDA_SETUP}"
conda activate ldsc39

python --version
python "${LDSC}" -h >/dev/null

for chromosome in $(seq 1 22); do
  if [[ ! -r "${REF_PREFIX}${chromosome}.l2.ldscore.gz" ]]; then
    echo "Missing LD-score chromosome ${chromosome}" >&2
    exit 1
  fi
  if [[ ! -r "${WEIGHT_PREFIX}${chromosome}.l2.ldscore.gz" ]]; then
    echo "Missing regression-weight chromosome ${chromosome}" >&2
    exit 1
  fi
done

mkdir -p "${OUTPUT_DIR}"
printf 'trait\texit_status\tconsole_log\tldsc_log\n' > "${RUN_STATUS}"

slugs=(asthma allergic_rhinitis bronchiectasis eosinophilic_disease bmi_irn nafld crp alt ast ggt hdl ldl triglycerides_fasting)
failures=0

for slug in "${slugs[@]}"; do
  input="${MUNGED_DIR}/${slug}.sumstats.gz"
  output_prefix="${OUTPUT_DIR}/${slug}"
  console_log="${output_prefix}.console.log"
  if [[ ! -r "${input}" ]]; then
    printf '%s\t127\t%s\t%s\n' "${slug}" "${console_log#${PROJECT_ROOT}/}" "${output_prefix#${PROJECT_ROOT}/}.log" >> "${RUN_STATUS}"
    failures=$((failures + 1))
    continue
  fi

  python "${LDSC}" \
    --h2 "${input}" \
    --ref-ld-chr "${REF_PREFIX}" \
    --w-ld-chr "${WEIGHT_PREFIX}" \
    --out "${output_prefix}" > "${console_log}" 2>&1
  code=$?
  printf '%s\t%s\t%s\t%s\n' "${slug}" "${code}" "${console_log#${PROJECT_ROOT}/}" "${output_prefix#${PROJECT_ROOT}/}.log" >> "${RUN_STATUS}"
  if [[ ${code} -ne 0 ]]; then
    failures=$((failures + 1))
  fi
done

echo "attempted_traits=${#slugs[@]}"
echo "failed_traits=${failures}"
exit "${failures}"
