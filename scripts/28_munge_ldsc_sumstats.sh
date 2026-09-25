#!/usr/bin/env bash
set -u
set -o pipefail

# Munge only traits that passed scripts/27_prepare_ldsc_sumstats.py QC.
# Run from any location; all analysis paths are resolved relative to the repo.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONDA_SETUP="${HOME}/miniforge3/etc/profile.d/conda.sh"
LDSC_MUNGE="${PROJECT_ROOT}/tools/ldsc/munge_sumstats.py"
HM3_ALLELES="${PROJECT_ROOT}/data/ldsc/prepared/w_hm3_merge_alleles.tsv"
PREMUNGE_QC="${PROJECT_ROOT}/results/ldsc/qc/ldsc_premunging_qc.tsv"
INPUT_DIR="${PROJECT_ROOT}/data/ldsc/prepared/raw_hm3"
OUTPUT_DIR="${PROJECT_ROOT}/data/ldsc/prepared/munged"
LOG_DIR="${PROJECT_ROOT}/results/ldsc/logs/munging"
ERROR_FILE="${PROJECT_ROOT}/results/ldsc/qc/ldsc_munging_errors.tsv"

if [[ ! -f "${CONDA_SETUP}" ]]; then
  echo "Missing Conda setup: ${CONDA_SETUP}" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "${CONDA_SETUP}"
conda activate ldsc39

python --version
python "${LDSC_MUNGE}" -h >/dev/null

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"
printf 'trait\tstage\terror\taction\n' > "${ERROR_FILE}"

slugs=(asthma allergic_rhinitis bronchiectasis eosinophilic_disease bmi_irn nafld crp alt ast ggt hdl ldl triglycerides_fasting)
traits=("Asthma" "Allergic rhinitis" "Bronchiectasis" "Eosinophilic disease" "BMI IRN" "NAFLD" "CRP" "ALT" "AST" "GGT" "HDL cholesterol" "LDL cholesterol" "Fasting triglycerides")
blocked=0
failed=0

for index in "${!slugs[@]}"; do
  slug="${slugs[$index]}"
  trait="${traits[$index]}"
  if [[ $# -gt 0 ]]; then
    selected=0
    for requested_slug in "$@"; do
      if [[ "${slug}" == "${requested_slug}" ]]; then
        selected=1
      fi
    done
    if [[ ${selected} -eq 0 ]]; then
      continue
    fi
  fi
  input="${INPUT_DIR}/${slug}_hm3.tsv"
  output_prefix="${OUTPUT_DIR}/${slug}"
  console_log="${LOG_DIR}/${slug}.console.log"
  ldsc_log="${LOG_DIR}/${slug}.log"
  status="$(awk -F '\t' -v wanted="${trait}" 'NR==1 {for(i=1;i<=NF;i++){if($i=="trait")r=i;if($i=="status")t=i}} NR>1 && $r==wanted {print $t}' "${PREMUNGE_QC}")"
  if [[ "${status}" != "ready" ]]; then
    printf '%s\tpremunging_qc\tstatus=%s\tskipped; review invalid or ambiguous rows\n' "${slug}" "${status:-missing}" >> "${ERROR_FILE}"
    blocked=$((blocked + 1))
    continue
  fi
  if [[ ! -f "${input}" ]]; then
    printf '%s\tinput\tmissing prepared file\tskipped\n' "${slug}" >> "${ERROR_FILE}"
    failed=$((failed + 1))
    continue
  fi

  python "${LDSC_MUNGE}" \
    --sumstats "${input}" \
    --out "${output_prefix}" \
    --snp SNP --a1 A1 --a2 A2 --p P --frq FRQ --N-col N \
    --signed-sumstats BETA,0 \
    --merge-alleles "${HM3_ALLELES}" > "${console_log}" 2>&1
  code=$?
  if [[ -f "${output_prefix}.log" ]]; then
    mv "${output_prefix}.log" "${ldsc_log}"
  fi
  if [[ ${code} -ne 0 ]]; then
    printf '%s\tmunge_sumstats\texit_code=%s\tinspect %s\n' "${slug}" "${code}" "${console_log#${PROJECT_ROOT}/}" >> "${ERROR_FILE}"
    failed=$((failed + 1))
    continue
  fi

  # The CBIIT ldsc39 merge implementation writes the full merge list and
  # represents filtered rows as blank A1/A2/Z/N fields. Retain only complete
  # rows so the final .sumstats.gz is a standard, fully usable LDSC table.
  filtered_tmp="${PROJECT_ROOT}/data/ldsc/prepared/tmp/${slug}.complete.sumstats.gz"
  gzip -dc "${output_prefix}.sumstats.gz" | \
    awk -F '\t' 'NR==1 || (NF==5 && $1!="" && $2!="" && $3!="" && $4!="" && $5!="")' | \
    gzip -c > "${filtered_tmp}"
  filter_code=$?
  if [[ ${filter_code} -ne 0 ]]; then
    printf '%s\tpost_munge_filter\texit_code=%s\tinspect raw munged output\n' "${slug}" "${filter_code}" >> "${ERROR_FILE}"
    failed=$((failed + 1))
    continue
  fi
  mv "${filtered_tmp}" "${output_prefix}.sumstats.gz"
done

echo "blocked_traits=${blocked}"
echo "failed_traits=${failed}"
if [[ ${blocked} -ne 0 || ${failed} -ne 0 ]]; then
  exit 2
fi
