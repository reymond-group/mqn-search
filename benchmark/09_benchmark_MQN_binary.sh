#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

MQN_BIN="${1:-${REPO_ROOT}/build/bin/mqn_new}"
INPUT_FILE="${2:-${REPO_ROOT}/datasets/ZINC_ChEMBL34_hi_conf_1M.smi}"
OUTPUT_FILE="${3:-${REPO_ROOT}/benchmark_out/ZINC_ChEMBL34_hi_conf_1M_MQN_benchmark.csv}"
LOG_FILE="${4:-${REPO_ROOT}/benchmark_out/09_benchmark_MQN_binary.log}"
TMPDIR="$(mktemp -d "${TMPDIR:-/tmp}/mqn_benchmark_XXXXXX")"
SMILES_ONLY_INPUT="${TMPDIR}/input_smiles_only.smi"

mkdir -p "$(dirname "${OUTPUT_FILE}")" "$(dirname "${LOG_FILE}")"

cleanup() {
    rm -rf "${TMPDIR}"
}
trap cleanup EXIT

if [ ! -x "${MQN_BIN}" ]; then
    echo "ERROR: MQN binary is missing or not executable: ${MQN_BIN}" >&2
    exit 1
fi

if [ ! -f "${INPUT_FILE}" ]; then
    echo "ERROR: input file not found: ${INPUT_FILE}" >&2
    exit 1
fi

START_EPOCH="$(date +%s)"
{
    echo "Started MQN binary benchmark"
    echo "binary=${MQN_BIN}"
    echo "input=${INPUT_FILE}"
    echo "output=${OUTPUT_FILE}"
    date
} > "${LOG_FILE}"

# The MQN binary expects one SMILES per line, so strip any trailing identifiers.
awk '{print $1}' "${INPUT_FILE}" > "${SMILES_ONLY_INPUT}"

{
    echo "prepared_input=${SMILES_ONLY_INPUT}"
    echo "prepared_records=$(wc -l < "${SMILES_ONLY_INPUT}")"
} >> "${LOG_FILE}"

"${MQN_BIN}" "${SMILES_ONLY_INPUT}" "${OUTPUT_FILE}"

END_EPOCH="$(date +%s)"
ELAPSED_SEC="$((END_EPOCH - START_EPOCH))"

{
    echo "Done processing"
    date
    echo "elapsed_seconds=${ELAPSED_SEC}"
} >> "${LOG_FILE}"
