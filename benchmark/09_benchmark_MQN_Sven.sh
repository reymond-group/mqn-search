#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

MQN_BIN="${1:-${REPO_ROOT}/build/bin/mqn_new}"
INPUT_FILE="${2:-${REPO_ROOT}/datasets/ZINC_ChEMBL34_hi_conf_1M.smi}"
OUTPUT_FILE="${3:-${REPO_ROOT}/benchmark_out/ZINC_ChEMBL34_hi_conf_1M_Sven-MQN_benchmark.csv}"
LOG_FILE="${4:-${REPO_ROOT}/benchmark_out/09_benchmark_MQN_Sven.log}"

mkdir -p "$(dirname "${OUTPUT_FILE}")" "$(dirname "${LOG_FILE}")"

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
    echo "Started Sven MQN benchmark"
    echo "binary=${MQN_BIN}"
    echo "input=${INPUT_FILE}"
    echo "output=${OUTPUT_FILE}"
    date
} > "${LOG_FILE}"

"${MQN_BIN}" "${INPUT_FILE}" "${OUTPUT_FILE}"

END_EPOCH="$(date +%s)"
ELAPSED_SEC="$((END_EPOCH - START_EPOCH))"

{
    echo "Done processing"
    date
    echo "elapsed_seconds=${ELAPSED_SEC}"
} >> "${LOG_FILE}"
