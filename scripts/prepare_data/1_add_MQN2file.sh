#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

INFILE="${REPO_ROOT}/datasets/ZINC_ChEMBL34_hi_conf_1M.smi"
OUTFILE="${REPO_ROOT}/datasets/ZINC_ChEMBL34_hi_conf_1M_MQN.csv"
MQN_BIN="${REPO_ROOT}/build/bin/mqn_new"
LOG="${REPO_ROOT}/datasets/1_addMQN.log"
TMPDIR="$(mktemp -d "${TMPDIR:-/tmp}/mqn_add_XXXXXX")"

if command -v nproc >/dev/null 2>&1; then
    DEFAULT_SPLITS="$(nproc)"
else
    DEFAULT_SPLITS="$(getconf _NPROCESSORS_ONLN)"
fi
SPLITS="${SPLITS:-${DEFAULT_SPLITS}}"

cleanup() {
    rm -rf "${TMPDIR}"
}
trap cleanup EXIT

echo "Starting processing file ${INFILE} on $(date)" > "${LOG}"
echo "Processing $(wc -l < "${INFILE}") SMILES" >> "${LOG}"

if [ ! -x "${MQN_BIN}" ]; then
    echo "ERROR: MQN binary is missing or not executable: ${MQN_BIN}" >> "${LOG}"
    exit 1
fi

# Split in chunks to parallize it
echo "Splitting the file into ${SPLITS} chunks under ${TMPDIR}" >> "${LOG}"
split -n "l/${SPLITS}" -d "${INFILE}" "${TMPDIR}/chunk_"

# Calculate all the MQN
pids=()
chunk_files=( "${TMPDIR}"/chunk_* )
for CHUNK in "${chunk_files[@]}"; do
    CHUNKOUT="${CHUNK}_OUT"
    # we need to split the file, since MQN_BIN is dumb
    cut -d' ' -f1 "${CHUNK}" > "${CHUNK}-tmp_smi"
    cut -d' ' -f2 "${CHUNK}" | sed '/^CHEMBL/! s/^/ZINC/' > "${CHUNK}-tmp_id"

    "${MQN_BIN}" "${CHUNK}-tmp_smi" "${CHUNKOUT}-tmp_mqn" &
    pids+=($!)
done

# Wait for all background processes to finish
for pid in "${pids[@]}"; do
    echo "checking PID: ${pid}"
    wait "$pid" || { echo "MQN calculation with PID $pid failed" >> "$LOG"; exit 1; }
done

echo "All chunks have finished calculating on $(date)" >> "${LOG}"

# Write the results for all the chunks
for CHUNK in "${chunk_files[@]}"; do
    # QC: check the consitency of the output here, and merge if all is good
    if [ "$(wc -l < "${CHUNK}_OUT-tmp_mqn")" -eq "$(wc -l < "${CHUNK}-tmp_smi")" ] && [ "$(wc -l < "${CHUNK}_OUT-tmp_mqn")" -eq "$(wc -l < "${CHUNK}-tmp_id")" ]; then    #  the sizes are OK, check the detail of the CHUNK/CHUNKOUT SMILES lists
        if [ -z "$(cut -d';' -f1 < "${CHUNK}_OUT-tmp_mqn" | diff - "${CHUNK}-tmp_smi")" ]; then
            echo "Output chunk ${CHUNK} OK. Merging the output." >> "${LOG}"
            # Format,merge and clean (take advantage of the process to fix the nonsense output format of MQN_BIN)
            paste -d',' <(cut -d';' -f1 < "${CHUNK}_OUT-tmp_mqn") "${CHUNK}-tmp_id" <(sed 's/^[^;]*;//; s/;/,/g; s/^/"[&/; s/$/&]"/' "${CHUNK}_OUT-tmp_mqn") >> "${CHUNK}_OUT"
            rm "${CHUNK}-tmp_smi" "${CHUNK}-tmp_id" "${CHUNK}_OUT-tmp_mqn"
        else
            echo "ERROR: mismatch between input and output SMILES lists. Will not merge output." >> "${LOG}"
        fi
    else
        echo "ERROR: input and output lengths not matching. Truncated output?" >> "${LOG}"
        echo "$(wc -l < "${CHUNK}_OUT-tmp_mqn") and $(wc -l < "${CHUNK}-tmp_smi") and $(wc -l < "${CHUNK}-tmp_id")" >> "${LOG}"
    fi
done

# Merge the chunks
echo "Merging the file from chunks" >> "${LOG}"
echo "SMILES,ID,MQN" > "$OUTFILE"
for CHUNKOUT in "${chunk_files[@]/%/_OUT}"; do
    echo "Merging ${CHUNKOUT}"
    cat "${CHUNKOUT}" >> "${OUTFILE}"
done
echo "Finished processing file ${OUTFILE} on $(date)" >> "${LOG}"
echo "Records in: $(wc -l < "${INFILE}")" >> "${LOG}"
echo "Records out: $(wc -l < "${OUTFILE}") header included." >> "${LOG}"
