
#!/bin/bash
# -*- coding: utf-8 -*-
# JCMM Feb 16, 2026
# PASS filter module
# Input:  INPUTVCF
# Output: PASSVCF

set -euo pipefail

echo "==================================================" | tee -a "${LOG}"
echo "[PASS FILTER MODULE] Started at $(date)" | tee -a "${LOG}"
echo "Input VCF: ${INPUTDIR}/${INPUTVCF}" | tee -a "${LOG}"
echo "==================================================" | tee -a "${LOG}"

# ------------------------------------------------------------------
# PASS filter
# ------------------------------------------------------------------

echo "[1/11] Filtering variants with FILTER=PASS..." | tee -a "${LOG}"

bcftools view \
    "${INPUTDIR}/${INPUTVCF}" \
    -f PASS \
    --threads "${THREADS}" \
    -Oz -o "${PASSVCF}" \
    2>&1 | tee -a "${LOG}"

# ------------------------------------------------------------------
# Index
# ------------------------------------------------------------------

bcftools index \
    --threads "${THREADS}" \
    -t "${PASSVCF}"

echo "PASS filter completed." | tee -a "${LOG}"
echo "Output PASSVCF: ${PASSVCF}" | tee -a "${LOG}"
echo "[PASS FILTER MODULE] Completed at $(date)" | tee -a "${LOG}"
echo "==================================================" | tee -a "${LOG}"