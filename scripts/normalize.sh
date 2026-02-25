#!/bin/bash
# -*- coding: utf-8 -*-
# JCMM Feb 16, 2026
# Normalization + genotype reduction module

set -euo pipefail

# ------------------------------------------------------------------
# Determine input VCF for normalization
# ------------------------------------------------------------------

if [[ "${INPUTTYPE}" == "site_only" || "${NUM_SAMPLES}" -eq 0 ]]; then
    echo "Site-only mode detected — using PASSVCF as normalization input." | tee -a "${LOG}"
    INPUT_FOR_NORM="${PASSVCF}"
else
    INPUT_FOR_NORM="${CLEANVCF}"
fi

echo "==================================================" | tee -a "${LOG}"
echo "[NORMALIZE MODULE] Started at $(date)" | tee -a "${LOG}"
echo "Input for normalization: ${INPUT_FOR_NORM}" | tee -a "${LOG}"
echo "==================================================" | tee -a "${LOG}"

# ------------------------------------------------------------------
# Safety check
# ------------------------------------------------------------------

if [[ ! -f "${INPUT_FOR_NORM}" ]]; then
    echo "ERROR: Normalization input not found: ${INPUT_FOR_NORM}" | tee -a "${LOG}"
    exit 1
fi

# ------------------------------------------------------------------
# Define outputs
# ------------------------------------------------------------------

export NORMVCF="${OUTPUTDIR}/${VCFBASENAME}.norm.${CHR}.vcf.gz"
export HTRHOMVCF="${TMPDIR}/${VCFBASENAME}.red.${CHR}.vcf.gz"
export PROTCOD="${TMPDIR}/${VCFBASENAME}.protcod.${CHR}.vcf.gz"

# ------------------------------------------------------------------
# Normalization
# ------------------------------------------------------------------

echo "[6/11] Running normalization..." | tee -a "${LOG}"

bcftools norm \
    -m -any \
    -f "${REFFASTA}" \
    --threads "${THREADS}" \
    -Oz \
    -o "${NORMVCF}" \
    "${INPUT_FOR_NORM}" \
    2>&1 | tee -a "${LOG}"

bcftools index --threads "${THREADS}" -t "${NORMVCF}"

echo "Normalization completed." | tee -a "${LOG}"

# ------------------------------------------------------------------
# Genotype-based reduction (remove hom-ref if GT exists)
# ------------------------------------------------------------------

echo "[7/11] Genotype-based reduction..." | tee -a "${LOG}"
echo "Input VCF: ${NORMVCF}" | tee -a "${LOG}"

NUM_SAMPLES_NORM=$(bcftools query -l "${NORMVCF}" | wc -l)

if [[ "${INPUTTYPE}" != "site_only" && "${NUM_SAMPLES_NORM}" -gt 0 ]]; then

    echo "Samples detected (${NUM_SAMPLES_NORM}). Checking GT field..." | tee -a "${LOG}"

    if bcftools view -h "${NORMVCF}" | grep -q '##FORMAT=<ID=GT'; then

        echo "GT field detected — removing hom-ref genotypes (0/0)." | tee -a "${LOG}"
        bcftools view "${NORMVCF}" \
            -i 'GT="het" || GT="homalt"' \
            --threads "${THREADS}" \
            -Oz -o "${HTRHOMVCF}"

        bcftools index --threads "${THREADS}" -t "${HTRHOMVCF}"

        echo "Genotype filtering completed." | tee -a "${LOG}"

    else
        echo "GT not defined in header. Skipping genotype filter." | tee -a "${LOG}"
        HTRHOMVCF="${NORMVCF}"
    fi

else
    echo "Site-only mode or no samples detected. Skipping genotype filter." | tee -a "${LOG}"
    HTRHOMVCF="${NORMVCF}"
fi

export HTRHOMVCF=${HTRHOMVCF}

# reduce for protein coding genes
echo "[6/11] filtering only protein coding genes..." | tee -a "${LOG}"

bcftools view -T "${REFCODING}" \
    -G "${HTRHOMVCF}" \
    --threads "${THREADS}" \
    -O z -o "${PROTCOD}"

bcftools index "${PROTCOD}"



echo "Final output for next step: ${PROTCOD}" | tee -a "${LOG}"
echo "[NORMALIZE MODULE] Completed at $(date)" | tee -a "${LOG}"
echo "==================================================" | tee -a "${LOG}"