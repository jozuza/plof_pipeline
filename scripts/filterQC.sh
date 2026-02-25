#!/bin/bash
# -*- coding: utf-8 -*-
# JCMM Feb 16, 2026
# Modular QC + filtering block
# Compatible with cohort and site-only VCFs

set -euo pipefail

echo "==================================================" | tee -a ${LOG}
echo "[FILTER MODULE] Started at $(date)" | tee -a ${LOG}
echo "Input VCF: ${PASSVCF}" | tee -a ${LOG}
echo "INPUTTYPE: ${INPUTTYPE}" | tee -a ${LOG}
echo "==================================================" | tee -a ${LOG}

# ------------------------------------------------------------------
# Define outputs
# ------------------------------------------------------------------

export CLEANVCF=${TMPDIR}/${VCFBASENAME}.clean.${CHR}.vcf.gz
export PREFIX=${TMPDIR}/${VCFBASENAME}.${CHR}

export HWE_FILE=${PREFIX}.hwe
export LMISS_FILE=${PREFIX}.lmiss
export IMISS_FILE=${PREFIX}.imiss

export LOW_HWE=${PREFIX}.low_hwe_sites.txt
export LOW_LMISS=${PREFIX}.low_qual_sites.txt
export LOW_SITES=${PREFIX}.low_sites_to_remove.txt
export LOW_INDV=${PREFIX}.low_qual_ind.txt

# ------------------------------------------------------------------
# Detect number of samples (defensive programming)
# ------------------------------------------------------------------

NUM_SAMPLES=$(bcftools query -l "${PASSVCF}" | wc -l)

if [ "${INPUTTYPE}" = "site_only" ] || [ "${NUM_SAMPLES}" -eq 0 ]; then

    echo "⚠️  Site-only or no samples detected (${NUM_SAMPLES}). Skipping genotype-based QC." | tee -a ${LOG}

    # Simply propagate PASSVCF → CLEANVCF
    bcftools view "${PASSVCF}" \
        --threads ${THREADS} \
        -Oz -o "${CLEANVCF}"

    bcftools index --threads ${THREADS} -t "${CLEANVCF}"

    echo "QC skipped. CLEANVCF ready." | tee -a ${LOG}
    echo "[FILTER MODULE] Completed at $(date)" | tee -a ${LOG}
    exit 0
fi

# ------------------------------------------------------------------
# COHORT MODE: Run QC
# ------------------------------------------------------------------

echo "🧬 Cohort VCF detected (${NUM_SAMPLES} samples). Running QC..." | tee -a ${LOG}

# ---------------- HWE ----------------
echo "[HWE] Running Hardy-Weinberg test..." | tee -a ${LOG}
vcftools --gzvcf "${PASSVCF}" --hardy --out "${PREFIX}"

if [ -f "${HWE_FILE}" ]; then
    awk 'NR>1 && $6 <= 1e-8 {print $1 "\t" $2}' "${HWE_FILE}" > "${LOW_HWE}"
else
    echo "No HWE file generated." | tee -a ${LOG}
    touch "${LOW_HWE}"
fi

echo "# Variants failing HWE (p ≤ 1e-8):" | tee -a ${LOG}
wc -l "${LOW_HWE}" | tee -a ${LOG}

# ---------------- Missingness (site) ----------------
echo "[Missingness - Site]" | tee -a ${LOG}
vcftools --gzvcf "${PASSVCF}" --missing-site --out "${PREFIX}"

if [ -f "${LMISS_FILE}" ]; then
    awk 'NR>1 && $6 >= 0.05 {print $1 "\t" $2}' "${LMISS_FILE}" > "${LOW_LMISS}"
else
    echo "No site missingness file generated." | tee -a ${LOG}
    touch "${LOW_LMISS}"
fi

echo "# Variants with missingness ≥ 5%:" | tee -a ${LOG}
wc -l "${LOW_LMISS}" | tee -a ${LOG}

# ---------------- Combine site exclusions ----------------
cat "${LOW_HWE}" "${LOW_LMISS}" 2>/dev/null | sort -u > "${LOW_SITES}"

echo "# Total sites to remove:" | tee -a ${LOG}
wc -l "${LOW_SITES}" | tee -a ${LOG}

# ---------------- Missingness (individual) ----------------
echo "[Missingness - Individual]" | tee -a ${LOG}
vcftools --gzvcf "${PASSVCF}" --missing-indv --out "${PREFIX}"

if [ -f "${IMISS_FILE}" ]; then
    awk 'NR>1 && $5 >= 0.05 {print $1}' "${IMISS_FILE}" > "${LOW_INDV}"
else
    echo "No individual missingness file generated." | tee -a ${LOG}
    touch "${LOW_INDV}"
fi

echo "# Individuals with missingness ≥ 5%:" | tee -a ${LOG}
wc -l "${LOW_INDV}" | tee -a ${LOG}

# ------------------------------------------------------------------
# Apply filtering
# ------------------------------------------------------------------

echo "[Filtering] Removing low-quality sites and individuals..." | tee -a ${LOG}

if [ -s "${LOW_INDV}" ] && [ -s "${LOW_SITES}" ]; then

    bcftools view \
        --samples-file ^"${LOW_INDV}" \
        -T ^"${LOW_SITES}" \
        "${PASSVCF}" \
        --threads ${THREADS} \
        -Oz -o "${CLEANVCF}"

elif [ -s "${LOW_INDV}" ]; then

    bcftools view \
        --samples-file ^"${LOW_INDV}" \
        "${PASSVCF}" \
        --threads ${THREADS} \
        -Oz -o "${CLEANVCF}"

elif [ -s "${LOW_SITES}" ]; then

    bcftools view \
        -T ^"${LOW_SITES}" \
        "${PASSVCF}" \
        --threads ${THREADS} \
        -Oz -o "${CLEANVCF}"

else

    echo "No exclusions detected. Passing VCF unchanged." | tee -a ${LOG}

    bcftools view \
        "${PASSVCF}" \
        --threads ${THREADS} \
        -Oz -o "${CLEANVCF}"
fi

bcftools index --threads ${THREADS} -t "${CLEANVCF}"

echo "[FILTER MODULE] Completed at $(date)" | tee -a ${LOG}
echo "Output CLEANVCF: ${CLEANVCF}" | tee -a ${LOG}
echo "==================================================" | tee -a ${LOG}