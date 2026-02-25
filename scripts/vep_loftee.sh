#!/bin/bash
# -*- coding: utf-8 -*-
# JCMM Feb 16, 2026
# VEP + LoFTEE annotation module

set -euo pipefail

# ------------------------------------------------------------------
# Determine input VCF for vep
# ------------------------------------------------------------------

export PROTCOD="${TMPDIR}/${VCFBASENAME}.protcod.${CHR}.vcf.gz"


# /data/results/vep_loftee/bravo.pub.norm.chr21.vcf.gz
echo "==================================================" | tee -a "${LOG}"
echo "[VEP MODULE] Started at $(date)" | tee -a "${LOG}"
echo "Input VCF: ${PROTCOD}" | tee -a "${LOG}"
echo "==================================================" | tee -a "${LOG}"

# ------------------------------------------------------------------
# Safety check
# ------------------------------------------------------------------

if [[ ! -f "${PROTCOD}" ]]; then
    echo "ERROR: Input VCF not found: ${PROTCOD}" | tee -a "${LOG}"
    exit 1
fi

# ------------------------------------------------------------------
# Define outputs
# ------------------------------------------------------------------
export VEPVCF="${TMPDIR}/${VCFBASENAME}.vep.${CHR}.vcf"
# export VEPTAB="${TMPDIR}/${VCFBASENAME}.vep.${CHR}.tab"
# export HCPOS="${TMPDIR}/${VCFBASENAME}.hcpos.${CHR}.tab"
export HCVCF="${TMPDIR}/${VCFBASENAME}.HC.${CHR}.vcf.gz"

# FINALVCF must already be defined upstream
if [[ -z "${FINALVCF:-}" ]]; then
    echo "ERROR: FINALVCF variable not defined." | tee -a "${LOG}"
    exit 1
fi

# ------------------------------------------------------------------
# 1. Shallow annotation (tab output)
# ------------------------------------------------------------------

echo "[8/11] Running shallow VEP annotation (LoFTEE)..." | tee -a "${LOG}"

/opt/vep/src/ensembl-vep/vep \
  --assembly GRCh38 \
  --cache \
  --offline \
  --format vcf \
  --vcf \
  --input_file "${PROTCOD}" \
  --output_file "${VEPVCF}" \
  --dir_cache /opt/vep/vep_cache \
  --dir_plugins /opt/vep_plugins/loftee \
  --plugin LoF,loftee_path:/opt/vep_plugins/loftee/,human_ancestor_fa:/opt/vep_plugins/loftee/human_ancestor.fa.gz,conservation_file:/opt/vep_plugins/loftee/loftee.sql \
  --pick \
  --minimal \
  --no_stats \
  --fork "${THREADS}" \
  --force_overwrite \
    2>&1 | tee -a "${LOG}"

echo "Shallow annotation completed." | tee -a "${LOG}"

# ------------------------------------------------------------------
# 2. Extract high-confidence (HC) LoF positions
# ------------------------------------------------------------------

# # echo "[9/11] Extracting high-confidence (HC) LoF variants..." | tee -a "${LOG}"
# bcftools view -i 'INFO/LoF="HC"' "${VEPVCF}" \
#     --threads "${THREADS}" \
#     -Oz -o "${HCVCF}"

# grep -v '^#' "${VEPTAB}" \
#     | grep 'HC' \
#     | awk '{print $1}' \
#     | sed 's/:/\t/' \
#     | sort -u \
#     | awk -F'[\t-]' '{
#         if (NF==2) {
#             print $1"\t"$2"\t"$2
#         } else {
#             print $1"\t"$2"\t"$3
#         }
#       }' > "${HCPOS}"

# if [[ ! -s "${HCPOS}" ]]; then
#     echo "WARNING: No high-confidence LoF variants detected." | tee -a "${LOG}"
# fi

# bcftools view \
#     -R "${HCPOS}" \
#     --threads "${THREADS}" \
#     -Oz \
#     -o "${HCVCF}" \
#     "${PROTCOD}"

# bcftools index --threads "${THREADS}" -t "${HCVCF}"

# echo "High-confidence subset VCF created: ${HCVCF}" | tee -a "${LOG}"

# ------------------------------------------------------------------
# 3. Full VEP annotation on HC variants
# ------------------------------------------------------------------

# echo "[10/11] Running full VEP annotation on HC LoF variants..." | tee -a "${LOG}"

# /opt/vep/src/ensembl-vep/vep \
#     --assembly GRCh38 \
#     --cache \
#     --offline \
#     --format vcf \
#     --vcf \
#     --input_file "${HCVCF}" \
#     --output_file "${FINALVCF}" \
#     --dir_cache /opt/vep/vep_cache \
#     --dir_plugins /opt/vep_plugins/loftee \
#     --plugin LoF,loftee_path:/opt/vep_plugins/loftee/,human_ancestor_fa:/opt/vep_plugins/loftee/human_ancestor.fa.gz,conservation_file:/opt/vep_plugins/loftee/loftee.sql \
#     --everything \
#     --fork "${THREADS}" \
#     --buffer_size 10000 \
#     --force_overwrite \
#     2>&1 | tee -a "${LOG}"

# echo "Full annotation completed." | tee -a "${LOG}"

# echo "[VEP MODULE] Completed at $(date)" | tee -a "${LOG}"
# echo "Final annotated VCF: ${FINALVCF}" | tee -a "${LOG}"
echo "==================================================" | tee -a "${LOG}"