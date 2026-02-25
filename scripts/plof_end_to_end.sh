#!/bin/bash
# -*- coding: utf-8 -*-
# JCMM Feb 16, 2026
# Main pipeline controller

set -euo pipefail

############################################################
# USER CONFIGURATION
############################################################

export INPUTTYPE="site_only"   # site_only | cohort
export VCFBASENAME="bravo.pub"
export EXTENSIONNAME="vcf.gz"
export BUCKET="gs://projects-usp/public-datasets/TOPMed/references"
export THREADS=7

############################################################
# FIXED PATHS
############################################################

export INPUTDIR="/data/vcf_files"
export TMPDIR="/data/temp"
export OUTPUTDIR="/data/results/vep_loftee"
export REFFASTA="/data/references/human_genome/Homo_sapiens_assembly38.fasta"
export REFCODING="/data/references/human_genome/coding_genes.bed"

mkdir -p "${INPUTDIR}" "${TMPDIR}" "${OUTPUTDIR}"

############################################################
# PREP CLOUD
############################################################

bash /scripts/prep_cloud.sh

############################################################
# CHROMOSOME LOOP
############################################################

for chr in {22..21}; do

    export CHR="chr${chr}"
    echo "=================================================="
    echo "🧬  Processing chromosome: ${CHR}"
    echo "=================================================="

    ########################################################
    # LOGGING
    ########################################################

    export LOG=${OUTPUTDIR}/${VCFBASENAME}.${CHR}.pipeline.log
    touch "${LOG}"

    echo "Pipeline started at $(date)" | tee -a "${LOG}"
    echo "INPUTTYPE: ${INPUTTYPE}" | tee -a "${LOG}"

    ########################################################
    # STEP 0 – DOWNLOAD
    ########################################################
    echo "==================================================" | tee -a ${LOG}
    export INPUTVCF="${CHR}.${VCFBASENAME}.${EXTENSIONNAME}"

    if [ -f "${INPUTDIR}/${INPUTVCF}" ]; then
        echo "⏭️  Step 0 (Download): Skipped." | tee -a "${LOG}"
    else
        echo "🚀  Step 0 (Download)" | tee -a "${LOG}"
        gsutil -m cp "${BUCKET}/${INPUTVCF}*"  ${INPUTDIR}
    fi
    echo "==================================================" | tee -a ${LOG}
    ########################################################
    # STEP 1 – PASS FILTER
    ########################################################

    export PASSVCF=${TMPDIR}/${VCFBASENAME}.pass.${CHR}.vcf.gz

    if [ -f "${PASSVCF}" ]; then
        echo "⏭️  Step 1 (PASS Filter): Skipped." | tee -a "${LOG}"
    else
        echo "🚀  Step 1 (PASS Filter)" | tee -a "${LOG}"
        bash /scripts/pass_filter.sh
    fi

    ########################################################
    # STEP 2A – QC FILTER
    ########################################################

    export CLEANVCF=${TMPDIR}/${VCFBASENAME}.clean.${CHR}.vcf.gz

    if [ -f "${CLEANVCF}" ]; then
        echo "⏭️  Step 2A (QC Filter): Skipped." | tee -a "${LOG}"
    else
        echo "🚀  Step 2A (QC Filter)" | tee -a "${LOG}"
        bash /scripts/filterQC.sh
    fi

    ########################################################
    # STEP 2B – NORMALIZATION
    ########################################################

    export NORMVCF=${OUTPUTDIR}/${VCFBASENAME}.norm.${CHR}.vcf.gz

    if [ -f "${NORMVCF}" ]; then
        echo "⏭️  Step 2B (Normalization): Skipped." | tee -a "${LOG}"
    else
        echo "🚀  Step 2B (Normalization)" | tee -a "${LOG}"
        bash /scripts/normalize.sh
    fi

    ########################################################
    # STEP 3 – VEP LOFTEE
    ########################################################

    export FINALVCF=${OUTPUTDIR}/${VCFBASENAME}.veploftee.${CHR}.vcf

    if [ -f "${FINALVCF}" ]; then
        echo "⏭️  Step 3 (VEP LOFTEE): Skipped." | tee -a "${LOG}"
    else
        echo "🚀  Step 3 (VEP LOFTEE)" | tee -a "${LOG}"
        bash /scripts/vep_loftee.sh
    fi

    echo "✅  Chromosome ${CHR} completed at $(date)" | tee -a "${LOG}"
    echo "==================================================" | tee -a "${LOG}"

done

echo "🎉  Pipeline completed successfully at $(date)"