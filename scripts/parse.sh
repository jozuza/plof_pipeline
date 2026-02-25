
#!/bin/bash
# -*- coding: utf-8 -*-
# JCMM Feb 16, 2026
# Parsing

############################################################
# USER CONFIGURATION
############################################################

export VCFBASENAME="bravo.pub"
export THREADS=7

############################################################
# FIXED PATHS
############################################################

export INPUTDIR="/data/results/vep_loftee"
export TMPDIR="/data/temp"
export OUTPUTDIR="/data/results/vep_loftee"

# 0. Generate a common header
echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tHet\tHom\t$(grep 'ID=CSQ' "${INPUTDIR}/${VCFBASENAME}.vep.chr20.vcf" | \
sed -e 's/.*Format: //;s/">//;s/|/\t/g')" > "${TEMPDIR}/${VCFBASENAME}.header.txt"

############################################################
# CHROMOSOME LOOP
############################################################

for chr in {22..1}; do

    export CHR="chr${chr}"
    echo "=================================================="
    echo "🧬  Processing chromosome: ${CHR}"
    echo "=================================================="
    echo "Pipeline started at $(date)"
    
    export VEPVCF="${INPUTDIR}/${VCFBASENAME}.vep.${CHR}.vcf"
    
    # 1. Filter for HC and LC 
    bcftools view --threads "${THREADS}" -i 'INFO/CSQ ~ "|HC|" || INFO/CSQ ~ "|LC|"' ${VEPVCF} > "${TEMPDIR}/${VCFBASENAME}.filt.${CHR}.vcf"
    
    # 2.  parse
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%Het\t%Hom\t%CSQ\n' "${TEMPDIR}/${VCFBASENAME}.filt.${CHR}.vcf" | sed 's/|/\t/g' > "${TEMPDIR}/${VCFBASENAME}.parse.${CHR}.tsv"

    # 3. Junta tudo na tabela final
    cat "${TEMPDIR}/${VCFBASENAME}.header.txt" "${TEMPDIR}/${VCFBASENAME}.parse.${CHR}.tsv" > "${OUTPUTDIR}/${VCFBASENAME}.vep.parsed.${CHR}.tsv"

    echo "✅  Chromosome ${CHR} completed at $(date)" 
    echo "=================================================="

done

echo "🎉  Pipeline completed successfully at $(date)"
