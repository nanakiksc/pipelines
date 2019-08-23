#!/usr/bin/env bash

# Usage:
# conda activate gatk
# ./germline_vc.sh <sample name>

SM=$1

BUILD=hg38
DATAD=${HOME}/data/${BUILD}
FASTA=${DATAD}/masked/${BUILD}.analysisSet.fa
EXOME=${HOME}/lab/pipelines/S07604514_hs_hg38/S07604514_Padded.bed
DBSNP=${DATAD}/dbsnp/All_20180418_chr.vcf.gz
MILLS=${DATAD}/gatk_resources/Mills_and_1000G_gold_standard.indels.${BUILD}.vcf.gz
PHASE=${DATAD}/gatk_resources/Homo_sapiens_assembly38.known_indels.vcf.gz

gatk HaplotypeCaller \
    --input 2_pre/${SM}.bam \
    --reference ${FASTA} \
    --intervals ${EXOME} \
    --output 3_var/${SM}.raw.vcf.gz

gatk CNNScoreVariants \
    --input 2_pre/${SM}.bam \
    --variant 3_var/${SM}.raw.vcf.gz \
    --reference ${FASTA} \
    --tensor-type read_tensor \
    --intervals ${EXOME} \
    --output 3_var/${SM}.cnn.vcf.gz

gatk FilterVariantTranches \
    --variant 3_var/${SM}.cnn.vcf.gz \
    --resource ${DBSNP} \
    --resource ${MILLS} \
    --resource ${PHASE} \
    --intervals ${EXOME} \
    --output 3_var/${SM}.vcf.gz
