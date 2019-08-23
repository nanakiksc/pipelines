#!/usr/bin/env bash

# Usage: ./anno.sh <sample_name>

SM=$1

BUILD=hg38
DATAD=${HOME}/data/${BUILD}
DNSFP=${DATAD}/dbnsfp/dbNSFP4.0a_variant.all.vcf.gz
ANNDB=${HOME}/utils/annovar/humandb/

mkdir -p 4_ann

bcftools view \
    --samples ${SM} \
    3_var/${SM}.vcf.gz | \
    bgzip > 3_var/${SM}.tumor.vcf.gz

tabix 3_var/${SM}.tumor.vcf.gz

bcftools annotate \
    --threads 2 \
    --annotations ${DNSFP} \
    --columns \
ID,\
INFO/SIFT_score,\
INFO/SIFT_pred,\
INFO/Polyphen2_HDIV_score,\
INFO/Polyphen2_HDIV_pred,\
INFO/Polyphen2_HVAR_score,\
INFO/Polyphen2_HVAR_pred,\
INFO/FATHMM_score,\
INFO/FATHMM_pred,\
INFO/gnomAD_exomes_AF,\
INFO/clinvar_clnsig,\
INFO/clinvar_trait \
    --output 3_var/${SM}.anno.vcf \
    3_var/${SM}.tumor.vcf.gz

rm 3_var/${SM}.tumor.vcf.gz 3_var/${SM}.tumor.vcf.gz.tbi

${HOME}/utils/annovar/convert2annovar.pl \
    --format vcf4old \
    --includeinfo \
    --outfile 4_ann/${SM}.ann \
    3_var/${SM}.anno.vcf

${HOME}/utils/annovar/annotate_variation.pl \
    --outfile 4_ann/${SM} \
    --dbtype refGene \
    --buildver ${BUILD} \
    --exonsort \
    --hgvs \
    4_ann/${SM}.ann \
    ${ANNDB}

rm 4_ann/${SM}.ann 4_ann/${SM}.log

${HOME}/lab/pipelines/format.py ${SM} > 4_ann/${SM}.xls

#rm ${SM}.variant_function ${SM}.exonic_variant_function
