#!/usr/bin/env bash

# Usage: ./vc_to.sh <tumor sample name>

TUMOR=$1

BUILD=hg38
DATAD=${HOME}/data/${BUILD}
FASTA=${DATAD}/masked/${BUILD}.analysisSet.fa
#FASTA=${HOME}/data/other_genomes/metagenomes/meta_hg38_mm10/meta_hg38_mm10.fa
DBSNP=${DATAD}/dbsnp/All_20180418_chr.vcf.gz
MILLS=${DATAD}/gatk_resources/Mills_and_1000G_gold_standard.indels.${BUILD}.vcf.gz
PHASE=${DATAD}/gatk_resources/Homo_sapiens_assembly38.known_indels.vcf.gz
EXOME=${HOME}/lab/pipelines/S07604514_hs_hg38/S07604514_Padded.bed
GRMLN=${DATAD}/axiom/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
GRMBA=${DATAD}/axiom/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.BIALLELIC.vcf.gz
#GRMLN=${DATAD}/gnomAD/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP_chr_header.vcf.gz
# Downloaded from: http://www.ensembl.org/info/docs/tools/vep/script/vep_example.html#gnomad

# TODO: use gnomAD instead of Axiom (currently no hg38 version available for gnomAD).

mkdir -p 3_var/filter_files

#############################################################
### Use custom af-of-alleles-not-in-resource only when    ###
### using germline-resource (with a value of 0.0000025    ###
### for gnomAD exomes and 0.00003125 for gnomAD genomes). ###
#############################################################
# GATK recommends 1 / (2 * database_size), which for gnomAD 2.0.1 exomes is:
# 1 / (2 * 170228) = 0.000002937237
# For Axiom:
# 1 / (2 * 1249) = 0.0004003203

gatk Mutect2 \
    --input 2_pre/${TUMOR}.bam \
    --tumor-sample ${TUMOR} \
    --reference ${FASTA} \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    --intervals ${EXOME} \
    --germline-resource ${GRMLN} \
    --af-of-alleles-not-in-resource 0.0004003203 \
    --output 3_var/${TUMOR}.raw.vcf.gz
#    --bam-output 3_var/${TUMOR}.HC.bam \
#    --panel-of-normals pon.vcf.gz \

if [ ! -f 3_var/filter_files/${TUMOR}.pileups.table ]
then
gatk GetPileupSummaries \
    --input 2_pre/${TUMOR}.bam \
    --variant ${GRMBA} \
    --intervals ${EXOME} \
    --output 3_var/filter_files/${TUMOR}.pileups.table
fi

if [ ! -f 3_var/filter_files/${TUMOR}.contamination ]
then
gatk CalculateContamination \
    --input 3_var/filter_files/${TUMOR}.pileups.table \
    --output 3_var/filter_files/${TUMOR}.contamination
fi

gatk FilterMutectCalls \
    --variant 3_var/${TUMOR}.raw.vcf.gz \
    --contamination-table 3_var/filter_files/${TUMOR}.contamination \
    --output 3_var/${TUMOR}.filtered.vcf.gz

rm 3_var/${TUMOR}.raw.vcf.gz 3_var/${TUMOR}.raw.vcf.gz.tbi Mutect2FilteringStats.tsv

if [ ! -f 3_var/filter_files/${TUMOR}.pre_adapter_detail_metrics ]
then
gatk CollectSequencingArtifactMetrics \
    --INPUT 2_pre/${TUMOR}.bam \
    --OUTPUT 3_var/filter_files/${TUMOR} \
    --REFERENCE_SEQUENCE ${FASTA}
fi

gatk FilterByOrientationBias \
    --variant 3_var/${TUMOR}.filtered.vcf.gz \
    --pre-adapter-detail-file 3_var/filter_files/${TUMOR}.pre_adapter_detail_metrics \
    --output 3_var/${TUMOR}.vcf.gz

rm 3_var/${TUMOR}.filtered.vcf.gz 3_var/${TUMOR}.filtered.vcf.gz.tbi 3_var/${TUMOR}.vcf.gz.summary

# gatk Mutect2 + CreateSomaticPanelOfNormals
