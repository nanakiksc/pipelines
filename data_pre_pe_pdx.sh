#!/usr/bin/env bash

# Usage: ./data_pre.sh <sample_name>

SM=$1

BWANC=4
BUILD=hg38
DATAD=${HOME}/data/${BUILD}
FASTA=${HOME}/data/other_genomes/metagenomes/meta_hg38_mm10/meta_hg38_mm10.fa
BWAIX=${HOME}/data/other_genomes/metagenomes/meta_hg38_mm10/bwa/meta_hg38_mm10
DBSNP=${DATAD}/dbsnp/All_20180418_chr.vcf.gz
MILLS=${DATAD}/gatk_resources/Mills_and_1000G_gold_standard.indels.${BUILD}.vcf.gz
PHASE=${DATAD}/gatk_resources/Homo_sapiens_assembly38.known_indels.vcf.gz

mkdir -p 2_pre

#${HOME}/utils/FastQC/fastqc \
#    ${SM} \
#    ${SM} \
#    ${TUMOR_FASTQ_R1} \
#    ${TUMOR_FASTQ_R2}

#cutadapt \
#    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
#    -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
#    -g CGACGCTCTTCCGATCT \
#    -G CGTGTGCTCTTCCGATCT \
#    -o 2_trim/${f}_1_trimmed.fastq.gz \
#    -p 2_trim/${f}_2_trimmed.fastq.gz \
#    1_raw/${f}_1.fastq.gz \
#    1_raw/${f}_2.fastq.gz

# Trim 3' TruSeq adapter.
#${HOME}/utils/Homer/bin/homerTools trim \
#    -5 GCTCTTCCGATCT \
#    -3 AGATCGGAAGAGCACACGTCT \
#    -mis 2 \
#    -minMatchLength 4 \
#    -min 20 \
#    ${SM}_1.fastq.gz \
#    ${SM}_2.fastq.gz

#rm ${SM}*.lengths

#${HOME}/perl5/bin/pairfq makepairs \
#    --forward ${SM}_1.fastq.gz.trimmed \
#    --reverse ${SM}_2.fastq.gz.trimmed \
#    --forw_paired ${SM}_1.sync.fastq.gz \
#    --rev_paired ${SM}_2.sync.fastq.gz \
#    --forw_unpaired ${SM}_1.unpaired.fastq.gz \
#    --rev_unpaired ${SM}_2.unpaired.fastq.gz \
#    --compress gzip

#rm ${SM}_*.fastq.gz.trimmed ${SM}_*.unpaired.fastq.gz

bwa mem -M -t ${BWANC} \
    -R '@RG\tID:1\tLB:Seq\tPL:illumina\tSM:'${SM}'\tPU:1.1.1' \
    ${BWAIX} \
    <(zcat 1_raw/${SM}_1.fastq.gz) \
    <(zcat 1_raw/${SM}_2.fastq.gz) | \
    samtools sort -o 2_pre/${SM}.raw.bam

#rm ${SM}_*.sync.fastq.gz

gatk MarkDuplicates \
    --INPUT 2_pre/${SM}.raw.bam \
    --METRICS_FILE 2_pre/${SM}.metrics.txt \
    --OUTPUT 2_pre/${SM}.dedup.bam

rm 2_pre/${SM}.raw.bam

gatk BaseRecalibrator \
    --input 2_pre/${SM}.dedup.bam \
    --reference ${FASTA} \
    --known-sites ${DBSNP} \
    --known-sites ${MILLS} \
    --known-sites ${PHASE} \
    --output 2_pre/${SM}.recal.table

gatk ApplyBQSR \
    --input 2_pre/${SM}.dedup.bam \
    --bqsr-recal-file 2_pre/${SM}.recal.table \
    --output 2_pre/${SM}.recal.bam

rm 2_pre/${SM}.dedup.bam 2_pre/${SM}.recal.table

samtools sort -o 2_pre/${SM}.bam 2_pre/${SM}.recal.bam

samtools index 2_pre/${SM}.bam 2_pre/${SM}.bai

rm 2_pre/${SM}.recal.bam 2_pre/${SM}.recal.bai
