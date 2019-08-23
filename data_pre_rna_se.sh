#!/usr/bin/env bash

# Usage: ./data_pre.sh <sample_name>

SM=$1

STARC=4
BUILD=hg38
DATAD=${HOME}/data/${BUILD}
STARI=${DATAD}/masked/STAR/
FASTA=${DATAD}/masked/GRCh38.primary_assembly.genome.fa
DBSNP=${DATAD}/dbsnp/All_20180418_chr.vcf.gz
MILLS=${DATAD}/gatk_resources/Mills_and_1000G_gold_standard.indels.${BUILD}.vcf.gz
PHASE=${DATAD}/gatk_resources/Homo_sapiens_assembly38.known_indels.vcf.gz

mkdir -p 2_pre

#${HOME}/utils/FastQC/fastqc \
#    ${SM} \
#    ${TUMOR_FASTQ_R1}

#cutadapt \
#    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
#    -g CGACGCTCTTCCGATCT \
#    -o 2_trim/${f}.trimmed.fastq.gz \
#    1_raw/${f}.fastq.gz

# Trim 3' TruSeq adapter.
#${HOME}/utils/Homer/bin/homerTools trim \
#    -5 GCTCTTCCGATCT \
#    -3 AGATCGGAAGAGCACACGTCT \
#    -mis 2 \
#    -minMatchLength 4 \
#    -min 20 \
#    ${SM}.fastq.gz

#rm ${SM}*.lengths

#${HOME}/utils/STAR/bin/Linux_x86_64/STAR \
#    --runThreadN ${STARC} \
#    --genomeDir ${STARI} \
#    --readFilesIn 1_raw/${SM}.fastq.gz \
#    --readFilesCommand zcat \
#    --outSAMtype BAM SortedByCoordinate \
#    --outStd BAM_SortedByCoordinate | \
#    samtools addreplacerg \
#    -r '@RG\tID:1\tLB:Seq\tPL:illumina\tSM:'${SM}'\tPU:1.1.1' \
#    -o 2_pre/${SM}.raw.bam -

#gatk MarkDuplicates \
#    --INPUT 2_pre/${SM}.raw.bam \
#    --METRICS_FILE 2_pre/${SM}.metrics.txt \
#    --OUTPUT 2_pre/${SM}.dedup.bam

#rm 2_pre/${SM}.raw.bam
#rm Log.out Log.progress.out Log.std.out SJ.out.tab

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

#rm 2_pre/${SM}.dedup.bam 2_pre/${SM}.recal.table

samtools sort -o 2_pre/${SM}.bam 2_pre/${SM}.recal.bam

samtools index 2_pre/${SM}.bam 2_pre/${SM}.bai

#rm 2_pre/${SM}.recal.bam 2_pre/${SM}.recal.bai