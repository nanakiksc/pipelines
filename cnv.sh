#!/usr/bin/env bash

# Usage: ./cnv.sh <tumor sample name>

TUMOR=$1

BUILD=hg38
DATAD=${HOME}/data/${BUILD}
EXOME=${HOME}/lab/pipelines/S07604514_hs_hg38/S07604514_Padded.bed
FASTA=${DATAD}/masked/${BUILD}.analysisSet.fa

mkdir -p 3_cnv

if [ ! -f 3_cnv/targets_C.preprocessed.interval_list ]
then
# Padding 150 bp. From default of 250, since intervals are already padded 100 bp.
gatk PreprocessIntervals \
    --intervals ${EXOME} \
    --reference ${FASTA} \
    --bin-length 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    --padding 150 \
    --output 3_cnv/targets_C.preprocessed.interval_list
fi

# TODO: AnnotateIntervals + DenoiseReadCounts to remove CG bias. Put it here.

gatk CollectFragmentCounts \
    --input ${TUMOR}.bam \
    --intervals targets_C.preprocessed.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output ${TUMOR}.counts.hdf5
