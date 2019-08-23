#!/usr/bin/env bash

# Usage: ./coverage.sh <file_name.bam>

BAM=$1
BED=${HOME}/lab/pipelines/S07604514_hs_hg38/S07604514_Padded.bed

echo ${BAM}

echo 'Average coverage:'
samtools depth -a \
    ${BAM} | \
    awk '{ n += 1; s += $3 } END { print s / n }'

echo 'Average coverage GENE (on target):'
samtools depth -a \
    -b ${BED} \
    ${BAM} | \
    awk '{ n += 1; s += $3 } END { print s / n }'

echo '% on target:'
ONT=$(samtools depth \
    -b ${BED} \
    ${BAM} | \
    awk '{ s += $3 } END { print s }')
ALL=$(samtools depth \
    ${BAM} | \
    awk '{ s += $3 } END { print s }')
echo 'scale = 6;' ${ONT} '/' ${ALL} | bc
