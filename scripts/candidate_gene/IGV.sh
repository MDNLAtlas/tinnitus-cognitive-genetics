#!/usr/bin/bash

# shellcheck disable=SC2164
# shellcheck disable=SC2086

cd ~/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/SV_BAMS
SAMPLE="LV29"

bam="$HOME/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/SV_BAMS/${SAMPLE}.recal.bam"
sub_bam="$HOME/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/MOCA/Q1/igv/${SAMPLE}.recal.bam"
samtools view -bh $bam chr2:178432025-178451105 > $sub_bam
samtools index $sub_bam

./tools/IGV_Linux_2.18.4/igv.sh

