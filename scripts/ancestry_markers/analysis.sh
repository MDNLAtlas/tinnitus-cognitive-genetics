#!/usr/bin/bash

# shellcheck disable=SC2164
# shellcheck disable=SC2086
# shellcheck disable=SC2045
# shellcheck disable=SC2012
# shellcheck disable=SC2013
# shellcheck disable=SC2126
# shellcheck disable=SC2126

PRJ_DIR="$HOME/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/MOCA/Q1"

mkdir -p $PRJ_DIR/ancestry_markers

### GBA LoF genes
cat $PRJ_DIR/GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.LoF.GBA.v3.csv| \
    sed -e 's/,/\t/g' -e 's/"//g'| \
    cut -f3| \
    sed '1d' > $PRJ_DIR/ancestry_markers/GBA_LoF_genes.txt

### Extract GBA LoF genes from MCI and nonMCI CSV
MCI_VCF="$HOME/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/MOCA/Q1/GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_Q1.tsv"

(head -1 $MCI_VCF; \
    awk -F'\t' 'NR==FNR {pat[$1]; next} ($13 in pat)' $PRJ_DIR/ancestry_markers/GBA_LoF_genes.txt \
    $MCI_VCF)| \
    cut -f1-14| \
    awk -F'\t' 'BEGIN{OFS="\t"} NR==1 {print; next} {print $1, $2, $3, $4, $5, $7/$5, $7, $11, $13, $14}' \
    > $PRJ_DIR/ancestry_markers/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_Q1.LoF_GBA_genes.tsv

nonMCI_VCF="$HOME/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/MOCA/Q1/GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_normal.tsv"

(head -1 $nonMCI_VCF; \
    awk -F'\t' 'NR==FNR {pat[$1]; next} ($13 in pat)' $PRJ_DIR/ancestry_markers/GBA_LoF_genes.txt \
    $nonMCI_VCF)| \
    cut -f1-14| \
    awk -F'\t' 'BEGIN{OFS="\t"} NR==1 {print; next} {print $1, $2, $3, $4, $5, $7/$5, $7, $11, $13, $14}' \
    > $PRJ_DIR/ancestry_markers/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_normal.LoF_GBA_genes.tsv

### GBA missense genes
cat $PRJ_DIR/GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.missense.GBA.v3.csv| \
    sed -e 's/,/\t/g' -e 's/"//g'| \
    cut -f3| \
    sed '1d' > $PRJ_DIR/ancestry_markers/GBA_missense_genes.txt

### Extract GBA LoF genes from MCI and nonMCI CSV
MCI_VCF="$HOME/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/MOCA/Q1/GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_Q1.tsv"

(head -1 $MCI_VCF; \
    awk -F'\t' 'NR==FNR {pat[$1]; next} ($13 in pat)' $PRJ_DIR/ancestry_markers/GBA_missense_genes.txt \
    $MCI_VCF)| \
    cut -f1-14| \
    awk -F'\t' 'BEGIN{OFS="\t"} NR==1 {print $1, $2, $3, $4, $5, $6, $7, $11, $13, $14; next} {print $1, $2, $3, $4, $5, $7/$5, $7, $11, $13, $14}' \
    > $PRJ_DIR/ancestry_markers/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_Q1.missense_GBA_genes.tsv

nonMCI_VCF="$HOME/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/MOCA/Q1/GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_normal.tsv"

(head -1 $nonMCI_VCF; \
    awk -F'\t' 'NR==FNR {pat[$1]; next} ($13 in pat)' $PRJ_DIR/ancestry_markers/GBA_missense_genes.txt \
    $nonMCI_VCF)| \
    cut -f1-14| \
    awk -F'\t' 'BEGIN{OFS="\t"} NR==1 {print $1, $2, $3, $4, $5, $6, $7, $11, $13, $14; next} {print $1, $2, $3, $4, $5, $7/$5, $7, $11, $13, $14}' \
    > $PRJ_DIR/ancestry_markers/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_normal.missense_GBA_genes.tsv
