#!/usr/bin/bash

# shellcheck disable=SC2013
# shellcheck disable=SC2164

cd "$HOME/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/MOCA/Q1/APOE_risk"

(head -1 cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.tsv; \
awk -F'\t' '{OFS="\t"} {if($40=="rs429358" || $40=="rs7412") print;}' cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.tsv) > APOE_GT_rs429358_rs7412.tsv

paste <(sed -n '1p' APOE_GT_rs429358_rs7412.tsv| cut -f40,55-| tr '\t' '\n') \
    <(sed -n '2p' APOE_GT_rs429358_rs7412.tsv| cut -f40,55-| tr '\t' '\n') \
    <(sed -n '3p' APOE_GT_rs429358_rs7412.tsv| cut -f40,55-| tr '\t' '\n')| head -n -1| sed -e 's/gnomADg/SAMPLE_ID/g' > APOE_GT_rs429358_rs7412_LF.txt

# Manually create UNITI_Leuven_clinical_data_THI_GUEF.csv from 
# 1. df_datapaper_clin.csv from /home/mdnl/Documents/Tam_14Mar25/UNITI dataset
# 2. Final_Export_Leuven_06-2025_clinical-blood.xlsx from /home/mdnl/Documents/Tam_14Mar25/UNITI dataset
# 3. LV04 TTDP4C from /home/mdnl/VibhasRDS/DATA_ANALYSIS_UNITI/UNITI/UNITIexportforAntonio-2023-11-10-105233.csv

paste <(echo -e "UNITI_PATIENT_ID\tSAMPLE_ID") <(cat UNITI_Leuven_clinical_data_THI_GUEF.csv| sed -e 's/,/\t/g'| head -1) > THI_GUEF_score.tsv
for sample in $(sed '1d' sampleID.tsv); do
    patientID=$(grep -w "$sample" UNITI_ID_patient_blood_sample_fixed.tsv| cut -f3)
    echo "Processing sample: $sample, Patient: $patientID"
    cat UNITI_Leuven_clinical_data_THI_GUEF.csv| sed -e 's/,/\t/g'| awk -v patientID="$patientID" -v sample="$sample" -F'\t' '{OFS="\t"} {if($2~patientID) print patientID"\t"sample"\t"$0}' >> THI_GUEF_score.tsv
done

# RRKRA
(head -1 cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.tsv; \
    awk -F'\t' '{OFS="\t"} {if($9~/splice_donor|splice_acceptor|frameshift|stop_gained|start_lost/ && $36=="HC" && $33>=20 && ($41!="" && $41<0.05) && ($44!="" && $44<0.05) && $11=="PRKRA") print;}' cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.tsv) > cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.LoF.PRKRA.tsv

paste <(sed -n '1p' cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.LoF.PRKRA.tsv| cut -f40,55-| tr '\t' '\n') \
    <(sed -n '2p' cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.LoF.PRKRA.tsv| cut -f40,55-| tr '\t' '\n') \
    <(sed -n '3p' cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.LoF.PRKRA.tsv| cut -f40,55-| tr '\t' '\n') \
    <(sed -n '4p' cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.LoF.PRKRA.tsv| cut -f40,55-| tr '\t' '\n')| head -n -1| sed -e 's/gnomADg/variantID/g' > cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.LoF.PRKRA.LF.tsv

(head -1 cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.tsv; \
    awk -F'\t' '{OFS="\t"} {if($9~/missense/ && $33>=20 && ($41!="" && $41<0.05) && ($44!="" && $44<0.05) && $11=="EIF4G1") print;}' cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.tsv) > cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.missense.EIF4G1.tsv

paste <(sed -n '1p' cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.missense.EIF4G1.tsv| cut -f40,55-| tr '\t' '\n') \
    <(sed -n '2p' cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.missense.EIF4G1.tsv| cut -f40,55-| tr '\t' '\n') \
    <(sed -n '3p' cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.missense.EIF4G1.tsv| cut -f40,55-| tr '\t' '\n') \
    <(sed -n '4p' cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.missense.EIF4G1.tsv| cut -f40,55-| tr '\t' '\n')| head -n -1| sed -e 's/gnomADg/variantID/g' > cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.missense.EIF4G1.LF.tsv

# GBA synaptic genes
(head -1 cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.tsv; \
    awk -F'\t' 'NR==FNR {pat[$1]; next} ($11 in pat)' ../GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCASub.genes.samples.GBA_genes.v3.txt \
    cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.tsv) > GBA_synaptic_GT.tsv

# PD genes
(head -1 cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.tsv; \
    awk -F'\t' 'NR==FNR {pat[$1]; next} ($11 in pat)' PD_genes.txt \
    cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.tsv) > PD_GT.tsv

# PD genes (missense and LoF)
(head -1 PD_GT.tsv; \
    awk -F'\t' '{OFS="\t"} {if($9~/splice_donor|splice_acceptor|frameshift|stop_gained|start_lost|missense/) print;}' PD_GT.tsv) > PD_GT_missense_LoF.tsv

# PD endophenotype genes (missense and LoF)
(head -1 cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.tsv; \
    awk -F'\t' '{OFS="\t"} {if ($9~/splice_donor|splice_acceptor|frameshift|stop_gained|start_lost|missense/ && $11~/GBA|TMEM175/) print;}' \
    cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.tsv) > PD_endophenotype_GT_missense_LoF.tsv

# Frontaltemporal dementia (missense and LoF)
(head -1 cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.tsv; \
    awk -F'\t' '{OFS="\t"} {if ($9~/splice_donor|splice_acceptor|frameshift|stop_gained|start_lost|missense/ && $11~/MAPT|GRN|TBK1|C9orf72|TREM2|BIN1/) print;}' \
    cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.tsv) > AD_FTD_LBD_GT_missense_LoF.tsv

# https://link.springer.com/article/10.1186/s40246-023-00499-z#MOESM1
# TREM2 p.(R47H) and p.(R62C)


