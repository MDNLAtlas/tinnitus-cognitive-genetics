#!/usr/bin/bash

# shellcheck disable=SC2086
# shellcheck disable=SC2043

conda activate variant_calling

# 1. GBA analysis for THI grades samples
for g in g{1..5}; do

    echo $g
    mkdir -p GBA_THI

    bcftools view -S GBA_THI/thi_${g}_samples.txt -Ou ../../UNITI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.vcf| \
    bcftools +fill-tags -Ou -- -t AC,AN,AF,AC_Hom,AC_Het | \
    bcftools view -i 'AC>0' -Oz -o GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.vcf.gz


    samples=$(cat GBA_THI/thi_${g}_samples.txt| tr '\n' '\t'| sed 's/\t$//g')
    
    n=$(wc -l GBA_THI/thi_${g}_samples.txt | awk '{print $1-1}')
    bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%INFO/AF\t%AC\t%AC_Het\t%AC_Hom\t%CSQ[\t%GT]\n' \
        -d -A tab -HH -p x \
        GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.vcf.gz| \
        sed -e 's/#CHROM/CHROM/g' -e "s/\(GT\t\)\{$n\}GT/$(echo -e "$samples")/g"| \
        awk '{OFS="\t"} NR==1 {for (i=10; i<=NF; i++) sub(/^x/, "", $i)} {print}' > GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.tsv

    (head -1 GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.tsv; \
        awk -F'\t' '{OFS="\t"} {if($11~/missense/ && $35>=20 && ($43!="." && $43<0.05) && ($46!="." && $46<0.05)) print;}' GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.tsv) > GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.missense.tsv

    (head -1 GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.tsv; \
        awk -F'\t' '{OFS="\t"} {if($11~/splice_donor|splice_acceptor|frameshift|stop_gained|start_lost/ && $38=="HC" && $35>=20 && ($43!="." && $43<0.05) && ($46!="." && $46<0.05)) print;}' GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.tsv) > GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.LoF.tsv
done

for g in mild severe; do

    echo $g
    mkdir -p GBA_THI

    bcftools view -S GBA_THI/thi_${g}_samples.txt -Ou ../../UNITI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.vcf| \
    bcftools +fill-tags -Ou -- -t AC,AN,AF,AC_Hom,AC_Het | \
    bcftools view -i 'AC>0' -Oz -o GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.vcf.gz


    samples=$(cat GBA_THI/thi_${g}_samples.txt| tr '\n' '\t'| sed 's/\t$//g')
    n=$(wc -l GBA_THI/thi_${g}_samples.txt | awk '{print $1-1}')
    bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%INFO/AF\t%AC\t%AC_Het\t%AC_Hom\t%CSQ[\t%GT]\n' \
        -d -A tab -HH -p x \
        GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.vcf.gz| \
        sed -e 's/#CHROM/CHROM/g' -e "s/\(GT\t\)\{$n\}GT/$(echo -e "$samples")/g"| \
        awk '{OFS="\t"} NR==1 {for (i=10; i<=NF; i++) sub(/^x/, "", $i)} {print}' > GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.tsv

    (head -1 GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.tsv; \
        awk -F'\t' '{OFS="\t"} {if($11~/missense/ && $35>=20 && ($43!="." && $43<0.05) && ($46!="." && $46<0.05)) print;}' GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.tsv) > GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.missense.tsv

    (head -1 GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.tsv; \
        awk -F'\t' '{OFS="\t"} {if($11~/splice_donor|splice_acceptor|frameshift|stop_gained|start_lost/ && $38=="HC" && $35>=20 && ($43!="." && $43<0.05) && ($46!="." && $46<0.05)) print;}' GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.tsv) > GBA_THI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.THI_${g}.LoF.tsv
done

# 2. GBA analysis for GUEF grades samples
for g in g{1..4}; do

    echo $g
    mkdir -p GBA_GUEF

    bcftools view -S GBA_GUEF/guef_${g}_samples.txt -Ou ../../UNITI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.vcf| \
    bcftools +fill-tags -Ou -- -t AC,AN,AF,AC_Hom,AC_Het | \
    bcftools view -i 'AC>0' -Oz -o GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.vcf.gz


    samples=$(cat GBA_GUEF/guef_${g}_samples.txt| tr '\n' '\t'| sed 's/\t$//g')
    n=$(wc -l GBA_GUEF/guef_${g}_samples.txt | awk '{print $1-1}')
    bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%INFO/AF\t%AC\t%AC_Het\t%AC_Hom\t%CSQ[\t%GT]\n' \
        -d -A tab -HH -p x \
        GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.vcf.gz| \
        sed -e 's/#CHROM/CHROM/g' -e "s/\(GT\t\)\{$n\}GT/$(echo -e "$samples")/g"| \
        awk '{OFS="\t"} NR==1 {for (i=10; i<=NF; i++) sub(/^x/, "", $i)} {print}' > GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.tsv

    (head -1 GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.tsv; \
        awk -F'\t' '{OFS="\t"} {if($11~/missense/ && $35>=20 && ($43!="." && $43<0.05) && ($46!="." && $46<0.05)) print;}' GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.tsv) > GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.missense.tsv

    (head -1 GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.tsv; \
        awk -F'\t' '{OFS="\t"} {if($11~/splice_donor|splice_acceptor|frameshift|stop_gained|start_lost/ && $38=="HC" && $35>=20 && ($43!="." && $43<0.05) && ($46!="." && $46<0.05)) print;}' GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.tsv) > GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.LoF.tsv
done

for g in severe; do

    echo $g
    mkdir -p GBA_GUEF

    bcftools view -S GBA_GUEF/guef_${g}_samples.txt -Ou ../../UNITI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.vcf| \
    bcftools +fill-tags -Ou -- -t AC,AN,AF,AC_Hom,AC_Het | \
    bcftools view -i 'AC>0' -Oz -o GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.vcf.gz


    samples=$(cat GBA_GUEF/guef_${g}_samples.txt| tr '\n' '\t'| sed 's/\t$//g')
    n=$(wc -l GBA_GUEF/guef_${g}_samples.txt | awk '{print $1-1}')
    bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%INFO/AF\t%AC\t%AC_Het\t%AC_Hom\t%CSQ[\t%GT]\n' \
        -d -A tab -HH -p x \
        GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.vcf.gz| \
        sed -e 's/#CHROM/CHROM/g' -e "s/\(GT\t\)\{$n\}GT/$(echo -e "$samples")/g"| \
        awk '{OFS="\t"} NR==1 {for (i=10; i<=NF; i++) sub(/^x/, "", $i)} {print}' > GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.tsv

    (head -1 GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.tsv; \
        awk -F'\t' '{OFS="\t"} {if($11~/missense/ && $35>=20 && ($43!="." && $43<0.05) && ($46!="." && $46<0.05)) print;}' GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.tsv) > GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.missense.tsv

    (head -1 GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.tsv; \
        awk -F'\t' '{OFS="\t"} {if($11~/splice_donor|splice_acceptor|frameshift|stop_gained|start_lost/ && $38=="HC" && $35>=20 && ($43!="." && $43<0.05) && ($46!="." && $46<0.05)) print;}' GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.tsv) > GBA_GUEF/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.GUEF_${g}.LoF.tsv
done

# 3. GBA analysis for PHQ9 grades samples
for g in g{1..4}; do

    echo $g
    mkdir -p GBA_PHQ9

    bcftools view -S GBA_PHQ9/phq9_${g}_samples.txt -Ou ../../UNITI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.vcf| \
    bcftools +fill-tags -Ou -- -t AC,AN,AF,AC_Hom,AC_Het | \
    bcftools view -i 'AC>0' -Oz -o GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.vcf.gz


    samples=$(cat GBA_PHQ9/phq9_${g}_samples.txt| tr '\n' '\t'| sed 's/\t$//g')
    n=$(wc -l GBA_PHQ9/phq9_${g}_samples.txt | awk '{print $1-1}')
    bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%INFO/AF\t%AC\t%AC_Het\t%AC_Hom\t%CSQ[\t%GT]\n' \
        -d -A tab -HH -p x \
        GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.vcf.gz| \
        sed -e 's/#CHROM/CHROM/g' -e "s/\(GT\t\)\{$n\}GT/$(echo -e "$samples")/g"| \
        awk '{OFS="\t"} NR==1 {for (i=10; i<=NF; i++) sub(/^x/, "", $i)} {print}' > GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.tsv

    (head -1 GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.tsv; \
        awk -F'\t' '{OFS="\t"} {if($11~/missense/ && $35>=20 && ($43!="." && $43<0.05) && ($46!="." && $46<0.05)) print;}' GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.tsv) > GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.missense.tsv

    (head -1 GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.tsv; \
        awk -F'\t' '{OFS="\t"} {if($11~/splice_donor|splice_acceptor|frameshift|stop_gained|start_lost/ && $38=="HC" && $35>=20 && ($43!="." && $43<0.05) && ($46!="." && $46<0.05)) print;}' GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.tsv) > GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.LoF.tsv
done

for g in severe; do

    echo $g
    mkdir -p GBA_PHQ9

    bcftools view -S GBA_PHQ9/phq9_${g}_samples.txt -Ou ../../UNITI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.vcf| \
    bcftools +fill-tags -Ou -- -t AC,AN,AF,AC_Hom,AC_Het | \
    bcftools view -i 'AC>0' -Oz -o GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.vcf.gz


    samples=$(cat GBA_PHQ9/phq9_${g}_samples.txt| tr '\n' '\t'| sed 's/\t$//g')
    n=$(wc -l GBA_PHQ9/phq9_${g}_samples.txt | awk '{print $1-1}')
    bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%INFO/AF\t%AC\t%AC_Het\t%AC_Hom\t%CSQ[\t%GT]\n' \
        -d -A tab -HH -p x \
        GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.vcf.gz| \
        sed -e 's/#CHROM/CHROM/g' -e "s/\(GT\t\)\{$n\}GT/$(echo -e "$samples")/g"| \
        awk '{OFS="\t"} NR==1 {for (i=10; i<=NF; i++) sub(/^x/, "", $i)} {print}' > GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.tsv

    (head -1 GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.tsv; \
        awk -F'\t' '{OFS="\t"} {if($11~/missense/ && $35>=20 && ($43!="." && $43<0.05) && ($46!="." && $46<0.05)) print;}' GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.tsv) > GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.missense.tsv

    (head -1 GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.tsv; \
        awk -F'\t' '{OFS="\t"} {if($11~/splice_donor|splice_acceptor|frameshift|stop_gained|start_lost/ && $38=="HC" && $35>=20 && ($43!="." && $43<0.05) && ($46!="." && $46<0.05)) print;}' GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.tsv) > GBA_PHQ9/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.PHQ9_${g}.LoF.tsv
done

# 3. GBA analysis for MOCA grades samples
for g in normal Q1; do

    echo $g
    mkdir -p GBA_revised

    bcftools view -S GBA_revised/moca_${g}_samples.txt -Ou ../../UNITI/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.vcf| \
    bcftools +fill-tags -Ou -- -t AC,AN,AF,AC_Hom,AC_Het | \
    bcftools view -i 'AC>0' -Oz -o GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_${g}.vcf.gz


    samples=$(cat GBA_revised/moca_${g}_samples.txt| tr '\n' '\t'| sed 's/\t$//g')
    n=$(wc -l GBA_revised/moca_${g}_samples.txt | awk '{print $1-1}')
    bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%INFO/AF\t%AC\t%AC_Het\t%AC_Hom\t%CSQ[\t%GT]\n' \
        -d -A tab -HH -p x \
        GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_${g}.vcf.gz| \
        sed -e 's/#CHROM/CHROM/g' -e "s/\(GT\t\)\{$n\}GT/$(echo -e "$samples")/g"| \
        awk '{OFS="\t"} NR==1 {for (i=10; i<=NF; i++) sub(/^x/, "", $i)} {print}' > GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_${g}.tsv

    (head -1 GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_${g}.tsv; \
        awk -F'\t' '{OFS="\t"} {if($11~/missense/ && $35>=20 && ($43!="." && $43<0.05) && ($46!="." && $46<0.05)) print;}' GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_${g}.tsv) > GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_${g}.missense.tsv

    (head -1 GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_${g}.tsv; \
        awk -F'\t' '{OFS="\t"} {if($11~/splice_donor|splice_acceptor|frameshift|stop_gained|start_lost/ && $38=="HC" && $35>=20 && ($43!="." && $43<0.05) && ($46!="." && $46<0.05)) print;}' GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_${g}.tsv) > GBA_revised/cohort_20230309.AB_GQ_DP.VQSR90.snp.recalibrated.PASS.annotated.MOCA_${g}.LoF.tsv
done