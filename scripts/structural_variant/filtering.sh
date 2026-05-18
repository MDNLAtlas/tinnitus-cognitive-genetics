#!/usr/bin/bash

# shellcheck disable=SC2164
# shellcheck disable=SC2086
# shellcheck disable=SC2045
# shellcheck disable=SC2012
# shellcheck disable=SC2013
# shellcheck disable=SC2126
# shellcheck disable=SC2126

cd $HOME/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/MOCA/Q1/

conda activate variant_calling

# Step 0.1: index compressed VCF files with tabix for cnvkit

for sample in $(ls cnvkit); do
    echo $sample
    cp cnvkit/$sample/${sample}.cnvcall.vcf cnvkit/$sample/${sample}.cnvcall.bk.vcf
    bgzip cnvkit/$sample/${sample}.cnvcall.bk.vcf
    tabix -p vcf cnvkit/$sample/${sample}.cnvcall.bk.vcf.gz
    mv cnvkit/$sample/${sample}.cnvcall.bk.vcf.gz cnvkit/$sample/${sample}.cnvcall.vcf.gz
    mv cnvkit/$sample/${sample}.cnvcall.bk.vcf.gz.tbi cnvkit/$sample/${sample}.cnvcall.vcf.gz.tbi
done

# Step 0.2: Create VCF status of each sample
samples=$(cat <(ls cnvkit| tr ' ' '\n') <(ls tiddit| tr ' ' '\n') <(ls manta| tr ' ' '\n')| sort| uniq)
echo -e "SampleID\tCNVkit\tManta\tTiddit" > script/sample_status.tsv

# Check each sample
for sample in $samples; do
    echo $sample
    # Check existence for each tool
    [[ -f cnvkit/$sample/${sample}.cnvcall.vcf.gz ]] && t1="Yes" || t1="No"
    [[ -f manta/$sample/${sample}.manta.diploid_sv.vcf.gz ]] && t2="Yes" || t2="No"
    [[ -f tiddit/$sample/${sample}.tiddit.vcf.gz ]] && t3="Yes" || t3="No"

    echo -e "${sample}\t${t1}\t${t2}\t${t3}" >> script/sample_status.tsv
done

grep -v "SampleTest" script/sample_status.tsv > script/sample_status_filtered.tsv
mv script/sample_status_filtered.tsv script/sample_status.tsv

# Step 0.3: generate uniq+short varID (ID1, ID2, ...) for each caller
# Help track the variants after filtering
for sample in $(sed '1d' script/sample_status.tsv| cut -f1); do
    echo $sample

    # CNVkit
    zcat cnvkit/$sample/${sample}.cnvcall.vcf.gz| \
    sed 's/.germline.call//g'| \
    awk 'BEGIN { OFS="\t" }
        /^#/ { print; next }
        { $3 = "ID" ++id; print }' | bgzip > cnvkit/$sample/${sample}.cnvcall.withID.vcf.gz
    tabix -p vcf cnvkit/$sample/${sample}.cnvcall.withID.vcf.gz

    # Manta
    zcat manta/$sample/${sample}.manta.diploid_sv.vcf.gz | \
    awk 'BEGIN { OFS="\t" }
        /^#/ { print; next }
        { $3 = "ID" ++id; print }' | bgzip > manta/$sample/${sample}.manta.diploid_sv.withID.vcf.gz
    tabix -p vcf manta/$sample/${sample}.manta.diploid_sv.withID.vcf.gz

    # Tiddit
    zcat tiddit/$sample/${sample}.tiddit.vcf.gz | \
    awk 'BEGIN { OFS="\t" }
        /^#/ { print; next }
        { $3 = "ID" ++id; print }' | bgzip > tiddit/$sample/${sample}.tiddit.withID.vcf.gz
    tabix -p vcf tiddit/$sample/${sample}.tiddit.withID.vcf.gz
done

# Step 1: Merge SVs from multiple callers using SVDB
# Refer to https://github.com/J35P312/SVDB
# Install svdb inside variant_calling conda environment
# git clone https://github.com/J35P312/SVDB.git
# cd SVDB
# pip install -e .
# svdb  --merge --help
while IFS=$'\t' read -r sample t1 t2 t3; do
    echo $sample
    mkdir -p merged/$sample
    if [[ "$t1" == "Yes" && "$t2" == "Yes" && "$t3" == "Yes" ]]; then
        echo "Processing $sample"
        svdb --merge \
            --vcf \
            cnvkit/$sample/${sample}.cnvcall.withID.vcf.gz:1 \
            manta/$sample/${sample}.manta.diploid_sv.withID.vcf.gz:2 \
            tiddit/$sample/${sample}.tiddit.withID.vcf.gz:3 \
            --priority 1,2,3 \
            --bnd_distance 2000 \
            --ins_distance 50 \
            --overlap 0.6 \
            --no_intra \
            --pass_only| bgzip > merged/$sample/${sample}.merged.vcf.gz
    fi
done < <(tail -n +2 script/sample_status.tsv)  # Skip header

# cp -r merged/ /home/mdnl/Documents/Tam_14Mar25/

# Kepp SVs detected by at least 2 callers
# dir="/home/mdnl/Documents/Tam_14Mar25/"
# for sample in `sed '1d' script/sample_status.tsv| cut -f1`; do
#     echo $sample
#     if [[ -f $dir/merged/$sample/${sample}.merged.vcf.gz ]]; then
#         echo "Filtering $sample"
#         zcat $dir/merged/$sample/${sample}.merged.vcf.gz | \
#           grep -e "^#" -e "SUPP_VEC=111" -e "SUPP_VEC=011" -e "SUPP_VEC=101" -e "SUPP_VEC=110" | \
#           bgzip > $dir/merged/$sample/${sample}.merged.flt.vcf.gz
#         tabix -p vcf $dir/merged/$sample/${sample}.merged.flt.vcf.gz
#     fi
# done

# Count no. of SVs in each merged VCF before and after filtering
echo -e "Sample\tnoOfSV_b\tnoOfSV_a" > script/sv_count.tsv
for sample in $(sed '1d' script/sample_status.tsv| cut -f1); do
    echo $sample

    echo "before filtering"
    b=$(zcat merged/$sample/${sample}.merged.vcf.gz | grep -v "^#" | wc -l)
    echo "after filtering"
    a=$(zcat merged/$sample/${sample}.merged.vcf.gz | grep -v "^#"| \
        grep -e "SUPP_VEC=111" -e "SUPP_VEC=011" -e "SUPP_VEC=101" -e "SUPP_VEC=110" | wc -l)
    echo -e "${sample}\t${b}\t${a}" >> script/sv_count.tsv
done
