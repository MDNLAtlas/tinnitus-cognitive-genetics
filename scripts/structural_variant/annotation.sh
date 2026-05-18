#!/usr/bin/bash

# shellcheck disable=SC2164
# shellcheck disable=SC2086
# shellcheck disable=SC2045
# shellcheck disable=SC2012
# shellcheck disable=SC2013
# shellcheck disable=SC2126
# shellcheck disable=SC2126

cd /home/mdnl/VibhasRDS/DATA_ANALYSIS_UNITI/MOCA/Q1/
# Step 2: Annotate merged VCF in each sample with VEP
# mkdir -p vep_db
# wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz
# wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz.tbi
# wget https://kircherlab.bihealth.org/download/CADD-SV/v1.1/1000G_phase3_SVs.tsv.gz
# wget https://kircherlab.bihealth.org/download/CADD-SV/v1.1/1000G_phase3_SVs.tsv.gz.tbi

# Step 2: Annotate merged SV VCF with annotSV
conda activate annotSV
# sample="SampleUG66"
# input="/home/mdnl/VibhasRDS/DATA_ANALYSIS_UNITI/MOCA/Q1/merged/$sample/${sample}.merged.flt.vcf.gz"
dir="/home/mdnl/Documents/Tam_14Mar25/"
annotDir="/home/mdnl/Documents/Tam_14Mar25/ref/AnnotSV_annotations"
for sample in $(sed '1d' script/sample_status.tsv| cut -f1); do
    echo $sample
    if [[ -f merged/$sample/${sample}.merged.vcf.gz ]]; then
        
        echo "Filtering $sample"
        mkdir -p $dir/merged/$sample
        cp merged/$sample/${sample}.merged.vcf.gz $dir/merged/$sample/${sample}.merged.vcf.gz
        zcat $dir/merged/$sample/${sample}.merged.vcf.gz | \
            grep -e "^#" -e "SUPP_VEC=111" -e "SUPP_VEC=011" -e "SUPP_VEC=101" -e "SUPP_VEC=110" | \
            bgzip > $dir/merged/$sample/${sample}.merged.flt.vcf.gz
        tabix -p vcf $dir/merged/$sample/${sample}.merged.flt.vcf.gz

        echo "Annotating $sample"
        input="$dir/merged/$sample/${sample}.merged.flt.vcf.gz"

        AnnotSV -SVinputFile $input \
            -annotationsDir $annotDir \
            -outputFile ${input/.vcf.gz/.annotSV.vcf.gz} \
            -genomeBuild GRCh38 
        
        mv $dir/merged/$sample/${sample}.merged.flt* merged/$sample/
    fi
done

# Step 2.1: Annotate CNV VCF with annotSV
# sample="SampleUG66"
# input="/home/mdnl/VibhasRDS/DATA_ANALYSIS_UNITI/MOCA/Q1/cnvkit/$sample/${sample}.cnvcall.withID.vcf.gz"
dir="/home/mdnl/Documents/Tam_14Mar25/"
annotDir="/home/mdnl/Documents/Tam_14Mar25/ref/AnnotSV_annotations"
for sample in $(sed '1d' script/sample_status.tsv| cut -f1); do
    echo $sample
    if [[ -f cnvkit/$sample/${sample}.cnvcall.withID.vcf.gz ]]; then
        
        echo "Filtering $sample"
        mkdir -p $dir/cnvkit/$sample
        cp cnvkit/$sample/${sample}.cnvcall.withID.vcf.gz $dir/cnvkit/$sample/${sample}.cnvcall.withID.vcf.gz

        echo "Annotating $sample"
        input="$dir/cnvkit/$sample/${sample}.cnvcall.withID.vcf.gz"

        AnnotSV -SVinputFile $input \
            -annotationsDir $annotDir \
            -outputFile ${input/.vcf.gz/.annotSV.vcf.gz} \
            -genomeBuild GRCh38 
        
        mv $dir/cnvkit/$sample/${sample}.cnvcall.withID.annotSV* merged/$sample/
    fi
done

# Step 2.2: Annotate SUPP_VEC=001 or SUPP_VEC=010 VCF with annotSV
# sample="SampleUG66"
# input="/home/mdnl/VibhasRDS/DATA_ANALYSIS_UNITI/MOCA/Q1/cnvkit/$sample/${sample}.cnvcall.withID.vcf.gz"
dir="/home/mdnl/Documents/Tam_14Mar25/"
annotDir="/home/mdnl/Documents/Tam_14Mar25/ref/AnnotSV_annotations"
for sample in $(sed '1d' script/sample_status.tsv| cut -f1); do
    echo $sample
    if [[ -f merged/$sample/${sample}.merged.vcf.gz ]]; then
        
        echo "Filtering $sample"
        mkdir -p $dir/merged/$sample
        cp merged/$sample/${sample}.merged.vcf.gz $dir/merged/$sample/${sample}.merged.vcf.gz

        zcat $dir/merged/$sample/${sample}.merged.vcf.gz | grep -e "^#" -e "SUPP_VEC=010"| \
            grep -e "^#" -e "PASS"| \
            bgzip > $dir/merged/$sample/${sample}.merged.flt.010.vcf.gz
        tabix -p vcf $dir/merged/$sample/${sample}.merged.flt.010.vcf.gz

        echo "Annotating $sample"
        input="$dir/merged/$sample/${sample}.merged.flt.010.vcf.gz"

        AnnotSV -SVinputFile $input \
            -annotationsDir $annotDir \
            -outputFile ${input/.vcf.gz/.annotSV.vcf.gz} \
            -genomeBuild GRCh38 
        
        mv $dir/merged/$sample/${sample}.merged.flt.010* merged/$sample/

        zcat $dir/merged/$sample/${sample}.merged.vcf.gz | grep -e "^#" -e "SUPP_VEC=001"| \
            grep -e "^#" -e "PASS"| \
            bgzip > $dir/merged/$sample/${sample}.merged.flt.001.vcf.gz
        tabix -p vcf $dir/merged/$sample/${sample}.merged.flt.001.vcf.gz

        echo "Annotating $sample"
        input="$dir/merged/$sample/${sample}.merged.flt.001.vcf.gz"

        AnnotSV -SVinputFile $input \
            -annotationsDir $annotDir \
            -outputFile ${input/.vcf.gz/.annotSV.vcf.gz} \
            -genomeBuild GRCh38 
        
        mv $dir/merged/$sample/${sample}.merged.flt.001* merged/$sample/
        rm $dir/merged/$sample/${sample}.merged.vcf.gz
    fi
done
