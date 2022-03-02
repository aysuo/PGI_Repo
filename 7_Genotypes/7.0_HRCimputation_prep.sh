#!/bin/bash

source $PGI_Repo/code/paths
pathGf=$1

if ! [[ -s $PGI_Repo/original_data/ref_data/HRC_imputation_qc/HRC-1000G-check-bim.pl ]]; then
    wget https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.10.zip -O $PGI_Repo/original_data/ref_data/HRC_imputation_qc/HRC-1000G-check-bim-v4.2.10.zip
    unzip $PGI_Repo/original_data/ref_data/HRC_imputation_qc/HRC-1000G-check-bim-v4.2.10.zip -d $PGI_Repo/original_data/ref_data/HRC_imputation_qc/
fi

if ! [[ -s $PGI_Repo/original_data/ref_data/HRC_imputation_qc/HRC.r1-1.GRCh37.wgs.mac5.sites.tab ]]; then
    wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz -O $PGI_Repo/original_data/ref_data/HRC_imputation_qc/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
    gunzip $PGI_Repo/original_data/ref_data/HRC_imputation_qc/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
fi


getFreq(){
    pathGf=$1
    plink2 --bfile ${pathGf} \
    --freq \
    --out ${pathGf}_frq
}

HRC_1000G_check_bim(){
    pathGf=$1

    perl $PGI_Repo/original_data/ref_data/HRC_imputation_qc/HRC-1000G-check-bim.pl \
        -b ${pathGf}.bim \
        -f ${pathGf}_frq.afreq \
        -r $PGI_Repo/original_data/ref_data/HRC_imputation_qc/HRC.r1-1.GRCh37.wgs.mac5.sites.tab \
        -h
    
    # "plink" calls plink1.7 on the server, change that to plink1.9
    sed -i 's/plink/plink1.9/g' Run-plink.sh
    sh Run-plink.sh | tee Run-plink.sh.log
}

convert2vcf(){
    pathGf=$1
    for chr in {1..22}; do
	    plink1.9 --bfile ${pathGf}-updated-chr$chr --recode vcf --out ${pathGf}_qc_chr$chr
    done
}

sortVcf(){
    pathGf=$1
    for chr in {1..22}; do
	    $vcfsort ${pathGf}_qc_chr$chr.vcf | bgzip -c > ${pathGf}_qc_chr${chr}_sorted.vcf.gz
    done
}

checkVcf(){
    pathGf=$1
    for chr in {1..22}; do
	    $python $checkvcf/checkVCF.py \
            -r $checkvcf/hs37d5.fa \
            -o ${pathGf}_qc_chr${chr}_vcfcheck ${pathGf}_qc_chr${chr}_sorted.vcf.gz
    done
}


################################################################################################

HRCqc(){
    pathGf=$1
    dirGf=$(echo $pathGf | rev | cut -d"/" -f1 --complement | rev)

    getFreq $pathGf
    HRC_1000G_check_bim $pathGf
    convert2vcf $pathGf
    sortVcf $pathGf
    checkVcf $pathGf
}

main(){
    HRCqc $pathGf
}

main