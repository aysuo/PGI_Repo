#!/bin/bash
refDir="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/ref_data/HRC_imputation_qc"
vcfsort="/disk/genetics/ukb/aokbay/bin/vcftools/src/perl/vcf-sort"
checkvcf="/disk/genetics/ukb/aokbay/bin/checkVCF"

pathGf=$1

getFreq(){
    pathGf=$1
    plink2 --bfile ${pathGf} \
    --freq \
    --out ${pathGf}_frq
}

HRC_1000G_check_bim(){
    pathGf=$1

    perl ${refDir}/HRC-1000G-check-bim.pl \
        -b ${pathGf}.bim \
        -f ${pathGf}_frq.afreq \
        -r ${refDir}/HRC.r1-1.GRCh37.wgs.mac5.sites.tab \
        -h
    
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
	    python2.7 $checkvcf/checkVCF.py \
            -r $checkvcf/hs37d5.fa \
            -o ${pathGf}_qc_chr${chr}_vcfcheck ${pathGf}_qc_chr${chr}_sorted.vcf.gz
    done
}


################################################################################################

HRCqc(){
    pathGf=$1
    dirGf=$(echo $pathGf | rev | cut -d"/" -f1 --complement | rev)
    cd $dirGf

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