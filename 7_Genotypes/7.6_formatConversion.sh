#!/bin/bash

source paths7


oxford2plink2(){
    gfIn=$1
    sampleIn=$2
    gfOut=$3

    for chr in {1..22}; do
        gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
        sampleInChr=$(echo "$sampleIn" | sed "s/\[1:22\]/${chr}/g")
        plink2 --gen $gfInChr \
            --sample $sampleInChr \
            --oxford-single-chr $chr \
            --make-pgen \
            --out ${gfOut}_chr$chr &
    done 
    wait
}


vcf2plink2(){
    gfIn=$1
    gfOut=$2
    snpidtype=$3

    if [[ $snpidtype == "ChrPosID" ]]; then
        for chr in {1..22}; 
        do
            gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
            plink2 --vcf $gfInChr \
                --make-pgen \
                --set-all-var-ids @:# \
                --double-id \
                --out ${gfOut}_chr$chr &
        done
    else
        for chr in {1..22};
        do
            gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
            plink2 --vcf $gfInChr \
                --double-id \
                --make-pgen \
                --out ${gfOut}_chr$chr &
        done
    fi

    wait
}

subsetHM3(){
    gfIn=$1
    snpidtype=$2
    gfOut=$3
    
    eval hm3='$'HM3_${snpidtype}

    if ! [[ $gfIn == *UKB* ]]; then
        sampleKeep=$(echo "$gfIn.psam" | sed "s/\[1:22\]/1/g")
    else
        shuf -n1000 $sample_UKB3 > $gfOut.sample
        sampleKeep=$gfOut.sample
    fi

    for chr in {1..22}; do
        gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")

        plink2 --pfile $gfInChr \
            --extract $hm3 \
            --keep $sampleKeep \
            --max-alleles 2 \
            --make-bed \
            --out ${gfOut}_chr$chr &
    done
    wait
    
    for chr in {1..22}; do
        echo ${gfOut}_chr$chr >> ${gfOut}_merge
    done

    plink1.9 --merge-list ${gfOut}_merge --make-bed --out ${gfOut}

    if [[ -f ${gfOut}-merge.missnp ]]; then
        for chr in {1..22}; do
            plink2 --bfile ${gfOut}_chr$chr --exclude ${gfOut}-merge.missnp --make-bed --out ${gfOut}_chr$chr.tmp &
        done
        wait

        rm ${gfOut}_merge
        for chr in {1..22}; do
            echo ${gfOut}_chr$chr.tmp >> ${gfOut}_merge
        done

        plink1.9 --merge-list ${gfOut}_merge --make-bed --out ${gfOut}
    fi

    rm -f ${gfOut}_merge ${gfOut}_chr* $gfOut.sample
}

rs2chrpos(){
    gfIn=$1
    gfOut=$2

    plink2 --bfile $gfIn \
        --set-all-var-ids @:# \
        --make-bed \
        --out $gfOut
}

kgp2rs(){
    gfIn=$1
    gfOut=$2

    awk -F"\t" '$5=="TRUE" && $6=="TRUE"{print $1,$2}' OFS="\t" $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv > $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv.both
    awk -F"\t" '$5=="TRUE" && $6=="FALSE"{print $1,$2}' OFS="\t" $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv > $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv.illumina
    awk -F"\t" '$5=="FALSE" && $6=="TRUE"{print $1,$2}' OFS="\t" $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv > $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv.gcc

    awk -F"\t" 'NR==FNR{a[$1]=$1;print;next}!($1 in a){print}' OFS="\t" $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv.both $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv.illumina > $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv.both.illumina
    awk -F"\t" 'NR==FNR{a[$1]=$1;print;next}!($1 in a){print}' OFS="\t" $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv.both.illumina $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv.gcc > $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv.preferred
    awk -F"\t" '!seen[$1]++{print}' OFS="\t" $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv.preferred > $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv.preferred.nodups
    rm $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv.both* $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv.illumina $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv.gcc $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv.preferred

    for chr in {1..22}; do
        plink2 --pfile ${gfIn}_chr$chr \
            --update-name $mainDir/original_data/ref_data/SNP_kgpID2rsID.csv.preferred.nodups 2 1 \
            --make-pgen \
            --out ${gfOut}_chr$chr
    done
}

#---------------------------------------------------------------------------------#

main () {
    for cohort in AH Dunedin EGCUT ELSA ERisk HRS2 HRS3 MCTFR Texas STRpsych STRtwge STRyatssstage UKB WLS; do
        mkdir -p $mainDir/derived_data/7_Genotypes/$cohort/plink2
        mkdir -p $mainDir/derived_data/7_Genotypes/$cohort/plink/HM3
    done

    vcf2plink2 "${p7_gf_1000G}" $mainDir/derived_data/7_Genotypes/1000G/plink2/1000Gph3 rsID
    subsetHM3  $mainDir/derived_data/7_Genotypes/1000G/plink2/1000Gph3_chr[1:22] rsID $mainDir/derived_data/7_Genotypes/1000G/plink/HM3/1000Gph3_HM3
    rs2chrpos $mainDir/derived_data/7_Genotypes/1000G/plink/HM3/1000Gph3_HM3 $mainDir/derived_data/7_Genotypes/1000G/plink/HM3/1000Gph3_HM3_chrpos

    for cohort in HRS2; do
        eval gf='$'p7_gf_${cohort} 
        eval sample='$'p7_sample_${cohort}
        oxford2plink2 "$gf" $sample $mainDir/derived_data/7_Genotypes/$cohort/plink2/$cohort
        kgp2rs $mainDir/derived_data/7_Genotypes/$cohort/plink2/$cohort $mainDir/derived_data/7_Genotypes/$cohort/plink2/$cohort.tmp
        rename ".tmp" "" $mainDir/derived_data/7_Genotypes/$cohort/plink2/*
        subsetHM3 $mainDir/derived_data/7_Genotypes/$cohort/plink2/${cohort}_chr[1:22] rsID $mainDir/derived_data/7_Genotypes/$cohort/plink/HM3/${cohort}_HM3
    done

    for cohort in WLS; do
        eval gf='$'p7_gf_${cohort} 
        eval sample='$'p7_sample_${cohort}
        oxford2plink2 "$gf" $sample $mainDir/derived_data/7_Genotypes/$cohort/plink2/$cohort
        subsetHM3 $mainDir/derived_data/7_Genotypes/$cohort/plink2/${cohort}_chr[1:22] rsID $mainDir/derived_data/7_Genotypes/$cohort/plink/HM3/${cohort}_HM3
    done

    for cohort in AH Dunedin EGCUT ELSA ERisk HRS3 MCTFR Texas STRpsych STRtwge STRyatssstage; do
        eval gf='$'p7_gf_${cohort} 
        vcf2plink2 "${gf}" $mainDir/derived_data/7_Genotypes/$cohort/plink2/$cohort ChrPosID
        subsetHM3 $mainDir/derived_data/7_Genotypes/$cohort/plink2/${cohort}_chr[1:22] ChrPosID $mainDir/derived_data/7_Genotypes/$cohort/plink/HM3/${cohort}_HM3
    done

    subsetHM3 $p7_gf_UKB rsID $mainDir/derived_data/7_Genotypes/$cohort/plink/HM3/${cohort}_HM3
    
}

main

