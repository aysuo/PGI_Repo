#!/bin/bash

source $PGI_Repo/code/paths

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
            if ! [[ -f ${gfOut}_chr$chr.pgen ]]
            then
                gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
                plink2 --vcf $gfInChr \
                    --make-pgen \
                    --set-all-var-ids @:# \
                    --double-id \
                    --out ${gfOut}_chr$chr &
            fi
        done
    else
        for chr in {1..22};
        do
            if ! [[ -f ${gfOut}_chr$chr.pgen ]]
            then
                gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
                plink2 --vcf $gfInChr \
                    --double-id \
                    --make-pgen \
                    --out ${gfOut}_chr$chr &
            fi
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
        shuf -n1000 $sample_orig_UKB3 > $gfOut.sample
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

    awk -F"\t" '$5=="TRUE" && $6=="TRUE"{print $1,$2}' OFS="\t" $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv > $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.both
    awk -F"\t" '$5=="TRUE" && $6=="FALSE"{print $1,$2}' OFS="\t" $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv > $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.illumina
    awk -F"\t" '$5=="FALSE" && $6=="TRUE"{print $1,$2}' OFS="\t" $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv > $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.gcc

    awk -F"\t" 'NR==FNR{a[$1]=$1;print;next}!($1 in a){print}' OFS="\t" $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.both $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.illumina > $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.both.illumina
    awk -F"\t" 'NR==FNR{a[$1]=$1;print;next}!($1 in a){print}' OFS="\t" $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.both.illumina $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.gcc > $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.preferred
    awk -F"\t" '!seen[$1]++{print}' OFS="\t" $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.preferred > $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.preferred.nodups
    rm $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.both* $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.illumina $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.gcc $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.preferred

    for chr in {1..22}; do
        plink2 --pfile ${gfIn}_chr$chr \
            --update-name $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.preferred.nodups 2 1 \
            --make-pgen \
            --out ${gfOut}_chr$chr
    done
}
