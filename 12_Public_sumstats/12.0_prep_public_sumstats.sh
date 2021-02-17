#!/bin/bash
dirCode="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/code/12_Public_sumstats"
dirData="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data"
dirOut="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/12_Public_sumstats"
LDgf="/disk/genetics/HRC/aokbay/LDgf/Full/perChr_rsid/HRC_geno05_mind01_maf001_hwe1e-10_rel025_nooutliers"

cd $dirOut
mkdir -p single
mkdir -p multi

# Clumping function
clump(){
    pheno=$1
    ss=$2
    version=$3

    mkdir -p $version/$pheno

    for chr in {1..22}; do
        plink1.9 --bfile ${LDgf}_chr$chr \
        --clump $ss \
        --clump-snp-field SNPID \
        --clump-p1 0.01 \
        --clump-p2 0.01 \
        --clump-field PVALUE \
        --clump-r2 0.1 \
        --clump-kb 1000000 \
        --out $version/$pheno/chr$chr &
    done
    wait
  
    awk '{print $3,$1,$4,$5}' OFS="\t" $version/$pheno/chr*.clumped > $version/$pheno/${pheno}_lead_SNPs_p1e-2.txt
    awk -F"\t" 'NR==FNR{a[$1]=$4;next}($2 in a && a[$2]<=5e-8) || FNR==1{print}' OFS="\t" $version/$pheno/${pheno}_lead_SNPs_p1e-2.txt $ss > $version/$pheno/${pheno}_${version}_p5e-8_sumstats.txt
    cat $version/$pheno/chr*.log > $version/$pheno/clump.log
    rm $version/$pheno/chr* $version/$pheno/clump.log
}


main(){
    ### Single-trait sumstats
    rm -f $dirCode/ss_clump_single $dirCode/ss_gwide_single
    for SSversion in clump gwide; do
        rm -f $dirCode/ss_clump_$SSversion
        while read pheno; do
            trait=$(awk -F"\t" -v pheno=$pheno '$1==pheno{print $2}' ../5_LDSC/singleMTAG/ER2_table_full.txt)
            ss=$(echo $dirData/4_MTAG_single/$pheno/${pheno}_trait_${trait}_formatted.txt)
            if ! [[ -f $ss ]]; then
                ss=$(echo $dirData/4_MTAG_single/$pheno/${pheno}_trait_formatted.txt)
            fi

            # Write version & sumstats into a file to keep record 
            echo -e "$pheno\t$ss" >> $dirCode/ss_clump_$SSversion

            if ! [[ -f single/$pheno/${pheno}_lead_SNPs_p1e-2.txt ]] && [[ $SSversion == "clump" ]]; then
                echo "Clumping single $pheno trait $trait .."
                clump $pheno $ss single > single/${pheno}.log
                echo "Clumping for $pheno trait $trait finished."
            elif ! [[ -f single/$pheno/${pheno}_single_gwide_sumstats.txt ]] && [[ $SSversion == "gwide" ]]; then
                echo "Copying genome-wide sumstats for single $pheno trait $trait .."
                mkdir -p single/$pheno
                cp $ss single/$pheno/${pheno}_single_gwide_sumstats.txt
                echo "Genome-wide sumstats for single $pheno trait $trait copied."
            fi
        done < $dirCode/version_${SSversion}_single
    done    

    ### Multi-trait sumstats
    while read pheno; do
        ss=$(echo $dirData/6_MTAG_multi/$pheno/${pheno}_trait_1_formatted.txt)

        if ! [[ -f multi/$pheno/${pheno}_lead_SNPs_p1e-2.txt ]]; then
            echo "Clumping $pheno - multi.."
            clump $pheno $ss multi > multi/${pheno}.log
            echo "Clumping for $pheno - multi finished."
        fi
    done < $dirCode/version_clump_multi
  
}

main