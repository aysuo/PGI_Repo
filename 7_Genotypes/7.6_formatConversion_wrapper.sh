#!/bin/bash

source $PGI_Repo/code/paths
source $PGI_Repo/code/7_Genotypes/7.6.0_formatConversion.sh


for cohort in AH Dunedin EGCUT ELSA ERisk HRS2 HRS3 MCTFR Texas STRpsych STRtwge STRyatssstage UKB WLS; do
    mkdir -p $PGI_Repo/derived_data/7_Genotypes/$cohort/plink2
    mkdir -p $PGI_Repo/derived_data/7_Genotypes/$cohort/plink/HM3
done

vcf2plink2 "${gf_orig_1000G}" $PGI_Repo/derived_data/7_Genotypes/1000G/plink2/1000Gph3 rsID
subsetHM3  $PGI_Repo/derived_data/7_Genotypes/1000G/plink2/1000Gph3_chr[1:22] rsID $PGI_Repo/derived_data/7_Genotypes/1000G/plink/HM3/1000Gph3_HM3
rs2chrpos $PGI_Repo/derived_data/7_Genotypes/1000G/plink/HM3/1000Gph3_HM3 $PGI_Repo/derived_data/7_Genotypes/1000G/plink/HM3/1000Gph3_HM3_chrpos

for cohort in HRS2; do
    eval gf='$'gf_orig_${cohort} 
    eval sample='$'sample_orig_${cohort}
    oxford2plink2 "$gf" $sample $PGI_Repo/derived_data/7_Genotypes/$cohort/plink2/$cohort
    kgp2rs $PGI_Repo/derived_data/7_Genotypes/$cohort/plink2/$cohort $PGI_Repo/derived_data/7_Genotypes/$cohort/plink2/$cohort.tmp
    rename ".tmp" "" $PGI_Repo/derived_data/7_Genotypes/$cohort/plink2/*
    subsetHM3 $PGI_Repo/derived_data/7_Genotypes/$cohort/plink2/${cohort}_chr[1:22] rsID $PGI_Repo/derived_data/7_Genotypes/$cohort/plink/HM3/${cohort}_HM3
done

for cohort in WLS; do
    eval gf='$'gf_orig_${cohort} 
    eval sample='$'sample_orig_${cohort}
    oxford2plink2 "$gf" $sample $PGI_Repo/derived_data/7_Genotypes/$cohort/plink2/$cohort
    subsetHM3 $PGI_Repo/derived_data/7_Genotypes/$cohort/plink2/${cohort}_chr[1:22] rsID $PGI_Repo/derived_data/7_Genotypes/$cohort/plink/HM3/${cohort}_HM3
done

for cohort in AH Dunedin EGCUT ELSA ERisk HRS3 MCTFR Texas STRpsych STRtwge STRyatssstage; do
    eval gf='$'gf_orig_${cohort} 
    vcf2plink2 "${gf}" $PGI_Repo/derived_data/7_Genotypes/$cohort/plink2/$cohort ChrPosID
    subsetHM3 $PGI_Repo/derived_data/7_Genotypes/$cohort/plink2/${cohort}_chr[1:22] ChrPosID $PGI_Repo/derived_data/7_Genotypes/$cohort/plink/HM3/${cohort}_HM3
done

subsetHM3 $gf_orig_UKB rsID $PGI_Repo/derived_data/7_Genotypes/$cohort/plink/HM3/${cohort}_HM3

