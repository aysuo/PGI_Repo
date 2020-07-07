#!/bin/bash

dirIn=/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/genotype_data/ERisk/genotyped
dirOut=/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes/ERisk/preImputationQC
dirCode=/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/code/7_Genotypes

cd $dirOut

awk 'NR==FNR{a[$1]=$2;next}($1 in a){$4=a[$1];print;next}{print}' $dirIn/Mothers.fam $dirIn/Kids.fam > $dirOut/Kids.fam

## Variant QC: Call rate> 0.98 , maf>0.01, hwe 1e-5
plink1.9 --bed $dirIn/Kids.bed \
--bim $dirIn/Kids.bim \
--fam $dirOut/Kids.fam \
--bmerge $dirIn/Mothers \
--make-bed \
--out ERisk

plink1.9 --bfile $dirOut/ERisk \
--geno 0.02 \
--maf 0.01 \
--hwe 1e-4 midp \
--autosome \
--make-bed \
--out ERisk_autosome_geno02_maf01_hwe1e-4

#--------------------------------------------------------------------#

## Get individuals with per chromosome subject-level missingness rate >0.05
for chr in {1..22}; do
    plink2 --bfile ERisk_autosome_geno02_maf01_hwe1e-4 \
    --chr $chr \
    --missing \
    --out ERisk_autosome_geno02_maf01_hwe1e-4_chr$chr &
done
wait

for chr in {1..22}; do
    sed 's/ \+/\t/g' ERisk_autosome_geno02_maf01_hwe1e-4_chr$chr.smiss | awk -F"\t" 'NR>1 && $5>0.05{print $1,$2}' OFS="\t" > chr${chr}_mind05.txt &
done
wait

cat chr* | sort | uniq > mind05_anychr.txt
rm chr*
# No one to be dropped

#--------------------------------------------------------------------#

## Sex check
plink1.9 --bfile ERisk \
--check-sex \
--out ERisk_sexcheck
# 6 problems

awk '$5!="OK" && NR>1{print $1,$2}' OFS="\t" ERisk_sexcheck.sexcheck > sex_problems.txt

#------------------------#

## Het/hom outliers
plink1.9 --bfile ERisk_autosome_geno02_maf01_hwe1e-4 \
--family \
--het \
--out ERisk_autosome_geno02_maf01_hwe1e-4_het

awk '$6>0.05 || $6<-0.05{print $1,$2}' OFS="\t" ERisk_autosome_geno02_maf01_hwe1e-4_het.het > het_hom_outliers.txt

#------------------------#

cat sex_problems.txt het_hom_outliers.txt > sexcheck_het_hom_outliers.txt

plink2 --bfile ERisk_autosome_geno02_maf01_hwe1e-4 \
--remove sexcheck_het_hom_outliers.txt \
--make-bed \
--out ERisk_autosome_geno02_maf01_hwe1e-4_sexcheck_hethom

#--------------------------------------------------------------------#

sh $dirCode/7.0_HRCimputation_prep.sh $dirOut/ERisk_autosome_geno02_maf01_hwe1e-4_sexcheck_hethom