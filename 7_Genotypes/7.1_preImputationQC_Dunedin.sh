#!/bin/bash

cd $gf_dir_genotyped_Dunedin/preImputationQC

## Variant QC: Call rate > 0.98 , maf>0.01, hwe 1e-4
plink2 --bfile $gf_orig_genotyped_Dunedin \
	--geno 0.02 \
	--maf 0.01 \
	--hwe 1e-4 midp \
    --autosome \
	--make-bed \
	--out Dunedin_autosome_geno02_maf01_hwe1e-4

#--------------------------------------------------------------------#

## Get individuals with per chromosome subject-level missingness rate > 0.05
for chr in {1..22}; do
    plink2 --bfile Dunedin_autosome_geno02_maf01_hwe1e-4 \
    --chr $chr \
    --missing \
    --out Dunedin_autosome_geno02_maf01_hwe1e-4_chr$chr &
done
wait

for chr in {1..22}; do
    sed 's/ \+/\t/g' Dunedin_autosome_geno02_maf01_hwe1e-4_chr$chr.smiss | awk -F"\t" 'NR>1 && $5>0.05{print $1,$2}' OFS="\t" > chr${chr}_mind05.txt &
done
wait

cat chr* | sort | uniq > mind05_anychr.txt
rm chr*
# No one to be dropped

#--------------------------------------------------------------------#

## Sex check
plink1.9 --bfile $gf_orig_genotyped_Dunedin \
--check-sex \
--out Dunedin_sexcheck
# No problems

#--------------------------------------------------------------------#

## Drop het/hom outliers
plink1.9 --bfile Dunedin_autosome_geno02_maf01_hwe1e-4 \
--het \
--out Dunedin_autosome_geno02_maf01_hwe1e-4_het

awk '$6>0.04 || $6<-0.04{print $1,$2}' OFS="\t" Dunedin_autosome_geno02_maf01_hwe1e-4_het.het > het_hom_outliers.txt

plink2 --bfile Dunedin_autosome_geno02_maf01_hwe1e-4 \
--remove het_hom_outliers.txt \
--make-bed \
--out Dunedin_autosome_geno02_maf01_hwe1e-4_hethom

#--------------------------------------------------------------------#

sh $PGI_Repo/code/7_Genotypes/7.0_HRCimputation_prep.sh Dunedin_autosome_geno02_maf01_hwe1e-4_hethom