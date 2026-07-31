#!/bin/bash

cd $PGI_Repo/derived_data/7_Genotypes/STRtwge/preImputationQC

## Variant QC: call rate> 0.98 , maf>0.01, hwe 1e-5
plink2 --bfile $PGI_Repo/original_data/genotype_data/STR/genotyped/twge/twge_scrambled \
--geno 0.02 \
--maf 0.01 \
--hwe 1e-5 midp \
--make-bed \
--out twge_geno02_maf01_hwe1e-5

#--------------------------------------------------------------------#

## Get individuals with per chromosome subject-level missingness rate >0.05
for chr in {1..22}
do
	plink2 --bfile twge_geno02_maf01_hwe1e-5 \
	--chr $chr \
	--missing \
	--out twge_geno02_maf01_hwe1e-5_chr$chr
done

for chr in {1..22}
do
	 sed 's/ \+/\t/g' twge_geno02_maf01_hwe1e-5_chr$chr.imiss | awk -F"\t" 'NR>1 && $7>0.05{print $2,$3}' OFS="\t" > chr${chr}_mind05.txt &
done
wait

cat chr* | sort | uniq > mind05_anychr.txt
rm twge_geno02_maf01_hwe1e-5_chr* chr*

#--------------------------------------------------------------------#

## Sex check
plink2 --bfile twge_geno02_maf01_hwe1e-5 \
	--check-sex \
	--out twge_geno02_maf01_hwe1e-5_sexcheck

sed 's/ \+/\t/g' twge_geno02_maf01_hwe1e-5_sexcheck.sexcheck |  awk -F"\t" '$6!="OK"{print $2,$3}' OFS="\t" > sex_mismatch.txt

#--------------------------------------------------------------------#

## Merge individuals failing sex check and per chr missingness and drop them
cat mind05_anychr.txt sex_mismatch.txt | sort | uniq > perchrmind_sexcheck_fail.txt
plink2 --bfile twge_geno02_maf01_hwe1e-5 \
	--remove perchrmind_sexcheck_fail.txt \
	--make-bed \
	--out twge_geno02_maf01_hwe1e-5_perchrmind_sexcheck

#--------------------------------------------------------------------#

## Drop het/hom outliers
plink2 --bfile twge_geno02_maf01_hwe1e-5_perchrmind_sexcheck \
	--autosome \
	--het \
	--out twge_geno02_maf01_hwe1e-5_perchrmind_sexcheck_het

sed 's/ \+/\t/g' twge_geno02_maf01_hwe1e-5_perchrmind_sexcheck_het.het | awk -F"\t" '$7>0.05 || $7<-0.03{print $2,$3}' > het_hom_outliers.txt

plink2 --bfile twge_geno02_maf01_hwe1e-5_perchrmind_sexcheck \
	--autosome \
	--remove het_hom_outliers.txt \
	--make-bed \
	--out twge_geno02_maf01_hwe1e-5_perchrmind_sexcheck_hethom

#--------------------------------------------------------------------#

## Drop ancestry outliers
N_outliers=1
touch ancestry_outliers.txt
while [ ${N_outliers} -gt 0 ]; do
	plink2 --bfile twge_geno02_maf01_hwe1e-5_perchrmind_sexcheck_hethom  \
	--remove ancestry_outliers.txt \
	--neighbour 1 3 \
	--out twge_geno02_maf01_hwe1e-5_perchrmind_sexcheck_hethom_IBS
	sed 's/ \+/\t/g' twge_geno02_maf01_hwe1e-5_perchrmind_sexcheck_hethom_IBS.nearest | awk -F"\t" '$6<-5{print $2,$3}' > tmp
	cat ancestry_outliers.txt tmp > tmp2
	mv tmp2 ancestry_outliers.txt
	N_outliers=$(wc -l tmp | cut -d" " -f1)
done
rm tmp*
## No outliers detected.


## Check for outliers in MDS plot
plink2 --bfile twge_geno02_maf01_hwe1e-5_perchrmind_sexcheck_hethom \
--within \
--cluster \
--mds-plot \
--out twge_geno02_maf01_hwe1e-5_perchrmind_sexcheck_hethom_mds
## No outliers

plink2 --bfile twge_geno02_maf01_hwe1e-5_perchrmind_sexcheck_hethom \
--freq \
--out twge_geno02_maf01_hwe1e-5_perchrmind_sexcheck_hethom_frq

#--------------------------------------------------------------------#

sh $PGI_Repo/code/7_Genotypes/7.0_HRCimputation_prep.sh $PGI_Repo/derived_data/7_Genotypes/STRtwge/preImputationQC/twge_geno02_maf01_hwe1e-5_perchrmind_sexcheck_hethom