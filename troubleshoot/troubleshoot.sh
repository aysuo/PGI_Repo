#!/bin/bash

## FIX MUNGE CHECK STATUS - doesn't work in first run

alias python=/homes/nber/aokbay/anaconda2/bin/python2.7
dirOut=/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/troubleshoot
hm3_snplist=/var/genetics/ukb/aokbay/reffiles/w_hm3.snplist
LDscores=/var/genetics/ukb/aokbay/bin/ldsc_old/eur_w_ld_chr/

allergyCat_old=/disk/genetics4/PGS/PGS_Repo/data/CLEANED_23andme/CLEANED.MAF_cat_allergy.dat_toQC.gz
allergyCat_new=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/3_QCd/23andMe/SEfilter/QC_AllergyCat-23andMe_2020_03_19/CLEANED.AllergyCat-23andMe
allergyDust_old=/disk/genetics4/PGS/PGS_Repo/data/CLEANED_23andme/CLEANED.MAF_dust_mites_allergy.dat_toQC.gz
allergyDust_new=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/3_QCd/23andMe/SEfilter/QC_AllergyDust-23andMe_2020_03_19/CLEANED.AllergyDust-23andMe
allergyPollen_old=/disk/genetics4/PGS/PGS_Repo/data/CLEANED_23andme/CLEANED.MAF_pollen_allergy.dat_toQC.gz
allergyPollen_new=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/3_QCd/23andMe/SEfilter/QC_AllergyPollen-23andMe_2020_03_19/CLEANED.AllergyPollen-23andMe
LDSC_old=/disk/genetics/ukb/aokbay/bin/ldsc_old
LDSC_new=/disk/genetics/ukb/aokbay/bin/ldsc

cd $dirOut

for pheno in allergyCat allergyDust allergyPollen; do 
#    eval ss='$'${pheno}_old
#    zcat $ss | sed 's/Z/ZZ/g' > ${pheno}_old.txt
    eval ${pheno}_old=${pheno}_old.txt
done


ldsc_munge_old () {
    ss=$1
    out=$2

	## Munge sumstats ##
	python ${LDSC_old}/munge_sumstats.py \
		--sumstats ${ss} \
		--out ${out} \
		--merge-alleles ${hm3_snplist}
}

ldsc_munge_new () {
    ss=$1
    out=$2

	## Munge sumstats ##
	python ${LDSC_new}/munge_sumstats.py \
		--sumstats ${ss} \
		--out ${out} \
		--merge-alleles ${hm3_snplist}
}

ldsc_h2_old () {
    ssMunged=$1
    out=$2

	## Get h2 and intercept ##
	python ${LDSC_old}/ldsc.py \
		--h2 ${ssMunged}  \
		--ref-ld-chr ${LDscores} \
		--w-ld-chr ${LDscores} \
		--out ${out}
}

ldsc_h2_new () {
    ssMunged=$1
    out=$2

	## Get h2 and intercept ##
	python ${LDSC_new}/ldsc.py \
		--h2 ${ssMunged}  \
		--ref-ld-chr ${LDscores} \
		--w-ld-chr ${LDscores} \
		--out ${out}
}


for pheno in allergyCat allergyDust allergyPollen; do 
    for ldsc_version in new; do
        for pheno_version in old new; do
            eval ss='$'${pheno}_${pheno_version}
            echo $ss
            ldsc_munge_${ldsc_version} $ss ${pheno}_ldsc${ldsc_version}_pheno${pheno_version}
            ldsc_h2_${ldsc_version} ${pheno}_ldsc${ldsc_version}_pheno${pheno_version}.sumstats.gz h2_${pheno}_ldsc${ldsc_version}_pheno${pheno_version}
        done
    done
done

awk 'NR==FNR{print;next}FNR>1{print}' OFS="\t" h2*.log > h2.txt




   




   