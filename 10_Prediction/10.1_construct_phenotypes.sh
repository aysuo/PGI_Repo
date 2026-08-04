#!/bin/bash
source $PGI_Repo/code/paths
cd $PGI_Repo/derived_data/10_Prediction

#-------------------------------------------------------------------------#

# WLS
stata -b $PGI_Repo/code/10_Prediction/10.1.0_save_WLS.do \
    $phen_dir_orig_WLS/wls_plg_13_08.dta \
    $phen_dir_WLS \
    $PGI_Repo/code/10_Prediction/10.1.0_save_WLS.do.log

Rscript $PGI_Repo/code/10_Prediction/10.1.1_construct_WLS_phenotypes.R \
    $phen_dir_WLS/tmp/WLS_renamed.csv \
    EUR \
    $phen_dir_WLS \
    $gf_dir_WLS/sampleQC/WLS_EUR_FID_IID.txt
    
    


#-------------------------------------------------------------------------#

# HRS
for ancestry in EUR EAS AFR AMR; do
    Rscript $PGI_Repo/code/10_Prediction/10.1.3_construct_HRS_phenotypes.R \
        $phen_dir_orig_HRS/randhrs1992_2018v2.csv \
        $phen_dir_HRS \
        $ancestry \
        NA \
        $gf_dir_HRS/sampleQC/HRS_${ancestry}_FID_IID.txt
done


# UKB
stata -b $PGI_Repo/code/10_Prediction/10.1.6_construct_UKB_phenotypes.do \
    $phen_dir_UKB/GWAS/tmp \
    $phen_dir_UKB/Prediction \
    $PGI_Repo/code/10_Prediction/10.1.6_construct_UKB3_phenotypes.do.log

for partition in 1 2 3
do
    sed -i 's/ $/ NA/g' $phen_dir_UKB/Prediction/UKB$partition/*.pheno
    sed -i 's/ $/ NA/g' $phen_dir_UKB/Prediction/UKB$partition/*.covar
done

