#!/bin/bash
source $PGI_Repo/code/paths
cd $PGI_Repo/derived_data/10_Prediction

#-------------------------------------------------------------------------#

# WLS
stata -b $PGI_Repo/code/10_Prediction/10.1.0_save_WLS.do \
    $WLS_pheno_data

Rscript $PGI_Repo/code/10_Prediction/10.1.1_construct_WLS_phenotypes.R

#-------------------------------------------------------------------------#

# HRS
Rscript $PGI_Repo/code/10_Prediction/10.1.2_construct_HRS_phenotypes.R \
    $PGI_Repo/original_data/prediction_phenotypes/HRS

stata -b $PGI_Repo/code/10_Prediction/10.1.3_construct_HRS_SWB_part1.do \
    $PGI_Repo/original_data/prediction_phenotypes/HRS

Rscript $PGI_Repo/code/10_Prediction/10.1.4_construct_HRS_SWB_part2.R

stata -b $PGI_Repo/code/10_Prediction/10.1.5_construct_HRS_CIDI.do \
    $PGI_Repo/original_data/prediction_phenotypes/HRS

# UKB
stata -b $PGI_Repo/code/10_Prediction/10.1.6_construct_UKB_phenotypes.do \
    $PGI_Repo/derived_data/1_UKB_GWAS/tmp/IDs_assignPartition_ordered.txt \
    $UKB_pheno_data_1 \
    $UKB_pheno_data_2 \
    $UKB_pheno_data_3 \
    $UKB_covar_data \
    $PGI_Repo/code/10_Prediction/10.1.6_construct_UKB3_phenotypes.do.log

for i in 1 2 3
do
    sed -i 's/ $/ NA/g' input/UKB$i/*.pheno
done