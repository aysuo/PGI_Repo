#!/bin/bash

cd $PGI_Repo
source $PGI_Repo/code/paths

# Get descriptives
Rscript code/10_Prediction/10.2.0_get_descriptives.R HRS $phen_dir_orig_HRS/randhrs1992_2018v2.csv $pc_dir_HRS $HRS_crosswalk > $PGI_Repo/code/10_Prediction/10.2.0_get_descriptives.HRS.log
Rscript code/10_Prediction/10.2.0_get_descriptives.R WLS $phen_dir_WLS/tmp/WLS_renamed.csv $pc_dir_WLS > $PGI_Repo/code/10_Prediction/10.2.0_get_descriptives.WLS.log
#Rscript code/10_Prediction/10.2.0_get_descriptives.R UKB3 $phen_dir_UKB/Prediction/tmp/pgs_repo.dta > $PGI_Repo/code/10_Prediction/10.2.0_get_descriptives.UKB3.log


# Run prediction
for cohort in HRS WLS
do
    eval pheno_dir='$'phen_dir_$cohort
    eval pgi_dir='$'pgi_dir_$cohort
    eval crosswalk='$'${cohort}_crosswalk
    eval pc_dir='$'$pc_dir_$cohort

    Rscript 10.2.1_predict_phenotypes.R \
        $cohort \
        EUR \
        "$pheno_dir" \
        "$pgi_dir" \
        "$crosswalk" \
        single_SBayesR,single_SBayesR_HM3 \
        single_SBayesR:single_SBayesR_HM3 \
        $PGI_Repo/derived_data/10_Prediction/output \
        1000 \
        "$pc_dir"/HRS_EUR_PCs.eigenvec
done

for cohort in HRS
do
    for ancestry in AFR EAS AMR; do
        eval pheno_dir='$'phen_dir_$cohort
        eval pgi_dir='$'pgi_dir_$cohort
        eval crosswalk='$'${cohort}_crosswalk
        eval pc_dir='$'$pc_dir_$cohort

        Rscript 10.2.1_predict_phenotypes.R \
            HRS \
            $ancestry \
            "$pheno_dir" \
            "$pgi_dir" \
            "$crosswalk" \
            single_SBayesR \
            NA \
            $PGI_Repo/derived_data/10_Prediction/output \
            1000 \
            "$pc_dir"/HRS_${ancestry}_PCs.eigenvec
    done
done


