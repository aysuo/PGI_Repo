#!/bin/bash

cd $PGI_Repo

# Get descriptives
Rscript code/10_Prediction/10.2.0_get_descriptives.R > $PGI_Repo/code/10_Prediction/10.2.0_get_descriptives.R.log

# Run prediction
for cohort in WLS UKB3 HRS2
do
    Rscript code/10_Prediction/10.2.1_predict_phenotypes.R $cohort 
done