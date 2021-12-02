#!/bin/bash

cd $PGI_Repo

# Get descriptives
Rscript code/10_Prediction/10.2.0_get_descriptives.R > $PGI_Repo/code/10_Prediction/10.2.0_get_descriptives.R.log

# Run prediction
for method in LDpred SBayesR
do
    Rscript code/10_Prediction/10.2.1_predict_phenotypes.R HRS2 $method
    Rscript code/10_Prediction/10.2.1_predict_phenotypes.R WLS $method
    Rscript code/10_Prediction/10.2.1_predict_phenotypes.R UKB3 $method
done