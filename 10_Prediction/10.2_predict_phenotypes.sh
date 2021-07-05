#!/bin/bash

cd $PGI_Repo/derived_data/10_Prediction

# Get descriptives
Rscript $PGI_Repo/code/10_Prediction/10.2.0_get_descriptives.R $PGI_Repo > $PGI_Repo/code/10_Prediction/10.2.0_get_descriptives.R.log

# Run prediction
Rscript $PGI_Repo/code/10_Prediction/10.2.1_predict_phenotypes.R HRS2 $PGI_Repo
Rscript $PGI_Repo/code/10_Prediction/10.2.1_predict_phenotypes.R WLS $PGI_Repo
Rscript $PGI_Repo/code/10_Prediction/10.2.1_predict_phenotypes.R UKB3 $PGI_Repo