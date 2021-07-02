#!/bin/bash

cd $mainDir/derived_data/10_Prediction

# Get descriptives
Rscript $mainDir/code/10_Prediction/10.2.0_get_descriptives.R $mainDir > $mainDir/code/10_Prediction/10.2.0_get_descriptives.R.log

# Run prediction
Rscript $mainDir/code/10_Prediction/10.2.1_predict_phenotypes.R HRS2 $mainDir
Rscript $mainDir/code/10_Prediction/10.2.1_predict_phenotypes.R WLS $mainDir
Rscript $mainDir/code/10_Prediction/10.2.1_predict_phenotypes.R UKB3 $mainDir