#!/bin/bash
cd /disk/genetics/PGS/Aysu/PGS_Repo_pipeline/code/10_Prediction

Rscript 10.2.1_predict_phenotypes.R HRS2 
Rscript 10.2.1_predict_phenotypes.R WLS 
Rscript 10.2.1_predict_phenotypes.R UKB3