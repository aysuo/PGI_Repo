#!/bin/bash

cd $PGI_Repo/derived_data/11_Figures/PGI_histograms

for cohort in AH Dunedin EGCUT ELSA ERisk HRS2 HRS3 MCTFR STRpsych STRtwge STRyatssstage Texas UKB1 UKB2 UKB3 WLS
do
    for score in single multi;
    do
        Rscript $PGI_Repo/code/11_Figures/11.0_PGI_histograms.R $cohort $score $PGI_Repo
    done
done