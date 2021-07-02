#!/bin/bash

cd $mainDir/derived_data/11_Figures/PGI_histograms

for cohort in AH Dunedin EGCUT ELSA ERisk HRS2 HRS3 MCTFR STRpsych STRtwge STRyatssstage Texas UKB1 UKB2 UKB3 WLS
do
    for score in single multi;
    do
        Rscript $mainDir/code/11_Figures/11.0_PGI_histograms.R $cohort $score $mainDir
    done
done