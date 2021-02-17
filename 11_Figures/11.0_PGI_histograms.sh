#!/bin/bash
dirCode="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/code"


for cohort in AH Dunedin EGCUT ELSA ERisk HRS2 HRS3 MCTFR STRpsych STRtwge STRyatssstage Texas UKB1 UKB2 UKB3 WLS
do
    for score in single multi;
    do
        Rscript $dirCode/11.0_PGI_histograms.R $cohort $score
    done
done