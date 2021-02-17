#!/bin/bash

dirCode="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/code/9_Scores"
dirIn="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/6_MTAG_multi"
dirOut="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/9_Scores/multi"


for cohort in AH Dunedin EGCUT ELSA ERisk HRS2 HRS3 MCTFR STRpsych STRtwge STRyatssstage Texas UKB1 UKB2 UKB3 WLS
do
    rm -f $dirCode/ss_multi_${cohort}
    # Get list of sumstats: Pheno name on first column (e.g. SWB-Okbay), file path on second
    for pheno in $(cat $dirCode/version_multi_$cohort); do
        if [[ -f $dirIn/$pheno/${pheno}_trait_formatted.txt ]]; then
            path="$dirIn/$pheno/${pheno}_trait_formatted.txt"
        else 
            path="$dirIn/$pheno/${pheno}_trait_1_formatted.txt"
        fi
        pheno=$(echo "$pheno-multi" | sed 's/[1-9]//g')
        echo $pheno $path >> $dirCode/ss_multi_${cohort}
    done
    sh $dirCode/9.0_PGS.sh multi $cohort $dirOut
done