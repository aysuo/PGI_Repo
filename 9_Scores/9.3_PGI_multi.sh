#!/bin/bash

for cohort in AH Dunedin EGCUT ELSA ERisk HRS2 HRS3 MCTFR STRpsych STRtwge STRyatssstage Texas UKB1 UKB2 UKB3 WLS
do
    rm -f $PGI_Repo/code/9_Scores/ss_multi_${cohort}
    # Get list of sumstats: Pheno name on first column (e.g. SWB-Okbay), file path on second
    for pheno in $(cat $PGI_Repo/code/9_Scores/version_multi_$cohort); do
        if [[ -f $PGI_Repo/derived_data/6_MTAG_multi/$pheno/${pheno}_trait_formatted.txt ]]; then
            path="$PGI_Repo/derived_data/6_MTAG_multi/$pheno/${pheno}_trait_formatted.txt"
        else 
            path="$PGI_Repo/derived_data/6_MTAG_multi/$pheno/${pheno}_trait_1_formatted.txt"
        fi
        pheno=$(echo "$pheno-multi" | sed 's/[1-9]//g')
        echo $pheno $path >> $PGI_Repo/code/9_Scores/ss_multi_${cohort}
    done
    sh $PGI_Repo/code/9_Scores/9.0_PGI.sh multi $cohort $PGI_Repo/derived_data/9_Scores/multi
done