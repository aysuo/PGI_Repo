#!/bin/bash

for cohort in AH Dunedin EGCUT ELSA ERisk HRS2 HRS3 MCTFR STRpsych STRtwge STRyatssstage Texas UKB1 UKB2 UKB3 WLS
do
    rm -f $mainDir/code/9_Scores/ss_multi_${cohort}
    # Get list of sumstats: Pheno name on first column (e.g. SWB-Okbay), file path on second
    for pheno in $(cat $mainDir/code/9_Scores/version_multi_$cohort); do
        if [[ -f $mainDir/derived_data/6_MTAG_multi/$pheno/${pheno}_trait_formatted.txt ]]; then
            path="$mainDir/derived_data/6_MTAG_multi/$pheno/${pheno}_trait_formatted.txt"
        else 
            path="$mainDir/derived_data/6_MTAG_multi/$pheno/${pheno}_trait_1_formatted.txt"
        fi
        pheno=$(echo "$pheno-multi" | sed 's/[1-9]//g')
        echo $pheno $path >> $mainDir/code/9_Scores/ss_multi_${cohort}
    done
    sh $mainDir/code/9_Scores/9.0_PGS.sh multi $cohort $mainDir/derived_data/9_Scores/multi
done