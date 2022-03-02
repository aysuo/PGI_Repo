#!/bin/bash

PGI(){
    cohort=$1
    method=$2

    rm -f $PGI_Repo/code/9_Scores/ss_single_${cohort}
    # Get list of sumstats: Pheno name on first column (e.g. SWB-Okbay), file path on second
    for pheno in $(cat $PGI_Repo/code/9_Scores/version_single_$cohort); do
        if ls $PGI_Repo/derived_data/4_MTAG_single/$pheno/*_trait_formatted*  1> /dev/null 2>&1
        then 
            path="$PGI_Repo/derived_data/4_MTAG_single/$pheno/${pheno}_trait_formatted.txt"
        else
            path="$PGI_Repo/derived_data/4_MTAG_single/$pheno/${pheno}_trait_1_formatted.txt"
        fi

        if [[ $method == "LDpred" ]]; then
            pheno=$(echo "$pheno-single" | sed 's/[1-9]//g')
        elif [[ $method == "SBayesR" ]]; then
            pheno=$(echo "$pheno-single")
            path=$(echo $path | sed 's/formatted/formatted_SBayesR/g')
        fi
        echo $pheno $path >> $PGI_Repo/code/9_Scores/ss_single_${cohort}
    done
    sh $PGI_Repo/code/9_Scores/9.0_PGS.sh single $cohort $PGI_Repo/derived_data/9_Scores/single $method
}


#for cohort in AH Dunedin EGCUT ELSA ERisk HRS2 HRS3 MCTFR STRpsych STRtwge STRyatssstage Texas UKB1 UKB2 UKB3 WLS
#do
#    PGI $cohort SBayesR
#done

for cohort in UKB1 UKB2 UKB3 HRS2 AH Dunedin EGCUT ELSA ERisk MCTFR STRpsych STRtwge STRyatssstage Texas WLS
do
    PGI $cohort SBayesR
done

