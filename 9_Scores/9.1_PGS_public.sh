#!/bin/bash

source paths9

# Get list of sumstats: Pheno name on first column (e.g. SWB-Okbay), file path on second
for cohort in Dunedin ERisk HRS2 WLS UKB3; 
do
    rm -f $mainDir/code/9_Scores/ss_public_${cohort}
	for pheno in $(cat $mainDir/code/9_Scores/version_public_$cohort); do 
			path=$(ls $mainDir/derived_data/3_QCd/public_scores/SEfilter/QC_${pheno}_*/*.gz)
			pheno=$(echo $path | rev | cut -d"/" -f1 | rev | cut -d"." -f2)
			echo $pheno $path >> $mainDir/code/9_Scores/ss_public_${cohort}
    done
	sh $mainDir/code/9_Scores/9.0_PGS.sh public $cohort $mainDir/derived_data/9_Scores/public
done