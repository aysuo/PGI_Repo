#!/bin/bash

# Get list of sumstats: Pheno name on first column (e.g. SWB-Okbay), file path on second
for cohort in Dunedin ERisk HRS2 WLS UKB3; 
do
    rm -f $PGI_Repo/code/9_Scores/ss_public_${cohort}
	for pheno in $(cat $PGI_Repo/code/9_Scores/version_public_$cohort); do 
			path=$(ls $PGI_Repo/derived_data/3_QCd/public_scores/SEfilter/QC_${pheno}_*/*.gz)
			pheno=$(echo $path | rev | cut -d"/" -f1 | rev | cut -d"." -f2)
			echo $pheno $path >> $PGI_Repo/code/9_Scores/ss_public_${cohort}
    done
	sh $PGI_Repo/code/9_Scores/9.0_PGS.sh public $cohort $PGI_Repo/derived_data/9_Scores/public
done