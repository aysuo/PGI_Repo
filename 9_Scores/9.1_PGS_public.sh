#!/bin/bash

dirCode="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code/9_Scores"
dirIn="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/3_QCd/public_scores/SEfilter"
dirOut="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/9_Scores/public"

cd /disk/genetics/PGS/Aysu/PGS_Repo_pipeline/code/9_Scores
# Get list of sumstats: Pheno name on first column (e.g. SWB-Okbay), file path on second
for cohort in HRS2 WLS; do
    rm -f $dirCode/ss_public_${cohort}
	#for pheno in $(cat $dirCode/version_public_$cohort); do 
	for pheno in COPD-Neale NEBmen-Neale FINSAT-Neale WORKSAT-Neale LONELY-Neale; do
			path=$(ls $dirIn/QC_${pheno}_*/*.gz)
			pheno=$(echo $path | rev | cut -d"/" -f1 | rev | cut -d"." -f2)
			echo $pheno $path >> $dirCode/ss_public_${cohort}
    done
done

#sh $dirCode/9.0_PGS.sh public WLS $dirOut
sh $dirCode/9.0_PGS.sh public HRS2 $dirOut