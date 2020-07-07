#!/bin/bash

dirCode="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/code/9_Scores"
dirIn="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/4_MTAG_single"
dirOut="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/9_Scores/single"

# Get list of sumstats: Pheno name on first column (e.g. SWB-Okbay), file path on second
for cohort in HRS2 WLS; do
    rm -f $dirCode/ss_single_${cohort}
	for pheno in $(cat $dirCode/MTAGversion_single_$cohort); do
        if [[ -f $dirIn/$pheno/${pheno}_trait_formatted.txt ]]; then
		    path="$dirIn/$pheno/${pheno}_trait_formatted.txt"
        else 
            path="$dirIn/$pheno/${pheno}_trait_1_formatted.txt"
        fi
		pheno=$(echo "$pheno-single" | sed 's/[1-9]//g')
		echo $pheno $path >> $dirCode/ss_single_${cohort}
    done
done

#sh $dirCode/9.0_PGS.sh single WLS $dirOut &
sh $dirCode/9.0_PGS.sh single HRS2 $dirOut &
wait
