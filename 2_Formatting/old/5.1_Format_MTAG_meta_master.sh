#!/bin/bash

dirIn="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/3_QCd"
dirOut="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/3.5_MTAG_meta"
dirCode="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code"

source $dirCode/2_Formatting/5.0_Format_MTAG_slave.sh

cd $dirOut

while read row; do
    pheno=$(echo $row | cut -d" " -f1)
    eval sumstats=$(echo $row | cut -d" " -f2 | sed 's_^_${dirIn}/_g' | sed 's_,_,${dirIn}/_g')

    if ! [[ $(find $pheno/${pheno}_mtag_meta_formatted.txt -type f -size +100000 2>/dev/null) ]]; then
        echo "Formatting MTAG output for $pheno.."
        format_MTAG $pheno $sumstats meta_freq &
    fi
done < $dirCode/3.5_MTAG_meta/MTAG_meta_filelist.txt
wait
echo "Formatting finished."