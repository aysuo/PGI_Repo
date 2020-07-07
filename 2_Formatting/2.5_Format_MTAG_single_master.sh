#!/bin/bash

dirIn="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/4_MTAG_single"
dirCode="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code"

source $dirCode/2_Formatting/2.5.1_Format_MTAG_slave.sh

cd $dirIn

i=0
for pheno in $dirIn/*; do
    format_MTAG $pheno 0 &
    let i+=1

    if [[ $i == 10 ]]; then
        wait
        i=0
    fi   
done
