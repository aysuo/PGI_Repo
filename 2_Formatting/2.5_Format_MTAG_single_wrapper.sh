#!/bin/bash

source $PGI_Repo/code/2_Formatting/2.5.1_Format_MTAG.sh

cd $PGI_Repo/derived_data/4_MTAG_single

i=0
for pheno in $PGI_Repo/derived_data/4_MTAG_single/*; do
    format_MTAG_SBayesR $pheno 0 &
    let i+=1

    # Run 10 at a time
    if [[ $i == 10 ]]; then
        wait
        i=0
    fi   
done
