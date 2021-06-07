#!/bin/bash

source paths2
source $mainDir/code/2_Formatting/2.5.1_Format_MTAG.sh

cd $mainDir/derived_data/6_MTAG_multi

i=0
for pheno in $mainDir/derived_data/6_MTAG_multi/*; do
    format_MTAG $pheno 1 &
    let i+=1

    # Run 10 at a time
    if [[ $i == 10 ]]; then
        wait
        i=0
    fi
done
