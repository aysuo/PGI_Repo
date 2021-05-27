#!/bin/bash

source paths2
source $p2_code/2.5.1_Format_MTAG.sh

cd $p2_MTAGsingle

i=0
for pheno in $p2_MTAGsingle/*; do
    format_MTAG $pheno 0 &
    let i+=1

    # Run 10 at a time
    if [[ $i == 10 ]]; then
        wait
        i=0
    fi   
done
