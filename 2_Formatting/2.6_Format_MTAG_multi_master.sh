#!/bin/bash

source paths2
source $p2_code/2.5.1_Format_MTAG.sh

cd $p2_MTAGmulti

i=0
for pheno in $p2_MTAGmulti/*; do
    format_MTAG $pheno 1 &
    let i+=1

    # Run 10 at a time
    if [[ $i == 10 ]]; then
        wait
        i=0
    fi
done
