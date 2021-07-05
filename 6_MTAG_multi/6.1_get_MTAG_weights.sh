#!/bin/bash

echo -e "MultiMTAG_code\tInput_files\tWeights" > $PGI_Repo/derived_data/6_MTAG_multi/MTAG_weights.txt
while read row
do
    pheno=$(echo $row | cut -d" " -f1)
    input_phenos=$(grep "^\-\-sumstats" $PGI_Repo/derived_data/6_MTAG_multi/$pheno/$pheno.log | sed "s/--sumstats //g ; s,$PGI_Repo/derived_data/4_MTAG_single/,,g" | awk -F"," '{for(i=1;i<=NF;i++) {split($i,a,"/"); printf a[1]","}; print ""}')  
    weights=$(grep -A1 "MTAG weight factors" $PGI_Repo/derived_data/6_MTAG_multi/$pheno/$pheno.log | tail -1 | cut -d"[" -f2 | sed 's/\]// ; s/ \+/,/g')
    echo -e "$pheno\t$input_phenos\t$weights" >> $PGI_Repo/derived_data/6_MTAG_multi/MTAG_weights.txt
done < $PGI_Repo/code/6_MTAG_multi/mtag_groups_all_versions.txt

sed -i 's/\t,/\t/g ; s/,\t/\t/g' $PGI_Repo/derived_data/6_MTAG_multi/MTAG_weights.txt