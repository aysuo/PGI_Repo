#!/bin/bash
dirCode=/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/code
dirData=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data


echo -e "MultiMTAG_code\tInput_files\tWeights" > $dirData/6_MTAG_multi/MTAG_weights.txt
while read row
do
    pheno=$(echo $row | cut -d" " -f1)
    input_phenos=$(grep "^\-\-sumstats" $dirData/6_MTAG_multi/$pheno/$pheno.log | sed "s/--sumstats //g ; s,$dirData/4_MTAG_single/,,g" | awk -F"," '{for(i=1;i<=NF;i++) {split($i,a,"/"); printf a[1]","}; print ""}')  
    weights=$(grep -A1 "MTAG weight factors" $dirData/6_MTAG_multi/$pheno/$pheno.log | tail -1 | cut -d"[" -f2 | sed 's/\]// ; s/ \+/,/g')
    echo -e "$pheno\t$input_phenos\t$weights" >> $dirData/6_MTAG_multi/MTAG_weights.txt
done < $dirCode/6_MTAG_multi/mtag_groups_all_versions.txt

sed -i 's/\t,/\t/g ; s/,\t/\t/g' $dirData/6_MTAG_multi/MTAG_weights.txt