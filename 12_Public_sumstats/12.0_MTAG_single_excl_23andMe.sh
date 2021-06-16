#!/bin/bash

source paths12
source $mainDir/code/4_MTAG_single/4.0_MTAG_single_functions.sh
source $mainDir/code/2_Formatting/2.5.1_Format_MTAG.sh

cd $mainDir/derived_data/12_Public_sumstats/single

checkStatusMTAG $mainDir/code/12_Public_sumstats/singleMTAG_excl_23andMe_input_filelist.txt

if [[ $status == 1 ]]; then
  while read row; do
    pheno=$(echo $row | cut -d" " -f1)
    sumstats=$(echo $row | cut -d" " -f2)
    eval sumstats=$sumstats

    mkdir -p ${pheno}

    MTAG_single $pheno $sumstats ${pheno}/${pheno} &
  done < $mainDir/code/12_Public_sumstats/singleMTAG_excl_23andMe_input_filelist.txt.rerun
  wait
fi

checkStatusMTAG $mainDir/code/12_Public_sumstats/singleMTAG_excl_23andMe_input_filelist.txt

while read row; do
  pheno=$(echo $row | cut -d" " -f1)
  format_MTAG $pheno 0 &
done < $mainDir/code/12_Public_sumstats/singleMTAG_excl_23andMe_input_filelist.txt
wait


