#!/bin/bash

source paths4 
source $mainDir/code/4_MTAG_single/4.0_MTAG_single_functions.sh

cd $mainDir/derived_data/4_MTAG_single

# Correct reverse-coded sumstats
for study in SWB-23andMe Risk-23andMe; do
  path=$(cut -f2 $mainDir/code/4_MTAG_single/singleMTAG_input_filelist.txt | awk -F"," -v study=$study '{for(i=1;i<=NF;i++) if ($i~study) print $i}'  | sort | uniq)
  eval path=$path
  # Unzip and rename original file with as *_revcoded
  unzipped=$(echo $path | sed 's/\.gz//g')
  if [[ $path == *.gz ]]; then
    gunzip $path
  fi 
  mv $unzipped ${unzipped}_revcoded
  # Reverse sign of effect 
  awk -F"\t" 'NR==1{print}NR>1{$8=-$8;print}' OFS="\t" ${unzipped}_revcoded > ${unzipped}
  gzip ${unzipped}
done


# Check which phenotypes aren't done
checkStatusMTAG $mainDir/code/4_MTAG_single/singleMTAG_input_filelist.txt

# Run MTAG for unfinished phenotypes, 20 at a time
if [[ $status == 1 ]]; then
  i=0
  while read row; do
    pheno=$(echo $row | cut -d" " -f1)
    sumstats=$(echo $row | cut -d" " -f2)
    eval sumstats=$sumstats

    mkdir -p $pheno

    MTAG_single $pheno $sumstats $pheno/$pheno &
    let i+=1

    if [[ $i == 20 ]]; then
      wait
      i=0
    fi
  done < $mainDir/code/4_MTAG_single/singleMTAG_input_filelist.txt.rerun
  wait
fi

# Check if everything ran successfully
checkStatusMTAG $mainDir/code/4_MTAG_single/singleMTAG_input_filelist.txt

