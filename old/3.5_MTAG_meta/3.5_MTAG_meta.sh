#!/bin/bash

dirIn="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/3_QCd"
dirOut="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/3.5_MTAG_meta"
dirCode="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code/"
mtag="/disk/genetics/ukb/aokbay/bin/mtag/mtag.py"

cd $dirOut

format_MTAG_meta(){
  pheno=$1
  ss=$2

  sumstats=$(echo $ss | sed 's/,/ /g')
  sumstats_unzipped=$(echo $sumstats | sed 's/\.gz//g')
  declare -a sumstats=$(echo "($sumstats)")
  declare -a sumstats_unzipped=$(echo "($sumstats_unzipped)")
  numSumstats=${#sumstats[@]}
  
  gunzip ${sumstats[0]}
  cut -f2,11  ${sumstats_unzipped[0]} > $dirOut/$pheno/sumN
  for (( i=1; i<$numSumstats; i++ )); do
    gunzip ${sumstats[$i]}
    awk -F"\t" 'NR==FNR{N[$1]=$2;next}($2 in N){print $2,$11+N[$2]}' OFS="\t" $dirOut/$pheno/sumN ${sumstats_unzipped[$i]} > $dirOut/$pheno/tmp
    mv $dirOut/$pheno/tmp $dirOut/$pheno/sumN
    gzip ${sumstats_unzipped[$i]} &
  done
  gzip ${sumstats_unzipped[0]}

  awk -F"\t" 'NR==FNR{N[$1]=$2;next}\
    FNR==1{print "cptid","SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","EFFECT","SE","PVALUE","N","IMPUTED","INFO","HWE_PVAL","CALLRATE"} \
    FNR>1{print $2":"$3,$1,$2,$3,$4,$5,$6,$7,$8,$9,N[$1],"NA","NA","NA","NA"}' OFS="\t" $dirOut/$pheno/sumN $dirOut/$pheno/${pheno}_mtag_meta.txt > $dirOut/$pheno/${pheno}_mtag_meta_formatted.txt

  rm $dirOut/$pheno/sumN
}


MTAG_meta(){
  pheno=$1
  ss=$2
  mkdir -p ${dirOut}/${pheno}

  echo " "
  echo "Running inverse variance weighted meta-analysis for $pheno.."
  echo " "
  python2.7 ${mtag} \
    --sumstats ${ss} \
    --out ${dirOut}/${pheno}/${pheno} \
    --snp_name SNPID \
    --chr_name CHR \
    --bpos_name POS \
    --a1_name EFFECT_ALLELE \
    --a2_name OTHER_ALLELE \
    --eaf_name EAF \
    --use_beta_se \
    --beta_name EFFECT \
    --se_name SE \
    --n_name N \
    --maf_min 0 \
    --perfect_gencov \
    --equal_h2 \
    --force 

    echo " "
    echo "Inverse variance weighted meta-analysis for $pheno finished."
    
}

checkStatus(){
  fileList=$1
  
  status=0
  rm -f $dirCode/3.5_MTAG_meta/MTAG_meta_rerun.txt
  while read row; do
    pheno=$(echo $row | cut -d" " -f1)
  
    if [[ $(find $dirOut/$pheno/${pheno}_mtag_meta.txt -type f -size +100000 2>/dev/null) ]]; then
      echo "MTAG meta-analysis for $pheno was successful, skipping.."
    else
      echo "MTAG meta-analysus for $pheno was not run before or has failed."
      grep $pheno $fileList >> $dirCode/3.5_MTAG_meta/MTAG_meta_rerun.txt
      status=1
    fi
  done < $fileList
}



main(){

  checkStatus $dirCode/3.5_MTAG_meta/MTAG_meta_filelist.txt

  if [[ $status == 1 ]]; then
    while read row; do
      pheno=$(echo $row | cut -d" " -f1)
      eval sumstats=$(echo $row | cut -d" " -f2 | sed 's_^_${dirIn}/_g' | sed 's_,_,${dirIn}/_g')
    
      MTAG_meta $pheno $sumstats &
  
     done < $dirCode/3.5_MTAG_meta/MTAG_meta_rerun.txt
  fi
  wait

  sh $dirCode/2_Formatting/5.1_Format_MTAG_meta_master.sh
}

main