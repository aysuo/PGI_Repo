#!/bin/bash

dirIn="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/"
dirOut="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/4_MTAG_single"
dirCode="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code/4_MTAG_single"
mtag="/disk/genetics/ukb/aokbay/bin/mtag/mtag.py"

cd $dirOut

# Correct reverse-coded sumstats
#for study in SWB-23andMe Risk-23andMe; do
#  path=$(cut -f2 $dirCode/singleMTAG_filelist.txt | awk -F"," -v study=$study '{for(i=1;i<=NF;i++) if ($i~study) print $i}'  | sort | uniq)
  
  # Unzip and rename original file with as *_revcoded
#  unzipped=$(echo $path | sed 's/\.gz//g')
#  if [[ $path == *.gz ]]; then
#    gunzip $dirIn/$path
#  fi 
#  mv $dirIn/$unzipped $dirIn/${unzipped}_revcoded
   # Reverse sign of effect 
#  awk -F"\t" 'NR==1{print}NR>1{$8=-$8;print}' OFS="\t" $dirIn/${unzipped}_revcoded > $dirIn/${unzipped}
#  gzip $dirIn/${unzipped}
#done


MTAG_single(){
  pheno=$1
  sumstats=$2
  out=$3

  echo " "
  echo "Running single-trait MTAG for $pheno.."
  echo " "
  python2.7 ${mtag} \
    --sumstats ${sumstats} \
    --out ${out} \
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
    --perfect_gencov 

    echo " "
    echo "Single-trait MTAG for $pheno finished."
    echo " "
}

checkStatus(){
    fileList=$1

    rm -f ${fileList}.error
    rm -f ${fileList}.rerun
    status=0

    while read row; do
      pheno=$(echo $row | cut -d" " -f1)
      
      if ! [[ $(ls ${dirOut}/${pheno}/${pheno}_trait* 2>/dev/null) ]]; then
        echo $pheno >> ${fileList}.error
        grep $pheno $fileList >> ${fileList}.rerun
        status=1
      fi
    done < $fileList
    
    if ! [[ -f ${fileList}.rerun ]]; then
        echo "Single-trait MTAG completed."
    else
        echo "Single-trait MTAG was not completed for some phenotypes:"
        cat ${fileList}.error
        echo ""
        echo "Errors are stored in ${fileList}.error"
    fi
}

main(){

  checkStatus $dirCode/singleMTAG_input_filelist.txt

  if [[ $status == 1 ]]; then
    i=0
    while read row; do
      pheno=$(echo $row | cut -d" " -f1)
      sumstats=$(echo $row | cut -d" " -f2)
      eval sumstats=$sumstats

      mkdir -p ${dirOut}/${pheno}

      MTAG_single $pheno $sumstats ${dirOut}/${pheno}/${pheno} &
      let i+=1

      if [[ $i == 20 ]]; then
        wait
        i=0
      fi
    done < $dirCode/singleMTAG_input_filelist.txt.rerun
    wait
  fi

  checkStatus $dirCode/singleMTAG_input_filelist.txt
}

main