#!/bin/bash
source $PGI_Repo/code/paths

MTAG_single(){
  pheno=$1
  sumstats=$2
  out=$3

  echo " "
  echo "Running single-trait MTAG for $pheno.."
  echo " "
  ${python} ${mtag} \
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

checkStatusMTAG(){
    fileList=$1

    rm -f ${fileList}.error
    rm -f ${fileList}.rerun
    status=0

    while read row; do
      pheno=$(echo $row | cut -d" " -f1)
      
      if ! [[ $(ls ${pheno}/${pheno}_trait* 2>/dev/null) ]]; then
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
