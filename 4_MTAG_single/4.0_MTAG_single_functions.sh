#!/bin/bash
source $PGI_Repo/code/paths

MTAG_single(){
  pheno=$1
  sumstats=$2
  out=$3
  force=$4
  nooverlap=$5
  metaformat=$6
  cptid=$7

  if [[ $force == 0 ]]
  then
    force=""
  else  
    force="--force"
  fi

  if [[ $nooverlap == 0 ]]
  then
    nooverlap=""
  else
    nooverlap="--no_overlap"
  fi

  if [[ $metaformat == 0 ]]
  then
    metaformat=""
  else
    metaformat="--meta_format"
  fi

  if [[ $cptid == 0 ]]
  then
    snpname="SNPID"
    ldrefpanel=""
  else
    snpname="cptid"
    ldrefpanel="--ld_ref_panel=${MTAGfolder}/ld_ref_panel/eur_w_ld_chr_cptid/"
  fi

  echo " "
  echo "Running single-trait MTAG for $pheno.."
  echo " "
  ${python} ${MTAGfolder}/mtag.py $nooverlap $force $metaformat $ldrefpanel \
    --sumstats ${sumstats} \
    --out ${out} \
    --snp_name ${snpname} \
    --chr_name CHR \
    --bpos_name POS \
    --a1_name EFFECT_ALLELE \
    --a2_name OTHER_ALLELE \
    --use_beta_se \
    --beta_name EFFECT \
    --se_name SE \
    --eaf_name EAF \
    --n_name N \
    --maf_min 0 \
    --perfect_gencov 

  # NEARSIGHTED6 cannot be analyzed by MTAG (omega not positive semi-definite). Just reformat the summary statistics to match MTAG output format.
  if [[ $pheno == "NEARSIGHTED6" ]]
  then
    awk -F"\t" 'BEGIN{print "SNP","CHR","BP","A1","A2","Z","N","FRQ","mtag_beta","mtag_se","mtag_z","mtag_pval"} \
      NR>1{print $2,$3,$4,$5,$6,$8/$9,$11,$7,$8,$9,$8/$9,$10}' ${sumstats} > ${pheno}_trait.txt 
  fi

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

