#!/bin/bash

dirIn="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/4_MTAG_single"
dirOut="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/6_MTAG_multi"
dirCode="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code/6_MTAG_multi"
mtag="/disk/genetics/ukb/aokbay/bin/mtag/mtag.py"
rgmatrix="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/5_LDSC/singleMTAG/rg_table.txt"

cd $dirOut

MTAG_multi(){
  pheno=$1
  sumstats=$2

  echo " "
  echo "Running multi-trait MTAG for $pheno.."
  echo " "
  python2.7 ${mtag} \
    --sumstats ${sumstats} \
    --out ${pheno}/${pheno} \
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
    --maf_min 0

    echo " "
    echo "Multi-trait MTAG for $pheno finished."
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
        awk -F"\t" -v pheno=$pheno '$1==pheno{print}' OFS="\t" $fileList >> ${fileList}.rerun
        status=1
      fi
    done < $fileList
    
    if ! [[ -f ${fileList}.rerun ]]; then
        echo "Multi-trait MTAG completed."
    else
        echo "Multi-trait MTAG was not completed for some phenotypes:"
        cat ${fileList}.error
        echo ""
        echo "Errors are stored in ${fileList}.error"
    fi
}

main(){
    awk -F"\t" 'NR==1{for (i=1;i<=NF;i++) pheno[i]=$i;next}\
        {group=""; for (i=2;i<=NF;i++) if ($i>=0.6 || $i<=-0.6) group=pheno[i]","group} \
        group != "" {print $1,group}' OFS="\t" $rgmatrix | sed 's/,$//g' > $dirCode/mtag_groups.txt

    checkStatus $dirCode/mtag_groups_all_versions.txt

    if [[ $status == 1 ]]; then
        j=0
        while read row; do
            pheno=$(echo $row | cut -d" " -f1)
            group=$(echo $row | cut -d" " -f2 | sed 's/,/ /g')

            declare -a group="($group)"
            
            for ((i=0;i<${#group[@]};i++)); do
                if [[ -f $dirIn/${group[$i]}/${group[$i]}_trait_1_formatted.txt ]]; then
                    sumstats[$i]=$dirIn/${group[$i]}/${group[$i]}_trait_1_formatted.txt
                else
                    sumstats[$i]=$dirIn/${group[$i]}/${group[$i]}_trait_formatted.txt
                fi
            done

            ss=$(echo ${sumstats[@]} | sed 's/ /,/g')
            mkdir -p ${dirOut}/${pheno}

            MTAG_multi $pheno $ss &
            let j+=1

            if [[ $j == 10 ]]; then
                wait
                j=0
            fi

            unset sumstats
        done < $dirCode/mtag_groups_all_versions.txt.rerun
        wait
    fi

  checkStatus $dirCode/mtag_groups_all_versions.txt
}

main