#!/bin/bash

source $PGI_Repo/code/paths
cd $PGI_Repo/derived_data/6_MTAG_multi

MTAG_multi(){
  pheno=$1
  sumstats=$2

  echo " "
  echo "Running multi-trait MTAG for $pheno.."
  echo " "
  ${python} ${mtag} \
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

checkStatusMultiMTAG(){
    fileList=$1

    rm -f ${fileList}.error
    rm -f ${fileList}.rerun
    status=0

    while read row; do
      pheno=$(echo $row | cut -d" " -f1)
      
      if ! [[ $(ls $PGI_Repo/derived_data/6_MTAG_multi/${pheno}/${pheno}_trait* 2>/dev/null) ]]; then
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
    # Get list of supplementary phenotypes for each phenotype using the rg table
    awk -F"\t" 'NR==1{for (i=1;i<=NF;i++) pheno[i]=$i;next}\
        {group=""; for (i=2;i<=NF;i++) if ($i>=0.6 || $i<=-0.6) group=pheno[i]","group} \
        group != "" {print $1,group}' OFS="\t" $PGI_Repo/derived_data/5_LDSC/singleMTAG/rg_table.txt | sed 's/,$//g' > $PGI_Repo/code/6_MTAG_multi/mtag_groups.txt

    # Have to get input files for different versions manually
    # Do that and write into "$PGI_Repo/code/6_MTAG_multi/mtag_groups_all_versions.txt"
    checkStatusMultiMTAG $PGI_Repo/code/6_MTAG_multi/mtag_groups_all_versions.txt

    if [[ $status == 1 ]]; then
        j=0
        while read row; do
            pheno=$(echo $row | cut -d" " -f1)
            group=$(echo $row | cut -d" " -f2 | sed 's/,/ /g')

            declare -a group="($group)"
            
            for ((i=0;i<${#group[@]};i++)); do
                if [[ -f $PGI_Repo/derived_data/4_MTAG_single/${group[$i]}/${group[$i]}_trait_1_formatted.txt ]]; then
                    sumstats[$i]=$PGI_Repo/derived_data/4_MTAG_single/${group[$i]}/${group[$i]}_trait_1_formatted.txt
                else
                    sumstats[$i]=$PGI_Repo/derived_data/4_MTAG_single/${group[$i]}/${group[$i]}_trait_formatted.txt
                fi
            done

            ss=$(echo ${sumstats[@]} | sed 's/ /,/g')
            mkdir -p $PGI_Repo/derived_data/6_MTAG_multi/${pheno}

            MTAG_multi $pheno $ss &
            let j+=1

            if [[ $j == 10 ]]; then
                wait
                j=0
            fi

            unset sumstats
        done < $PGI_Repo/code/6_MTAG_multi/mtag_groups_all_versions.txt.rerun
        wait
    fi

  checkStatusMultiMTAG $PGI_Repo/code/6_MTAG_multi/mtag_groups_all_versions.txt
}

main