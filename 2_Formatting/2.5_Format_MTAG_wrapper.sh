#!/bin/bash

source $PGI_Repo/code/2_Formatting/2.5.1_Format_MTAG.sh
source $PGI_Repo/code/paths
cd $PGI_Repo/derived_data/4_MTAG_single

#phenos_skip_Nfilter=()

j=0
for pheno in $PGI_Repo/derived_data/4_MTAG_single/*; do
    pheno_name=$(basename "$pheno")
    if [[ " ${phenos_skip_Nfilter[*]} " =~ " $pheno_name " ]]
    then
        format_MTAG_SBayesR "$pheno" 0 1 &
    else
        format_MTAG_SBayesR "$pheno" 0 0 &
    fi
    let j+=1

    # Run 10 at a time
    if [[ $j == 10 ]]; then
        wait
        j=0
    fi   
done


# ## PRCA, ALZ1,3,5: Convert cptid to rsid
for pheno in $PGI_Repo/derived_data/4_MTAG_single/ALZ1 $PGI_Repo/derived_data/4_MTAG_single/ALZ5 $PGI_Repo/derived_data/4_MTAG_single/ALZ3 $PGI_Repo/derived_data/4_MTAG_single/PRCA* 
do
  for file in $pheno/*_trait_*_formatted_SBayesR.txt
   do
       mv $file ${file}_ChrPosID
       file=$(echo $file | sed 's/_ChrPosID//g' )
       awk -F"\t" 'NR==FNR{a[$2]=$1;next}FNR==1{FS=" ";print;next}($1 in a){FS=" "; $1=a[$1] ; print}' $HRC_rsid2chrpos_map ${file}_ChrPosID > ${file} &
   done
   wait
done
