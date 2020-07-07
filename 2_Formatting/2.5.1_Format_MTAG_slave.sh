#!/bin/bash

format_MTAG(){
    dirIn=$1
    multi=$2

    pheno=$(echo $dirIn | rev | cut -d"/" -f1 | rev)
    echo "Formatting MTAG output for $pheno.."

    if [[ $multi == 0 ]]; then
        sumstats=$(ls $dirIn/*trait*)
    else
        sumstats=$(ls $dirIn/*trait_1.txt)
    fi
    
    out=$(echo $sumstats | sed 's/\.txt/_formatted\.txt/g')
    declare -a sumstats=$(echo "($sumstats)")
    declare -a out=$(echo "($out)")

    numSumstats=${#sumstats[@]}

    cut -f1,7  ${sumstats[0]} > $dirIn/sumN
    for (( i=1; i<$numSumstats; i++ )); do
        awk -F"\t" 'NR==FNR{N[$1]=$2;next}($1 in N){print $1,$7+N[$1]}' OFS="\t" $dirIn/sumN ${sumstats[$i]} > $dirIn/tmp
        mv $dirIn/tmp $dirIn/sumN
    done

    for (( i=0; i<$numSumstats; i++ )); do
        out=$(echo ${sumstats[$i]} | sed 's/\.txt/_formatted\.txt/g')
        gwasN=$(grep -A $(($i+1)) "GWAS equiv. (max) N" $pheno/$pheno.log | tail -1  | awk '{print $NF}')
        echo $gwasN
        awk -F"\t" -v gwasN=$gwasN 'NR==FNR{N[$1]=$2;next} \
            FNR==1 { for (i=1; i<=NF; i++) { ix[$i] = i }; print "cptid","SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","EFFECT","SE","PVALUE","Ntot","N" }\
            FNR>1 {print $ix["CHR"]":"$ix["BP"],$ix["SNP"],$ix["CHR"],$ix["BP"],$ix["A1"],$ix["A2"],$ix["FRQ"],$ix["mtag_beta"],$ix["mtag_se"],$ix["mtag_pval"],N[$ix["SNP"]],gwasN}' OFS="\t" $dirIn/sumN ${sumstats[$i]} > ${out[$i]} &
    done 
    wait
    
    rm $dirIn/sumN
    echo "Formatting for $pheno finished."
}


