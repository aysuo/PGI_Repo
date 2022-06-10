#!/bin/bash

format_MTAG_LDpred(){
    dirIn=$1
    multi=$2

    # Get phenotype name from path (e.g. dirIn=/XXX/YYY/NEURO1 -> pheno=NEURO1)
    pheno=$(echo $dirIn | rev | cut -d"/" -f1 | rev)
    echo "Formatting MTAG output for $pheno.."

    # If the input is single-trait MTAG output, will format all output files
    if [[ $multi == 0 ]]; then
        sumstats=$(ls $dirIn/*trait* | sed '/SBayesR/d')
    else
    # If it's multi-trait MTAG output, will format only the target phenotype 
        sumstats=$(ls $dirIn/*trait_1.txt)
    fi
    
    # Formatted output will be named XXX_formatted.txt
    out=$(echo $sumstats | sed 's/\.txt/_formatted\.txt/g')
    declare -a sumstats=$(echo "($sumstats)")
    declare -a out=$(echo "($out)")

    numSumstats=${#sumstats[@]}

    # Get total N by summing up the sample sizes across MTAG output files
    cut -f1,7  ${sumstats[0]} > $dirIn/sumN
    for (( i=1; i<$numSumstats; i++ )); do
        awk -F"\t" 'NR==FNR{N[$1]=$2;next}($1 in N){print $1,$7+N[$1]}' OFS="\t" $dirIn/sumN ${sumstats[$i]} > $dirIn/tmp
        mv $dirIn/tmp $dirIn/sumN
    done

    for (( i=0; i<$numSumstats; i++ )); do
        # Get GWAS-equivalent N from MTAG log file
        gwasN=$(grep -A $(($i+1)) "GWAS equiv. (max) N" $pheno/$pheno.log | tail -1  | awk '{print $NF}')
        # Format file, adding in total N and GWAS-eq N
        awk -F"\t" -v gwasN=$gwasN 'NR==FNR{N[$1]=$2;next} \
            FNR==1 { for (i=1; i<=NF; i++) { ix[$i] = i }; print "cptid","SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","EFFECT","SE","PVALUE","Ntot","N" }\
            FNR>1 {print $ix["CHR"]":"$ix["BP"],$ix["SNP"],$ix["CHR"],$ix["BP"],$ix["A1"],$ix["A2"],$ix["FRQ"],$ix["mtag_beta"],$ix["mtag_se"],$ix["mtag_pval"],N[$ix["SNP"]],gwasN}' OFS="\t" $dirIn/sumN ${sumstats[$i]} > ${out[$i]} &
    done 
    wait
    
    rm $dirIn/sumN
    echo "Formatting for $pheno finished."
}




    
# TO DO: Make sure LDSC accepts the SBayesR column names (may need to alter script 5.0 to specify col names with flags)
format_MTAG_SBayesR(){
    dirIn=$1
    multi=$2

    # Get phenotype name from path (e.g. dirIn=/XXX/YYY/NEURO1 -> pheno=NEURO1)
    pheno=$(echo $dirIn | rev | cut -d"/" -f1 | rev)
    echo "Formatting MTAG output for $pheno.. (SBayesR format)"

    # If the input is single-trait MTAG output, will format all output files
    if [[ $multi == 0 ]]; then
        sumstats=$(ls $dirIn/*trait* | sed '/formatted/d')
    else
    # If it's multi-trait MTAG output, will format only the target phenotype 
        sumstats=$(ls $dirIn/*trait_1.txt)
    fi
    
    # Formatted output will be named XXX_formatted_SBayesR.txt
    out=$(echo $sumstats | sed 's/\.txt/_formatted_SBayesR\.txt/g')
    declare -a sumstats=$(echo "($sumstats)")
    declare -a out=$(echo "($out)")

    numSumstats=${#sumstats[@]}
    Nsnps=$(($(wc -l ${sumstats[0]} | cut -d" " -f1)-1))
    medianSNP=$((($Nsnps/2)+1))

    for (( i=0; i<$numSumstats; i++ )); do
        # Get GWAS-equivalent N from MTAG log file
        gwasN=$(grep -A $(($i+1)) "GWAS equiv. (max) N" $dirIn/$pheno.log | tail -1  | awk '{print $NF}')
        
        # Get median effective N 
        awk -F"\t" 'NR==1 { for (i=1; i<=NF; i++) { ix[$i] = i } } \
            NR>1 {print 1/(2*$ix["mtag_se"]*$ix["mtag_se"]*$ix["FRQ"]*(1-$ix["FRQ"]))}' OFS="\t" ${sumstats[$i]} | sort -g > $dirIn/effN_$i
        medianN=$(awk -v M=${medianSNP} 'NR==M{print}' $dirIn/effN_$i )
        maxN=$(tail -1 $dirIn/effN_$i | cut -d" " -f2)
        minN=$(head -1 $dirIn/effN_$i | cut -d" " -f2)
        
        # Format file, adding in total N and GWAS-eq N
        awk -F"\t" -v gwasN=$gwasN -v medianN=$medianN 'NR==1 { for (i=1; i<=NF; i++) { ix[$i] = i }; print "SNP","A1","A2","freq","b","se","p","N" }\
            NR>1 {effN=1/(2*$ix["mtag_se"]*$ix["mtag_se"]*$ix["FRQ"]*(1-$ix["FRQ"]))}\
            NR>1 && effN>0.8*medianN {print $ix["SNP"],$ix["A1"],$ix["A2"],$ix["FRQ"],$ix["mtag_beta"],$ix["mtag_se"],$ix["mtag_pval"],gwasN}' ${sumstats[$i]} > ${out[$i]} 
        NsnpsOut=$(wc -l ${out[$i]} | cut -d" " -f1)
     
        echo "$pheno - $i: median effN = $medianN, min effN = $minN, max effN = $maxN, GWAS-eq N = $gwasN, $(($Nsnps-$NsnpsOut)) SNPs dropped." >> $dirIn/../SBayesR_Nfiltering.txt
    done
    rm $dirIn/effN*
    echo "Formatting for $pheno finished."
}