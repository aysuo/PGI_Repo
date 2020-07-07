format_MTAG(){
    pheno=$1
    ss=$2
    freq=$3

    sumstats=$(echo $ss | sed 's/,/ /g')
    sumstats_unzipped=$(echo $sumstats | sed 's/\.gz//g')
    declare -a sumstats=$(echo "($sumstats)")
    declare -a sumstats_unzipped=$(echo "($sumstats_unzipped)")
    numSumstats=${#sumstats[@]}
  
    #for (( i=0; i<$numSumstats; i++ )); do
    #    if [[ ${sumstats[$i]} == *.gz ]]; then
    #        gunzip ${sumstats[$i]}
    #    fi
    #done

    cut -f2,11  ${sumstats_unzipped[0]} > $dirOut/$pheno/sumN
    for (( i=1; i<$numSumstats; i++ )); do
        awk -F"\t" 'NR==FNR{N[$1]=$2;next}($2 in N){print $2,$11+N[$2]}' OFS="\t" $dirOut/$pheno/sumN ${sumstats_unzipped[$i]} > $dirOut/$pheno/tmp
        mv $dirOut/$pheno/tmp $dirOut/$pheno/sumN
    done
    
    #for (( i=0; i<$(($numSumstats-1)); i++ )); do
    #    if [[ ${sumstats[$i]} == *.gz ]]; then
    #        gzip ${sumstats_unzipped[$i]} &
    #    fi
    #done

    awk -F"\t" -v freq=$freq 'NR==FNR{N[$1]=$2;next} \
    FNR==1 { for (i=1; i<=NF; i++) { ix[$i] = i }; print "cptid","SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","EFFECT","SE","PVALUE","N","Neff" }\
    FNR>1 && $ix[freq]>0 && $ix[freq]<1{print $ix["CHR"]":"$ix["BP"],$ix["SNP"],$ix["CHR"],$ix["BP"],$ix["A1"],$ix["A2"],$ix[freq],$ix["mtag_beta"],$ix["mtag_se"],$ix["mtag_pval"],N[$ix["SNP"]],1/(($ix["mtag_se"]^2)*2*$ix[freq]*(1-$ix[freq]))}' OFS="\t" $dirOut/$pheno/sumN $dirOut/$pheno/${pheno}_mtag_meta.txt > $dirOut/$pheno/${pheno}_mtag_meta_formatted.txt

  rm $dirOut/$pheno/sumN
  wait
}
