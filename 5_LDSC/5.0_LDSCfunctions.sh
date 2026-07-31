#!/bin/bash

source $PGI_Repo/code/paths

ldsc_munge () {
    ss=$1
    snplist=$2
    out=$3

    eval snps='$'${snplist}_rsID

	## Munge sumstats ##
	$python ${LDSC}/munge_sumstats.py \
		--sumstats ${ss} \
		--out ${out} \
		--merge-alleles ${snps}
}

ldsc_h2 () {
    pheno=$1
    ssMunged=$2
    snplist=$3
    out=$4

    eval LDscores='$'${snplist}_LDscores

    phenos_highchi2=("BL_CHOL","BL_LDL","BL_HDL","BL_nonHDL")
	pheno=$(echo $pheno | cut -d"-" -f1)
    
    ## Get h2 and intercept ##
    if [[ ${phenos_highchi2[@]} =~ $pheno  ]]
    then
        chi2cutoff=1000
    else
        chi2cutoff=inf
    fi

	$python ${LDSC}/ldsc.py \
		--h2 ${ssMunged}  \
		--ref-ld-chr ${LDscores} \
		--w-ld-chr ${LDscores} \
        --two-step ${chi2cutoff} \
		--out ${out}
}


ldsc_rg () {
    ssPairMunged=$1
    snplist=$2
    out=$3

    eval LDscores='$'${snplist}_LDscores

    phenos_highchi2=("BL_CHOL","BL_LDL","BL_HDL","BL_nonHDL")
	pheno1=$(echo ${ssPairMunged} | cut -d"," -f1 | cut -d"/" -f1)
    pheno2=$(echo ${ssPairMunged} | cut -d"," -f2 | cut -d"/" -f1)
    
    if [[ ${phenos_highchi2[@]} =~ $pheno1 ]] || [[ ${phenos_highchi2[@]} =~ $pheno2 ]]
    then
        chi2cutoff=1000
    else
        chi2cutoff=inf
    fi

    ## Get rg 
	$python ${LDSC}/ldsc.py \
	    --rg ${ssPairMunged}  \
	    --ref-ld-chr ${LDscores} \
	    --w-ld-chr ${LDscores} \
        --two-step ${chi2cutoff} \
	    --out ${out}
}

###############################################################################

munge () {
    fileListMunge=$1
    snplist=$2

    i=1
    while read row; do
        pheno=$(echo ${row} | cut -d" " -f1)
        ss=$(echo ${row} | cut -d" " -f2)
        eval ss=$ss

        if [[ $ss == *.gz ]]
        then 
            gunzip $ss
            ss=$(echo $ss | sed 's/\.gz//g')
        fi

        GWAS=$(echo $ss | rev | cut -d"/" -f1 | rev | sed 's/CLEANED\.//g' | sed 's/_formatted//g' | sed 's/_SBayesR//g' | sed 's/\.txt//g')

        mkdir -p $pheno/munged
        echo "Munging sumstats for ${GWAS}.."
        ldsc_munge $ss ${snplist} ${pheno}/munged/${GWAS}  &
        
        let i+=1

        if [[ $i == 10 ]]; then
            wait
            i=0
        fi
    
    done  < $fileListMunge
    wait
}


h2 () {
    fileListh2=$1
    snplist=$2

    i=1
    while read row; do
        pheno=$(echo ${row} | cut -d" " -f1)
        GWAS=$(echo ${row} | cut -d" " -f2)

        mkdir -p ${pheno}/h2
        echo "Estimating h2 for ${GWAS}.."
        ldsc_h2 ${pheno} ${pheno}/munged/${GWAS}.sumstats.gz ${snplist} ${pheno}/h2/h2_${GWAS} &
        
        let i+=1

        if [[ $i == 10 ]]; then
            wait
            i=0
        fi
    done < $fileListh2
    wait
}


rg () { 
    fileListrg=$1
    snplist=$2

    i=1
    while read row; do
        pheno1=$(echo ${row} | cut -d" " -f1 | cut -d"," -f1)
        pheno2=$(echo ${row} | cut -d" " -f1 | cut -d"," -f2)
        GWAS1=$(echo ${row} | cut -d" " -f2 | cut -d"," -f1)
        GWAS2=$(echo ${row} | cut -d" " -f2 | cut -d"," -f2)
        mkdir -p ${pheno1}/rg

        echo "Estimating rg between ${GWAS1} and ${GWAS2}.."
        ldsc_rg ${pheno1}/munged/${GWAS1}.sumstats.gz,${pheno2}/munged/${GWAS2}.sumstats.gz \
            ${snplist} \
            ${pheno1}/rg/rg_${GWAS1}_${GWAS2} &
        let i+=1

        if [[ $i == 10 ]]; then
            wait
            i=0
        fi
    done < $fileListrg
    wait
}

###############################################################################
# Check for which phenotypes in file list munging/h2/rg hasn't finished
checkStatusLDSC () {
    fileList=$1
    analysis=$2

    echo "Checking status of ${analysis} for ${fileList}.."

    rm -f ${fileList}.${analysis}.error
    rm -f ${fileList}.${analysis}.rerun
    
    phenoList=$(cut -f1 $fileList| sed 's/[0-9]$//g' | sed 's/[0-9]$//g' |  sort | uniq)
        
    for pheno in $phenoList; do    
        sumstats=$(grep ^$pheno[1-9] $fileList | cut -f2 | sed 's/,/\n/g' | sort | uniq)
        GWASs=$(echo $sumstats | awk '{for (i=1;i<=NF;i++) {N=split($i,a,"/"); print a[N]}}' | sed 's/CLEANED\.//g'|  sed 's/_formatted//g' | sed 's/_SBayesR//g' | sed 's/\.txt//g' | sed 's/\.gz//g' )

        for GWAS in $GWASs; do
            case $analysis in 
                munge)
                    if ! [[ -f ${pheno}/munged/${GWAS}.sumstats.gz ]] || [[ $(grep "ERROR converting summary statistics:" ${pheno}/munged/${GWAS}.log) ]]; then
                        echo "Munging for $GWAS was either not run before or it failed.." >> ${fileList}.munge.error
                        ss=$(echo $sumstats | awk -v GWAS=$GWAS '{for(i=1;i<=NF;i++) if ($i~GWAS"_") print $i}')
                        echo -e "$pheno\t$ss" >> ${fileList}.munge.rerun
                    fi
                    ;;
                h2)
                    if ! [[ -f ${pheno}/h2/h2_${GWAS}.log ]] || ! [[ $(grep "Total Observed scale h2" ${pheno}/h2/h2_${GWAS}.log) ]]; then
                        echo "LDSC h2 estimation for $GWAS was either not run before or it failed.." >> ${fileList}.h2.error
                        echo -e "$pheno\t$GWAS" >> ${fileList}.h2.rerun
                    fi
                    ;;
                rg) 
                    GWASs2=$(echo $GWASs | cut --complement -d" " -f1)

                    for GWAS2 in $GWASs2; do
                        if ! [[ $GWAS == $GWAS2 ]]; then
                            if ! [[ -f ${pheno}/rg/rg_${GWAS}_${GWAS2}.log ]] || ! [[ $(grep "Summary of Genetic Correlation Results" ${pheno}/rg/rg_${GWAS}_${GWAS2}.log) ]]; then
                                echo "LDSC rg estimation between $GWAS and $GWAS2 was either not run before or it failed.." >> ${fileList}.rg.error
                                echo -e "$pheno,$pheno\t$GWAS,$GWAS2" >> ${fileList}.rg.rerun
                            fi
                        fi
                        GWASs=${GWASs2}
                    done
                    ;;
            esac
        done
    done

    if [[ $analysis == "rg_meta" ]]; then
        phenos=$(cut -f1 ${fileList})
        GWASs=$(cut -f2 ${fileList})
        declare -a phenos=$(echo "($phenos)")
        declare -a GWASs=$(echo "($GWASs)")
        
        N_phenos=${#phenos[@]}

        for (( i=$((${N_phenos}-1)); i>=0; i-- )); do
            for (( j=$((${N_phenos}-1)); j>$i; j-- )); do
                if ! [[ -f ${phenos[$i]}/rg/rg_${GWASs[$i]}_${GWASs[$j]}.log ]] || ! [[ $(grep "Summary of Genetic Correlation Results" ${phenos[$i]}/rg/rg_${GWASs[$i]}_${GWASs[$j]}.log) ]]; then
                    echo "LDSC rg estimation between ${phenos[$i]} and ${phenos[$j]} was either not run before or it failed.." >> ${fileList}.rg_meta.error
                    echo -e "${phenos[$i]},${phenos[$j]}\t${GWASs[$i]},${GWASs[$j]}" >> ${fileList}.rg_meta.rerun
                fi
            done
        done
    fi

    status=0
    if ! [[ -f ${fileList}.${analysis}.rerun ]]; then
        echo "${analysis} of sumstats is complete."
    else
        echo "${analysis} of sumstats has finished but there were errors:"
        cat ${fileList}.${analysis}.error
        echo ""
        echo "Errors are stored in ${fileList}.${analysis}.error"
        status=1
    fi
}


###############################################################################

# FORMAT RESULTS FOR SINGLE-MTAG INPUT FILES (h2 of each input GWAS and rg between input GWASs, needed for QC)

# Write LDSC h^2 results for single-MTAG input files into a table
LDSC_h2_stats () {
    fileList=$1
    phenoList=$(cut -f1 $fileList| sed 's/[0-9]$//g' | sed 's/[0-9]$//g' | sort | uniq)

    for pheno in ${phenoList}; do
        echo -e "File\tSNPs\th2\tSE\tLambda_GC\tMeanChi2\tIntercept\tSE" > $pheno/h2_${pheno}.txt
	    for h2log in $pheno/h2/*.log; do
            GWAS=$(echo ${h2log} | cut -d"." -f1 | sed 's/h2_//g')
		    SNPs=$(grep "After merging with regression SNP LD" ${h2log} | cut -d" " -f7)
		    h2=$(grep "Total Observed scale h2" ${h2log} | cut -d":" -f2 | cut -d" " -f2)
		    h2_SE=$(grep "Total Observed scale h2" ${h2log} | cut -d":" -f2 | cut -d" " -f3)
		    Lambda_GC=$(grep "Lambda GC" ${h2log} | cut -d":" -f2)
		    MeanChi2=$(grep "Mean Chi^2" ${h2log} | cut -d":" -f2)
		    Intercept=$(grep "Intercept" ${h2log} | cut -d":" -f2 | cut -d" " -f2)
		    Intercept_SE=$(grep "Intercept" ${h2log} | cut -d":" -f2 | cut -d" " -f3 | sed 's/(//g' | sed 's/)//g')  
        
		    echo -e "${GWAS}\t${SNPs}\t${h2}\t${h2_SE}\t${Lambda_GC}\t${MeanChi2}\t${Intercept}\t${Intercept_SE}" >> $pheno/h2_${pheno}.txt
        done
    done
}

# Write LDSC rg between single-MTAG input files for each phenotype into a table
LDSC_rg_stats () {
    fileList=$1
    phenoList=$(cut -f1 $fileList| sed 's/[0-9]$//g' | sed 's/[0-9]$//g' | sort | uniq)

    for pheno in $phenoList; do
        rm -f $pheno/rg_${pheno}.txt
        sumstats=$(grep ^$pheno[1-9] $fileList | cut -f2 | sed 's/,/\n/g' | sort | uniq)
        GWASs=$(echo $sumstats | awk '{for (i=1;i<=NF;i++) {N=split($i,a,"/"); print a[N]}}' | sed 's/CLEANED\.//g'| sed 's/_mtag_meta_formatted\.txt//g')
        declare -a GWASs=$(echo "($GWASs)")
        N_GWASs=${#GWASs[@]}
        
        for (( i=$((${N_GWASs}-1)); i>=0; i-- )); do
            echo -e -n "\t"${GWASs[$i]} >> $pheno/rg_${pheno}.txt
        done
        
        for (( i=$((${N_GWASs}-1)); i>=0; i-- )); do
	        echo -e -n "\n"${GWASs[$i]} >> $pheno/rg_${pheno}.txt
		    for (( j=$((${N_GWASs}-1)); j>$i; j-- )); do
			    rg=$(grep -A 3 "Summary of Genetic Correlation Results" $pheno/rg/rg_${GWASs[$i]}_${GWASs[$j]}.log | sed -n '3p' | awk '{print $3,"("$4")"}')
                echo -e -n "\t"$rg >> $pheno/rg_${pheno}.txt
            done
	    done
	done
}

GWAS_h2_table () {
    fileList=$1
    phenoList=$(cut -f1 $fileList| sed 's/[0-9]$//g' | sed 's/[0-9]$//g' | sort | uniq)

    rm h2_table.txt
    echo -e "\tFile\tSNPs\th2\tSE\tLambda_GC\tMeanChi2\tIntercept\tSE" > h2_table.txt
    for pheno in ${phenoList}; do
        echo -e -n ${pheno} >> h2_table.txt
	    for h2log in $pheno/h2/*.log; do
            GWAS=$(echo ${h2log} | sed "s,${pheno}/h2/h2_${pheno}-,," | cut -d"." -f1)
		    SNPs=$(grep "After merging with regression SNP LD" ${h2log} | cut -d" " -f7)
		    h2=$(grep "Total Observed scale h2" ${h2log} | cut -d":" -f2 | cut -d" " -f2)
		    h2_SE=$(grep "Total Observed scale h2" ${h2log} | cut -d":" -f2 | cut -d" " -f3)
		    Lambda_GC=$(grep "Lambda GC" ${h2log} | cut -d":" -f2)
		    MeanChi2=$(grep "Mean Chi^2" ${h2log} | cut -d":" -f2)
		    Intercept=$(grep "Intercept" ${h2log} | cut -d":" -f2 | cut -d" " -f2)
		    Intercept_SE=$(grep "Intercept" ${h2log} | cut -d":" -f2 | cut -d" " -f3 | sed 's/(//g' | sed 's/)//g')  
        
		    echo -e "\t${GWAS}\t${SNPs}\t${h2}\t${h2_SE}\t${Lambda_GC}\t${MeanChi2}\t${Intercept}\t${Intercept_SE}" >> h2_table.txt
        done
    done
}

GWAS_rg_table () {
    fileList=$1
    phenoList=$(cut -f1 $fileList| sed 's/[0-9]$//g' | sed 's/[0-9]$//g' | sort | uniq)

    echo -e "\tGWAS1\tGWAS2\trg\trg_SE\tgcov_int\tgcov_int_se" > rg_table.txt
    for pheno in $phenoList; do
        sumstats=$(grep ^$pheno[1-9] $fileList | cut -f2 | sed 's/,/\n/g' | sort | uniq)
        GWASs=$(echo $sumstats | awk '{for (i=1;i<=NF;i++) {N=split($i,a,"/"); print a[N]}}' | sed 's/CLEANED\.//g'| sed 's/_mtag_meta_formatted\.txt//g')
        declare -a GWASs=$(echo "($GWASs)")
        N_GWASs=${#GWASs[@]}
        
        if [[ ${N_GWASs} == 1 ]]
        then
            echo -e ${pheno} >> rg_table.txt
        else
            echo -e -n ${pheno} >> rg_table.txt
        fi
        
        for (( i=$((${N_GWASs}-1)); i>=0; i-- )); do
            GWAS1=$(echo ${GWASs[$i]} | sed "s/${pheno}-//g")
		    for (( j=$((${N_GWASs}-1)); j>$i; j-- )); do
                GWAS2=$(echo ${GWASs[$j]} | sed "s/${pheno}-//g")
			    rg=$(grep -A 3 "Summary of Genetic Correlation Results" $pheno/rg/rg_${GWASs[$i]}_${GWASs[$j]}.log | sed -n '3p' | awk '{print $3}')
                rg_SE=$(grep -A 3 "Summary of Genetic Correlation Results" $pheno/rg/rg_${GWASs[$i]}_${GWASs[$j]}.log | sed -n '3p' | awk '{print $4}')
                gcov=$(grep -A 3 "Summary of Genetic Correlation Results" $pheno/rg/rg_${GWASs[$i]}_${GWASs[$j]}.log | sed -n '3p' | awk '{print $11}')
                gcov_SE=$(grep -A 3 "Summary of Genetic Correlation Results" $pheno/rg/rg_${GWASs[$i]}_${GWASs[$j]}.log | sed -n '3p' | awk '{print $12}')
                #echo -e "\t"${GWASs[$i]}"\t"${GWASs[$j]}"\t"$rg"\t"$rg_SE"\t"$gcov"\t"$gcov_SE >> rg_table.txt
                echo -e "\t"${GWAS1}"\t"${GWAS2}"\t"$rg"\t"$rg_SE"\t"$gcov"\t"$gcov_SE >> rg_table.txt
            done
	    done
	done
}

###############################################################################

# Calculate E(R^2) based on largest sample size (PHENO1) and largest h^2 MTAG output
# Write results into table
ER2_table() {
    fileList=$1
    dirIn=$2
    MTAGtype=$3

    phenoList=$(cut -f1 $fileList| sed 's/[0-9]$//g' | sed 's/[0-9]$//g' | sort | uniq)

    if [[ $MTAGtype == "single" ]]; then
        rm -f $PGI_Repo/code/5_LDSC/rg_meta_filelist
    fi

    echo -e "GWAS\t#SNPs MTAG\tMeanChi2\t#SNPs ldsc\th2\tSE\tGWAS equivalent N\tE(R2)" > ER2_table.txt

    for pheno in ${phenoList}; do
        numSumstats=$(ls $dirIn/${pheno}1/*_formatted_SBayesR.txt | wc -l)
        maxh2=0

        for (( i=1; i<=$numSumstats; i++ )); do        
            if [[ $MTAGtype == "single"  && $numSumstats == 1 ]]; then
                eval tag="trait"
            else
                eval tag="trait_$i"
            fi    

            h2log=${pheno}/h2/h2_${pheno}1_${tag}.log
            h2=$(grep "Total Observed scale h2" $h2log | cut -d":" -f2 | cut -d" " -f2)
            
            if [[ $h2 > $maxh2 ]]; then
                maxh2=$h2
                maxh2tag=$tag
                maxh2index=$i
            fi
        done

        gwasN=$(grep -A $maxh2index "GWAS equiv. (max) N" $dirIn/${pheno}1/${pheno}1.log | tail -1  | awk '{print $NF}')
        SNPsMTAG=$(grep -A $maxh2index "GWAS equiv. (max) N" $dirIn/${pheno}1/${pheno}1.log | tail -1  | awk '{print $3}')
        MeanChi2=$(grep -A $maxh2index "GWAS equiv. (max) N" $dirIn/${pheno}1/${pheno}1.log | tail -1  | awk '{print $7}')
        GWAS=${pheno}1_${maxh2tag}
		SNPsLDSC=$(grep "After merging with regression SNP LD" ${pheno}/h2/h2_${pheno}1_${maxh2tag}.log | cut -d" " -f7)
		h2_SE=$(grep "Total Observed scale h2" ${pheno}/h2/h2_${pheno}1_${maxh2tag}.log | cut -d":" -f2 | cut -d" " -f3 | sed 's/(//g' | sed 's/)//g')
        ER2=$(awk -v h2=$maxh2 -v N=$gwasN 'BEGIN{print h2/(1+(60000/(h2*N)))}')

        echo -e "${pheno}\t${SNPsMTAG}\t${MeanChi2}\t${SNPsLDSC}\t${maxh2}\t${h2_SE}\t${gwasN}\t${ER2}" >> ER2_table.txt

        # Write the trait number with largest h^2 into a file (will use those files to estimate pairwise rg's)
        if [[ $MTAGtype == "single" ]]; then
            echo -e "${pheno}\t${GWAS}" >> $PGI_Repo/code/5_LDSC/rg_meta_filelist
        fi
    done  
}

# Table for E(R^2) - observed R^2 comparison
# Also writes which MTAG output (trait_#) was used (i.e. which had largest h^2)
ER2_table_full() {
    fileList=$1
    dirIn=$2
    MTAGtype=$3

    phenoList=$(cut -f1 $fileList)

    echo -e "GWAS\tTrait Nr\t#SNPs MTAG\tMeanChi2\t#SNPs ldsc\th2\tSE\tGWAS equivalent N\tE(R2)" > ER2_table_full.txt


    for pheno in ${phenoList}; do
        phenoDir=$(echo $pheno | sed 's/[0-9]$//g' | sed 's/[0-9]$//g' )
        numSumstats=$(ls $dirIn/${pheno}/*_formatted_SBayesR.txt | wc -l)
        maxh2=0

        for (( i=1; i<=$numSumstats; i++ )); do        
            if [[ $MTAGtype == "single"  && $numSumstats == 1 ]]; then
                eval tag="trait"
            else
                eval tag="trait_$i"
            fi    

            h2log=${phenoDir}/h2/h2_${pheno}_${tag}.log
            h2=$(grep "Total Observed scale h2" $h2log | cut -d":" -f2 | cut -d" " -f2)
            
            if [[ $h2 > $maxh2 ]]; then
                maxh2=$h2
                maxh2tag=$tag
                maxh2index=$i
            fi
        done

        gwasN=$(grep -A $maxh2index "GWAS equiv. (max) N" $dirIn/${pheno}/${pheno}.log | tail -1  | awk '{print $NF}')
        SNPsMTAG=$(grep -A $maxh2index "GWAS equiv. (max) N" $dirIn/${pheno}/${pheno}.log | tail -1  | awk '{print $3}')
        MeanChi2=$(grep -A $maxh2index "GWAS equiv. (max) N" $dirIn/${pheno}/${pheno}.log | tail -1  | awk '{print $7}')
        GWAS=${pheno}_${maxh2tag}
		SNPsLDSC=$(grep "After merging with regression SNP LD" ${phenoDir}/h2/h2_${pheno}_${maxh2tag}.log | cut -d" " -f7)
		h2_SE=$(grep "Total Observed scale h2" ${phenoDir}/h2/h2_${pheno}_${maxh2tag}.log | cut -d":" -f2 | cut -d" " -f3 | sed 's/(//g' | sed 's/)//g')
        ER2=$(awk -v h2=$maxh2 -v N=$gwasN 'BEGIN{print h2/(1+(60000/(h2*N)))}')

        echo -e "${pheno}\t$maxh2index\t${SNPsMTAG}\t${MeanChi2}\t${SNPsLDSC}\t${maxh2}\t${h2_SE}\t${gwasN}\t${ER2}" >> ER2_table_full.txt
    done  
}


# Write rg between single-MTAG output for different phenotypes into a table
# Lines that are commented out are to make a table with SE's (decided to include no SEs in the table for the paper because it gets too busy)
rg_table(){
    fileList=$1

    phenos=$(cut -f1 ${fileList})
    GWASs=$(cut -f2 ${fileList})
    declare -a phenos=$(echo "($phenos)")
    declare -a GWASs=$(echo "($GWASs)")
    N_phenos=${#phenos[@]}
    echo $N_phenos

    for (( i=0; i<${N_phenos} ; i++ )); do
        echo "Getting rg's for ${phenos[$i]}"
        eval rg_${phenos[$i]}_${phenos[$i]}=1
        for (( j=$(($i+1)); j<${N_phenos}; j++ )); do
            echo "Evaluating rg between ${phenos[$i]} and ${phenos[$j]}"
            eval rg_${phenos[$i]}_${phenos[$j]}=$(grep -A 3 "Summary of Genetic Correlation Results" ${phenos[$i]}/rg/rg_${GWASs[$i]}_${GWASs[$j]}.log | sed -n '3p' | awk '{print $3}')
            #eval rg_se_${phenos[$i]}_${phenos[$j]}="$(grep -A 3 "Summary of Genetic Correlation Results" ${phenos[$i]}/rg/rg_${GWASs[$i]}_${GWASs[$j]}.log | sed -n '3p' | awk '{print $3,"("$4")"}')"
            eval rg_${phenos[$j]}_${phenos[$i]}='$'rg_${phenos[$i]}_${phenos[$j]}
            #eval rg_se_${phenos[$j]}_${phenos[$i]}="'$'rg_se_${phenos[$i]}_${phenos[$j]}"
        done
    done

    echo -e "\t${phenos[@]}" | sed 's/ /\t/g' > rg_table.txt
    #echo -e "\t${phenos[@]}" > rg_se_table.txt
    for (( i=0; i<${N_phenos} ; i++ )); do
        echo -e -n "${phenos[$i]}" >> rg_table.txt
        #echo -e -n "${phenos[$i]}" >> rg_se_table.txt
        for (( j=0; j<${N_phenos}; j++ )); do
            eval rg='$'rg_${phenos[$i]}_${phenos[$j]}
            echo -e -n "\t$rg" >> rg_table.txt
            #echo -e -n "\trg_se_${phenos[$i]}_${phenos[$j]}" >> rg_se_table.txt
        done
        echo -e -n "\n" >> rg_table.txt
        #echo -e "\n" >> rg_se_table.txt
    done    
}
