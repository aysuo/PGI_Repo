#!/bin/bash

source paths5

ldsc_munge () {
    ss=$1
    out=$2

	## Munge sumstats ##
	$python ${LDSC}/munge_sumstats.py \
		--sumstats ${ss} \
		--out ${out} \
		--merge-alleles ${hm3snps}
}

ldsc_h2 () {
    ssMunged=$1
    out=$2

	## Get h2 and intercept ##
	$python ${LDSC}/ldsc.py \
		--h2 ${ssMunged}  \
		--ref-ld-chr ${LDscores} \
		--w-ld-chr ${LDscores} \
        --two-step inf \
		--out ${out}
}

ldsc_rg () {
    ssPairMunged=$1
    out=$2

    ## Get rg 
	$python ${LDSC}/ldsc.py \
	    --rg ${ssPairMunged}  \
	    --ref-ld-chr ${LDscores} \
	    --w-ld-chr ${LDscores} \
        --two-step inf \
	    --out ${out}
}

###############################################################################

munge () {
    fileListMunge=$1

    i=1
    while read row; do
        pheno=$(echo ${row} | cut -d" " -f1)
        ss=$(echo ${row} | cut -d" " -f2)
        eval ss=$ss
        study=$(echo $ss | rev | cut -d"/" -f1 | rev | sed 's/CLEANED\.//g' | sed 's/_formatted\.txt//g')

        mkdir -p $pheno/munged
        echo "Munging sumstats for ${study}.."
        ldsc_munge $ss ${pheno}/munged/${study} &
        
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

    i=1
    while read row; do
        pheno=$(echo ${row} | cut -d" " -f1)
        study=$(echo ${row} | cut -d" " -f2)

        mkdir -p ${pheno}/h2
        echo "Estimating h2 for ${study}.."
        ldsc_h2 ${pheno}/munged/${study}.sumstats.gz ${pheno}/h2/h2_${study} &
        
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

    i=1
    while read row; do
        pheno1=$(echo ${row} | cut -d" " -f1 | cut -d"," -f1)
        pheno2=$(echo ${row} | cut -d" " -f1 | cut -d"," -f2)
        study1=$(echo ${row} | cut -d" " -f2 | cut -d"," -f1)
        study2=$(echo ${row} | cut -d" " -f2 | cut -d"," -f2)
        mkdir -p ${pheno1}/rg

        echo "Estimating rg between ${study1} and ${study2}.."
        ldsc_rg ${pheno1}/munged/${study1}.sumstats.gz,${pheno2}/munged/${study2}.sumstats.gz \
            ${pheno1}/rg/rg_${study1}_${study2} &
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
    
    phenoList=$(cut -f1 $fileList| sed 's/[1-9]//g' | sort | uniq)
    
    for pheno in $phenoList; do    
        sumstats=$(grep ^$pheno[1-9] $fileList | cut -f2 | sed 's/,/\n/g' | sort | uniq)
        studies=$(echo $sumstats | awk '{for (i=1;i<=NF;i++) {N=split($i,a,"/"); print a[N]}}' | sed 's/CLEANED\.//g'| sed 's/_formatted\.txt//g')
        
        for study in $studies; do
            case $analysis in 
                munge)
                    if ! [[ -f ${pheno}/munged/${study}.sumstats.gz ]] || [[ $(grep "ERROR converting summary statistics:" ${pheno}/munged/${study}.log) ]]; then
                        echo "Munging for $study was either not run before or it failed.." >> ${fileList}.munge.error
                        ss=$(echo $sumstats | awk -v study=$study '{for(i=1;i<=NF;i++) if ($i~study) print $i}')
                        echo -e "$pheno\t$ss" >> ${fileList}.munge.rerun
                    fi
                    ;;
                h2)
                    if ! [[ -f ${pheno}/h2/h2_${study}.log ]] || ! [[ $(grep "Total Observed scale h2" ${pheno}/h2/h2_${study}.log) ]]; then
                        echo "LDSC h2 estimation for $study was either not run before or it failed.." >> ${fileList}.h2.error
                        echo -e "$pheno\t$study" >> ${fileList}.h2.rerun
                    fi
                    ;;
                rg) 
                    studies2=$(echo $studies | cut --complement -d" " -f1)

                    for study2 in $studies2; do
                        if ! [[ $study == $study2 ]]; then
                            if ! [[ -f ${pheno}/rg/rg_${study}_${study2}.log ]] || ! [[ $(grep "Summary of Genetic Correlation Results" ${pheno}/rg/rg_${study}_${study2}.log) ]]; then
                                echo "LDSC rg estimation between $study and $study2 was either not run before or it failed.." >> ${fileList}.rg.error
                                echo -e "$pheno,$pheno\t$study,$study2" >> ${fileList}.rg.rerun
                            fi
                        fi
                        studies=${studies2}
                    done
                    ;;
            esac
        done
    done

    if [[ $analysis == "rg_meta" ]]; then
        phenos=$(cut -f1 ${fileList})
        studies=$(cut -f2 ${fileList})
        declare -a phenos=$(echo "($phenos)")
        declare -a studies=$(echo "($studies)")
        
        N_phenos=${#phenos[@]}

        for (( i=$((${N_phenos}-1)); i>=0; i-- )); do
            for (( j=$((${N_phenos}-1)); j>$i; j-- )); do
                if ! [[ -f ${phenos[$i]}/rg/rg_${studies[$i]}_${studies[$j]}.log ]] || ! [[ $(grep "Summary of Genetic Correlation Results" ${phenos[$i]}/rg/rg_${studies[$i]}_${studies[$j]}.log) ]]; then
                    echo "LDSC rg estimation between ${phenos[$i]} and ${phenos[$j]} was either not run before or it failed.." >> ${fileList}.rg_meta.error
                    echo -e "${phenos[$i]},${phenos[$j]}\t${studies[$i]},${studies[$j]}" >> ${fileList}.rg_meta.rerun
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

# FORMAT RESULTS FOR SINGLE-MTAG INPUT FILES (h2 of each input ss and rg between input ss, needed for QC)

# Write LDSC h^2 results for single-MTAG input files into a table
LDSC_h2_stats () {
    fileList=$1
    phenoList=$(cut -f1 $fileList| sed 's/[1-9]//g' | sort | uniq)

    for pheno in ${phenoList}; do
        echo -e "File\tSNPs\th2\tSE\tLambda_GC\tMeanChi2\tIntercept\tSE" > $pheno/h2_${pheno}.txt
	    for h2log in $pheno/h2/*.log; do
            study=$(echo ${h2log} | cut -d"." -f1 | sed 's/h2_//g')
		    SNPs=$(grep "After merging with regression SNP LD" ${h2log} | cut -d" " -f7)
		    h2=$(grep "Total Observed scale h2" ${h2log} | cut -d":" -f2 | cut -d" " -f2)
		    h2_SE=$(grep "Total Observed scale h2" ${h2log} | cut -d":" -f2 | cut -d" " -f3)
		    Lambda_GC=$(grep "Lambda GC" ${h2log} | cut -d":" -f2)
		    MeanChi2=$(grep "Mean Chi^2" ${h2log} | cut -d":" -f2)
		    Intercept=$(grep "Intercept" ${h2log} | cut -d":" -f2 | cut -d" " -f2)
		    Intercept_SE=$(grep "Intercept" ${h2log} | cut -d":" -f2 | cut -d" " -f3 | sed 's/(//g' | sed 's/)//g')  
        
		    echo -e "${study}\t${SNPs}\t${h2}\t${h2_SE}\t${Lambda_GC}\t${MeanChi2}\t${Intercept}\t${Intercept_SE}" >> $pheno/h2_${pheno}.txt
        done
    done
}

# Write LDSC rg between single-MTAG input files for each phenotype into a table
LDSC_rg_stats () {
    fileList=$1
    phenoList=$(cut -f1 $fileList| sed 's/[1-9]//g' | sort | uniq)

    for pheno in $phenoList; do
        rm $pheno/rg_${pheno}.txt
        sumstats=$(grep ^$pheno[1-9] $fileList | cut -f2 | sed 's/,/\n/g' | sort | uniq)
        studies=$(echo $sumstats | awk '{for (i=1;i<=NF;i++) {N=split($i,a,"/"); print a[N]}}' | sed 's/CLEANED\.//g'| sed 's/_mtag_meta_formatted\.txt//g')
        declare -a studies=$(echo "($studies)")
        N_studies=${#studies[@]}
        
        for (( i=$((${N_studies}-1)); i>=0; i-- )); do
            echo -e -n "\t"${studies[$i]} >> $pheno/rg_${pheno}.txt
        done
        
        for (( i=$((${N_studies}-1)); i>=0; i-- )); do
	        echo -e -n "\n"${studies[$i]} >> $pheno/rg_${pheno}.txt
		    for (( j=$((${N_studies}-1)); j>$i; j-- )); do
			    rg=$(grep -A 3 "Summary of Genetic Correlation Results" $pheno/rg/rg_${studies[$i]}_${studies[$j]}.log | sed -n '3p' | awk '{print $3,"("$4")"}')
                echo -e -n "\t"$rg >> $pheno/rg_${pheno}.txt
            done
	    done
	done
}

###############################################################################

# FORMAT RESULTS FOR SINGLE-/MULTI-MTAG OUTPUT FILES

# Calculate E(R^2) based on largest sample size (PHENO1) and largest h^2 MTAG output
# Write results into table
ER2_table() {
    fileList=$1
    dirIn=$2
    single=$3

    phenoList=$(cut -f1 $fileList| sed 's/[1-9]//g' | sort | uniq)

    if [[ $single == 1 ]]; then
        rm -f $dirCode/5_LDSC/rg_meta_filelist
    fi

    echo -e "Study\t#SNPs MTAG\tMeanChi2\t#SNPs ldsc\th2\tSE\tGWAS equivalent N\tE(R2)" > ER2_table.txt

    for pheno in ${phenoList}; do
        numSumstats=$(ls $dirIn/${pheno}1/*_formatted.txt | wc -l)
        maxh2=0

        for (( i=1; i<=$numSumstats; i++ )); do        
            if [[ $single == 1  && $numSumstats == 1 ]]; then
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
        study=${pheno}1_${maxh2tag}
		SNPsLDSC=$(grep "After merging with regression SNP LD" ${pheno}/h2/h2_${pheno}1_${maxh2tag}.log | cut -d" " -f7)
		h2_SE=$(grep "Total Observed scale h2" ${pheno}/h2/h2_${pheno}1_${maxh2tag}.log | cut -d":" -f2 | cut -d" " -f3 | sed 's/(//g' | sed 's/)//g')
        ER2=$(awk -v h2=$maxh2 -v N=$gwasN 'BEGIN{print h2/(1+(60000/(h2*N)))}')

        echo -e "${pheno}\t${SNPsMTAG}\t${MeanChi2}\t${SNPsLDSC}\t${maxh2}\t${h2_SE}\t${gwasN}\t${ER2}" >> ER2_table.txt

        # Write the trait number with largest h^2 into a file (will use those files to estimate pairwise rg's)
        if [[ $single == 1 ]]; then
            echo -e "${pheno}\t${study}" >> $dirCode/5_LDSC/rg_meta_filelist
        fi
    done  
}

# Table for E(R^2) - observed R^2 comparison
# Also writes which MTAG output (trait_#) was used (i.e. which had largest h^2)
ER2_table_full() {
    fileList=$1
    dirIn=$2
    single=$3

    phenoList=$(cut -f1 $fileList)

    echo -e "Study\tTrait Nr\t#SNPs MTAG\tMeanChi2\t#SNPs ldsc\th2\tSE\tGWAS equivalent N\tE(R2)" > ER2_table_full.txt


    for pheno in ${phenoList}; do
        phenoDir=$(echo $pheno |  sed 's/[1-9]//g')
        numSumstats=$(ls $dirIn/${pheno}/*_formatted.txt | wc -l)
        maxh2=0

        for (( i=1; i<=$numSumstats; i++ )); do        
            if [[ $single == 1  && $numSumstats == 1 ]]; then
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
        study=${pheno}_${maxh2tag}
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
    studies=$(cut -f2 ${fileList})
    declare -a phenos=$(echo "($phenos)")
    declare -a studies=$(echo "($studies)")
    N_phenos=${#phenos[@]}
    echo $N_phenos

    for (( i=0; i<${N_phenos} ; i++ )); do
        echo "Getting rg's for ${phenos[$i]}"
        eval rg_${phenos[$i]}_${phenos[$i]}=1
        for (( j=$(($i+1)); j<${N_phenos}; j++ )); do
            echo "Evaluating rg between ${phenos[$i]} and ${phenos[$j]}"
            eval rg_${phenos[$i]}_${phenos[$j]}=$(grep -A 3 "Summary of Genetic Correlation Results" ${phenos[$i]}/rg/rg_${studies[$i]}_${studies[$j]}.log | sed -n '3p' | awk '{print $3}')
            #eval rg_se_${phenos[$i]}_${phenos[$j]}="$(grep -A 3 "Summary of Genetic Correlation Results" ${phenos[$i]}/rg/rg_${studies[$i]}_${studies[$j]}.log | sed -n '3p' | awk '{print $3,"("$4")"}')"
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
