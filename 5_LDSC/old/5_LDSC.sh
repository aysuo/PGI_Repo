#!/bin/bash

## FIX MUNGE CHECK STATUS - doesn't work in first run

alias python=/homes/nber/aokbay/anaconda2/bin/python2.7
dirCode=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code/5_LDSC
dirOut_study=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/5_LDSC/study_level
dirIn_study=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/3_QCd
dirOut_meta=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/5_LDSC/meta_level
hm3_snplist=/var/genetics/ukb/aokbay/reffiles/w_hm3.snplist
LDscores=/var/genetics/ukb/aokbay/bin/ldsc_old/eur_w_ld_chr/
LDSC=/disk/genetics/ukb/aokbay/bin/ldsc



ldsc_munge () {
    ss=$1
    out=$2

	## Munge sumstats ##
	python ${LDSC}/munge_sumstats.py \
		--sumstats ${ss} \
		--out ${out} \
		--merge-alleles ${hm3_snplist}
}


ldsc_h2 () {
    ssMunged=$1
    out=$2

	## Get h2 and intercept ##
	python ${LDSC}/ldsc.py \
		--h2 ${ssMunged}  \
		--ref-ld-chr ${LDscores} \
		--w-ld-chr ${LDscores} \
		--out ${out}
}

ldsc_rg () {
    ssPairMunged=$1
    out=$2

    ## Get rg 
	python ${LDSC}/ldsc.py \
	--rg ${ssPairMunged}  \
	--ref-ld-chr ${LDscores} \
	--w-ld-chr ${LDscores} \
	--out ${out}
}

munge () {
    fileListMunge=$1
    level=$2
    eval dirOut='$'dirOut_${level}
    eval dirIn='$'dirIn_${level}
    cd $dirOut

    i=1
    while read row; do
        pheno=$(echo ${row} | cut -d" " -f1)
        ss=$(echo ${row} | cut -d" " -f2)
        study=$(echo $ss | rev | cut -d"/" -f1 | rev | cut -d"." -f2)

        echo "Munging sumstats for ${study}.."
        ldsc_munge $dirIn/$ss ${pheno}/munged/${study} &
        
        let i+=1

        if [[ $i == 10 ]]; then
            wait
        fi
         
    done  < $fileListMunge
    wait
}

h2 () {
    fileListh2=$1
    level=$2
    eval dirOut='$'dirOut_${level}
    eval dirIn='$'dirIn_${level}
    cd $dirOut

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
        fi
    done < $fileListh2
    wait
}

rg () { 
    fileListrg=$1
    level=$2
    eval dirOut='$'dirOut_${level}
    eval dirIn='$'dirIn_${level}
    cd $dirOut

    i=1
    while read row; do
        pheno=$(echo ${row} | cut -d" " -f1)
        study1=$(echo ${row} | cut -d" " -f2 | cut -d"," -f1)
        study2=$(echo ${row} | cut -d" " -f2 | cut -d"," -f2)
        mkdir -p ${pheno}/rg

        echo "Estimating rg between ${study1} and ${study2}.."
        ldsc_rg ${pheno}/munged/${study1}.sumstats.gz,${pheno}/munged/${study2}.sumstats.gz \
            ${pheno}/rg/rg_${study1}_${study2} &
        let i+=1

        if [[ $i == 10 ]]; then
            wait
        fi
    done < $fileListrg
    wait
}


checkStatus () {
    fileList=$1
    level=$2
    analysis=$3
    eval dirOut='$'dirOut_${level}
    eval dirIn='$'dirIn_${level}
    cd $dirOut

    echo "Checking status of ${level}-level ${analysis} for ${fileList}.."

    rm -f $dirCode/error_${analysis}_${level}.txt
    rm -f $dirCode/rerun_${analysis}_${level}.txt
    
    phenoList=$(cut -f1 $fileList| sed 's/[1-9]//g' | sort | uniq)
    
    for pheno in $phenoList; do    
        sumstats=$(grep ^$pheno[1-9] $fileList | cut -f2 | sed 's/,/\n/g' | sort | uniq)
        studies=$(echo $sumstats | awk '{for (i=1;i<=NF;i++) {N=split($i,a,"/"); print a[N]}}' | sed 's/CLEANED\.//g'| sed 's/\.gz//g')
        
        for study in $studies; do
            case $analysis in 
                munge)
                    if ! [[ -f ${pheno}/munged/${study}.sumstats.gz ]] || [[ $(grep "ERROR converting summary statistics:" ${pheno}/munged/${study}.log) ]]; then
                        echo "Error: Munging for $study was unsuccessful.." >> $dirCode/error_munge_$level.txt
                        ss=$(echo $sumstats | awk -v study=$study '{for(i=1;i<=NF;i++) if ($i~study) print $i}')
                        echo -e "$pheno\t$ss" >> $dirCode/rerun_munge_$level.txt
                    fi
                    ;;
                h2)
                    if ! [[ -f ${pheno}/h2/h2_${study}.log ]] || ! [[ $(grep "Total Observed scale h2" ${pheno}/h2/h2_${study}.log) ]]; then
                        echo "Error: LDSC h2 estimation for $study  was unsuccessful.." >> $dirCode/error_h2_$level.txt
                        echo -e "$pheno\t$study" >> $dirCode/rerun_h2_$level.txt
                    fi
                    ;;
                rg) 
                    studies2=$(echo $studies | cut --complement -d" " -f1)

                    for study2 in $studies2; do
                        if ! [[ $study == $study2 ]]; then
                            if ! [[ -f ${pheno}/rg/rg_${study}_${study2}.log ]] || ! [[ $(grep "Summary of Genetic Correlation Results" ${pheno}/rg/rg_${study}_${study2}.log) ]]; then
                                echo "Error: LDSC rg estimation between $study and $study2 was unsuccessful .." >> $dirCode/error_rg_$level.txt
                                echo -e "$pheno\t$study,$study2" >> $dirCode/rerun_rg_$level.txt
                            fi
                        fi
                        studies=${studies2}
                    done
                    ;;
            esac
        done
    done

    status=0
    if ! [[ -f $dirCode/error_${analysis}_${level}.txt ]]; then
        echo "${analysis} of sumstats is complete."
    else
        echo "${analysis} of sumstats has finished but there were errors:"
        cat $dirCode/error_${analysis}_${level}.txt
        echo ""
        echo "Errors are stored in $dirCode/error_${analysis}_${level}.txt"
        status=1
    fi
}


LDSC() {
    fileList=$1
    level=$2

    for analysis in munge h2 rg; do
        checkStatus $fileList $level $analysis

        if [[ $status == 1 ]]; then
            echo "Rerunning ${level}-level analysis for the unfinished files in $fileList.."
            $analysis $dirCode/rerun_${analysis}_${level}.txt $level
            checkStatus $fileList $level $analysis
        fi

    done
}

    


# input: cohort list
LDSC_h2_stats () {
    fileList=$1
    level=$2
    eval dirOut='$'dirOut_${level}
    cd $dirOut

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



# input: cohort list
LDSC_rg_stats () {
    fileList=$1
    level=$2
    eval dirOut='$'dirOut_${level}
    cd $dirOut

    phenoList=$(cut -f1 $fileList| sed 's/[1-9]//g' | sort | uniq)

    for pheno in $phenoList; do
        rm $pheno/rg_${pheno}.txt
        sumstats=$(grep ^$pheno[1-9] $fileList | cut -f2 | sed 's/,/\n/g' | sort | uniq)
        studies=$(echo "$sumstats" | awk '{for (i=1;i<=NF;i++) {N=split($i,a,"/"); print a[N]}}' | sed 's/CLEANED\.//g'| sed 's/\.gz//g')
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

########################################

#LDSC $dirCode/singleMTAG_filelist.txt study
#LDSC_h2_stats $dirCode/singleMTAG_filelist.txt study
LDSC_rg_stats $dirCode/singleMTAG_filelist.txt study



   




   