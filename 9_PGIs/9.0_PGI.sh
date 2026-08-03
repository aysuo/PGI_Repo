#!/bin/bash

source $PGI_Repo/code/paths

score=$1
cohort=$2
method=$3
gen=$4
snps=$5
snpidtype=$6
valgf=$7
sample=$8
ancestry=$9
valbim=${10}
phased=${11}
onlyweights=${12}
sBayesRsnps=${13}

eval pgi_dir='$'pgi_dir_${cohort}
eval gf_dir='$'gf_dir_${cohort}
mkdir -p ${pgi_dir}/${score}_${method}/tmp
pgi_dir=${pgi_dir}/${score}_${method}

cd ${PGI_Repo}/derived_data/9_Scores/${score}_${method}
mkdir -p logs tmp weights

imputeNrobust=("CPD1-single" "CPD4-single" "CPD5-single" "CPD6-single" "EVERSMOKE1-single" "EVERSMOKE4-single" "EVERSMOKE5-single" "EVERSMOKE6-single")

##############################################################
############## Define input files and parameters #############
##############################################################

########################## LDpred ############################
P=1
########################## SBayesR ###########################
if [[ $sBayesRsnps == "NA" || $sBayesRsnps == "2.9m" ]]; then
	LDmatrices=$SBayesR_LDmatrices
elif [[ $sBayesRsnps == "HM3" ]]; then
	LDmatrices=$SBayesR_LDmatrices_HM3
fi

rm tmp/SBayesR_LDmatrices
for chr in {1..22}; do
	echo $LDmatrices | sed "s/\[1:22\]/$chr/" >> tmp/SBayesR_LDmatrices
done
pi=0.95,0.02,0.02,0.01
gamma=0.0,0.01,0.1,1
chainLength=10000
burnIn=2000

##############################################################

LDpred(){
	fileList=$1

	eval LDgf='$'HRC_LDgf_${snpidtype}

	if [[ $snpidtype == "rs" ]]; then
		snpid="SNPID"; else
		snpid="cptid"
	fi
	
	i=0
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		ssPath=$(echo $row | cut -d" " -f2)

		rm -f coord/${cohort}_${pheno}.coord
		nohup bash $PGI_Repo/code/9_Scores/9.0.1_LDpred.sh \
			--sumstats=$ssPath \
			--snpid=$snpid \
			--out=${cohort}_${pheno} \
			--LDgf=$LDgf \
			--Valbim=$valbim \
			--P=$P > logs/ldpred_${pheno}_${cohort}.log &
			
		let i+=1
		
		if [[ $i == 6 ]]; then
			wait
			i=0
		fi

	done < $fileList
	wait
}

##############################################################


SBayesR(){
	fileList=$1
	sBayesRsnps=$2

	i=0
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		ssPath=$(echo $row | cut -d" " -f2)	
			
		if [[ ${imputeNrobust[@]} =~ $pheno  ]]
		then
			imputeNrobustFlag="--impute-n --robust"
		else
			imputeNrobustFlag=""
		fi

		if [[ $sBayesRsnps == "NA" || $sBayesRsnps == "2.9m" ]]; then 
			out=${pheno}_weights_SBayesR
		elif [[ $sBayesRsnps == "HM3" ]]; then
			out=${pheno}_weights_SBayesR_HM3
		fi

		nohup $gctb --sbayes R $imputeNrobustFlag \
			--mldm tmp/SBayesR_LDmatrices \
			--seed 123 \
			--exclude-mhc \
			--pi $pi \
			--gamma $gamma \
			--gwas-summary $ssPath \
			--chain-length $chainLength \
			--burn-in $burnIn \
			--out-freq 100 \
			--out tmp/${out} 2>&1 | tee "logs/${out}.log" &
			let i+=1
	
		if [[ $i == 3 ]]; then
			wait
			i=0
		fi

	done < $fileList
	wait

	i=0
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		awk '{$12=$3":"$4;print}' OFS="\t" tmp/${out}.snpRes > weights/${out}.txt &
		let i+=1

		if [[ $i == 10 ]]; then
			wait
			i=0
		fi
	done < $fileList
	wait
}

##############################################################


checkStatusPGI(){
	step=$1
	sBayesRsnps=$2

	echo "Checking status.."
	
	rm -f $PGI_Repo/code/9_Scores/${cohort}_${ancestry}_${gen}_${score}_${step}_rerun
	
	status=1
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		phenoNoNum=$(echo $pheno | sed 's/[0-9]-/-/g' | sed 's/[0-9]-/-/g')
		
		case $step in
			weight)
				case $method in 
					LDpred)
						if ! ls pickled/${cohort}_${pheno}_*.pkl.gz 1> /dev/null 2>&1; then
							grep $pheno $PGI_Repo/code/9_Scores/ss_single_${cohort}_${gen}_${ancestry} >> $PGI_Repo/code/9_Scores/${cohort}_${ancestry}_${gen}_${score}_${step}_rerun
							echo "LDpred for $pheno in $cohort has not been run yet or was unsuccessful."
							status=0
						fi
						;;
					SBayesR)
						if [[ $sBayesRsnps == "NA" || $sBayesRsnps == "2.9m" ]]; then 
							weights=weights/${pheno}_weights_SBayesR.txt
						elif [[ $sBayesRsnps == "HM3" ]]; then
							weights=weights/${pheno}_weights_SBayesR_HM3.txt
						fi

                        if ! [[ $(find ${weights} -type f -size +100 2>/dev/null) ]]; then
							grep $pheno $PGI_Repo/code/9_Scores/ss_single_${cohort}_${gen}_${ancestry} >> $PGI_Repo/code/9_Scores/${cohort}_${ancestry}_${gen}_${score}_${step}_rerun
							echo "SBayesR for $pheno (snps: $sBayesRsnps) in $cohort has not been run yet or was unsuccessful."
							status=0
						fi
						;;
				esac
				;;
			PGI)	
				if [[ $gen == "parent" ]]; then
					parental_tag="_parental"
				elif [[ $gen == "proband" ]]; then
					parental_tag=""
				fi	

				case $method in 
					LDpred)
						if ! [[ $(find ${pgi_dir}/PGI_${cohort}_${ancestry}${cohort_tag}_${pheno}_LDpred_p*.txt -type f -size +5 2>/dev/null) ]]; then 
							grep $pheno $PGI_Repo/code/9_Scores/ss_single_${cohort}_${gen}_${ancestry} >> $PGI_Repo/code/9_Scores/${cohort}_${ancestry}_${gen}_${score}_${step}_rerun
							echo "makePGI (LDpred) for $pheno in $cohort $gen has not been run yet or was unsuccessful."
							status=0
						fi
						;;
					SBayesR)
						if [[ $snps == "NA" ]]; then
							dirToCheck=${pgi_dir}
						else
							dirToCheck=${pgi_dir}/restrictedSNPs
						fi

						if [[ $sBayesRsnps == "NA" || $sBayesRsnps == "2.9m" ]]; then 
							PGI=${dirToCheck}/PGI_${cohort}_${ancestry}${parental_tag}_${phenoNoNum}_SBayesR.txt
						elif [[ $sBayesRsnps == "HM3" ]]; then
							PGI=${dirToCheck}/PGI_${cohort}_${ancestry}${parental_tag}_${phenoNoNum}_SBayesR_HM3.txt
						fi

						if ! [[ $(find ${dirToCheck}/PGI_${cohort}_${ancestry}${parental_tag}_${phenoNoNum}_SBayesR.* -type f -size +5 2>/dev/null) ]]; then 
							grep $pheno $PGI_Repo/code/9_Scores/ss_single_${cohort}_${gen}_${ancestry} >> $PGI_Repo/code/9_Scores/${cohort}_${ancestry}_${gen}_${score}_${step}_rerun
							echo "makePGI (SBayesR) for $pheno (snps: $sBayesRsnps) in $cohort $gen has not been run yet or was unsuccessful."
							status=0
						fi
						;;
				esac
				;;
		esac
	done < $PGI_Repo/code/9_Scores/ss_single_${cohort}_${gen}_${ancestry}

	if [[ -f $PGI_Repo/code/9_Scores/${cohort}_${ancestry}_${gen}_${score}_${step}_rerun ]]; then
		mv $PGI_Repo/code/9_Scores/${cohort}_${ancestry}_${gen}_${score}_${step}_rerun $PGI_Repo/code/9_Scores/ss_single_${cohort}_${gen}_${ancestry}_${step}
	fi
}


makePGI(){
	fileList=$1
	sBayesRsnps=$2

	i=0
	while read row; do

		pheno=$(echo $row | cut -d" " -f1)
		phenoNoNum=$(echo $pheno | sed 's/[0-9]-/-/g' | sed 's/[0-9]-/-/g' )
		
		if [[ $method == "LDpred" ]]; then
			weights="weights/${cohort}_${pheno}_weights_LDpred_p*.txt"
			if [[ $snpidtype ==  "rs" || $snpidtype == "ChrPosIDhg19" ]]; then
					cols="3,4,7"
			elif [[ $snpidtype == "ChrPosIDhg38" ]]; then		
				if [[ ! -f $weights.hg38 ]]; then
					echo "Lifting over LDpred weights for $pheno to hg38..."
					awk 'NR>1{gsub(/chrom_/,"",$1); print "chr"$1,$2-1,$2,$0}' OFS="\t" $weights > $weights.lift.bed
					$liftover -bedPlus=3 $weights.lift.bed $chain_hg19toHg38 $weights.lifted $weights.unlifted
					awk -F"\t" 'BEGIN{OFS="\t"; print "chrom_hg19","pos_hg19","sid_hg19","nt1","nt2","raw_beta","ldpred_beta","ChrPosID_hg38"} \
						{gsub(/chr/,"",$1); $11=$1":"$3; print $4,$5,$6,$7,$8,$9,$10,$11}' OFS="\t" $weights.lifted > $weights.hg38
					rm $weights.lift.bed $weights.lifted
				fi
				weights=$weights.hg38
				cols="8,4,7"
			fi		
		elif [[ $method == "SBayesR" ]]; then	
			if [[ $sBayesRsnps == "NA" || $sBayesRsnps == "2.9m" ]]; then 
				weights=weights/${pheno}_weights_SBayesR.txt
				out=PGI_${cohort}_${ancestry}_${phenoNoNum}_${method}
			elif [[ $sBayesRsnps == "HM3" ]]; then
				weights=weights/${pheno}_weights_SBayesR_HM3.txt
				out=PGI_${cohort}_${ancestry}_${phenoNoNum}_${method}_HM3
			fi

			if [[ $snpidtype ==  "rs" ]]; then
				cols="2,5,8"
			elif [[ $snpidtype == "ChrPosIDhg19" ]]; then
				cols="12,5,8"
			elif [[ $snpidtype == "ChrPosIDhg38" ]]; then
				if [[ ! -f $weights.hg38 ]]; then
					echo "Lifting over SBayesR weights for $pheno to hg38..."
					awk -F"\t" 'NR>1{print "chr"$3,$4-1,$4,$0}' OFS="\t" $weights > $weights.lift.bed
						$liftover -bedPlus=3 $weights.lift.bed $chain_hg19toHg38 $weights.lifted $weights.unlifted
					awk -F"\t" 'BEGIN{OFS="\t"; print "Id","Name","Chrom_hg19","Position_hg19","A1","A2","A1Frq","A1Effect","SE","PIP","LastSampleEff","ChrPosID_hg19","ChrPosID_hg38"} \
						{gsub(/chr/,"",$1); $16=$1":"$3; print $4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' OFS="\t" $weights.lifted > $weights.hg38
					rm $weights.lift.bed $weights.lifted
				fi
				weights=$weights.hg38
				cols="13,5,8"
			fi

		fi

		if [[ $snps != "NA" ]]; then
			mkdir -p ${pgi_dir}/restrictedSNPs
			bash $PGI_Repo/code/9_Scores/9.0.2_make_PGI.sh \
				--weights=${weights} \
				--weightCols=${cols} \
				--valgf=${valgf} \
				--sampleKeep=${sample} \
				--extract=${snps} \
				--out=${pgi_dir}/restrictedSNPs/${out} &
		else
			bash $PGI_Repo/code/9_Scores/9.0.2_make_PGI.sh \
				--weights=${weights} \
				--weightCols=${cols} \
				--valgf=${valgf} \
				--sampleKeep=${sample} \
				--out=${pgi_dir}/${out} &
		fi

		let i+=1
		
		if [[ $i == 1 ]]; then
			wait
			i=0
		fi
	done < $fileList
	wait
}


makefPGI(){
    fileList=$1

	source $snipar_venv

    # Get paths
	if [[ $cohort == "UKB1" ]]
		then
			gf_plink=$gf_dir/degree1_relatives_plink/UKB_degree1_chr@
			gf_parental=$gf_dir/parental/chr_@
		else
			gf_bgen=$gf_dir/parental/input/bgen/${cohort}_chr@
			gf_plink=$gf_dir/parental/input/plink/${cohort}_chr@
			gf_parental=$gf_dir/parental/output/${cohort}_parental_chr@
	fi

	if [[ $cohort == "STRtwge" ]]
		then
			am_adj="--no_am_adj"
	fi

    i=0
    while read row; do
        pheno=$(echo $row | cut -d" " -f1 | tr -d '#')
        phenoNoNum=$(echo $pheno | sed 's/[0-9]-/-/g' | sed 's/[0-9]-/-/g' )
        
        if [[ $method == "LDpred" ]]
			then
				a1="nt1"
				a2="nt2"
				beta="ldpred_beta"
				if [[ $snpidtype ==  "rs" || $snpidtype == "ChrPosIDhg19" ]]
					then
						weights="weights/${cohort}_${pheno}_weights_LDpred_p*.txt"
						snpid="sid"
				elif [[ $snpidtype == "ChrPosIDhg38" ]]
					then
						weights="weights/${cohort}_${pheno}_weights_LDpred_p*.txt.hg38"
						snpid="ChrPosID_hg38"
				fi
		elif [[ $method == "SBayesR" ]]
			then
				a1="A1"
                a2="A2"
                beta="A1Effect"
				if [[ $snpidtype ==  "rs" ]]
					then
						weights="weights/${pheno}_weights_SBayesR.txt"
						snpid="Name"
				elif [[ $snpidtype == "ChrPosIDhg19" ]]
					then
						weights="weights/${pheno}_weights_SBayesR.txt"
						snpid="Chrom:Position"
				elif [[ $snpidtype == "ChrPosIDhg38" ]]
					then
						weights="weights/${pheno}_weights_SBayesR.txt.hg38"
						snpid="ChrPosID_hg38"
				fi
		fi

        #Calculate PGI for parents with imputed genotypes: Use --no-am-adj for STRtwge
		if [[ $phased == 1 ]]
			then
				pgs.py \
				${pgi_dir}/PGS_${cohort}_parental_${phenoNoNum}_${method} \
				--bgen $gf_bgen \
				--imp $gf_parental \
				--beta_col $beta \
				--SNP $snpid \
				--A1 $a1 \
				--A2 $a2 \
				--weights $weights \
				$am_adj
			else
				pgs.py \
				${pgi_dir}/PGS_${cohort}_parental_${phenoNoNum}_${method} \
				--bed $gf_plink \
				--imp $gf_parental \
				--beta_col $beta \
				--SNP $snpid \
				--A1 $a1 \
				--A2 $a2 \
				--weights $weights \
				$am_adj
		fi


        let i+=1

        if [[ $i == 1 ]]; then
            wait
            i=0
        fi
    done < $fileList

    wait
}




main(){
	echo "----------------------------------------------------------------------"
	echo -n "PGI ($method) on $cohort-$gen started on "
	date
	echo ""
	start=$(date +%s)

	pass=1
	status=0

	while [[ $status == 0 ]]; do
		echo ""
		echo "$method weight step - pass $pass.."
		
		checkStatusPGI weight
		status=$status

		if [[ $status == 0 ]]; then
			$method $PGI_Repo/code/9_Scores/ss_single_${cohort}_${gen}_${ancestry}_weight
			pass=$(($pass+1))
		fi

		checkStatusPGI weight
		pass=$(($pass+1))
		if [[ $pass > 1 ]]; then
			break
		fi
	done

	if [[ $status == 0  ]]; then
		echo "$method weight step cannot be completed for $cohort $gen $score scores. Check for errors in input files."
	fi

	if [[ $onlyweights == 1 ]];
	then
		echo "Skipping PGI step for $cohort as no genetic data available."
		exit 0
	fi

	status=0
	pass=1

	while [[ $status == 0 ]]; do
		echo ""
		echo "makePGI $pass.."
		
		checkStatusPGI PGI
		status=$status
	
		if [[ $status == 0 && $gen == "proband" ]]; then
			makePGI $PGI_Repo/code/9_Scores/ss_single_${cohort}_${gen}_${ancestry}_PGI
			pass=$(($pass+1))
		elif [[ $status == 0 && $gen == "parent" ]]; then
			makefPGI $PGI_Repo/code/9_Scores/ss_single_${cohort}_${gen}_${ancestry}_PGI 
			pass=$(($pass+1))
		fi

		checkStatusPGI PGI

		if [[ $pass > 3 ]]; then
			break
		fi
	done

	if [[ $status == 0  ]]; then
		echo "$method - PGI stage cannot be completed for $cohort $gen $score. Check for errors in input files."
	fi

	
	rm -f $PGI_Repo/code/9_Scores/ss_single_${cohort}_${gen}_${ancestry}

	echo ""
	echo "--------------------------------------------------"
	echo ""

	echo ""
	echo -n "Finished getting PGIs for $cohort $gen on "
	date

	end=$(date +%s)
	echo "Analysis took $(( ($end - $start)/60 )) minutes."
	echo "-----------------------------------------------------" 
	echo ""

}

main