#!/bin/bash

source $PGI_Repo/code/paths

score=$1
cohort=$2
method=$3

eval pgi_dir='$'pgi_dir_${cohort}
mkdir -p ${pgi_dir}/${score}_${method}/tmp
pgi_dir=${pgi_dir}/${score}_${method}

cd ${PGI_Repo}/derived_data/9_Scores/${score}_${method}
mkdir -p logs tmp weights

##############################################################
############## Define input files and parameters #############
##############################################################

eval gf_dir='$'gf_dir_${cohort}

case $cohort in
	"UKB1" | "UKB2" | "UKB3")
		valgf=$gf_plink2_UKB
		part=$(echo $cohort | sed 's/UKB//g')
		sample="$PGI_Repo/derived_data/1_UKB_GWAS/partitions/UKB_part${part}_eid.txt"
		snpidtype="rs"
		;;
	*)		
		valgf="${gf_dir}/plink2/${cohort}_chr[1:22]"
		sample="${gf_dir}/sampleQC/${cohort}_EUR_FID_IID.txt"
		case $i in
			"HRS2" | "WLS") 
				snpidtype="rs" 
				;;
			*) 
				snpidtype="chrpos" 
				;;
		esac
		;;
esac	

########################## LDpred ############################

case $cohort in
	"UKB1" | "UKB2" | "UKB3")
		valbim="${gf_dir}/plink/HM3/UKB_HM3"
		;;
	*)		
		valbim="${gf_dir}/plink/HM3/${cohort}_HM3"
		;;
esac	
P=1


########################## SBayesR ###########################

rm tmp/SBayesR_LDmatrices
for chr in {1..22}; do
	echo $SBayesR_LDmatrices | sed "s/\[1:22\]/$chr/" >> tmp/SBayesR_LDmatrices
done
pi=0.95,0.02,0.02,0.01
gamma=0.0,0.01,0.1,1
chainLength=10000
burnIn=2000

##############################################################
##############################################################
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

	i=0
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		ssPath=$(echo $row | cut -d" " -f2)	
	
		nohup $gctb --sbayes R \
    		--mldm tmp/SBayesR_LDmatrices \
        	--exclude-mhc \
     		--seed 123 \
     		--pi $pi \
     		--gamma $gamma \
     		--gwas-summary $ssPath \
     		--chain-length $chainLength \
     		--burn-in $burnIn \
     		--out-freq 100 \
     		--out tmp/${pheno}_weights_SBayesR 2>&1 | tee "logs/${pheno}_weights_SBayesR.log" &
	
		let i+=1
	
		if [[ $i == 1 ]]; then
			wait
			i=0
		fi

	done < $fileList
	wait

	i=0
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		awk '{$12=$3":"$4;print}' OFS="\t" tmp/${pheno}_weights_SBayesR.snpRes > weights/${pheno}_weights_SBayesR.txt &
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

	echo "Checking status.."
	
	rm -f $PGI_Repo/code/9_Scores/${cohort}_${score}_${step}_rerun
	
	status=1
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		phenoNoNum=$(echo $pheno | sed 's/[1-9]//g')
		
		case $step in
			weight)
				case $method in 
					LDpred)
						if ! ls pickled/${cohort}_${pheno}_*.pkl.gz 1> /dev/null 2>&1; then
							grep $pheno $PGI_Repo/code/9_Scores/ss_${score}_${cohort} >> $PGI_Repo/code/9_Scores/${cohort}_${score}_${step}_rerun
							echo "LDpred for $pheno in $cohort has not been run yet or was unsuccessful."
							status=0
						fi
						;;
					SBayesR)
                        if ! [[ $(find weights/${pheno}_weights_SBayesR.txt -type f -size +100 2>/dev/null) ]]; then
							grep $pheno $PGI_Repo/code/9_Scores/ss_${score}_${cohort} >> $PGI_Repo/code/9_Scores/${cohort}_${score}_${step}_rerun
							echo "SBayesR for $pheno in $cohort has not been run yet or was unsuccessful."
							status=0
						fi
						;;
				esac
				;;
			PGI)
				case $method in 
					LDpred)
						if ! [[ $(find ${pgi_dir}/PGS_${cohort}_${pheno}_LDpred_p*.txt -type f -size +20 2>/dev/null) ]]; then 
							grep $pheno $PGI_Repo/code/9_Scores/ss_${score}_${cohort} >> $PGI_Repo/code/9_Scores/${cohort}_${score}_${step}_rerun
							echo "makePGI (LDpred) for $pheno in $cohort has not been run yet or was unsuccessful."
							status=0
						fi
						;;
					SBayesR)
						if ! [[ $(find ${pgi_dir}/PGS_${cohort}_${phenoNoNum}_SBayesR.txt -type f -size +20 2>/dev/null) ]]; then 
							grep $pheno $PGI_Repo/code/9_Scores/ss_${score}_${cohort} >> $PGI_Repo/code/9_Scores/${cohort}_${score}_${step}_rerun
							echo "makePGI (SBayesR) for $pheno in $cohort has not been run yet or was unsuccessful."
							status=0
						fi
						;;
				esac
				;;
		esac


	done < $PGI_Repo/code/9_Scores/ss_${score}_${cohort}

	if [[ -f $PGI_Repo/code/9_Scores/${cohort}_${score}_${step}_rerun ]]; then
		mv $PGI_Repo/code/9_Scores/${cohort}_${score}_${step}_rerun $PGI_Repo/code/9_Scores/ss_${score}_${cohort}_${step}
	fi
}


makePGI(){
	fileList=$1

	i=0
	while read row; do

		pheno=$(echo $row | cut -d" " -f1)
		phenoNoNum=$(echo $pheno | sed 's/[1-9]//g')
		
		if [[ $method == "LDpred" ]]
		then
			cols="3,4,7"
			weights="weights/${cohort}_${pheno}_weights_LDpred_p*.txt"	
		elif [[ $method == "SBayesR" ]]
			then
				weights="weights/${pheno}_weights_SBayesR.txt"
				if [[ $snpidtype == "rs" ]]
        			then
            			cols="2,5,8"
       				else
            			cols="12,5,8"
    			fi
		fi

		bash $PGI_Repo/code/9_Scores/9.0.2_make_PGI.sh \
			--weight=${weights} \
			--weightCols=${cols} \
			--valgf=${valgf} \
			--sampleKeep=${sample} \
			--out=${pgi_dir}/PGS_${cohort}_${phenoNoNum}_${method} &

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
	echo -n "PGI ($method) on $cohort started on "
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
			$method $PGI_Repo/code/9_Scores/ss_${score}_${cohort}_weight
			pass=$(($pass+1))
		fi

		checkStatusPGI weight

		if [[ $pass > 1 ]]; then
			break
		fi
	done

	if [[ $status == 0  ]]; then
		echo "$method weight step cannot be completed for $cohort $score scores. Check for errors in input files."
	fi

	status=0
	pass=1

	while [[ $status == 0 ]]; do
		echo ""
		echo "makePGI $pass.."
		
		checkStatusPGI PGI
		status=$status

		if [[ $status == 0 ]]; then
			makePGI $PGI_Repo/code/9_Scores/ss_${score}_${cohort}_PGI
			pass=$(($pass+1))
		fi

		checkStatusPGI PGI

		if [[ $pass > 3 ]]; then
			break
		fi
	done

	if [[ $status == 0  ]]; then
		echo "$method - PGI stage cannot be completed for $cohort $score scores. Check for errors in input files."
	fi

	
	rm -f $PGI_Repo/code/9_Scores/ss_${score}_${cohort}

	echo ""
	echo "--------------------------------------------------"
	echo ""

	echo ""
	echo -n "Finished getting PGS for $cohort on "
	date

	end=$(date +%s)
	echo "Analysis took $(( ($end - $start)/60 )) minutes."
	echo "-----------------------------------------------------" 
	echo ""

}

main