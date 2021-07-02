#!/bin/bash

source $mainDir/code/paths

score=$1
cohort=$2
dirOut=$3

cd $dirOut

##############################################################
########## Define LDpred input files and parameters ##########
##############################################################

for i in AH Dunedin EGCUT ERisk ELSA HRS2 HRS3 MCTFR STRpsych STRtwge STRyatssstage Texas WLS; do
	declare valbim_${i}="$mainDir/derived_data/7_Genotypes/${i}/plink/HM3/${i}_HM3"
	declare valgf_${i}="$mainDir/derived_data/7_Genotypes/${i}/plink2/${i}_chr[1:22]"
	declare sample_${i}="$mainDir/derived_data/7_Genotypes/${i}/sampleQC/${i}_EUR_FID_IID.txt"
	case $i in
		"HRS2" | "WLS") 
			declare snpidtype_${i}="rs" 
			;;
		*) 
			declare snpidtype_${i}="chrpos" 
			;;
	esac
done

for part in 1 2 3; do
	declare valbim_UKB${part}="$mainDir/derived_data/7_Genotypes/UKB/plink/HM3/UKB_HM3"
	declare valgf_UKB${part}=$gf_plink2_UKB
	declare sample_UKB${part}="$mainDir/derived_data/1_UKB_GWAS/partitions/UKB_part${part}_eid.txt"
	declare snpidtype_UKB${part}="rs"
done

##----------------------------------------------------------##

P=1

##############################################################

LDpred(){
	fileList=$1
	cohort=$2
	mkdir -p logs

	eval snpidtype='$'snpidtype_${cohort}
	eval LDgf='$'HRC_LDgf_${snpidtype}
	eval valbim='$'valbim_${cohort}

	if [[ $snpidtype == "rs" ]]; then
		snpid="SNPID"; else
		snpid="cptid"
	fi
	
	i=0
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		ssPath=$(echo $row | cut -d" " -f2)

		rm -f coord/${cohort}_${pheno}.coord
		nohup bash $mainDir/code/9_Scores/9.0.1_LDpred.sh \
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

checkStatusPGI(){
	fileList=$1
	cohort=$2
	step=$3

	echo "Checking status.."
	
	rm -f $mainDir/code/9_Scores/${cohort}_${score}_${step}_rerun
	
	status=1
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		
		case $step in
			LDpred)
				if ! ls pickled/${cohort}_${pheno}_*.pkl.gz 1> /dev/null 2>&1; then
					grep $pheno $mainDir/code/9_Scores/ss_${score}_${cohort} >> $mainDir/code/9_Scores/${cohort}_${score}_${step}_rerun
					echo "LDpred for $score $pheno score for $cohort was unsuccessful."
					status=0
				fi
				;;
			makePGS)
				if ! [[ $(find scores/PGS_${cohort}_${pheno}*.txt -type f -size +100 2>/dev/null) ]]; then 
					grep $pheno $mainDir/code/9_Scores/ss_${score}_${cohort} >> $mainDir/code/9_Scores/${cohort}_${score}_${step}_rerun
					echo "makePGS for $score $pheno score for $cohort was unsuccessful."
					status=0
				fi
				;;
		esac

	done < $mainDir/code/9_Scores/ss_${score}_${cohort}

	if [[ -f $mainDir/code/9_Scores/${cohort}_${score}_${step}_rerun ]]; then
		mv $mainDir/code/9_Scores/${cohort}_${score}_${step}_rerun $mainDir/code/9_Scores/ss_${score}_${cohort}_${step}
	fi
}


makePGS(){
	fileList=$1
	cohort=$2
	
	eval valgf='$'valgf_${cohort}
	eval sample='$'sample_${cohort}

	i=0
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		for weight in weights/${cohort}_${pheno}_weights_LDpred_p*.txt; do
			p=$(echo $weight | sed "s,weights/${cohort}_${pheno}_weights_LDpred_p,,g" | sed 's/\.txt//g')
			bash $mainDir/code/9_Scores/9.0.2_make_PGS.sh \
				--weight=weights/${cohort}_${pheno}_weights_LDpred_p*.txt \
				--weightCols=3,4,7 \
				--valgf=${valgf} \
				--sampleKeep=${sample} \
				--out=${cohort}_${pheno}_LDpred_p$P &
		done
		let i+=1
		
		if [[ $i == 5 ]]; then
			wait
			i=0
		fi
	done < $fileList
	wait
}


PGS(){
	score=$1
	cohort=$2
	
	echo "----------------------------------------------------------------------"
	echo -n "PGS on $cohort started on "
	date
	echo ""
	start=$(date +%s)

	pass=1
	status=0

	while [[ $status == 0 ]]; do
		echo ""
		echo "LDpred pass $pass.."
		
		checkStatusPGI $score $cohort LDpred
		status=$status

		if [[ $status == 0 ]]; then
			LDpred $mainDir/code/9_Scores/ss_${score}_${cohort}_LDpred $cohort
			pass=$(($pass+1))
		fi

		if [[ $pass > 1 ]]; then
			break
		fi
	done

	if [[ $status == 0  ]]; then
		echo "LDpred stage cannot be completed for $cohort $score scores. Check for errors in input files."
	fi

	status=0
	pass=1

	while [[ $status == 0 ]]; do
		echo ""
		echo "makePGS $pass.."
		
		checkStatusPGI $score $cohort makePGS
		status=$status

		if [[ $status == 0 ]]; then
			makePGS $mainDir/code/9_Scores/ss_${score}_${cohort}_makePGS $cohort
			pass=$(($pass+1))
		fi

		if [[ $pass > 3 ]]; then
			break
		fi
	done

	if [[ $status == 0  ]]; then
		echo "makePGS stage cannot be completed for $cohort $score scores. Check for errors in input files."
	fi

	
	rm -f $mainDir/code/9_Scores/ss_${score}_${cohort}

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


PGS $score $cohort