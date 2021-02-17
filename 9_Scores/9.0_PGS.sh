#!/bin/bash

dirCode="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/code/9_Scores"
dirGeno="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes"
score=$1
cohort=$2
dirOut=$3


##############################################################
########## Define LDpred input files and parameters ##########
##############################################################
LDgf_rs="/disk/genetics/HRC/aokbay/LDgf/HM3/HRC_HM3_geno02_mind02_rel025_nooutliers"
LDgf_chrpos="/disk/genetics/HRC/aokbay/LDgf/HM3/HRC_HM3_geno02_mind02_rel025_nooutliers_ChrPosID"

for i in "HRS3" "HRS2" "WLS" "Dunedin" "ERisk" "AH" "STRpsych" "STRtwge" "STRyatssstage" "Texas" "ELSA" "EGCUT" "MCTFR"; do
	declare valbim_${i}="${dirGeno}/${i}/plink/HM3/${i}_HM3"
	declare valgf_${i}="${dirGeno}/${i}/plink2/${i}_chr[1:22]"
	declare sample_${i}="${dirGeno}/${i}/sampleQC/${i}_EUR_FID_IID.txt"
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
	declare valbim_UKB${part}="${dirGeno}/UKB/plink/HM3/UKB_HM3"
	declare valgf_UKB${part}="/disk/genetics2/ukb/orig/UKBv3/imputed_plink2_HM3/ukb_imp_chr[1:22]_v3_HM3_nodup"
	declare sample_UKB${part}="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/1_UKB_GWAS/partitions/UKB_part${part}_eid.txt"
	declare snpidtype_UKB${part}="rs"
done

P=1


# 1:score (public / single / multi), 2:cohort 
LDpred(){
	fileList=$1
	cohort=$2
	dirOut=$3
	mkdir -p $dirOut/logs

	eval snpidtype='$'snpidtype_${cohort}
	eval LDgf='$'LDgf_${snpidtype}
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
		nohup bash $dirCode/9.0.1_LDpred.sh \
			--sumstats=$ssPath \
			--snpid=$snpid \
			--out=${cohort}_${pheno} \
			--LDgf=$LDgf \
			--Valbim=$valbim \
			--P=$P > $dirOut/logs/ldpred_${pheno}_${cohort}.log &
			
		let i+=1
		
		if [[ $i == 6 ]]; then
			wait
			i=0
		fi

	done < $fileList
	wait
}

checkStatus(){
	fileList=$1
	cohort=$2
	dirOut=$3
	step=$4

	echo "Checking status.."
	
	rm -f $dirCode/${cohort}_${score}_${step}_rerun
	
	status=1
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		
		case $step in
			LDpred)
				if ! ls $dirOut/pickled/${cohort}_${pheno}_*.pkl.gz 1> /dev/null 2>&1; then
					grep $pheno $dirCode/ss_${score}_${cohort} >> $dirCode/${cohort}_${score}_${step}_rerun
					echo "LDpred for $score $pheno score for $cohort was unsuccessful."
					status=0
				fi
				;;
			makePGS)
				if ! [[ $(find $dirOut/scores/PGS_${cohort}_${pheno}*.txt -type f -size +100 2>/dev/null) ]]; then 
					#! ls $dirOut/scores/PGS_${cohort}_${pheno}*.txt 1> /dev/null 2>&1; then
					grep $pheno $dirCode/ss_${score}_${cohort} >> $dirCode/${cohort}_${score}_${step}_rerun
					echo "makePGS for $score $pheno score for $cohort was unsuccessful."
					status=0
				fi
				;;
		esac

	done < $dirCode/ss_${score}_${cohort}

	if [[ -f $dirCode/${cohort}_${score}_${step}_rerun ]]; then
		mv $dirCode/${cohort}_${score}_${step}_rerun $dirCode/ss_${score}_${cohort}_${step}
	fi
}


makePGS(){
	fileList=$1
	cohort=$2
	dirOut=$3
	
	eval valgf='$'valgf_${cohort}
	eval sample='$'sample_${cohort}

	i=0
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		for weight in weights/${cohort}_${pheno}_weights_LDpred_p*.txt; do
			p=$(echo $weight | sed "s,weights/${cohort}_${pheno}_weights_LDpred_p,,g" | sed 's/\.txt//g')
			bash $dirCode/9.0.2_make_PGS.sh \
				--weight=weights/${cohort}_${pheno}_weights_LDpred_p*.txt \
				--weightCols=3,4,7 \
				--valgf=${valgf} \
				--sampleKeep=${sample} \
				--out=${cohort}_${pheno}_LDpred_p$P &
		done
		let i+=1
		
		if [[ $i == 1 ]]; then
			wait
			i=0
		fi
	done < $fileList
	wait
}


PGS(){
	score=$1
	cohort=$2
	dirOut=$3
	cd $dirOut
	
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
		
		checkStatus $score $cohort $dirOut LDpred
		status=$status

		if [[ $status == 0 ]]; then
			LDpred $dirCode/ss_${score}_${cohort}_LDpred $cohort $dirOut
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
		
		checkStatus $score $cohort $dirOut makePGS
		status=$status

		if [[ $status == 0 ]]; then
			makePGS $dirCode/ss_${score}_${cohort}_makePGS $cohort $dirOut
			pass=$(($pass+1))
		fi

		if [[ $pass > 3 ]]; then
			break
		fi
	done

	if [[ $status == 0  ]]; then
		echo "makePGS stage cannot be completed for $cohort $score scores. Check for errors in input files."
	fi

	
	rm -f $dirCode/ss_${score}_${cohort}

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


PGS $score $cohort $dirOut