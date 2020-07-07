#!/bin/bash

dirCode="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code/9_Scores"
score=$1
cohort=$2
dirOut=$3

##############################################################
########## Define LDpred input files and parameters ##########
##############################################################
LDgf_rs="/disk/genetics/HRC/aokbay/LDgf/HM3/HRC_HM3_geno02_mind02_rel025_nooutliers"
LDgf_chrpos="/disk/genetics/HRC/aokbay/LDgf/HM3/HRC_HM3_geno02_mind02_rel025_nooutliers_ChrPosID"

valbim_HRS2="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes/HRS2/plink/HM3/HRS2_HM3"
#valbim_AH="/disk/genetics4/projects/EA4/derived_data/PGS/reffiles/AH_HM3"
valbim_WLS="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes/WLS/plink/HM3/WLS_HM3"
#valbim_STR_PSYCH="/disk/genetics4/projects/EA4/derived_data/PGS/reffiles/STR-PSYCH_HM3"
#valbim_STR_TWGE="/disk/genetics4/projects/EA4/derived_data/PGS/reffiles/STR-TWGE_HM3"
#valbim_STR_YATSSSTAGE="/disk/genetics4/projects/EA4/derived_data/PGS/reffiles/STR-YATSSSTAGE_HM3"
#valbim_UKB=$LDgf_rs

valgf_HRS2="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes/HRS2/plink2/HRS2_chr[1:22]"
#valgf_AH="/disk/genetics/dbgap/data/derived/addhealth_HRCimputation/gen/AddHealth_HRC_chr[1:22].gen.gz"
valgf_WLS="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes/WLS/plink2/WLS_chr[1:22]"
#valgf_STR_PSYCH="/disk/genetics/PGS/PGS_Repo/data/GENOTYPES/STR/STR_Salty/imputed/gen/STR_PSYCH_HRC_chr[1:22].gen.gz"
#valgf_STR_TWGE="/disk/genetics/PGS/PGS_Repo/data/GENOTYPES/STR/STR_Twingene/imputed/gen/STR_TWGE_HRC_chr[1:22].gen.gz"
#valgf_STR_YATSSSTAGE="/disk/genetics/PGS/PGS_Repo/data/GENOTYPES/STR/STR_YATSS_STAGE/imputed/gen/STR_YATSS-STAGE_HRC_chr[1:22].gen.gz"
#valgf_UKB="/disk/genetics2/ukb/orig/UKBv3/imputed_plink_HM3/ukb_imp_chr[1:22]_v3_HM3_nodup"

sample_HRS2="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes/HRS2/sampleQC/HRS2_EUR_FID_IID.txt"
#sample_AH="/disk/genetics2/dbgap/data/derived/addhealth_HRCimputation/gen/AddHealth_HRC_chr1.samples"
sample_WLS="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes/WLS/sampleQC/WLS_EUR_FID_IID.txt"
#sample_STR_PSYCH="/disk/genetics4/PGS/PGS_Repo/data/GENOTYPES/STR/STR_Salty/imputed/gen/STR_PSYCH_HRC_chr1.sample"
#sample_STR_TWGE="/disk/genetics4/PGS/PGS_Repo/data/GENOTYPES/STR/STR_Twingene/imputed/gen/STR_TWGE_HRC_chr1.sample"
#sample_STR_YATSSSTAGE="/disk/genetics4/PGS/PGS_Repo/data/GENOTYPES/STR/STR_YATSS_STAGE/imputed/gen/STR_YATSS-STAGE_HRC_chr1.sample"
#sample_UKB="/disk/genetics2/ukb/orig/UKBv2/linking/ukb_imp_v2.sample"

snpidtype_HRS2=rs
#snpidtype_AH=chrpos
snpidtype_WLS=rs
#snpidtype_STR-PSYCH=chrpos
#snpidtype_STR-TWGE=chrpos
#snpidtype_STR-YATSSSTAGE=chrpos
#snpidtype_UKB=rs

P=1


# 1:score (public / single / multi), 2:cohort 
LDpred(){
	score=$1
	cohort=$2
	dirOut=$3

	mkdir -p $dirOut/logs

	eval snpidtype='$'snpidtype_${cohort}
	eval LDgf='$'LDgf_${snpidtype}
	eval valbim='$'valbim_${cohort}
	
	i=0
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		ssPath=$(echo $row | cut -d" " -f2)

		rm -f coord/${cohort}_${pheno}.coord
		nohup bash $dirCode/9.0.1_LDpred.sh \
			--sumstats=$ssPath \
			--out=${cohort}_${pheno} \
			--LDgf=$LDgf \
			--Valbim=$valbim \
			--P=$P > $dirOut/logs/ldpred_${pheno}_${cohort}.log &
			
		let i+=1
		
		if [[ $i == 6 ]]; then
			wait
			i=0
		fi

	done < $dirCode/ss_${score}_${cohort}
	wait
}

checkStatus(){
	score=$1
	cohort=$2
	dirOut=$3
	step=$4

	echo "Checking status.."
	
	rm -f $dirCode/${cohort}_${score}_rerun
	
	status=1
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		
		if [[ $step == "LDpred" ]] && ! ls $dirOut/pickled/${cohort}_${pheno}_*.pkl.gz 1> /dev/null 2>&1; then
			grep $pheno $dirCode/ss_${score}_${cohort} >> $dirCode/${cohort}_${score}_rerun
			echo "LDpred for $score $pheno score for $cohort was unsuccessful."
			status=0
		elif [[ $step == "makePGS" ]] && ! ls $dirOut/scores/PGS_${cohort}_${pheno}* 1> /dev/null 2>&1; then
			grep $pheno $dirCode/ss_${score}_${cohort} >> $dirCode/${cohort}_${score}_rerun
			echo "makePGS for $score $pheno score for $cohort was unsuccessful."
			status=0
		fi
	done < $dirCode/ss_${score}_${cohort}

	if [[ -f $dirCode/${cohort}_${score}_rerun ]]; then
		mv $dirCode/${cohort}_${score}_rerun $dirCode/ss_${score}_${cohort}
	fi
}


makePGS(){
	score=$1
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
		
		if [[ $i == 6 ]]; then
			wait
			i=0
		fi

	done < $dirCode/ss_${score}_${cohort}
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
			LDpred $score $cohort $dirOut
			pass=$(($pass+1))
		fi

		if [[ $pass > 3 ]]; then
			break
		fi
	done

	if [[ $status == 0  ]]; then
		echo "LDpred stage cannot be completed for $cohort $score scores. Check for errors in input files."
	fi

	status=0
	checkStatus $score $cohort $dirOut makePGS
	status=$status

	if [[ $status == 0 ]]; then
		makePGS $score $cohort $dirOut
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