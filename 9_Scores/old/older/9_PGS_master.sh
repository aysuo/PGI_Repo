#!/bin/bash

dirCode="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code/9_Scores"
dirIn_public="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/3_QCd/public_scores/SEfilter"
dirOut_public="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/9_Scores/public"

#dirIn_single=
#dirOut_single
#dirIn_multi=
#dirOut_multi=



##############################################################
########## Define LDpred input files and parameters ##########
##############################################################
LDgf_rs="/var/genetics2/HRC/aokbay/LDgf/HM3/HRC_HM3_geno02_mind02_rel025_nooutliers"
LDgf_chrpos="/disk/genetics2/HRC/aokbay/LDgf/HM3/HRC_HM3_geno02_mind02_rel025_nooutliers_ChrPosID"

valbim_HRS="/disk/genetics4/projects/EA4/derived_data/PGS/reffiles/HRS_HM3"
valbim_AH="/disk/genetics4/projects/EA4/derived_data/PGS/reffiles/AH_HM3"
valbim_WLS="/disk/genetics4/projects/EA4/derived_data/PGS/reffiles/WLS_HM3"
valbim_STR_PSYCH="/disk/genetics4/projects/EA4/derived_data/PGS/reffiles/STR-PSYCH_HM3"
valbim_STR_TWGE="/disk/genetics4/projects/EA4/derived_data/PGS/reffiles/STR-TWGE_HM3"
valbim_STR_YATSSSTAGE="/disk/genetics4/projects/EA4/derived_data/PGS/reffiles/STR-YATSSSTAGE_HM3"
valbim_UKB=$LDgf_rs

valgf_HRS="/disk/genetics/dbgap/data/HRS_Jun19_2018/HRS/files/62567/PhenoGenotypeFiles/RootStudyConsentSet_phs000428.CIDR_Aging_Omni1.v2.p2.c1.NPR/GenotypeFiles/phg000515.v1.HRS_phase123_imputation.genotype-imputed-data.c1/imputed/HRS2_chr[1:22].gprobs.gz"
valgf_AH="/disk/genetics/dbgap/data/derived/addhealth_HRCimputation/gen/AddHealth_HRC_chr[1:22].gen.gz"
valgf_WLS="/disk/genetics/WLS_DBGAP/derived/gen/WLS_chr[1:22].gen.gz"
valgf_STR_PSYCH="/disk/genetics/PGS/PGS_Repo/data/GENOTYPES/STR/STR_Salty/imputed/gen/STR_PSYCH_HRC_chr[1:22].gen.gz"
valgf_STR_TWGE="/disk/genetics/PGS/PGS_Repo/data/GENOTYPES/STR/STR_Twingene/imputed/gen/STR_TWGE_HRC_chr[1:22].gen.gz"
valgf_STR_YATSSSTAGE="/disk/genetics/PGS/PGS_Repo/data/GENOTYPES/STR/STR_YATSS_STAGE/imputed/gen/STR_YATSS-STAGE_HRC_chr[1:22].gen.gz"
valgf_UKB="/disk/genetics2/ukb/orig/UKBv3/imputed_plink_HM3/ukb_imp_chr[1:22]_v3_HM3_nodup"

sample_HRS="/disk/genetics/dbgap/data/HRS_Jun19_2018/HRS/files/62567/PhenoGenotypeFiles/RootStudyConsentSet_phs000428.CIDR_Aging_Omni1.v2.p2.c1.NPR/GenotypeFiles/phg000515.v1.HRS_phase123_imputation.genotype-imputed-data.c1/phased/HRS2_chr1.sample.gz"
sample_AH="/disk/genetics2/dbgap/data/derived/addhealth_HRCimputation/gen/AddHealth_HRC_chr1.samples"
sample_WLS="/disk/genetics3/WLS_DBGAP/derived/gen/WLS_chr1.sample"
sample_STR_PSYCH="/disk/genetics4/PGS/PGS_Repo/data/GENOTYPES/STR/STR_Salty/imputed/gen/STR_PSYCH_HRC_chr1.sample"
sample_STR_TWGE="/disk/genetics4/PGS/PGS_Repo/data/GENOTYPES/STR/STR_Twingene/imputed/gen/STR_TWGE_HRC_chr1.sample"
sample_STR_YATSSSTAGE="/disk/genetics4/PGS/PGS_Repo/data/GENOTYPES/STR/STR_YATSS_STAGE/imputed/gen/STR_YATSS-STAGE_HRC_chr1.sample"
sample_UKB="/disk/genetics2/ukb/orig/UKBv2/linking/ukb_imp_v2.sample"

genoFormat_HRS=dosage
genoFormat_AH=dosage
genoFormat_WLS=dosage
genoFormat_STR-PSYCH=dosage
genoFormat_STR-TWGE=dosage
genoFormat_STR-YATSSSTAGE=dosage
genoFormat_UKB=hardcall

valgfFormat_HRS=gen
valgfFormat_AH=gen
valgfFormat_WLS=gen
valgfFormat_STR-PSYCH=gen
valgfFormat_STR-TWGE=gen
valgfFormat_STR-YATSSSTAGE=gen
valgfFormat_UKB=plink2Chr

snpidtype_HRS=rs
snpidtype_AH=chrpos
snpidtype_WLS=rs
snpidtype_STR-PSYCH=chrpos
snpidtype_STR-TWGE=chrpos
snpidtype_STR-YATSSSTAGE=chrpos
snpidtype_UKB=rs

rsid_public=SNPID
chrposid_public=cptid
chr_public=CHR
bp_public=BP
effall_public=EFFECT_ALLELE
altall_public=OTHER_ALLELE
eaf_public=EAF
zscore_public=NA
effect_public=EFFECT
efftype_public=LINREG
se_public=SE
pval_public=PVALUE
info_public=INFO      
N_public=N

P=1


# 1:score (public / single / multi), 2:cohort 
LDpred(){
	score=$1
	cohort=$2

	mkdir -p $dirOut/logs

	eval rsid='$'rsid_${score} 
	eval chrposid='$'chrposid_${score} 
	eval chr='$'chr_${score}
	eval bp='$'bp_${score}
	eval effall='$'effall_${score}
	eval altall='$'altall_${score}
	eval eaf='$'eaf_${score}
	eval effect='$'effect_${score}
	eval zscore='$'zscore_${score}
	eval efftype='$'efftype_${score}
	eval se='$'se_${score}
	eval pval='$'pval_${score}
	eval info='$'info_${score}
	eval N='$'N_${score}

	eval snpidtype='$'snpidtype_${cohort}
	eval LDgf='$'LDgf_${snpidtype}
	eval valbim='$'valbim_${cohort}


	i=0
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		ssPath=$(echo $row | cut -d" " -f2)

		if ! [[ -f sumstats/sumstats_${pheno}_${snpidtype}.txt ]]; then
			bash $dirCode/9.1_format_sumstats.sh \
			--snpidtype=$snpidtype \
			--rsid=$rsid \
			--chrposid=$chrposid \
			--chr=$chr \
			--bp=$bp \
			--effall=$effall \
			--altall=$altall \
			--eaf=$eaf \
			--effect=$effect \
			--zscore=$zscore \
			--efftype=$efftype \
			--se=$se \
			--pval=$pval \
			--N=$N \
			--sumstats=$ssPath \
			--out=${pheno}_${snpidtype} > $dirOut/logs/format_${pheno}_${snpidtype}.log
		else
			echo "Sumstats in LDpred format for ${pheno} ${snpidtype} already exists."
		fi
	
		rm -f coord/${cohort}_${pheno}.coord
		nohup bash $dirCode/9.2_LDpred.sh \
			--efftype=$efftype \
			--sumstats=sumstats/sumstats_${pheno}_${snpidtype}.txt \
			--out=${cohort}_${pheno} \
			--LDgf=$LDgf \
			--Valbim=$valbim \
			--P=$P > $dirOut/logs/ldpred_${pheno}_${cohort}.log &
			
		let i+=1
		
		if [[ $i == 3 ]]; then
			wait
			i=0
		fi

	done < $dirCode/ss_${score}_${cohort}
	wait
}

checkStatus(){
	score=$1
	cohort=$2

	echo "Checking status.."
	
	rm -f $dirCode/${cohort}_${score}_rerun
	
	status=1
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		
		if ! ls $dirOut/pickled/${cohort}_${pheno}_*.pkl.gz 1> /dev/null 2>&1; then
			grep $pheno $dirCode/ss_${score}_${cohort} >> $dirCode/${cohort}_${score}_rerun
			echo "$score $pheno score for $cohort was unsuccessful."
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
	
	eval valgf='$'valgf_${cohort}
	eval valgfFormat='$'valgfFormat_${cohort}
	eval genoFormat='$'genoFormat_${cohort}
	eval sample='$'sample_${cohort}

	i=0
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)

		for weight in weights/${cohort}_${pheno}_weights_LDpred_p*.txt; do
			p=$(echo $weight | sed "s,weights/${cohort}_${pheno}_weights_LDpred_p,,g" | sed 's/\.txt//g')
			bash $dirCode/9.3_make_PGS.sh \
			--weight=weights/${cohort}_${pheno}_weights_LDpred_p*.txt \
			--weightCols=3,4,7 \
			--valgf=${valgf} \
			--sample=${sample} \
			--valgfFormat=${valgfFormat} \
			--genoFormat=${genoFormat} \
			--out=${cohort}_${pheno}_LDpred_p$P &
		done
		let i+=1
		
		if [[ $i == 5 ]]; then
			wait
			i=0
		fi

	done < $dirCode/ss_${score}_${cohort}
}


PGS(){
	score=$1
	cohort=$2
	
	eval dirOut='$'dirOut_${score}
	eval dirIn='$'dirIn_${score}
	cd $dirOut
	
	echo "----------------------------------------------------------------------"
	echo -n "PGS on $cohort started on "
	date
	echo ""
	start=$(date +%s)


	
	# Get list of sumstats for public scores: Pheno name on first column (e.g. SWB-Okbay), file path on second
	rm -f $dirCode/ss_${score}_${cohort}
	for dir in ${dirIn}/QC*; do
		if [[ -d $dir ]]; then 
			path=$(ls $dir/*.gz)
			pheno=$(echo $path | rev | cut -d"/" -f1 | rev | cut -d"." -f2)
			echo $pheno $path >> $dirCode/ss_${score}_${cohort}
		fi
	done

	pass=1
	status=0

	while [[ $status == 0 ]]; do
		echo ""
		echo "LDpred pass $pass.."
		
		checkStatus $score $cohort
		status=$status
		echo status=$status

		if [[ $status == 0 ]]; then
			LDpred $score $cohort
			pass=$(($pass+1))
		fi


		echo status=$status
		if [[ $pass > 3 ]]; then
			break
		fi

	done

	if [[ $status == 0  ]]; then
		echo "LDpred stage cannot be completed for $cohort $score scores. Check for errors in input files."
	else
		makePGS $score $cohort
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

main(){
	PGS public WLS
	#PGS public HRS
}

main

