#!/bin/bash

cd /disk/genetics4/projects/EA4/derived_data/PGS

dirCode="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code/9_Scores"
dirIn_public="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/3_QCd/public_scores/SEfilter"
dirOut_public="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/9_Scores/public"
#dirIn_single=
#dirOut_single
#dirIn_multi=
#dirOut_multi=


# Get list of sumstats for public scores: Pheno name on first column (e.g. SWB-Okbay), file path on second
rm -f $dirCode/ss_public
for dir in ${dirIn_public}/QC*; do
	if [[ -d $dir ]]; then 
		path=$(ls $dir/*.gz)
		pheno=$(echo $path | rev | cut -d"/" -f1 | rev | cut -d"." -f2)
		echo $pheno $path >> $dirCode/ss_public
	fi
done



##############################################################
########## Define LDpred input files and parameters ##########
##############################################################
LDgf_rs=/var/genetics2/HRC/aokbay/LDgf/HM3/HRC_HM3_geno02_mind02_rel025_nooutliers
LDgf_chrpos=/disk/genetics2/HRC/aokbay/LDgf/HM3/HM3_geno02_mind02_rel025_nooutliers_ChrPosID

valbim_HRS=/disk/genetics4/PGS/PGS_Repo/data/GENOTYPES/HRS_prelim/HRS_prelim
valbim_AH=/disk/genetics2/dbgap/data/derived/addhealth_HRCimputation/plink/AddHealth_HRC_chr1-22
valbim_WLS=/disk/genetics3/WLS_DBGAP/derived/plink/WLS_chr1-22
#valbim_STR-PSYCH=
#valbim_STR-TWGE=
#valbim_STR-YATSSSTAGE=
#valbim_UKB=

valgf_HRS=/disk/genetics/dbgap/aokbay/HRS/1kG_dosage/HRS_chr[1:22].gen.gz
valgf_AH=/disk/genetics/dbgap/data/derived/addhealth_HRCimputation/gen/AddHealth_HRC_chr[1:22].gen.gz
valgf_WLS=/disk/genetics/WLS_DBGAP/derived/gen/WLS_chr[1:22].gen.gz
valgf_STR-PSYCH=/disk/genetics/PGS/PGS_Repo/data/GENOTYPES/STR/STR_Salty/imputed/gen/STR_PSYCH_HRC_chr[1:22].gen.gz
valgf_STR-TWGE=/disk/genetics/PGS/PGS_Repo/data/GENOTYPES/STR/STR_Twingene/imputed/gen/STR_TWGE_HRC_chr[1:22].gen.gz
valgf_STR-YATSSSTAGE=/disk/genetics/PGS/PGS_Repo/data/GENOTYPES/STR/STR_YATSS_STAGE/imputed/gen/STR_YATSS-STAGE_HRC_chr[1:22].gen.gz
valgf_UKB=/disk/genetics2/ukb/orig/UKBv3/imputed_data/ukb_imp_chr[1:22]_v3.bgen

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
valgfFormat_UKB=bgen

snpidtype_HRS=rs
snpidtype_AH=chrpos
snpidtype_WLS=rs
snpidtype_STR-PSYCH=rs
snpidtype_STR-TWGE=rs
snpidtype_STR-YATSSSTAGE=rs
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
PGS(){

	score=$1
	cohort=$2

	eval pathOut='$'dirOut_${score}

	mkdir -p $pathOut/logs
	cd $pathOut

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
	eval valgf='$'valgf_${cohort}
	eval valgfFormat='$'valgfFormat_${cohort}
	eval genoFormat='$'genoFormat_${cohort}
	eval valbim='$'valbim_${cohort}


	i=0
	while read row; do
		pheno=$(echo $row | cut -d" " -f1)
		ssPath=$(echo $row | cut -d" " -f2)

		if ! [[ -f sumstats/${pheno}_${snpidtype}.txt ]]; then
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
			--out=${pheno}_${snpidtype} > $pathOut/logs/format_${pheno}_${snpidtype}.log
		fi


		nohup bash $dirCode/9.2_LDpred.sh \
			--efftype=$efftype \
			--sumstats=sumstats/sumstats_${pheno}_${snpidtype}.txt \
			--out=${cohort}_${pheno} \
			--LDgf=$LDgf \
			--Valbim=$valbim \
			--P=$P > $pathOut/logs/ldpred_${pheno}_${cohort}.log &
			
		let i+=1
		
		if [[ $i == 3 ]]; then
			wait
			i=0
		fi

		for weight in weights/${cohort}_${pheno}_weights_LDpred_p*.txt; do
			p=$(echo $weight | sed "s,weights/${cohort}_${pheno}_weights_LDpred_p,,g" | sed 's/\.txt//g')
			bash $dirCode/9.3_make_PGS.sh \
			--weight=weights/${cohort}_${pheno}_weights_LDpred_p*.txt \
			--weightCols=3,4,7 \
			--valgf=${valgf} \
			--valgfFormat=${valgfFormat} \
			--genoFormat=${genoFormat} \
			--out=${cohort}_${pheno}_LDpred_p$P &
		done

	done < $dirCode/ss_${score}

}

main()
{
	PGS public HRS
	PGS public WLS 
}

