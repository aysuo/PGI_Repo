#!/bin/bash

dirIn="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/1_UKB_GWAS/output"
dirLogs="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/1_UKB_GWAS/logs"
dirOut="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/UKB"
dirCode="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code/2_Formatting"
pathCallrate="/disk/genetics2/ukb/orig/UKBv2/snpstats/ukb_imp_v2.lmiss"
pathHWE="/disk/genetics2/ukb/orig/UKBv2/snpstats/ukb_imp_v2.hwe"

cd $dirOut

format_UKB(){
	pheno=$1
	part=$2
	
	# Get N from log files
	N=$(grep -a "Number of indivs with no missing phenotype(s) to use:" $dirLogs/$pheno/UKB_${pheno}_part${part}_BOLTLMM.log | cut -d" " -f10)

	# Add N and callrate
	awk -F"\t" -v N=$N 'NR==FNR{a[$2]=1-$5;next} \
	FNR==1{print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE"}
	($1 in a && FNR>1) {print $1,$2,$3,$5,$6,$7,$9,$10,$11,$8,N,"1",a[$1];next}
	(FNR>1) {print $1,$2,$3,$5,$6,$7,$9,$10,$11,$8,N,"1","NA"}' OFS="\t" $pathCallrate $dirIn/$pheno/UKB_${pheno}_part${part}_BOLTLMM > ${pheno}-UKB${part}_tmp

	# Add HWE p-value
	awk -F"\t" 'NR==FNR{a[$2]=$9;next} \
	FNR==1{print $0, "HWE_PVAL"} \
	($1 in a && FNR>1) {print $0,a[$1];next} \
	(FNR>1) {print $0,"NA"}' OFS="\t" $pathHWE ${pheno}-UKB${part}_tmp > ${pheno}-UKB${part}.txt

	rm ${pheno}-UKB${part}_tmp
}

main() {
	#phenos=$(ls -d $dirIn/* | rev | cut -d"/" -f1 | rev | sed '/^old/d')
	for part in 1 2 3; do
		for pheno in $phenos; do
			echo "Formatting $pheno part $part"
			format_UKB $pheno $part &
		done
		wait
	done
}

main



