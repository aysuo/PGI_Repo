#!/bin/bash

source $PGI_Repo/code/paths

cd $PGI_Repo/derived_data/2_Formatted/UKB


format_UKB(){
	inputGWAS=$1
	inputLOG=$2
	output=$3
	
	# Get N from log files
	N=$(grep -a "Number of indivs with no missing phenotype(s) to use:" $inputLOG | cut -d" " -f10)

	# Add N and callrate
	awk -F"\t" -v N=$N 'NR==FNR{a[$2]=1-$5;next} \
	FNR==1{print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE"}
	($1 in a && FNR>1) {print $1,$2,$3,$5,$6,$7,$9,$10,$11,$8,N,"1",a[$1];next}
	(FNR>1) {print $1,$2,$3,$5,$6,$7,$9,$10,$11,$8,N,"1","NA"}' OFS="\t" $UKB_callrate $inputGWAS > ${output}_tmp

	# Add HWE p-value
	awk -F"\t" 'NR==FNR{a[$2]=$9;next} \
	FNR==1{print $0, "HWE_PVAL"} \
	($1 in a && FNR>1) {print $0,a[$1];next} \
	(FNR>1) {print $0,"NA"}' OFS="\t" $UKB_hwe ${output}_tmp > ${output}

	rm ${output}_tmp

	gzip ${output}
}


main() {
	phenos=$(ls -d $PGI_Repo/derived_data/1_UKB_GWAS/output/* | rev | cut -d"/" -f1 | rev | sed '/^old/d')

	for pheno in $phenos; do
	if ! [[ ${pheno} == "ANOREX" ]]
	then
		for part in 1 2 3; do
			if [[ -f ${pheno}-UKB${part}.txt ]]
			then
				echo "Formatted file for $pheno part $part already exists, skipping"
			else
				echo "Formatting $pheno part $part"
				format_UKB \
					$PGI_Repo/derived_data/1_UKB_GWAS/output/$pheno/UKB_${pheno}_part${part}_BOLTLMM \
					$PGI_Repo/derived_data/1_UKB_GWAS/logs/$pheno/UKB_${pheno}_part${part}_BOLTLMM.log \
					${pheno}-UKB${part}.txt &
			fi
		done
		wait
	else
		echo "Formatting $pheno full sample"
		format_UKB \
			$PGI_Repo/derived_data/1_UKB_GWAS/output/$pheno/UKB_${pheno}_BOLTLMM \
			$PGI_Repo/derived_data/1_UKB_GWAS/logs/$pheno/UKB_${pheno}_BOLTLMM.log \
			${pheno}-UKB.txt &
	fi
	done
}

main



