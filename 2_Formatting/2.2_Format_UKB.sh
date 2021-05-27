#!/bin/bash

source paths2

cd $p2_UKBout

format_UKB(){
	pheno=$1
	part=$2
	
	# Get N from log files
	N=$(grep -a "Number of indivs with no missing phenotype(s) to use:" $p2_UKBlog/$pheno/UKB_${pheno}_part${part}_BOLTLMM.log | cut -d" " -f10)

	# Add N and callrate
	awk -F"\t" -v N=$N 'NR==FNR{a[$2]=1-$5;next} \
	FNR==1{print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE"}
	($1 in a && FNR>1) {print $1,$2,$3,$5,$6,$7,$9,$10,$11,$8,N,"1",a[$1];next}
	(FNR>1) {print $1,$2,$3,$5,$6,$7,$9,$10,$11,$8,N,"1","NA"}' OFS="\t" $p2_UKBcallrate $p2_UKBin/$pheno/UKB_${pheno}_part${part}_BOLTLMM > ${pheno}-UKB${part}_tmp

	# Add HWE p-value
	awk -F"\t" 'NR==FNR{a[$2]=$9;next} \
	FNR==1{print $0, "HWE_PVAL"} \
	($1 in a && FNR>1) {print $0,a[$1];next} \
	(FNR>1) {print $0,"NA"}' OFS="\t" $p2_UKBhwe ${pheno}-UKB${part}_tmp > ${pheno}-UKB${part}.txt

	rm ${pheno}-UKB${part}_tmp
}

main() {
	phenos=$(ls -d $p2_UKBin/* | rev | cut -d"/" -f1 | rev | sed '/^old/d')
	for part in 1 2 3; do
		for pheno in $phenos; do
			echo "Formatting $pheno part $part"
			format_UKB $pheno $part &
		done
		wait
	done
}

main



