#!/bin/bash

pathOut=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/3_QCd/public
pathIn=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public
pathCode=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code/3_QC

ls -1 $pathIn/*.txt > $pathCode/1.1_public_sumstats.txt
sed -n '1,14p' $pathCode/1.1_public_sumstats.txt > $pathCode/1.1_public_sumstats_batch1.txt
sed -n '15,27p' $pathCode/1.1_public_sumstats.txt > $pathCode/1.1_public_sumstats_batch2.txt

for batch in 1 2; do
	while read file; do
		tag=$(echo $file | rev | cut -f1 -d"/" | cut -d"." -f2 | rev )
		nohup sh $pathCode/0.1_EasyQC.sh \
		--fileIn $file \
		--sep TAB \
		--miss NA \
		--pathOut $pathOut \
		--tag $tag \
		--SNPID SNPID \
		--SNPIDtype rs \
		--EA EFFECT_ALLELE \
		--OA OTHER_ALLELE \
		--EAF EAF \
		--EFFECT BETA \
		--SE SE \
		--P P \
		--N N \
		--EFFECTtype BETA \
		--CHR CHR \
		--BP BP \
		--INFO INFO \
		--cutoff_MAF 0.01 \
		--cutoff_INFO 0.7 \
		--XY 0 \
		--INDEL 0 \
		--cptref "/disk/genetics/ukb/aokbay/reffiles/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.rsid_map" \
		--afref "/disk/genetics/ukb/aokbay/reffiles/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.cptid.maf001.gz" \
		--snpStd 0 \
		--SDy 1 > $pathOut/QC_${tag}_$(date +"%Y_%m_%d").log &
	done < $pathCode/1.1_public_sumstats_batch$batch.txt
	wait

	while read file; do
		tag=$(echo $file | rev | cut -f1 -d"/" | cut -d"." -f2 | rev )
		Rscript $pathCode/0.2_QCplots.R \
		--file $pathOut/QC_${tag}_$(date +"%Y_%m_%d") \
		--sdy 1 &
	done < $pathCode/1.1_public_sumstats_batch$batch.txt
	wait
done


#mkdir $pathOut/Results
#rm $pathOut/Results/EasyQC_rep.txt
#ls -d $pathOut/*/ | sed '/Results/d' | sed '/OLD/d'> public_sumstats_dirs
#while read dir; do 
#	cat ${dir}*.rep >> $pathOut/Results/EasyQC_rep.txt
#	cp ${dir}*.png $pathOut/Results
#done < public_sumstats_dirs
#rm public_sumstats_dirs


################################################################################################3

## Get estimated SD of phenotype from regressions on 1/sqrt(2*N*MAF*(1-MAF)) and apply SE ratio filter
mkdir $pathOut/SEfilter

for batch in 1 2; do
	while read file; do
		tag=$(echo $file | rev | cut -f1 -d"/" | cut -d"." -f2 | rev )
		SD=$(sed -n '/^SE_reg/p' $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/CLEANED.*.log  | awk '{print $2}')
		nohup sh $pathCode/0.1_EasyQC.sh \
		--fileIn $file \
		--sep TAB \
		--miss NA \
		--pathOut $pathOut/SEfilter \
		--tag $tag \
		--SNPID SNPID \
		--SNPIDtype rs \
		--EA EFFECT_ALLELE \
		--OA OTHER_ALLELE \
		--EAF EAF \
		--EFFECT BETA \
		--SE SE \
		--P P \
		--N N \
		--EFFECTtype BETA \
		--CHR CHR \
		--BP BP \
		--INFO INFO \
		--cutoff_MAF 0.01 \
		--cutoff_INFO 0.7 \
		--XY 0 \
		--INDEL 0 \
		--cptref "/disk/genetics/ukb/aokbay/reffiles/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.rsid_map" \
		--afref "/disk/genetics/ukb/aokbay/reffiles/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.cptid.maf001.gz" \
		--snpStd 0 \
		--SDy $SD \
		--cutoff_SE 2 > $pathOut/SEfilter/QC_${tag}_$(date +"%Y_%m_%d").log &
	done < $pathCode/1.1_public_sumstats_batch$batch.txt
	wait

	while read file; do
		tag=$(echo $file | rev | cut -f1 -d"/" | cut -d"." -f2 | rev )
		SD=$(sed -n '/^SE_reg/p' $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/CLEANED.*.log  | awk '{print $2}')
		Rscript $pathCode/0.2_QCplots.R \
		--file $pathOut/SEfilter/QC_${tag}_$(date +"%Y_%m_%d") \
		--sdy $SD &
	done < $pathCode/1.1_public_sumstats_batch$batch.txt
	wait
done

mkdir $pathOut/SEfilter/Results
rm $pathOut/SEfilter/Results/EasyQC_rep.txt
ls -d $pathOut/SEfilter/*/ | sed '/Results/d' | sed '/OLD/d'> public_sumstats_SEfilter_dirs
while read dir; do 
	cat ${dir}*.rep >> $pathOut/SEfilter/Results/EasyQC_rep.txt
	cp ${dir}*.png $pathOut/SEfilter/Results
done < public_sumstats_SEfilter_dirs
rm public_sumstats_SEfilter_dirs



