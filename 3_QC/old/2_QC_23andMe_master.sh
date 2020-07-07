#!/bin/bash

pathOut=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/3_QCd/23andMe
pathIn=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/23andMe
pathCode=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code/3_QC


easyQC(){
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
		--IMPUTED IMPUTED \
		--HWE HWE_PVAL \
		--CALLRATE CALLRATE \
		--cutoff_HWE 1e-20 \
		--cutoff_CALLRATE 0.95 \
		--cutoff_MAF 0.01 \
		--cutoff_INFO 0.7 \
		--XY 0 \
		--INDEL 0 \
		--cptref "/disk/genetics/ukb/aokbay/reffiles/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.rsid_map" \
		--afref "/disk/genetics/ukb/aokbay/reffiles/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.cptid.maf001.gz" \
		--snpStd 0 \
		--SDy 1 > $pathOut/QC_${tag}_$(date +"%Y_%m_%d").log &
	done < $1
	wait
}

qcPlots(){
	while read file; do
		tag=$(echo $file | rev | cut -f1 -d"/" | cut -d"." -f2 | rev )
		Rscript $pathCode/0.2_QCplots.R \
		--file $pathOut/QC_${tag}_$(date +"%Y_%m_%d") \
		--sdy 1 &
	done < $1
	wait
}

# $1: filelist, $2: EasyQC or QCplots, $3: Output tag (e.g. batch number)
checkStatus(){
	while read file; do
		tag=$(echo $file | rev | cut -f1 -d"/" | cut -d"." -f2 | rev )
		Nplots_EasyQC="3"
		Nplots_QCplots="5"
		if [[ -f $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/*.png ]];then
			Nplots=$(ls $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/*png | wc -l)
		else
			Nplots=0
		fi
		rm $pathCode/${2}_${3}_status
		if [[ $Nplots < $Nplots_${2} ]]; then
			echo "$tag $2 FALSE" >> $pathCode/${2}_${3}_status
			rerun=1
		else
			echo "$tag $2 TRUE" >> $pathCode/${2}_${3}_status
		fi
	done < $1
	return $rerun
}

rerun(){
	while read file; do
		tag=$(echo $file | rev | cut -f1 -d"/" | cut -d"." -f2 | rev )
		status=$(echo $row | cut -d" " -f3)
	
		if [[ $status == "FALSE" ]]; then
			echo "Rerunning $tag $1 .."
			grep $tag $pathCode/2.1_23andMe_sumstats_batch${2}.txt >> ${1}_${2}_retry
			if [[ $1 == "EasyQC" ]]; then
				rm -r $pathOut/QC_${tag}_$(date +"%Y_%m_%d")
			fi
		fi
	done < ${1}_${2}_status

	$analysis ${1}_${2}_retry
}


main(){
	#ls -1 $pathIn/23andMe* > $pathCode/2.1_23andMe_sumstats.txt
	## Add file tags manually as first field (file names are not self-explanatory)

	sed -n '1,10p' $pathCode/2.1_23andMe_sumstats.txt > $pathCode/2.1_23andMe_sumstats_batch1.txt
	sed -n '11,20p' $pathCode/2.1_23andMe_sumstats.txt > $pathCode/2.1_23andMe_sumstats_batch2.txt
	sed -n '21,31p' $pathCode/2.1_23andMe_sumstats.txt > $pathCode/2.1_23andMe_sumstats_batch3.txt
	sed -n '32,43p' $pathCode/2.1_23andMe_sumstats.txt > $pathCode/2.1_23andMe_sumstats_batch4.txt

	for batch in 1; do
		echo "-----------------------------------------------------"
		echo -n "Batch $batch started on "
		date
		echo ""
		start=$(date +%s)

		for analysis in easyQC qcPlots; do
			echo "-----------------------------------------------------"
			echo "Running $analysis on batch $batch.."
			echo ""
			$analysis $pathCode/2.1_23andMe_sumstats_batch$batch.txt 
			echo ""
			echo "$analysis on batch $batch finished."
			echo "Checking status.."
			checkStatus $pathCode/2.1_23andMe_sumstats_batch$batch.txt $analysis $batch
			if [[ $ == 1 ]]; then
				rerun ${analysis} ${batch}
			else
				echo "All files successfully finished."
			fi
		done
		echo "-----------------------------------------------------"		
		echo -n "Batch $batch finished running on "
		date
		end=$(date +%s)
		echo "Analysis took $(( ($end - $start)/60 )) minutes."
		echo "-----------------------------------------------------" 
	done

	echo "Script finished."
	rm $pathCode/2.1_23andMe_sumstats_batch*.txt 
}








## Check which ones are not finished, remove those QC dirs for rerun
## To do: Automate this to write list to 2.1_23andMe_sumstats_retry.txt and rerun



mkdir $pathOut/Results
rm $pathOut/Results/EasyQC_rep.txt
ls -d $pathOut/*/ | sed '/Results/d' > $pathCode/23andMe_sumstats_dirs
while read dir; do 
	cat ${dir}*.rep >> $pathOut/Results/EasyQC_rep.txt
	cp ${dir}*.png $pathOut/Results
done < $pathCode/23andMe_sumstats_dirs
rm $pathCode/23andMe_sumstats_dirs

