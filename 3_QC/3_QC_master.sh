#!/bin/bash

pathOut=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/3_QCd
pathIn=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted
pathCode=/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code/3_QC

cd $pathOut

# $1: filelist, $2: cohort $3: SE filter
easyQC(){
	fileList=$1
	cohort=$2
	SEfilter=$3

	echo "-----------------------------------------------------"
	echo "Running EasyQC on $cohort with SE filter = $SEfilter."
	echo ""

	i=0
	while read file; do
		tag=$(echo $file | rev | cut -f1 -d"/" | cut -d"." -f3 | rev )

		echo "Processing $tag.."
		if [[ $SEfilter == 0 ]]; then
			cutoff_SE=50
			SDy=1
			outDir=$pathOut/$cohort
		else
			mkdir -p $pathOut/$cohort/SEfilter
			cutoff_SE=2
			SDy=$(sed -n '/^SE_reg/p' $pathOut/$cohort/QC_${tag}_*/CLEANED.*.log  | awk '{print $2}')
			outDir=$pathOut/$cohort/SEfilter
		fi

		rm -r -f $outDir/QC_${tag}_*
		
		nohup sh $pathCode/3.1_EasyQC.sh \
		--fileIn $file \
		--sep TAB \
		--miss NA \
		--pathOut $outDir \
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
		--cutoff_SE $cutoff_SE \
		--SDy $SDy > $outDir/QC_${tag}_$(date +"%Y_%m_%d").log &

		let i+=1
		
		if [[ $i == 10 ]]; then
			wait
			i=0
		fi

	done < $fileList
	wait

	echo "Done."
	echo ""
	echo "-----------------------------------------------------"
}


qcPlots(){
	fileList=$1
	cohort=$2
	SEfilter=$3

	echo "-----------------------------------------------------"
	echo "Getting QC plots  $cohort for with SE filter = $SEfilter."
	echo ""

	i=0
	while read file; do
		tag=$(echo $file | rev | cut -f1 -d"/" | cut -d"." -f3 | rev )
		if [[ $SEfilter == "1" ]]; then
			outDir=$pathOut/$cohort/SEfilter/
		else
			outDir=$pathOut/$cohort
		fi

		file=$(echo $outDir/QC_${tag}_* | cut -d" " -f1)

		echo "Processing $tag.."
		Rscript $pathCode/3.2_QCplots.R \
		--file $file \
		--sdy 1 >> $outDir/QC_${tag}_$(date +"%Y_%m_%d").log &

		let i+=1
		
		if [[ $i == 10 ]]; then
			wait
			i=0
		fi

	done < $fileList
	wait

	echo "Done."
	echo ""
	echo "-----------------------------------------------------"
}

# $1: filelist, $2: cohort, $3:analysis, $4:SEfilter
checkStatus(){
	fileList=$1
	cohort=$2
	analysis=$3
	SEfilter=$4

	echo "Checking status.."
	
	rm -f $pathCode/${cohort}_${analysis}_status
	rm -f $pathCode/${cohort}_${analysis}_rerun
	
	status=1
	while read file; do
		tag=$(echo $file | rev | cut -f1 -d"/" | cut -d"." -f3 | rev )
		targetN_easyQC="3"
		targetN_qcPlots="5"
		eval targetN='$'{targetN_${analysis}}

		if [[ $SEfilter == 0 ]]; then
			inputDir="$pathOut/$cohort"
		else
			inputDir="$pathOut/$cohort/SEfilter"
		fi

		if ls $inputDir/QC_${tag}_*/*png 1> /dev/null 2>&1; then
			Nplots=$(ls $inputDir/QC_${tag}_*/*png | wc -l)
		else
			Nplots=0
		fi
		
		if [[ $Nplots < $targetN ]]; then
			grep $tag $fileList >> $pathCode/${cohort}_${analysis}_rerun
			echo "$analysis for $cohort $tag was unsuccessful."
			status=0
		fi

	done < $fileList

	if [[ -f $pathCode/${cohort}_${analysis}_rerun ]]; then
		mv $pathCode/${cohort}_${analysis}_rerun ${fileList}
	fi
}

#1 cohort
QC(){
	cohort=$1
	inputDir=$pathIn/$cohort

	echo "----------------------------------------------------------------------"
	echo -n "QC on $cohort started at"
	date
	echo ""
	start=$(date +%s)


	for SEfilter in 0 1; do
		ls -d $pathIn/$cohort/*.gz > $pathCode/${cohort}_inputFileList
		pass=1
		status=0

		while [[ $status == 0 ]]; do
			echo "--------------------------------------------------"
			echo "EasyQC pass $pass.."

			checkStatus $pathCode/${cohort}_inputFileList $cohort easyQC $SEfilter
			status=$status

			if [[ $status == 0 ]]; then
				easyQC $pathCode/${cohort}_inputFileList $cohort $SEfilter
				pass=$(($pass+1))
			fi

			if [[ $pass > 1 ]]; then
				break
			fi

		done

		if [[ $status == 0  ]]; then
			echo "easyQC could not be completed for $cohort - SE filter=$SEfilter . Check for errors in input files."
		fi
		
		echo ""
		echo "--------------------------------------------------"
		echo ""

		ls -d $pathIn/$cohort/*.gz > $pathCode/${cohort}_inputFileList
		pass=1
		status=0

		while [[ $status == 0 ]]; do
			echo "qcPlots pass $pass :"
			
			checkStatus $pathCode/${cohort}_inputFileList $cohort qcPlots $SEfilter
			status=$status
		
			if [[ $status == 0 ]]; then
				qcPlots $pathCode/${cohort}_inputFileList $cohort $SEfilter
				pass=$(($pass+1))
			fi

			if [[ $pass > 3 ]]; then
				break
			fi
		done

		if [[ $status == 0  ]]; then
			echo "qcPlots could not be completed for $cohort - SE filter=$SEfilter . Check for errors in input files."
		fi

		rm $pathCode/${cohort}_inputFileList

		if [[ $SEfilter == 0 ]]; then
			outDir=$pathOut/$cohort
		else
			outDir=$pathOut/$cohort/SEfilter
		fi

		# Copy QC plots and *.rep files into Results directory 
		rm -r -f $outDir/Results
		mkdir $outDir/Results
		ls -d $outDir/QC_* > $pathCode/${cohort}_SE${SEfilter}_QCdirs
		while read dir; do
			if [[ -d $dir ]]; then
				cat ${dir}/*.rep >> $outDir/Results/${cohort}_SE${SEfilter}_easyQC_rep.txt
				cp ${dir}/*.png $outDir/Results
			fi
		done < $pathCode/${cohort}_SE${SEfilter}_QCdirs
		rm $pathCode/${cohort}_SE${SEfilter}_QCdirs


	done	
	echo ""
	echo -n "Finished running QC on $cohort at"
	date

	end=$(date +%s)
	echo "Analysis took $(( ($end - $start)/60 )) minutes."
	echo "-----------------------------------------------------" 
	echo ""

}

main(){
	#QC 23andMe 
	#QC public
	QC public_scores
	#QC UKB
}

main