#!/bin/bash

source $PGI_Repo/code/paths

cd $PGI_Repo/derived_data/3_QCd

# EasyQC SNP filtering
easyQC(){
	fileList=$1
	cohort=$2
	SEfilter=$3

	echo "-----------------------------------------------------"
	echo "Running EasyQC on $cohort with SE filter = $SEfilter."
	echo ""

	i=0
	while read file; do
		pheno=$(echo $file | rev | cut -f1 -d"/" | cut -d"." -f3 | rev )

		echo "Processing $pheno.."
		if [[ $SEfilter == 0 ]]; then
			cutoff_SE=100
			SDy=1
			outDir=$PGI_Repo/derived_data/3_QCd/$cohort
		else
			mkdir -p $PGI_Repo/derived_data/3_QCd/$cohort/SEfilter
			cutoff_SE=2
			SDy=$(sed -n '/^SE_reg/p' $PGI_Repo/derived_data/3_QCd/$cohort/QC_${pheno}_*/CLEANED.*.log  | awk '{print $2}')
			outDir=$PGI_Repo/derived_data/3_QCd/$cohort/SEfilter
		fi

		rm -r -f $outDir/QC_${pheno}_*
		
		nohup sh $PGI_Repo/code/3_QC/3.1_EasyQC.sh \
		--fileIn $file \
		--sep TAB \
		--miss NA \
		--pathOut $outDir \
		--tag $pheno \
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
		--cptref $HRC_EasyQC_cptref \
		--afref $HRC_EasyQC_afref \
		--snpStd 0 \
		--cutoff_SE $cutoff_SE \
		--SDy $SDy > $outDir/QC_${pheno}_$(date +"%Y_%m_%d").log &

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

# Predicted vs reported SE plot, per chromosome SE ratio plots
qcPlots(){
	fileList=$1
	cohort=$2
	SEfilter=$3

	echo "-----------------------------------------------------"
	echo "Getting QC plots  $cohort for with SE filter = $SEfilter."
	echo ""

	i=0
	while read file; do
		pheno=$(echo $file | rev | cut -f1 -d"/" | cut -d"." -f3 | rev )
		if [[ $SEfilter == "1" ]]; then
			outDir=$PGI_Repo/derived_data/3_QCd/$cohort/SEfilter/
		else
			outDir=$PGI_Repo/derived_data/3_QCd/$cohort
		fi

		file=$(echo $outDir/QC_${pheno}_* | cut -d" " -f1)

		echo "Processing $pheno.."
		Rscript $PGI_Repo/code/3_QC/3.2_QCplots.R \
		--file $file \
		--ref $HRC_qcplotsRef \
		--sdy 1 >> $outDir/QC_${pheno}_$(date +"%Y_%m_%d").log &

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

# Function to check if QC steps finished successfully (looks at number of plots produced)
# EasyQC step should result in 3 plots, qcPlots in another 3
# Writes unsuccessful ones into $fileList to be rerun  
checkStatus(){
	fileList=$1
	cohort=$2
	analysis=$3
	SEfilter=$4

	echo "Checking status.."
	
	rm -f $PGI_Repo/code/3_QC/${cohort}_${analysis}_status
	rm -f $PGI_Repo/code/3_QC/${cohort}_${analysis}_rerun
	
	status=1
	while read file; do
		pheno=$(echo $file | rev | cut -f1 -d"/" | cut -d"." -f3 | rev )
		targetN_easyQC="3"
		targetN_qcPlots="6"
		eval targetN='$'{targetN_${analysis}}

		if [[ $SEfilter == 0 ]]; then
			inputDir="$PGI_Repo/derived_data/3_QCd/$cohort"
		else
			inputDir="$PGI_Repo/derived_data/3_QCd/$cohort/SEfilter"
		fi

		if ls $inputDir/QC_${pheno}_*/*png 1> /dev/null 2>&1; then
			Nplots=$(ls $inputDir/QC_${pheno}_*/*png | wc -l)
		else
			Nplots=0
		fi
		
		if [[ $Nplots < $targetN ]]; then
			grep $pheno $fileList >> $PGI_Repo/code/3_QC/${cohort}_${analysis}_rerun
			echo "$analysis for $cohort $pheno was unsuccessful."
			status=0
		fi

	done < $fileList

	if [[ -f $PGI_Repo/code/3_QC/${cohort}_${analysis}_rerun ]]; then
		mv $PGI_Repo/code/3_QC/${cohort}_${analysis}_rerun ${fileList}
	fi
}


QC(){
	cohort=$1
	inputDir=$PGI_Repo/derived_data/2_Formatted/$cohort

	echo "----------------------------------------------------------------------"
	echo -n "QC on $cohort started at"
	date
	echo ""
	start=$(date +%s)

	# First run without SE filter, use those results to estimate SDy with qcPlots()
	# Then rerun with SE filter, using estimated SDy to calculate predicted SEs
	for SEfilter in 0 1; do
		ls -d $PGI_Repo/derived_data/2_Formatted/$cohort/*.gz > $PGI_Repo/code/3_QC/${cohort}_inputFileList
		pass=1
		status=0

		while [[ $status == 0 ]]; do
			echo "--------------------------------------------------"
			echo "EasyQC pass $pass.."

			checkStatus $PGI_Repo/code/3_QC/${cohort}_inputFileList $cohort easyQC $SEfilter
			status=$status

			if [[ $status == 0 ]]; then
				easyQC $PGI_Repo/code/3_QC/${cohort}_inputFileList $cohort $SEfilter
				pass=$(($pass+1))
			fi

			if [[ $pass > 3 ]]; then
				break
			fi

		done

		if [[ $status == 0  ]]; then
			echo "easyQC could not be completed for $cohort - SE filter=$SEfilter . Check for errors in input files."
		fi
		
		echo ""
		echo "--------------------------------------------------"
		echo ""

		ls -d $PGI_Repo/derived_data/2_Formatted/$cohort/*.gz > $PGI_Repo/code/3_QC/${cohort}_inputFileList
		pass=1
		status=0

		while [[ $status == 0 ]]; do
			echo "qcPlots pass $pass :"
			
			checkStatus $PGI_Repo/code/3_QC/${cohort}_inputFileList $cohort qcPlots $SEfilter
			status=$status
		
			if [[ $status == 0 ]]; then
				qcPlots $PGI_Repo/code/3_QC/${cohort}_inputFileList $cohort $SEfilter
				pass=$(($pass+1))
			fi

			if [[ $pass > 3 ]]; then
				break
			fi
		done

		if [[ $status == 0  ]]; then
			echo "qcPlots could not be completed for $cohort - SE filter=$SEfilter . Check for errors in input files."
		fi

		rm $PGI_Repo/code/3_QC/${cohort}_inputFileList

		if [[ $SEfilter == 0 ]]; then
			outDir=$PGI_Repo/derived_data/3_QCd/$cohort
		else
			outDir=$PGI_Repo/derived_data/3_QCd/$cohort/SEfilter
		fi

		# Copy QC plots and *.rep files into Results directory 
		rm -r -f $outDir/Results
		mkdir $outDir/Results
		ls -d $outDir/QC_* > $PGI_Repo/code/3_QC/${cohort}_SE${SEfilter}_QCdirs
		while read dir; do
			if [[ -d $dir ]]; then
				cat ${dir}/*.rep >> $outDir/Results/${cohort}_SE${SEfilter}_easyQC_rep.txt
				cp ${dir}/*.png $outDir/Results
			fi
		done < $PGI_Repo/code/3_QC/${cohort}_SE${SEfilter}_QCdirs
		rm $PGI_Repo/code/3_QC/${cohort}_SE${SEfilter}_QCdirs


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
	#QC UKB
	QC public
	#QC public_scores
}

main