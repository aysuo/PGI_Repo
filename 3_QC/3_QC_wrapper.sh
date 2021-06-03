#!/bin/bash

source paths3

cd $p3_QC
export R_LIBS=$p3_Rlib/:$R_LIBS

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
			outDir=$p3_QC/$cohort
		else
			mkdir -p $p3_QC/$cohort/SEfilter
			cutoff_SE=2
			SDy=$(sed -n '/^SE_reg/p' $p3_QC/$cohort/QC_${pheno}_*/CLEANED.*.log  | awk '{print $2}')
			outDir=$p3_QC/$cohort/SEfilter
		fi

		rm -r -f $outDir/QC_${pheno}_*
		
		nohup sh $p3_code/3.1_EasyQC.sh \
		--fileIn $file \
		--sep TAB \
		--miss NA \
		--p3_QC $outDir \
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
		--cptref $p3_cptref \
		--afref $p3_afref \
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
			outDir=$p3_QC/$cohort/SEfilter/
		else
			outDir=$p3_QC/$cohort
		fi

		file=$(echo $outDir/QC_${pheno}_* | cut -d" " -f1)

		echo "Processing $pheno.."
		Rscript $p3_code/3.2_QCplots.R \
		--file $file \
		--ref $p3_qcplotsRef \
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
	
	rm -f $p3_code/${cohort}_${analysis}_status
	rm -f $p3_code/${cohort}_${analysis}_rerun
	
	status=1
	while read file; do
		pheno=$(echo $file | rev | cut -f1 -d"/" | cut -d"." -f3 | rev )
		targetN_easyQC="3"
		targetN_qcPlots="6"
		eval targetN='$'{targetN_${analysis}}

		if [[ $SEfilter == 0 ]]; then
			inputDir="$p3_QC/$cohort"
		else
			inputDir="$p3_QC/$cohort/SEfilter"
		fi

		if ls $inputDir/QC_${pheno}_*/*png 1> /dev/null 2>&1; then
			Nplots=$(ls $inputDir/QC_${pheno}_*/*png | wc -l)
		else
			Nplots=0
		fi
		
		if [[ $Nplots < $targetN ]]; then
			grep $pheno $fileList >> $p3_code/${cohort}_${analysis}_rerun
			echo "$analysis for $cohort $pheno was unsuccessful."
			status=0
		fi

	done < $fileList

	if [[ -f $p3_code/${cohort}_${analysis}_rerun ]]; then
		mv $p3_code/${cohort}_${analysis}_rerun ${fileList}
	fi
}


QC(){
	cohort=$1
	inputDir=$p3_formatted/$cohort

	echo "----------------------------------------------------------------------"
	echo -n "QC on $cohort started at"
	date
	echo ""
	start=$(date +%s)

	# First run without SE filter, use those results to estimate SDy with qcPlots()
	# Then rerun with SE filter, using estimated SDy to calculate predicted SEs
	for SEfilter in 0 1; do
		ls -d $p3_formatted/$cohort/*.gz > $p3_code/${cohort}_inputFileList
		pass=1
		status=0

		while [[ $status == 0 ]]; do
			echo "--------------------------------------------------"
			echo "EasyQC pass $pass.."

			checkStatus $p3_code/${cohort}_inputFileList $cohort easyQC $SEfilter
			status=$status

			if [[ $status == 0 ]]; then
				easyQC $p3_code/${cohort}_inputFileList $cohort $SEfilter
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

		ls -d $p3_formatted/$cohort/*.gz > $p3_code/${cohort}_inputFileList
		pass=1
		status=0

		while [[ $status == 0 ]]; do
			echo "qcPlots pass $pass :"
			
			checkStatus $p3_code/${cohort}_inputFileList $cohort qcPlots $SEfilter
			status=$status
		
			if [[ $status == 0 ]]; then
				qcPlots $p3_code/${cohort}_inputFileList $cohort $SEfilter
				pass=$(($pass+1))
			fi

			if [[ $pass > 3 ]]; then
				break
			fi
		done

		if [[ $status == 0  ]]; then
			echo "qcPlots could not be completed for $cohort - SE filter=$SEfilter . Check for errors in input files."
		fi

		rm $p3_code/${cohort}_inputFileList

		if [[ $SEfilter == 0 ]]; then
			outDir=$p3_QC/$cohort
		else
			outDir=$p3_QC/$cohort/SEfilter
		fi

		# Copy QC plots and *.rep files into Results directory 
		rm -r -f $outDir/Results
		mkdir $outDir/Results
		ls -d $outDir/QC_* > $p3_code/${cohort}_SE${SEfilter}_QCdirs
		while read dir; do
			if [[ -d $dir ]]; then
				cat ${dir}/*.rep >> $outDir/Results/${cohort}_SE${SEfilter}_easyQC_rep.txt
				cp ${dir}/*.png $outDir/Results
			fi
		done < $p3_code/${cohort}_SE${SEfilter}_QCdirs
		rm $p3_code/${cohort}_SE${SEfilter}_QCdirs


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
	QC 23andMe 
	QC UKB
	QC public
	QC public_scores
}

main