#!/bin/bash
source $PGI_Repo/code/paths
export R_LIBS=$HOME/R/x86_64-redhat-linux-gnu-library/4.3/:$R_LIBS

# Get correlation between paternal-maternal / parental and proband PGIs, and snipar-proband and regular proband PGIs
for cohort in AH ALSPAC ERisk GS GSOEP MCTFR MCS STRtwge STRgsa STRpsych Texas WLS UKB1
do
    eval pgiDir='$'pgi_dir_${cohort}/single_SBayesR
    Rscript $PGI_Repo/code/13_Parental_imputation/13.2.0_PGIcorrelations.R $cohort $pgiDir $PGI_Repo/derived_data/13_Parental_imputation/QC/correlations
done

# Get list of SNPs in hdf5 filesc
for cohort in AH ALSPAC ERisk GS GSOEP MCTFR MCS STRtwge STRgsa STRpsych Texas WLS UKB1
do
    if [[ $cohort == "UKB1" ]]
    then 
        eval gfDir='$'gf_dir_${cohort}/parental
    else
        eval gfDir='$'gf_dir_${cohort}/parental/output
    fi
    mkdir -p $PGI_Repo/derived_data/13_Parental_imputation/QC/SNPlists/$cohort
    Rscript $PGI_Repo/code/13_Parental_imputation/13.2.1_get_imputed_SNPs.R $cohort $gfDir $PGI_Repo/derived_data/13_Parental_imputation/QC/SNPlists/$cohort
done

# Get number of SNPs overlapping with weights
for cohort in AH ALSPAC ERisk GS GSOEP MCTFR MCS STRtwge STRgsa STRpsych Texas WLS UKB1
do
    rm -f $PGI_Repo/derived_data/13_Parental_imputation/QC/SNPlists/${cohort}_weights_overlapping.snps
    if [[ $cohort == "UKB1" ]]
    then
        snpidcol=2
    elif [[ $cohort == "STRgsa" || $cohort == "MCS" ]]
    then
        snpidcol=13
    else 
        snpidcol=12
    fi

    while read pheno 
    do
        if [[ $cohort == "STRgsa" || $cohort == "MCS" ]]
        then
            weights="$PGI_Repo/derived_data/9_Scores/single_SBayesR/weights/${pheno}-single_weights_SBayesR.txt.hg38"
        else
            weights="$PGI_Repo/derived_data/9_Scores/single_SBayesR/weights/${pheno}-single_weights_SBayesR.txt"
        fi
        Nsnps=$(awk -v snpid=$snpidcol 'NR==FNR{a[$snpid]=$snpid;next}($2 in a){x++}END{print x}' $weights $PGI_Repo/derived_data/13_Parental_imputation/QC/SNPlists/$cohort/$cohort.bim)
        echo -e "$pheno\t$Nsnps" >> $PGI_Repo/derived_data/13_Parental_imputation/QC/SNPlists/${cohort}_weights_overlapping.snps
    done < $PGI_Repo/code/9_Scores/version_single_$cohort
done

