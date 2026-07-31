#!/bin/bash

source $PGI_Repo/code/5_LDSC/5.0_LDSCfunctions.sh

LDSC_mtag() {
    fileList=$1
    dirIn=$2
    MTAGtype=$3 # single or multi
    SNPlist=$4  # HM3 or SBayesR_2pt8M
    ER2table=$5 # ER2table = "full" if making table for the E(R2) vs R2 analysis, "normal" otherwise

    dirOut=$PGI_Repo/derived_data/5_LDSC/${MTAGtype}MTAG/$SNPlist

    cd $dirOut
        
    for analysis in munge h2; do
        checkStatusLDSC $fileList $analysis
    
        if [[ $status == 1 ]]; then
            echo "Rerunning LDSC for the unfinished files in $fileList.."
            $analysis ${fileList}.${analysis}.rerun $SNPlist 
            checkStatusLDSC $fileList $analysis
        fi
    done

    if [[ $ER2table == "normal" ]]; then
        ER2_table $fileList $dirIn $MTAGtype
    else
        ER2_table_full $fileList $dirIn $MTAGtype
    fi

    if [[ $ER2table == "normal" ]] && [[ $MTAGtype == "single" ]]; then
        checkStatusLDSC $PGI_Repo/code/5_LDSC/rg_meta_filelist rg_meta
        rg $PGI_Repo/code/5_LDSC/rg_meta_filelist.rg_meta.rerun $SNPlist
        rg_table $PGI_Repo/code/5_LDSC/rg_meta_filelist
    fi 
}

###############################################################################

main(){
    # Estimate h^2, rg for single-trait MTAG output, write E(R^2) and rg tables.
    rm -f $PGI_Repo/code/5_LDSC/singleMTAG_output_filelist.txt
    for phenodir in $PGI_Repo/derived_data/4_MTAG_single/*; do
        # Use largest sample size version (v1, e.g. NEURO1) for each phenotype
        if [[ $phenodir == *1 ]]; then
            pheno=$(echo $phenodir | rev | cut -d"/" -f1 | rev)
            ss=$(echo $phenodir/*trait*_formatted_SBayesR.txt | sed 's/ /,/g')
            echo -e "$pheno\t$ss" >> $PGI_Repo/code/5_LDSC/singleMTAG_output_filelist.txt
        fi
    done

    LDSC_mtag $PGI_Repo/code/5_LDSC/singleMTAG_output_filelist.txt $PGI_Repo/derived_data/4_MTAG_single single HM3 normal
    LDSC_mtag $PGI_Repo/code/5_LDSC/singleMTAG_output_filelist.txt $PGI_Repo/derived_data/4_MTAG_single single SBayesR_2pt8M normal

    # E(R^2) - observed R^2 comparison table for single-trait PGIs
    # Calculate E(R^2) based on all GWAS used to make PGIs for validation cohorts and largest h^2 MTAG output
    rm -f $PGI_Repo/code/5_LDSC/singleMTAG_output_filelist.txt
    versions=$(cat $PGI_Repo/code/9_Scores/version_single_* | sort | uniq)
    for pheno in $versions; do
        if [[ $pheno != *1 ]]; then
            ss=$(echo $PGI_Repo/derived_data/4_MTAG_single/$pheno/*trait*_formatted_SBayesR*.txt | sed 's/ /,/g')
            echo -e "$pheno\t$ss" >> $PGI_Repo/code/5_LDSC/singleMTAG_output_filelist.txt
        fi
    done

    LDSC_mtag $PGI_Repo/code/5_LDSC/singleMTAG_output_filelist.txt $PGI_Repo/derived_data/4_MTAG_single single HM3 full
    LDSC_mtag $PGI_Repo/code/5_LDSC/singleMTAG_output_filelist.txt $PGI_Repo/derived_data/4_MTAG_single single SBayesR_2pt8M full

    rm -f $PGI_Repo/code/5_LDSC/singleMTAG_output_filelist.txt
}

main

