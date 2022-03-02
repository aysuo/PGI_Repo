#!/bin/bash

source $PGI_Repo/code/5_LDSC/5.0_LDSCfunctions.sh

LDSC_mtag() {
    fileList=$1
    dirIn=$2
    dirOut=$3
    # single = 1 if single-trait MTAG, 0 if multi-trait
    single=$4
    # ER2table = "full" if making table for the E(R2) vs R2 analysis, "normal" otherwise
    ER2table=$5

    cd $dirOut
        
    for analysis in munge h2; do
        checkStatusLDSC $fileList $analysis
    
        if [[ $status == 1 ]]; then
            echo "Rerunning LDSC for the unfinished files in $fileList.."
            $analysis ${fileList}.${analysis}.rerun
            checkStatusLDSC $fileList $analysis
        fi
    done

    if [[ $ER2table == "normal" ]]; then
        ER2_table $fileList $dirIn $single
    else
        ER2_table_full $fileList $dirIn $single
    fi

    if [[ $ER2table == "normal" ]] && [[ $single == 1 ]]; then
        checkStatusLDSC $PGI_Repo/code/5_LDSC/rg_meta_filelist rg_meta
        rg $PGI_Repo/code/5_LDSC/rg_meta_filelist.rg_meta.rerun
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
            ss=$(echo $phenodir/*trait*_formatted*.txt | sed 's/ /,/g')
            echo -e "$pheno\t$ss" >> $PGI_Repo/code/5_LDSC/singleMTAG_output_filelist.txt
        fi
    done

    LDSC_mtag $PGI_Repo/code/5_LDSC/singleMTAG_output_filelist.txt $PGI_Repo/derived_data/4_MTAG_single $PGI_Repo/derived_data/5_LDSC/singleMTAG 1 normal

    # E(R^2) - observed R^2 comparison table for single-trait PGIs
    # Calculate E(R^2) based on all GWAS used to make PGIs for validation cohorts and largest h^2 MTAG output
    versions=$(cat $PGI_Repo/code/9_Scores/version_single_* | sort | uniq)
    for pheno in $versions; do
        if [[ $pheno != *1 ]]; then
            ss=$(echo $PGI_Repo/derived_data/4_MTAG_single/$pheno/*trait*_formatted*.txt | sed 's/ /,/g')
            echo -e "$pheno\t$ss" >> $PGI_Repo/code/5_LDSC/singleMTAG_output_filelist.txt
        fi
    done

    LDSC_mtag $PGI_Repo/code/5_LDSC/singleMTAG_output_filelist.txt $PGI_Repo/derived_data/4_MTAG_single $PGI_Repo/derived_data/5_LDSC/singleMTAG 1 full

    # Estimate h^2 for multi-trait MTAG output, write E(R^2) table
    rm -f $PGI_Repo/code/5_LDSC/multiMTAG_output_filelist.txt
    for phenodir in $PGI_Repo/derived_data/6_MTAG_multi/*; do
        if [[ $phenodir == *1 ]]; then
            pheno=$(echo $phenodir | rev | cut -d"/" -f1 | rev)
            ss=$(echo $phenodir/*trait_1_formatted.txt)
            echo -e "$pheno\t$ss" >> $PGI_Repo/code/5_LDSC/multiMTAG_output_filelist.txt
        fi
    done

    LDSC_mtag $PGI_Repo/code/5_LDSC/multiMTAG_output_filelist.txt $PGI_Repo/derived_data/6_MTAG_multi $PGI_Repo/derived_data/5_LDSC/multiMTAG 0 normal

    # E(R^2) - observed R^2 comparison table for multi-trait PGIs
    versions=$(cat $PGI_Repo/code/9_Scores/version_multi_* | sort | uniq)
    for pheno in $versions; do
        if [[ $pheno != *1 ]]; then
            ss=$(echo $PGI_Repo/derived_data/6_MTAG_multi/$pheno/*trait_1_formatted.txt)
            echo -e "$pheno\t$ss" >> $PGI_Repo/code/5_LDSC/multiMTAG_output_filelist.txt
        fi
    done

    LDSC_mtag $PGI_Repo/code/5_LDSC/multiMTAG_output_filelist.txt $PGI_Repo/derived_data/6_MTAG_multi $PGI_Repo/derived_data/5_LDSC/multiMTAG 0 full
}

main

