#!/bin/bash

source $PGI_Repo/code/5_LDSC/5.0_LDSCfunctions.sh
cd $PGI_Repo/derived_data/5_LDSC/study_level

LDSC_study() {
    fileList=$1
    SNPlist=$2

    for analysis in munge h2 rg; do
        checkStatusLDSC $fileList $analysis

        if [[ $status == 1 ]]; then
            echo "Rerunning study-level analysis for the unfinished files in $fileList.."
            $analysis ${fileList}.${analysis}.rerun $SNPlist
            checkStatusLDSC $fileList $analysis
        fi
    done
}

###############################################################################

main(){
    LDSC_study $PGI_Repo/code/4_MTAG_single/singleMTAG_input_filelist.txt HM3
    LDSC_h2_stats $PGI_Repo/code/4_MTAG_single/singleMTAG_input_filelist.txt
    LDSC_rg_stats $PGI_Repo/code/4_MTAG_single/singleMTAG_input_filelist.txt
    study_h2_table $PGI_Repo/code/4_MTAG_single/singleMTAG_input_filelist.txt
    study_rg_table $PGI_Repo/code/4_MTAG_single/singleMTAG_input_filelist.txt
}

main
