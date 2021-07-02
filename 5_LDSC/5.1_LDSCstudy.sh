#!/bin/bash

source $mainDir/code/5_LDSC/5.0_LDSCfunctions.sh

LDSC_study() {
    fileList=$1
    
    cd $mainDir/derived_data/5_LDSC/study_level

    for analysis in munge h2 rg; do
        checkStatusLDSC $fileList $analysis

        if [[ $status == 1 ]]; then
            echo "Rerunning study-level analysis for the unfinished files in $fileList.."
            $analysis ${fileList}.${analysis}.rerun
            checkStatusLDSC $fileList $analysis
        fi

    done

    LDSC_h2_stats $fileList
    LDSC_rg_stats $fileList
}

###############################################################################

main(){
    LDSC_study $mainDir/code/4_MTAG_single/singleMTAG_input_filelist.txt
}

main


   




   