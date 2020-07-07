#!/bin/bash
dirCode=/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/code
dirIn=/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data
dirOut=/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/5_LDSC/study_level

###############################################################################

source $dirCode/5_LDSC/5.0_LDSCfunctions.sh

LDSC_study() {
    fileList=$1
    
    cd $dirOut

    for analysis in munge h2 rg; do
        checkStatus $fileList $analysis

        if [[ $status == 1 ]]; then
            echo "Rerunning study-level analysis for the unfinished files in $fileList.."
            $analysis ${fileList}.${analysis}.rerun
            checkStatus $fileList $analysis
        fi

    done

    LDSC_h2_stats $fileList
    LDSC_rg_stats $fileList
}

###############################################################################

main(){
    LDSC_study $dirCode/4_MTAG_single/singleMTAG_input_filelist.txt
}

main


   




   