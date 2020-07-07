#!/bin/bash
dirCode=/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/code
dirData=/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data

###############################################################################

source $dirCode/5_LDSC/5.0_LDSCfunctions.sh

LDSC_mtag() {
    fileList=$1
    dirIn=$2
    dirOut=$3
    single=$4

    cd $dirOut
        
    for analysis in munge h2; do
        checkStatus $fileList $analysis
    
        if [[ $status == 1 ]]; then
            echo "Rerunning LDSC for the unfinished files in $fileList.."
            $analysis ${fileList}.${analysis}.rerun
            checkStatus $fileList $analysis
        fi
    done

    ER2_table $fileList $dirIn $single

    if [[ $single == 1 ]]; then
        checkStatus $dirCode/5_LDSC/rg_meta_filelist rg_meta
        rg $dirCode/5_LDSC/rg_meta_filelist.rg_meta.rerun
        rg_table $dirCode/5_LDSC/rg_meta_filelist
    fi
       
}

###############################################################################

main(){
    # Single-trait MTAG
    #rm -f $dirCode/5_LDSC/singleMTAG_output_filelist.txt
    #for phenodir in ${dirData}/4_MTAG_single/*; do
    #    if [[ $phenodir == *1 ]]; then
    #        pheno=$(echo $phenodir | rev | cut -d"/" -f1 | rev)
    #        ss=$(echo $phenodir/*trait*_formatted*.txt | sed 's/ /,/g')
    #        echo -e "$pheno\t$ss" >> $dirCode/5_LDSC/singleMTAG_output_filelist.txt
    #    fi
    #done

    #LDSC_mtag $dirCode/5_LDSC/singleMTAG_output_filelist.txt $dirData/4_MTAG_single $dirData/5_LDSC/singleMTAG 1

    rm -f $dirCode/5_LDSC/multiMTAG_output_filelist.txt
    for phenodir in ${dirData}/6_MTAG_multi/*; do
        if [[ $phenodir == *1 ]]; then
            pheno=$(echo $phenodir | rev | cut -d"/" -f1 | rev)
            ss=$(echo $phenodir/*trait_1_formatted.txt)
            echo -e "$pheno\t$ss" >> $dirCode/5_LDSC/multiMTAG_output_filelist.txt
        fi
    done

    LDSC_mtag $dirCode/5_LDSC/multiMTAG_output_filelist.txt $dirData/6_MTAG_multi $dirData/5_LDSC/multiMTAG 0
}

main

