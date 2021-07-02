#!/bin/bash

source $mainDir/code/paths

cd $mainDir/derived_data/1_UKB_GWAS
mkdir -p output
mkdir -p logs

## pheno column name
pheno_col="resid"

## QC filters and bolt options
bgen_MAF=0.01 # MAF filter for bgen genotype data
bgen_INFO=0.7 # INFO filter for bgen genotype data
num_Threads=7 # number of cores to use


run_BOLTLMM(){
  pheno=$1
  mkdir -p output/$pheno
  mkdir -p logs/$pheno
  
  echo "----------------------------------------------------------------------"

  echo "Running GWAS for phenotype: ${pheno}"

  echo -n "GWAS started at"
  date
  echo ""
  start=$(date +%s)

  
  for j in 1 2 3 # for partitions 1 through 3
  do
    #### set file names
    pheno_file="input/UKB_${pheno}_part${j}.pheno"

    nohup ${bolt}/bolt --lmm \
      --LDscoresFile=${bolt}/tables/LDSCORE.1000G_EUR.tab.gz \
      --fam=${UKB_model_fam} \
      --bed=${UKB_model_bed} \
      --bim=${UKB_model_bim} \
      --bgenFile=${UKB_bgen} \
      --sampleFile=${UKB_bgen_sample} \
      --numThreads=${num_Threads} \
      --bgenMinMAF=${bgen_MAF} \
      --bgenMinINFO=${bgen_INFO} \
      --phenoFile=${pheno_file} \
      --phenoCol=${pheno_col} \
      --statsFile=output/${pheno}/UKB_${pheno}_part${j}_BOLTLMM_plink \
      --statsFileBgenSnps=output/${pheno}/UKB_${pheno}_part${j}_BOLTLMM > logs/${pheno}/UKB_${pheno}_part${j}_BOLTLMM.log &
  done
  wait

  echo ""
  echo -n "Script finished at"
  date
  echo ""
  
  end=$(date +%s)
  echo "Running GWAS for $pheno took $(( ($end - $start)/60 )) minutes."
  echo "-----------------------------------------------------" 
  echo ""
}  
  

main(){
  # Get list of phenotypes
  pheno_list=$(ls input/*.pheno | rev | cut -d"/" -f1 | rev | sed 's/UKB_//g' | sed 's/_part[1-3]\.pheno//g' | sort | uniq)

  for pheno in $pheno_list; do
    run_BOLTLMM $pheno
  done
}

main


########################################################
########################################################

#### End of script