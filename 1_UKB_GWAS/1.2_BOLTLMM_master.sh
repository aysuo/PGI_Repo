#!/bin/bash

########################################################
####### SCRIPT FOR RUNNING GWAS WITH BOLT LMM ##########
########################################################

#### Requirements
  ## BOLT LMM software (http://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/)
  ## Phenotype data in the form of [FID IID resid]
    # ID format should be the same across phenotype and genotype data
    # Missing phenotype should have -9 in resid column
  ## Genotype data in plink format that contain model SNPs (~500K) to be used in estimating GRM
  ## Imputed genotype data in the format of your choice to be used for association analysis

#### Documentation
  ## For detailed documentation on running bolt (especially with UKB data), refer to wiki (https://github.com/omeed-maghzian/ssgac/wiki/GWAS-in-UKB2-with-BOLT-LMM)




########################################################
###################### USER INPUT ######################
########################################################

#### Fill in the quotes below.

    ## Working directory
    WD="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/1_UKB_GWAS"
    cd $WD
    mkdir -p output
    mkdir -p logs

    ## path to bolt-lmm package
    bolt_path="/disk/genetics2/pub/software/bolt/BOLT-LMM_v2.3"

    ## path to genotype data in plink format that contain model SNPs
    genotype_plink_fam="/disk/genetics2/ukb/orig/UKBv2/bolt_modelsnps/ukb2_maf01_info60_rsq30_nodup_recoded_pruned_chr1.fam"
    genotype_plink_bed="/disk/genetics2/ukb/orig/UKBv2/bolt_modelsnps/ukb2_maf01_info60_rsq30_nodup_recoded_pruned_chr{1:22}.bed"
    genotype_plink_bim="/disk/genetics2/ukb/orig/UKBv2/bolt_modelsnps/ukb2_maf01_info60_rsq30_nodup_recoded_pruned_chr{1:22}.bim"

    ## path to imputed genotype data in the bgen/sample format (BGEN 1.2 format)
    genotype_bgen="/disk/genetics2/ukb/orig/UKBv2/decrypted/ukb_imp_chr{1:22}_v2.bgen"
    genotype_sample="/disk/genetics2/ukb/orig/UKBv2/linking/ukb_imp_v2.sample"

    ## path to phenotype data
    pheno_col="resid"

    ## QC filters and bolt options
    bgen_MAF=0.01 # MAF filter for bgen genotype data
    bgen_INFO=0.7 # INFO filter for bgen genotype data
    num_Threads=7 # number of cores to use


########################################################
################## END OF USER INPUT ###################
########################################################


# 1: pheno
run_BOLTLMM(){
  pheno=$1
  mkdir -p ${WD}/output/$pheno
  mkdir -p ${WD}/logs/$pheno
  
  echo "----------------------------------------------------------------------"

  echo "Running GWAS for phenotype: ${pheno}"

  echo -n "GWAS started at"
  date
  echo ""
  start=$(date +%s)

  
  for j in `seq 1 3` # for partitions 1 through 3
  do
    #### set file names
    pheno_file="input/UKB_${pheno}_part${j}.pheno"

    nohup ${bolt_path}/bolt --lmm \
      --LDscoresFile=${bolt_path}/tables/LDSCORE.1000G_EUR.tab.gz \
      --fam=${genotype_plink_fam} \
      --bed=${genotype_plink_bed} \
      --bim=${genotype_plink_bim} \
      --bgenFile=${genotype_bgen} \
      --sampleFile=${genotype_sample} \
      --numThreads=${num_Threads} \
      --bgenMinMAF=${bgen_MAF} \
      --bgenMinINFO=${bgen_INFO} \
      --phenoFile=${pheno_file} \
      --phenoCol=${pheno_col} \
      --statsFile=${WD}/output/${pheno}/UKB_${pheno}_part${j}_BOLTLMM_plink \
      --statsFileBgenSnps=${WD}/output/${pheno}/UKB_${pheno}_part${j}_BOLTLMM > ${WD}/logs/${pheno}/UKB_${pheno}_part${j}_BOLTLMM.log &
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
  pheno_list=$(ls $WD/input/*.pheno | rev | cut -d"/" -f1 | rev | sed 's/UKB_//g' | sed 's/_part[1-3]\.pheno//g' | sort | uniq)

  for pheno in $pheno_list; do
    run_BOLTLMM $pheno
  done
}

main


########################################################
########################################################

#### End of script