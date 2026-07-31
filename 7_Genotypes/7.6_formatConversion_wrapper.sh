#!/bin/bash

source $PGI_Repo/code/paths
source $PGI_Repo/code/7_Genotypes/7.6.0_formatConversion.sh

cohorts_ChrPosIDhg19=(ALSPAC AH1kG AH Dunedin EGCUT ELSA ERisk GS HRS MCTFR MIDUSomni10 MIDUSomni11 NCDS Texas PSID STRtwge STRpsych)
cohorts_ChrPosIDhg38=(MCTFR BCS70 STRgsa MCS NSHD NCDS)
cohorts_rsID=()
cohorts_vcf=(AH1kG AH Dunedin EGCUT ELSA ERisk GS HRS MCTFR MIDUSomni10 MIDUSomni11 NCDS Texas PSID STRtwge STRpsych MCTFR BCS70 STRgsa MCS NSHD)
cohorts_bgen=(ALSPAC)


for cohort in AH1kG AH ALSPAC BCS70 Dunedin EGCUT ELSA ERisk GS HRS MCS MCTFR MIDUS MIDUSomni10 MIDUSomni11 Texas PSID STRpsych STRtwge STRgsa WLS; do
do    
    eval gf_out='$'gf_dir_${cohort} 
    eval gf_orig='$'gf_orig_${cohort} 
    eval sample='$'sample_orig_${cohort}
    mkdir -p ${gf_out}/plink2
    mkdir -p ${gf_out}/plink/HM3

    if [[ " ${cohorts_ChrPosIDhg19[@]} " =~ " ${cohort} " ]]; then
        snpidtype=ChrPosIDhg19
    elif [[ " ${cohorts_ChrPosIDhg38[@]} " =~ " ${cohort} " ]]; then
        snpidtype=ChrPosIDhg38
    elif [[ " ${cohorts_rsID[@]} " =~ " ${cohort} " ]]; then
        snpidtype=rsID
    fi

    eval hm3='$'HM3_${snpidtype}

    if [[ " ${cohorts_vcf[@]} " =~ " ${cohort} " ]]; then
        vcf2plink2 "${gf_orig}" ${gf_out}/plink2/$cohort $snpidtype "$sample"
    elif [[ " ${cohorts_bgen[@]} " =~ " ${cohort} " ]]; then
        bgen2plink2 "${gf_orig}" ${gf_out}/plink2/$cohort $snpidtype "$sample"
    fi

    subsetSNPs ${gf_out}/plink2/${cohort}_chr[1:22] $hm3 ${gf_out}/plink/HM3/${cohort}_HM3
    mergePlink ${gf_out}/plink/HM3/${cohort}_HM3_chr[1:22] ${gf_out}/plink/HM3/${cohort}_HM3
done

# Process 1000G reference panel
vcf2plink2 "${gf_orig_1000G}" $PGI_Repo/derived_data/7_Genotypes/1000G/plink2/1000Gph3 rsID
subsetSNPs  $PGI_Repo/derived_data/7_Genotypes/1000G/plink2/1000Gph3_chr[1:22] $HM3_rsID $PGI_Repo/derived_data/7_Genotypes/1000G/plink/HM3/1000Gph3_HM3
mergePlink $PGI_Repo/derived_data/7_Genotypes/1000G/plink/HM3/1000Gph3_HM3_chr[1:22] $PGI_Repo/derived_data/7_Genotypes/1000G/plink/HM3/1000Gph3_HM3
rs2chrpos $PGI_Repo/derived_data/7_Genotypes/1000G/plink/HM3/1000Gph3_HM3 $PGI_Repo/derived_data/7_Genotypes/1000G/plink/HM3/1000Gph3_HM3_chrpos
liftOver $PGI_Repo/derived_data/7_Genotypes/1000G/plink/HM3/1000Gph3_HM3_chrpos $PGI_Repo/derived_data/7_Genotypes/1000G/plink/HM3/1000Gph3_HM3_chrpos_hg38 $chain_hg19toHg38

# NCDS and MIDUS have special format conversion scripts because they are large and require special handling. Run those scripts separately.
sh $PGI_Repo/code/7_Genotypes/7.6.0.0_formatConversion_NCDS.sh
sh $PGI_Repo/code/7_Genotypes/7.6.0.1_formatConversion_MIDUS.sh

## UKB relatives (from before the switch to RAP)
# UKB relatives - plink1 format
# echo "
# library(rhdf5)
# pedigree <- h5read('$gf_dir_UKB1/parental/chr_1.hdf5', '/pedigree') 
# pedigree <- as.data.frame(t(pedigree))
# write.table(pedigree, '$gf_dir_UKB1/degree1_relatives_plink/UKBrel_pedigree.txt' ,row.names=F, quote=F)" > $PGI_Repo/code/7_Genotypes/pedigree.R
# Rscript $PGI_Repo/code/7_Genotypes/pedigree.R

# # Get list of individuals in hdf5
# awk 'NR>2 {print $2, $2} $5=="True" {print $3, $3} $6=="True" {print $4, $4}' $gf_dir_UKB1/degree1_relatives_plink/UKBrel_pedigree.txt | sort | uniq > $gf_dir_UKB1/degree1_relatives_plink/UKB_degree1.sample
# rm $PGI_Repo/code/7_Genotypes/pedigree.R $gf_dir_UKB1/degree1_relatives_plink/UKBrel_pedigree.txt

# # Convert to plink1
# bgen2plink "$UKBv3_bgen" $gf_dir_UKB1/degree1_relatives_plink/UKB_degree1 rsID $UKBv3_bgen_sample_eid $gf_dir_UKB1/degree1_relatives_plink/UKB_degree1.sample

# # Restrict to SNPs in hdf5
# cut -d" " -f2 $PGI_Repo/derived_data/13_Parental_imputation/QC/SNPlists/UKB1/UKB1.bim > $gf_dir_UKB1/degree1_relatives_plink/SNPlist
# subsetSNPs $gf_dir_UKB1/degree1_relatives_plink/UKB_degree1_chr[1:22] $gf_dir_UKB1/degree1_relatives_plink/SNPlist $gf_dir_UKB1/degree1_relatives_plink/UKB_degree1_tmp
# rm $gf_dir_UKB1/degree1_relatives_plink/UKB_degree1_chr*
# rename degree1_tmp degree1 $gf_dir_UKB1/degree1_relatives_plink/*
# rm $gf_dir_UKB1/degree1_relatives_plink/SNPlist $gf_dir_UKB1/degree1_relatives_plink/UKBrel_pedigree.txt