#!/bin/bash

source $PGI_Repo/code/paths
source $PGI_Repo/code/8_PCs/8.0_getPCs.sh

# High LD regions obtained from 
# https://github.com/cran/plinkQC/blob/master/inst/extdata/high-LD-regions-hg19-GRCh37.txt
# https://github.com/cran/plinkQC/blob/master/inst/extdata/high-LD-regions-hg38-GRCh38.txt
# https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)#cite_note-3
# Price et al. (2008) Long-Range LD Can Confound Genome Scans in Admixed Populations. Am. J. Hum. Genet. 86, 127-147
# Weale M. (2010) Quality Control for Genome-Wide Association Studies from Michael R. Barnes and Gerome Breen (eds.), Genetic Variation: Methods and Protocols, Methods in Molecular Biology, vol. 628, DOI 10.1007/978-1-60327-367-1_19, © Springer Science+Business Media, LLC 2010
# Anderson, Carl A., et al. "Data quality control in genetic case-control association studies." Nature protocols 5.9 (2010): 1564-1573.

# Cohorts that do not have imputation accuracy in plink2 files
cohorts_preFilterInfo=(ALSPAC NCDS MIDUS PSID)

cohorts_infoColName_info=(ALSPAC MIDUS)
cohorts_infoColName_R2=(AH AH1kG BCS70 Dunedin ELSA ERisk EstBB GS GSOEP HRS MCS MCTFR NSHD PSID STRpsych STRtwge Texas WLS WLS1kG)
cohorts_infoColName_minInfo=(NCDS)
cohorts_ChrPosIDhg38=(BCS70 MCTFR MCS NCDS NSHD STRgsa)
cohorts_ChrPosIDhg19=(AH AH1kG ALSPAC BCS70 Dunedin ELSA ERisk EstBB GS GSOEP HRS MCS MCTFR MIDUS MIDUSomni10 MIDUSomni11 NCDS NSHD PSID STRpsych STRtwge Texas WLS WLS1kG)


infoCutoff=0.7
mafCutoff=0.01
pruneWindow=1000
pruneShift=5
pruneR2=0.1
relCutoff=0.05


getPCs() {

    cohort=$1
    ancestry=$2

    eval gf_dir='$'gf_dir_${cohort}
    eval pc_dir='$'pc_dir_${cohort}

    cd $pc_dir

    if [[ " ${cohorts_infoColName_info[@]} " =~ " ${cohort} "  ]]
    then
        infoCol=info
    elif [[ " ${cohorts_infoColName_R2[@]} " =~ " ${cohort} "  ]]
    then
        infoCol=R2
    elif [[ " ${cohorts_infoColName_minInfo[@]} " =~ " ${cohort} "  ]]
    then
        infoCol=minInfo
    fi

    if [[ " ${cohorts_preFilterInfo[@]} " =~ " ${cohort} "  ]]
    then
        echo "Cohort $cohort does not have imputation accuracy in plink2 files. Filtering SNPs based on info score from imputation summary statistics."
        filterInfo $cohort $infoCutoff $infoCol ${cohort}_infofiltered.snps
        pre_filterInfo="yes"
    else
        pre_filterInfo="no"
    fi

    if [[ $cohort == "PSID" ]]
    then
        awk '{split($1,a,":"); print a[1]":"a[2]}' ${cohort}_infofiltered.snps > ${cohort}_infofiltered.snps.tmp
        mv ${cohort}_infofiltered.snps.tmp ${cohort}_infofiltered.snps
    fi

    if [[ " ${cohorts_ChrPosIDhg38[@]} " =~ " ${cohort} "  ]]
    then
        exclusion_regions=${PGI_Repo}/code/8_PCs/high-LD-regions-hg38-GRCh38.txt
    elif [[ " ${cohorts_ChrPosIDhg19[@]} " =~ " ${cohort} "  ]]
    then
        exclusion_regions=${PGI_Repo}/code/8_PCs/high-LD-regions-hg19-GRCh37.txt
    fi
        
    prune ${cohort} \
        ${exclusion_regions} \
        ${pre_filterInfo} \
        ${infoCol} \
        ${infoCutoff} \
        ${mafCutoff} \
        ${pruneWindow} \
        ${pruneShift} \
        ${pruneR2} \
        ${gf_dir}/sampleQC/${cohort}_${ancestry}_FID_IID.txt \
        ${cohort}_${ancestry}_maf_info_filtered_highLDexcluded

    subset_unrelated ${cohort}_${ancestry}_maf_info_filtered_highLDexcluded_pruned ${relCutoff}
    PCs ${cohort}_${ancestry}_maf_info_filtered_highLDexcluded_pruned ${cohort}_${ancestry}_PCs
}

EURcohorts=(AH ALSPAC BCS70 Dunedin ELSA ERisk GS GSOEP HRS MCS MCTFR MIDUS NSHD STRgsa STRpsych STRtwge Texas WLS)
AFRcohorts=(HRS AH1kG MCS MIDUS PSID)
EAScohorts=(HRS AH1kG MCTFR MIDUS)
SAScohorts=(BCS70 MCS)
AMRcohorts=(HRS AH1kG MCTFR PSID)

for ancestry in AFR EAS AMR SAS 
do
    eval "cohorts=(\${${ancestry}cohorts[@]})"
    for cohort in "${cohorts[@]}" 
    do
        echo "==================================="
        echo COHORT: $cohort, ANCESTRY: $ancestry

        eval pc_dir='$'pc_dir_${cohort}

        if [[ -f ${pc_dir}/${cohort}_${ancestry}_PCs.eigenvec ]]
        then
            echo "PCs already exist for $cohort, skipping..."
            continue
        fi
        getPCs $cohort $ancestry > $PGI_Repo/code/8_PCs/8_getPCs_${cohort}_${ancestry}.log
    done
done