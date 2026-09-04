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
cohorts_infoFromStats=(ALSPAC NCDS MIDUS PSID)

cohorts_infoColName_info=(ALSPAC MIDUS)
cohorts_infoColName_R2=(AH AH1kG BCS70 Dunedin ELSA ERisk EstBB GS GSOEP HRS MCS MCTFR NSHD PSID STRpsych STRtwge Texas WLS WLS1kG)
cohorts_infoColName_minInfo=(NCDS)
cohorts_ChrPosIDhg38=(BCS70 MCTFR MCS NCDS NSHD STRgsa)
cohorts_ChrPosIDhg19=(1000G AH AH1kG ALSPAC BCS70 Dunedin ELSA ERisk EstBB GS GSOEP HRS MCS MCTFR MIDUS MIDUSomni10 MIDUSomni11 NCDS NSHD PSID STRpsych STRtwge Texas WLS WLS1kG)


infoCutoff=0.7
mafCutoff=0.01
pruneWindow=1000
pruneShift=5
pruneR2=0.1
relCutoff=0.05



getPCs() {

    cohort=$1
    ancestry=$2

    echo "==================================="
    echo "Getting PCs for cohort $cohort, ancestry $ancestry".

    eval gf_dir='$'gf_dir_${cohort}
    eval pc_dir='$'pc_dir_${cohort}

    if [[ " ${cohorts_infoColName_info[@]} " =~ " ${cohort} "  ]]
    then
        echo "Info column name is 'info'."
        infoCol=info
    elif [[ " ${cohorts_infoColName_R2[@]} " =~ " ${cohort} "  ]]
    then
        echo "Info column name is 'R2'."
        infoCol=R2
    elif [[ " ${cohorts_infoColName_minInfo[@]} " =~ " ${cohort} "  ]]
    then
        echo "Info column name is 'minInfo'."
        infoCol=minInfo
    fi

    if [[ " ${cohorts_infoFromStats[@]} " =~ " ${cohort} "  ]]
    then
        echo "Cohort $cohort does not have imputation accuracy in plink2 files. Obtaining list of SNPs with info>$infoCutoff from imputation summary statistics."
        getInfoFilteredSNPs $cohort $infoCutoff $infoCol ${pc_dir}/${cohort}_infofiltered.snps
        infoFromStats="yes"
    else
        infoFromStats="no"
    fi

    if [[ $cohort == "PSID" ]]
    then
        awk '{split($1,a,":"); print a[1]":"a[2]}' $pc_dir/${cohort}_infofiltered.snps > $pc_dir/${cohort}_infofiltered.snps.tmp
        mv $pc_dir/${cohort}_infofiltered.snps.tmp $pc_dir/${cohort}_infofiltered.snps
    fi

    if [[ " ${cohorts_ChrPosIDhg38[@]} " =~ " ${cohort} "  ]]
    then
        echo "Using high-LD regions for hg38."
        exclusion_regions=${PGI_Repo}/code/8_PCs/high-LD-regions-hg38-GRCh38.txt
        snpidtype="ChrPosIDhg38"
    elif [[ " ${cohorts_ChrPosIDhg19[@]} " =~ " ${cohort} "  ]]
    then
        echo "Using high-LD regions for hg19."
        exclusion_regions=${PGI_Repo}/code/8_PCs/high-LD-regions-hg19-GRCh37.txt
        snpidtype="ChrPosIDhg19"
    fi

    if [[ $(wc -l < ${gf_dir}/sampleQC/${cohort}_${ancestry}_FID_IID.txt ) -gt 500 ]]
    then
        echo "Filtering and pruning genotype data.."
        prune ${cohort} \
            ${exclusion_regions} \
            ${infoFromStats} \
            ${infoCol} \
            ${infoCutoff} \
            ${mafCutoff} \
            ${pruneWindow} \
            ${pruneShift} \
            ${pruneR2} \
            ${gf_dir}/sampleQC/${cohort}_${ancestry}_FID_IID.txt \
            ${pc_dir}/${cohort}_${ancestry}_maf_info_filtered_highLDexcluded

        echo "Obtain unrelated individuals.."
        subset_unrelated ${pc_dir}/${cohort}_${ancestry}_maf_info_filtered_highLDexcluded_pruned ${relCutoff}
        N=$(grep -c "unrelated" ${pc_dir}/${cohort}_${ancestry}_maf_info_filtered_highLDexcluded_pruned.clusters)
    else
        N=$(wc -l < ${gf_dir}/sampleQC/${cohort}_${ancestry}_FID_IID.txt )
    fi

    if [[ $N -lt 500 ]]
    then
        echo "Cohort $cohort has less than 500 unrelated individuals (N=$N) of $ancestry ancestry. Going to estimate PC weights using 1000 Genomes."
        rm -f ${pc_dir}/${cohort}_${ancestry}_maf_info_filtered_highLDexcluded_pruned* 

        if ! [[ -f ${gf_dir_1000G}/sampleQC/1000Gph3_${ancestry}_FID_IID.txt ]]
        then 
            awk -F"\t" -v anc=$ancestry '$6==anc{print $1}' $pop1000G > ${gf_dir_1000G}/sampleQC/1000Gph3_${ancestry}_FID_IID.txt
        fi

        # if [[ ! -f ${gf_dir_1000G}/PCs/1000Gph3_${ancestry}_${snpidtype}.eigenvec ]]
        # then
            if [[ ! -f ${gf_dir_1000G}/plink/1000Gph3_${ancestry}_maf_filtered_highLDexcluded_rsID_pruned.bed ]]
            then
                echo "Filtering and pruning 1000 Genomes genotype data.."
                prune 1000G \
                    ${exclusion_regions} \
                    no \
                    NA \
                    NA \
                    ${mafCutoff} \
                    ${pruneWindow} \
                    ${pruneShift} \
                    ${pruneR2} \
                    ${gf_dir_1000G}/sampleQC/1000Gph3_${ancestry}_FID_IID.txt \
                    ${gf_dir_1000G}/plink/1000Gph3_${ancestry}_maf_filtered_highLDexcluded_rsID
            fi
            
            if [[ $snpidtype == "ChrPosIDhg19" || $snpidtype == "ChrPosIDhg38" ]] && [[ ! -f  ${gf_dir_1000G}/plink/1000Gph3_${ancestry}_maf_filtered_highLDexcluded_ChrPosIDhg19_pruned.bed ]]
            then
                echo "Converting 1000 Genomes genotype data snpid's to Chr:Pos format (hg19)"
                rs2chrpos ${gf_dir_1000G}/plink/1000Gph3_${ancestry}_maf_filtered_highLDexcluded_rsID_pruned ${gf_dir_1000G}/plink/1000Gph3_${ancestry}_maf_filtered_highLDexcluded_ChrPosIDhg19_pruned
            fi

            if [[ $snpidtype == "ChrPosIDhg38" && ! -f  ${gf_dir_1000G}/plink/1000Gph3_${ancestry}_maf_filtered_highLDexcluded_ChrPosIDhg38_pruned.bed ]]
            then
                echo "Lifting coordinates from hg19 to hg38 and updating Chr:Pos format snpid's accordingly."
                liftOver ${gf_dir_1000G}/plink/1000Gph3_${ancestry}_maf_filtered_highLDexcluded_ChrPosIDhg19_pruned ${gf_dir_1000G}/plink/1000Gph3_${ancestry}_maf_filtered_highLDexcluded_ChrPosIDhg38_pruned $chain_hg19toHg38
            fi
        # fi

        echo "Obtaining PCs using 1000 Genomes data.."
        PC_project ${gf_dir_1000G}/plink/1000Gph3_${ancestry}_maf_filtered_highLDexcluded_${snpidtype}_pruned \
            ${gf_dir_1000G}/PCs/1000Gph3_${ancestry}_${snpidtype} \
            ${cohort} \
            ${gf_dir}/sampleQC/${cohort}_${ancestry}_FID_IID.txt  \
            ${pc_dir}/${cohort}_${ancestry}_PCs
        echo "$cohort $ancestry 1000G" >> PCinfo.txt
    else
        echo "Obtaining PCs.."
        PCs ${pc_dir}/${cohort}_${ancestry}_maf_info_filtered_highLDexcluded_pruned ${pc_dir}/${cohort}_${ancestry}_PCs 
        echo "$cohort $ancestry own" >> PCinfo.txt
    fi

    
}

EURcohorts=(AH ALSPAC BCS70 ELSA ERisk GS GSOEP HRS MCS MCTFR MIDUS NSHD PSID STRgsa STRpsych STRtwge Texas WLS) #Dunedin EstBB TwinLife UKB FinnGen
AFRcohorts=(HRS AH1kG MCS MIDUS PSID)
EAScohorts=(HRS AH1kG MCTFR MIDUS)
SAScohorts=(BCS70 MCS)
AMRcohorts=(HRS AH1kG MCTFR PSID)

for ancestry in EUR EAS AFR AMR SAS 
do
    eval "cohorts=(\${${ancestry}cohorts[@]})"
    for cohort in "${cohorts[@]}" 
    do
        echo "==================================="
        echo COHORT: $cohort, ANCESTRY: $ancestry

        eval pc_dir='$'pc_dir_${cohort}

        if [[ -f ${pc_dir}/${cohort}_${ancestry}_PCs.eigenvec ||  -f ${pc_dir}/${cohort}_${ancestry}_PCs.sscore ]]
        then
            echo "PCs already exist for $cohort, skipping..."
            continue
        fi
        getPCs $cohort $ancestry > $PGI_Repo/code/8_PCs/8_getPCs_${cohort}_${ancestry}.log 
    done
done