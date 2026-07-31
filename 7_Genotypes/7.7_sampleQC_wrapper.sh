#!/bin/bash

source $PGI_Repo/code/paths
source $PGI_Repo/code/7_Genotypes/7.7.0_sampleQC.sh

cohortshg38=("BCS70" "MCTFR" "MCS" "NCDS" "NSHD" "STRgsa")

# Admixture population number (K)
k=5

# Mahalanobis delta cutoffs
delta_EUR=2
delta_AFR=2
delta_EAS=2
delta_SAS=5
delta_AMR=5

# Admixture ancestry component filter for Mahalanobis
admix_threshold=0.5

echo -e "cohort\tancestry\tN" > $PGI_Repo/derived_data/7_Genotypes/sampleQC_summary.txt

for cohort in WLS1kG AH1kG ALSPAC BCS70 ELSA ERisk GS GSOEP HRS MCS MCTFR MIDUS NCDS NSHD PSID Texas STRpsych STRgsa STRtwge 
do
    echo "##########################################################################"
    echo "COHORT: $cohort"
    echo ""

    if [[ ${cohortshg38[@]} =~ $cohort  ]]
    then
        snpidtype="ChrPosIDhg38"
    else
        snpidtype="ChrPosIDhg19"
    fi 
    eval gf_dir='$'gf_dir_${cohort}
    PCA $cohort "NA" $snpidtype ${cohort}_1kG_HM3_PCs
    plotPCs $cohort $gf_dir/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec $gf_dir/sampleQC/${cohort}_PCA.pdf

    if ! [[ -f $gf_dir/sampleQC/admixture/${cohort}_1kG_ancestry_proportions.txt ]]
    then 
        admixture $cohort $k
    fi
    
    for ancestry in EUR AFR EAS SAS AMR
    do
        if [[ $ancestry == "EUR" ]]
        then
            extractAncestry $cohort 5 $ancestry
        else
            extractAncestryMahalanobis $cohort $ancestry $delta_EUR $delta_AFR $delta_EAS $delta_SAS $delta_AMR $admix_threshold
        fi
        plotPCs $cohort $gf_dir/sampleQC/${cohort}_${ancestry}_1kG_HM3_PCs.eigenvec $gf_dir/sampleQC/${cohort}_${ancestry}_PCA.pdf
        
        N=$(wc -l < $gf_dir/sampleQC/${cohort}_${ancestry}_FID_IID.txt)
        echo -e "$cohort\t$ancestry\t$N" >> $PGI_Repo/derived_data/7_Genotypes/sampleQC_summary.txt
    done
    echo ""
done