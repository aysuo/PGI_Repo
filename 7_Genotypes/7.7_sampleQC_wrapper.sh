#!/bin/bash

source $PGI_Repo/code/paths
source $PGI_Repo/code/7.7.0_sample_QC.sh


for cohort in AH Dunedin EGCUT ELSA ERisk HRS2 HRS3 MCS MCTFR NCDS Texas STRpsych STRtwge STRyatssstage WLS; do
    if [[ $cohort == "HRS2" || $cohort == "WLS" ]]; then
        snpidtype="rs"
    else
        snpidtype="chrpos"
    fi 
    eval gf_dir='$'gf_dir_${cohort}
    PCA $cohort "NA" $snpidtype ${cohort}_1kG_HM3_PCs
    plotPCs $cohort $gf_dir/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec $gf_dir/sampleQC/${cohort}_PCA.pdf
    extractAncestry $cohort 5 EUR
    plotPCs $cohort $gf_dir/sampleQC/${cohort}_EUR_1kG_HM3_PCs.eigenvec $gf_dir/sampleQC/${cohort}_EUR_PCA.pdf
done



