

#!/bin/bash

source paths2

cd $mainDir/pgs


mergePGI(){
    cohort=$1
    public=$2

    # Include public PGI in the package if public=1
    if [[ $public == 1 ]]
        then
            PGStypes="single multi public"
            $mainDir/pgs="$mainDir/pgs/withPublicPGI"
        else
            PGStypes="single multi"
    fi

    # Merge all PGI available for cohort 
    tmp=$(ls -1 $mainDir/derived_data/9_Scores/single/scores/PGS_${cohort}_* | head -1)
    cut -f1,2 $tmp > $mainDir/pgs/tmp_$cohort 
    for PGStype in $PGStypes; do
        for PGS in $(ls -1 $mainDir/derived_data/9_Scores/$PGStype/scores/PGS_${cohort}_*); do 
            awk -F"\t" 'NR==FNR{pgs[$2]=$5;next}{print $0,pgs[$2]}' OFS="\t" $PGS $mainDir/pgs/tmp_$cohort > $mainDir/pgs/tmp2_$cohort
            mv $mainDir/pgs/tmp2_$cohort $mainDir/pgs/tmp_$cohort
        done
    done

    # Rename public PGI column headers as $PHENO-public
    if [[ $public == 1 ]]
        then
            awk -F"\t" 'NR==1{  for(i=3;i<=NF;i++) { if(! ($i ~ "single" || $i ~ "multi" ) ) gsub(/-\w+/,"-public",$i) } } {print}' OFS="\t" $mainDir/pgs/tmp_$cohort > $mainDir/pgs/${cohort}_PGI.txt
        else
            mv $mainDir/pgs/tmp_$cohort $mainDir/pgs/${cohort}_PGI.txt
    fi

    # For all cohorts except for UKB, merge in PCs (UKB PCs are available from UKB)
    if ! [[ "$cohort" = UKB* ]]
        then
            awk '(NR==1){print $0,"PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20";next} \
                NR==FNR{a[$2]=$0;next} \
                {print a[$2],$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' OFS="\t" $mainDir/pgs/${cohort}_PGI.txt $mainDir/derived_data/8_PCs/$cohort/${cohort}_PCs.eigenvec > $mainDir/pgs/${cohort}_PGI_20PCs.txt
            rm $mainDir/pgs/${cohort}_PGI.txt
    fi

    # Rename PGI columns (remove _LDpred_p1, cohort name, change PGS to PGI)
    sed -i 's/_LDpred_p1//g' $mainDir/pgs/${cohort}_*.txt
    sed -i "s/${cohort}_//g" $mainDir/pgs/${cohort}_*.txt
    sed -i "s/PGS/PGI/g" $mainDir/pgs/${cohort}_*.txt
}


fixIDs(){
    cohort=$1
    # If FID and IID have the format FID_IID FID_IID, split them
    awk -F"\t" 'NR==1{print;next} \
        {split($1,a,/_/) ; $1=a[1]; $2=a[2] ; print}' OFS="\t" $mainDir/pgs/${cohort}_PGI_20PCs.txt > $mainDir/pgs/tmp_$cohort
    mv $mainDir/pgs/tmp_$cohort $mainDir/pgs/${cohort}_PGI_20PCs.txt
}


# Merge PGI & PCs for the validation cohorts (includes public PGI for the comparison analyses)
for cohort in Dunedin ERisk HRS2 UKB3 WLS; do
     mergePGI $cohort 1
done

# Merge PGI & PCs for all cohorts
for cohort in AH Dunedin EGCUT ELSA ERisk HRS3 MCTFR STRpsych STRtwge STRyatssstage Texas UKB1 UKB2 UKB3 WLS; do
    mergePGI $cohort 0
done

# Fix FID-IID if necessary
for cohort in AH Dunedin ELSA ERisk; do
    fixIDs $cohort
done

# Merge with ReadMe, Supp Tables, User Guide and zip
cd v1.0 
for cohort in AH Dunedin EGCUT ELSA ERisk MCTFR STRpsych STRtwge STRyatssstage Texas UKB1 UKB2 UKB3 WLS; do
    mkdir -p $cohort
    cp ../${cohort}*.txt $cohort/${cohort}_PGIrepo_v1.0.txt

    if [[ $cohort = UKB* ]]; then
        sed -i '/PC1/d' $cohort/ReadMe.txt
    fi

    cp UserGuide.pdf $cohort/UserGuide_v1.0.pdf
    cp SupplementaryTables.xlsx $cohort/SupplementaryTables.xlsx
    zip -r ${cohort}_PGIrepo_v1.0.zip $cohort
done

# Prepare HRS manually because ReadMe is different (contains info about PC shuffling) (also rename HRS3 as HRS)  
mkdir -p HRS
cp ../HRS3*.txt HRS/HRS_PGIrepo_v1.0.txt
cp UserGuide.pdf HRS/UserGuide_v1.0.pdf
cp SupplementaryTables.xlsx HRS/SupplementaryTables.xlsx
zip -r HRS_PGIrepo_v1.0.zip HRS