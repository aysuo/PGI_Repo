#!/bin/bash

source $PGI_Repo/code/paths


# mergePGI_v1(){
#     cohort=$1
#     public=$2

#     # Include public PGI in the package if public=1
#     if [[ $public == 1 ]]
#         then
#             PGItypes="single multi public"
#             $PGI_Repo/pgs="$PGI_Repo/pgs/withPublicPGI"
#         else
#             PGItypes="single multi"
#     fi

#     # Merge all PGI available for cohort 
#     tmp=$(ls -1 $PGI_Repo/derived_data/9_Scores/single/scores/PGI_${cohort}_* | head -1)
#     cut -f1,2 $tmp > $PGI_Repo/pgs/tmp_$cohort 
#     for PGItype in $PGItypes; do
#         for PGI in $(ls -1 $PGI_Repo/derived_data/9_Scores/$PGItype/scores/PGI_${cohort}_*); do 
#             awk -F"\t" 'NR==FNR{pgs[$2]=$5;next}{print $0,pgs[$2]}' OFS="\t" $PGI $PGI_Repo/pgs/tmp_$cohort > $PGI_Repo/pgs/tmp2_$cohort
#             mv $PGI_Repo/pgs/tmp2_$cohort $PGI_Repo/pgs/tmp_$cohort
#         done
#     done

#     # Rename public PGI column headers as $PHENO-public
#     if [[ $public == 1 ]]
#         then
#             awk -F"\t" 'NR==1{  for(i=3;i<=NF;i++) { if(! ($i ~ "single" || $i ~ "multi" ) ) gsub(/-\w+/,"-public",$i) } } {print}' OFS="\t" $PGI_Repo/pgs/tmp_$cohort > $PGI_Repo/pgs/${cohort}_PGI.txt
#         else
#             mv $PGI_Repo/pgs/tmp_$cohort $PGI_Repo/pgs/${cohort}_PGI.txt
#     fi

#     # For all cohorts except for UKB, merge in PCs (UKB PCs are available from UKB)
#     if ! [[ "$cohort" = UKB* ]]
#         then
#             awk '(NR==1){print $0,"PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20";next} \
#                 NR==FNR{a[$2]=$0;next} \
#                 {print a[$2],$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' OFS="\t" $PGI_Repo/pgs/${cohort}_PGI.txt $PGI_Repo/derived_data/8_PCs/$cohort/${cohort}_PCs.eigenvec > $PGI_Repo/pgs/${cohort}_PGI_20PCs.txt
#             rm $PGI_Repo/pgs/${cohort}_PGI.txt
#     fi

#     # Rename PGI columns (remove _LDpred_p1, cohort name, change PGI to PGI)
#     sed -i 's/_LDpred_p1//g' $PGI_Repo/pgs/${cohort}_*.txt
#     sed -i "s/${cohort}_//g" $PGI_Repo/pgs/${cohort}_*.txt
#     sed -i "s/PGI/PGI/g" $PGI_Repo/pgs/${cohort}_*.txt
# }

mergePGI_v2(){
    cohort=$1
    inputDir=$2
    outputDir=$3
    
    echo "Merging PGIs for cohort: $cohort"
    
    tmp=$(ls -1 $inputDir/PGI_${cohort}*_SBayesR.txt | head -1)
    cut -f1-2 $tmp > $outputDir/tmp/tmp
    
    for file in $inputDir/PGI_${cohort}*_SBayesR.txt; do 
        # Extract phenotype from file name
        echo "Processing file: $file"
        pheno=$(echo $file | rev | cut -d"/" -f1 | rev | cut -d"-" -f1 | sed "s/PGI_${cohort}_//g")
                
        awk -v P=$pheno -F"\t" 'NR==1{$5="PGI_"P}NR==FNR{pgi[$2]=$5;next}{print $0,pgi[$2]}' OFS="\t" $file $outputDir/tmp/tmp > $outputDir/tmp/tmp2
        mv $outputDir/tmp/tmp2 $outputDir/tmp/tmp
        echo "Processed phenotype: $pheno"
    done
    mv $outputDir/tmp/tmp $outputDir/tmp/pgi.txt
}

mergePCs(){
    cohort=$1
    outputDir=$2

    eval pc_dir='$'pc_dir_${cohort}

    echo "Merging PCs for cohort: $cohort"

    # If PCs are available, merge them
    if [[ -f $pc_dir/${cohort}_PCs.eigenvec ]]; then
        awk '(NR==1){print $0,"PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20";next} \
            NR==FNR{a[$2]=$0;next} \
            {print a[$2],$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' OFS="\t" $outputDir/tmp/pgi.txt $pc_dir/${cohort}_PCs.eigenvec > $outputDir/${cohort}_PGIrepo_v2.0.txt
    else
        echo "No PCs available for cohort: $cohort"
        cp $outputDir/tmp/pgi.txt $outputDir/${cohort}_PGIrepo_v2.0.txt
    fi
}

mergeParental(){
    cohort=$1
    inputDir=$2
    outputDir=$3

    echo "Merging parental PGIs for cohort: $cohort"

    tmp=$(ls -1 $inputDir/PGI_${cohort}_parental_* | head -1)
    echo "Adding FATHER_ID and MOTHER_ID as 3rd and 4th columns to the main PGI file"

    awk 'NR==FNR{FATHERID[$2]=$3;MOTHERID[$2]=$4;next}
        {FS="\t"; if ($2 in FATHERID) {line=$1 OFS $2 OFS FATHERID[$2] OFS MOTHERID[$2]} else {line=$1 OFS $2 OFS "NA" OFS "NA"};
        for(i=3;i<=NF;i++){line=line OFS $i}; print line}' OFS="\t" $tmp $outputDir/tmp/pgi.txt > $outputDir/tmp/tmp_parental

    echo "Adding parental PGIs to the main PGI file"
    
    for file in $inputDir/PGI_${cohort}_parental_*; do 
        echo "Processing parental file: $file"
        pheno=$(echo $file | rev | cut -d"/" -f1 | rev | cut -d"-" -f1 | sed "s/PGI_${cohort}_parental_//g")
        
        # Two different versions of snipar have been used to make the parental PGIs and the old version doesn't output MOTHER_ID and FATHER_ID, so we cannot refer to columns by number
        if [[ $(grep "paternal" $file) ]]; then
            awk -v P=$pheno '
                NR==FNR      {for (i=1; i<=NF; i++){ ix[$i] = i }}
                NR==1        {$ix["proband"]="PGI_"P"_proband";$ix["paternal"]="PGI_"P"_paternal";$ix["maternal"]="PGI_"P"_maternal"}
                NR==FNR      {pgi[$2]=$ix["proband"] OFS $ix["paternal"] OFS $ix["maternal"];next}
                {FS="\t"; if ($2 in pgi) {print $0,pgi[$2]} else {print $0,"NA","NA","NA"}}' OFS="\t" $file $outputDir/tmp/tmp_parental > $outputDir/tmp/tmp2_parental    
        else
            awk -v P=$pheno '
                NR==FNR      {for (i=1; i<=NF; i++){ ix[$i] = i }}
                NR==1        {$ix["proband"]="PGI_"P"_proband";$ix["parental"]="PGI_"P"_parental"}
                NR==FNR      {pgi[$2]=$ix["proband"] OFS $ix["parental"];next}
                {FS="\t"; if ($2 in pgi) {print $0,pgi[$2]} else {print $0,"NA","NA"}}' OFS="\t" $file $outputDir/tmp/tmp_parental > $outputDir/tmp/tmp2_parental
        fi 

        mv $outputDir/tmp/tmp2_parental $outputDir/tmp/tmp_parental
        echo "Processed phenotype: $pheno"
    done

    mv $outputDir/tmp/tmp_parental $outputDir/tmp/pgi.txt
}            

fixIDs(){
    cohort=$1
    outputDir=$2
    format=$3

    echo "Fixing FID and IID for cohort: $cohort"

    if [[ $(grep "FATHER_ID" $outputDir/${cohort}_PGIrepo_v2.0.txt) ]]; then
        awk -F"\t" -v format=$format 'NR==1{print;next} \
            {split($1,a,/_/) ; split($3,b,/_/); split($4,c,/_/);
            if (format=="fid_iid") {$1=a[1]; $2=a[2] ; if ($3!="NA") {$3=b[2]}; if ($4!="NA") $4=c[2]};
            if (format=="x_iid_x_iid") {$1=a[1]"_"a[2]; $2=a[1]"_"a[2] ; if ($3!="NA") {$3=b[1]"_"b[2]}; if ($4!="NA") $4=c[1]"_"c[2]};
            print}' OFS="\t"  $outputDir/${cohort}_PGIrepo_v2.0.txt > $outputDir/tmp/tmp
    else
        awk -F"\t" -v format=$format 'NR==1{print;next} \
            {split($1,a,/_/) ; 
            if (format=="fid_iid") {$1=a[1]; $2=a[2]};
            if (format=="x_iid_x_iid") {$1=a[1]"_"a[2]; $2=a[1]"_"a[2]} ; 
            print}' OFS="\t" $outputDir/${cohort}_PGIrepo_v2.0.txt > $outputDir/tmp/tmp
    fi

    mv $outputDir/tmp/tmp $outputDir/${cohort}_PGIrepo_v2.0.txt
}


######################## VERSION 1.0 #########################

# # Merge PGI & PCs for the validation cohorts (includes public PGI for the comparison analyses)
# for cohort in Dunedin ERisk HRS2 UKB3 WLS; do
#      mergePGI $cohort 1
# done

# # Merge PGI & PCs for all cohorts
# for cohort in AH Dunedin EGCUT ELSA ERisk HRS3 MCTFR STRpsych STRtwge STRyatssstage Texas UKB1 UKB2 UKB3 WLS; do
#     mergePGI $cohort 0
# done

# # Fix FID-IID if necessary
# for cohort in AH Dunedin ELSA ERisk; do
#     fixIDs $cohort
# done

# # Merge with ReadMe, Supp Tables, User Guide and zip
# cd v1.0 
# for cohort in AH Dunedin EGCUT ELSA ERisk MCTFR STRpsych STRtwge STRyatssstage Texas UKB1 UKB2 UKB3 WLS; do
#     mkdir -p $cohort
#     cp ../${cohort}*.txt $cohort/${cohort}_PGIrepo_v1.0.txt

#     if [[ $cohort = UKB* ]]; then
#         sed -i '/PC1/d' $cohort/ReadMe.txt
#     fi

#     cp UserGuide.pdf $cohort/UserGuide_v1.0.pdf
#     cp SupplementaryTables.xlsx $cohort/SupplementaryTables.xlsx
#     zip -r ${cohort}_PGIrepo_v1.0.zip $cohort
# done

# # Prepare HRS manually because ReadMe is different (contains info about PC shuffling) (also rename HRS3 as HRS)  
# mkdir -p HRS
# cp ../HRS3*.txt HRS/HRS_PGIrepo_v1.0.txt
# cp UserGuide.pdf HRS/UserGuide_v1.0.pdf
# cp SupplementaryTables.xlsx HRS/SupplementaryTables.xlsx
# zip -r HRS_PGIrepo_v1.0.zip HRS

################### VERSION 1.1 (EA4 update) ###################

# cd v1.1
# for cohort in EGCUT STRpsych STRtwge STRyatssstage Texas UKB1 UKB2 UKB3 WLS; do
#     mkdir $cohort
#     awk -F"\t" 'NR==FNR{EApgi[$2]=$5;next}FNR==1{print;next}($2 in EApgi){$16=EApgi[$2];print}' OFS="\t" $PGI_Repo/derived_data/9_Scores/single_SBayesR/scores/PGI_${cohort}_EA-single_SBayesR.txt ../v1.0/${cohort}/${cohort}_PGIrepo_v1.0.txt > ${cohort}/${cohort}_PGIrepo_v1.1.txt &
# done

# for cohort in AH Dunedin ELSA ERisk; do
#     mkdir $cohort
#     awk -F"\t" 'NR==1{print;next} \
#         {split($1,a,/_/) ; $1=a[1]; $2=a[2] ; print}' OFS="\t" $PGI_Repo/derived_data/9_Scores/single_SBayesR/scores/PGI_${cohort}_EA-single_SBayesR.txt > ${cohort}/ea4_${cohort}
#     awk -F"\t" 'NR==FNR{EApgi[$2]=$5;next}FNR==1{print;next}($2 in EApgi){$16=EApgi[$2];print}' OFS="\t" ${cohort}/ea4_${cohort} ../v1.0/${cohort}/${cohort}_PGIrepo_v1.0.txt > ${cohort}/${cohort}_PGIrepo_v1.1.txt
#     rm ${cohort}/ea4_${cohort}
# done

# for cohort in EGCUT STRpsych STRtwge STRyatssstage Texas UKB1 UKB2 UKB3 WLS AH Dunedin ELSA ERisk; do
#     cp ReadMe.txt $cohort/ReadMe.txt
#     # Edit ReadMe's manually, add the N of new EA PGI, etc
#     mv $cohort/ReadMe.txt $cohort/ReadMe_${cohort}_PGIrepo_v1.1.txt
# done

# for cohort in EGCUT STRpsych STRtwge STRyatssstage Texas UKB1 UKB2 UKB3 WLS AH Dunedin ELSA ERisk; do
#     zip -r ${cohort}_PGIrepo_v1.1.zip $cohort
# done


######################### VERSION 2.0 ##########################
family_cohorts="AH ALSPAC ERisk GS GSOEP MCS MCTFR PSID STRtwge STRgsa STRpsych Texas WLS"
fid_iid_format_cohorts="AH Dunedin ELSA ERisk GS GSOEP MIDUS WLS"
x_iid_x_iid_format_cohorts="BCS70 MCS NCDS"


for cohort in AH ALSPAC BCS70 Dunedin ELSA ERisk GS GSOEP HRS MCTFR MCS MIDUS NCDS NSHD PSID STRtwge STRpsych STRgsa Texas WLS; do
    echo "==================================="
    echo "Processing cohort: $cohort"
    
    eval pgiDir='$'pgi_dir_${cohort}
    dirIn=$pgiDir/single_SBayesR
    dirOut=$pgiDir/release/v2.0
    mkdir -p $dirOut/tmp

    # Merge PGI files
    mergePGI_v2 $cohort $dirIn $dirOut

    # Merge parental PGIs if family cohort
    if [[ $family_cohorts == *"$cohort"* ]]; then
        mergeParental $cohort $dirIn $dirOut
    fi
    
    # Merge PCs 
    mergePCs $cohort $dirOut
    
    # Fix FID and IID if necessary
    if [[ $fid_iid_format_cohorts == *"$cohort"* ]]; then
        fixIDs $cohort $dirOut fid_iid
    elif [[ $x_iid_x_iid_format_cohorts == *"$cohort"* ]]; then
        fixIDs $cohort $dirOut x_iid_x_iid
    fi    
    
    # Clean-up
    rm -rf $dirOut/tmp

    echo "Cohort $cohort processed successfully."
done  

####################

# Add batch to STR
Merge STR subcohorts
awk -F"\t" 'NR==1{print $0,"BATCH";next} 
    NR==FNR{print $0,"twingene";next}
    FNR>1{print $0,"psych"}' OFS="\t" $pgi_dir_STRtwge/release/v2.0/STRtwge_PGIrepo_v2.0.txt $pgi_dir_STRpsych/release/v2.0/STRpsych_PGIrepo_v2.0.txt > $pgi_dir_STRpsych/release/v2.0/tmp

awk -F"\t" 'NR==FNR{print $0;next} 
    FNR>1{print $0,"gsa"}' OFS="\t" $pgi_dir_STRpsych/release/v2.0/tmp $pgi_dir_STRgsa/release/v2.0/STRgsa_PGIrepo_v2.0.txt > $pgi_dir_STRpsych/release/v2.0/STR_PGIrepo_v2.0.txt
rm $pgi_dir_STRpsych/release/v2.0/tmp
cp $pgi_dir_STRpsych/release/v2.0/STR_PGIrepo_v2.0.txt $pgi_dir_STRgsa/release/v2.0/STR_PGIrepo_v2.0.txt

####################

# Add batch to MIDUS 
awk -F"\t" 'NR==FNR{split($1,a,/_/); $1=a[1]; b[$1]=$1; next}
            FNR==1{print $0,"BATCH";next}
            ($1 in b){print $0,"omni10";next}
            {print $0, "omni11"}' OFS="\t" $gf_dir_MIDUSomni10/plink2/MIDUSomni10_chr1.nodup.psam $pgi_dir_MIDUS/release/v2.0/MIDUS_PGIrepo_v2.0.txt > $pgi_dir_MIDUS/release/v2.0/tmp

mv $pgi_dir_MIDUS/release/v2.0/tmp $pgi_dir_MIDUS/release/v2.0/MIDUS_PGIrepo_v2.0.txt

####################

# Add batch to NCDS
awk 'NR==FNR{FS=","; gsub(/"/,"",$1); gsub(/"/,"",$2); a[$1]=$2; next} 
    FNR==1{FS="\t"; print $0,"BATCH";next}
    {FS="\t"; if ($1 in a) {print $0,a[$1]}}' OFS="\t" $array_indicator_NCDS $pgi_dir_NCDS/release/v2.0/NCDS_PGIrepo_v2.0.txt > $pgi_dir_NCDS/release/v2.0/tmp
mv $pgi_dir_NCDS/release/v2.0/tmp $pgi_dir_NCDS/release/v2.0/NCDS_PGIrepo_v2.0.txt

####################

# Shuffle HRS PCs
sh $PGI_Repo/code/2_Formatting/2.7_Shuffle_HRS_PCs.sh

# Check if there are NaNs in the PGI files, column count, row count
echo -e "Cohort\tNcol\tNrow\tNaN" > $PGI_Repo/code/2_Formatting/2.7_PGI_stats.log
for cohort in AH ALSPAC BCS70 Dunedin ELSA ERisk GS GSOEP HRS MCTFR MCS MIDUS NCDS STRtwge STRpsych STRgsa Texas UKB1 UKB2 UKB3 WLS; do
    eval pgiDir='$'pgi_dir_${cohort}

    if [[ $(grep "nan" $pgiDir/release/v2.0/${cohort}_PGIrepo_v2.0.txt) ]]; then
        NaN=1
    else
        NaN=0
    fi
    awk -F"\t" -v C=$cohort -v NaN=$NaN 'END{print C, NF, NR, NaN}' OFS="\t" $pgiDir/release/v2.0/${cohort}_PGIrepo_v2.0.txt >> $PGI_Repo/code/2_Formatting/2.7_PGI_stats.log
done

for cohort in AH ALSPAC BCS70 Dunedin ELSA ERisk GS GSOEP HRS MCTFR MCS MIDUS NCDS NSHD PSID STRtwge STRpsych STRgsa Texas WLS; do
    eval pgiDir='$'pgi_dir_${cohort}
    cd $pgiDir/release/v2.0/
    zip SSGAC_PGI_Repository_v2_${cohort}.zip SSGAC*
    #zip --encrypt SSGAC_PGI_Repository_v2_${cohort}.zip SSGAC*
done

