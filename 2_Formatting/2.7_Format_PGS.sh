dirPGS="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/9_Scores"
dirOut="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/pgs"
dirPC="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/8_PCs"

cd $dirOut


mergePGS(){
    cohort=$1
    public=$2

    if [[ $public == 1 ]]
        then
            PGStypes="single multi public"
            dirOut="$dirOut/withPublicPGI"
        else
            PGStypes="single multi"
    fi

    tmp=$(ls -1 $dirPGS/single/scores/PGS_${cohort}_* | head -1)
    cut -f1,2 $tmp > $dirOut/tmp_$cohort 

    for PGStype in $PGStypes; do
        for PGS in $(ls -1 $dirPGS/$PGStype/scores/PGS_${cohort}_*); do 
            awk -F"\t" 'NR==FNR{pgs[$2]=$5;next}{print $0,pgs[$2]}' OFS="\t" $PGS $dirOut/tmp_$cohort > $dirOut/tmp2_$cohort
            mv $dirOut/tmp2_$cohort $dirOut/tmp_$cohort
        done
    done

    if [[ $public == 1 ]]
        then
            awk -F"\t" 'NR==1{  for(i=3;i<=NF;i++) { if(! ($i ~ "single" || $i ~ "multi" ) ) gsub(/-\w+/,"-public",$i) } } {print}' OFS="\t" $dirOut/tmp_$cohort > $dirOut/${cohort}_PGI.txt
        else
            mv $dirOut/tmp_$cohort $dirOut/${cohort}_PGI.txt
    fi

    if ! [[ "$cohort" = UKB* ]]
        then
            awk '(NR==1){print $0,"PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20";next} \
                NR==FNR{a[$2]=$0;next} \
                {print a[$2],$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' OFS="\t" $dirOut/${cohort}_PGI.txt $dirPC/$cohort/${cohort}_PCs.eigenvec > $dirOut/${cohort}_PGI_20PCs.txt
            rm $dirOut/${cohort}_PGI.txt
    fi

    sed -i 's/_LDpred_p1//g' $dirOut/${cohort}_*.txt
    sed -i "s/${cohort}_//g" $dirOut/${cohort}_*.txt
    sed -i "s/PGS/PGI/g" $dirOut/${cohort}_*.txt
}


fixIDs(){
    cohort=$1
    awk -F"\t" 'NR==1{print;next} \
        {split($1,a,/_/) ; $1=a[1]; $2=a[2] ; print}' OFS="\t" $dirOut/${cohort}_PGI_20PCs.txt > $dirOut/tmp_$cohort
    mv $dirOut/tmp_$cohort $dirOut/${cohort}_PGI_20PCs.txt
}


for cohort in Dunedin ERisk HRS2 UKB3 WLS; do
     mergePGS $cohort 1
done

for cohort in AH Dunedin EGCUT ELSA ERisk HRS3 MCTFR STRpsych STRtwge STRyatssstage Texas UKB1 UKB2 UKB3 WLS; do
    mergePGS $cohort 0
done

for cohort in AH Dunedin ELSA ERisk; do
    fixIDs $cohort
done

cd v1.0 
for cohort in AH Dunedin EGCUT ELSA ERisk HRS3 MCTFR STRpsych STRtwge STRyatssstage Texas UKB1 UKB2 UKB3 WLS; do
    mkdir -p $cohort
    cp ../${cohort}*.txt $cohort/${cohort}_PGIrepo_v1.0.txt
    sed  "s/<cohort>/$cohort/g" ReadMe.txt > $cohort/ReadMe.txt
    cp UserGuide.pdf $cohort/UserGuide_v1.0.pdf
    cp SupplementaryTables.xlsx $cohort/SupplementaryTables.xlsx
    zip -r ${cohort}_PGIrepo_v1.0.zip $cohort
done
    