dirOut="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/12_h2"
dirGeno="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes"
dirPheno="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/10_Prediction/input"
dirPC="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/8_PCs"
HRS2crosswalk="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/prediction_phenotypes/HRS/HRS_GENOTYPEV2_XREF.txt"

for cohort in AH Dunedin EGCUT ERisk ELSA HRS2 HRS3 MCTFR STRpsych STRtwge STRyatssstage Texas WLS; do
    declare gf_$cohort="$dirGeno/$cohort/plink/HM3/${cohort}_HM3"
    declare sample_${cohort}="${dirGeno}/${cohort}/sampleQC/${cohort}_EUR_FID_IID.txt"
done

makeGRM_HM3(){
    cohort=$1
    eval gf='$'gf_${cohort}
    eval sample='$'sample_${cohort}
    mkdir -p $dirOut/$cohort/grm

    gcta64 --bfile $gf \
        --keep $sample \
        --maf 0.01 \
        --make-grm \
        --out $dirOut/$cohort/grm/${cohort}_hm3_maf01 \
        --thread-num 10

    gcta64 --grm $dirOut/$cohort/grm/${cohort}_hm3_maf01 \
        --grm-cutoff 0.025 \
        --make-grm \
        --out $dirOut/$cohort/grm/${cohort}_hm3_maf01_rel025

    rm $dirOut/$cohort/grm/${cohort}_hm3_maf01.*
}

estimate_h2(){
    cohort=$1
    pheno=$2
    mkdir -p $dirOut/$cohort/h2 

    gcta64 --grm $dirOut/$cohort/grm/${cohort}_hm3_maf01_rel025 \
        --pheno "$dirPheno/$cohort/tmp_$pheno.pheno" \
        --reml \
        --qcovar $dirPC/$cohort/${cohort}_10PCs.eigenvec \
        --out $dirOut/$cohort/h2/${cohort}_${pheno} \
        --thread-num 10
}

h2_table(){
    cohort=$1
    
    echo -e "Phenotype\th2\tSE\tPval\tN" > $dirOut/$cohort/h2_$cohort.txt
    for phenoPath in $dirPheno/$cohort/*; do
        pheno=$(echo $phenoPath | rev | cut -f1 -d"/" | rev | cut -f1 -d".")
        h2=$(grep "V(G)/Vp" $dirOut/$cohort/h2/${cohort}_${pheno}.hsq | cut -f2)
        SE=$(grep "V(G)/Vp" $dirOut/$cohort/h2/${cohort}_${pheno}.hsq | cut -f3)
        Pval=$(grep "Pval" $dirOut/$cohort/h2/${cohort}_${pheno}.hsq | cut -f2)
        n=$(grep "^n" $dirOut/$cohort/h2/${cohort}_${pheno}.hsq | cut -f2)


        echo -e "$pheno\t$h2\t$SE\t$Pval\t$n" >> $dirOut/$cohort/h2_$cohort.txt
    done
}

for cohort in WLS; do
    i=0
    #makeGRM_HM3 $cohort

    cut -f1-12 -d" " $dirPC/$cohort/${cohort}_PCs.eigenvec > $dirPC/$cohort/${cohort}_10PCs.eigenvec

    for phenoPath in $dirPheno/$cohort/*; do
        pheno=$(echo $phenoPath | rev | cut -f1 -d"/" | rev | cut -f1 -d".")
        if [[ $cohort == "WLS" ]]; then
            awk -F"," '{print $1,$1,$3}' OFS="\t" $phenoPath > $dirPheno/$cohort/tmp_$pheno.pheno
        elif [[ $cohort == "HRS2" ]]; then
            eval gf='$'gf_${cohort}
            if [[ $(grep hhidpn $dirPheno/$cohort/$pheno.pheno) ]]; then
                awk -F"," 'NR==FNR{IID[$5]=$3;next}($1 in IID){print IID[$1],$2}' $HRS2crosswalk $phenoPath > $dirPheno/$cohort/temp_$pheno.pheno
            else
                awk -F"," 'NR==FNR{IID[$1FS$2]=$3;next}($1FS$2 in IID){print IID[$1FS$2],$3}' $HRS2crosswalk $phenoPath > $dirPheno/$cohort/temp_$pheno.pheno
            fi
            awk 'NR==FNR{FID[$2]=$1;next}($1 in FID){print FID[$1],$1,$2}' ${gf}.fam $dirPheno/$cohort/temp_$pheno.pheno > $dirPheno/$cohort/tmp_$pheno.pheno
            rm $dirPheno/$cohort/temp_$pheno.pheno
        fi
        
        estimate_h2 $cohort $pheno
        
        let i+=1
    
        if [[ $i == 10 ]]; then 
            wait
            i=0
        fi 
    done
    wait
    rm $dirPheno/$cohort/tmp_*
    
    h2_table $cohort
done

