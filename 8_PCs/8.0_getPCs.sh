#!/bin/bash

dirCode="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/code/8_PCs"
dirOut="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/8_PCs"

gf_WLS="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes/WLS/plink2/WLS_chr[1:22]"
sample_WLS="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes/WLS/sampleQC/WLS_EUR_FID_IID.txt"

gf_HRS2="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes/HRS2/plink2/HRS2_chr[1:22]"
sample_HRS2="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes/HRS2/sampleQC/HRS2_EUR_FID_IID.txt"


info_HRS2="/disk/genetics/dbgap/data/HRS_Jun19_2018/HRS/files/62567/PhenoGenotypeFiles/RootStudyConsentSet_phs000428.CIDR_Aging_Omni1.v2.p2.c1.NPR/GenotypeFiles/phg000515.v1.HRS_phase123_imputation.genotype-imputed-data.c1/metrics/HRS2_chr[1:22].metrics.gz"
info_WLS="/disk/genetics/WLS_DBGAP/derived/gen/WLS_chr[1:22].info"


filterInfo(){
    cohort=$1

    mkdir -p $dirOut/$cohort
    eval info='$'info_${cohort}

    rm -f $dirOut/${cohort}/${cohort}_info70.snps
    for chr in {1..22}; do
        infoChr=$(echo "$info" | sed "s/\[1:22\]/$chr/g")
        if [[ $infoChr == *.gz ]]; then 
            zcat $infoChr | awk 'NR==1{for(i=1;i<=NF;i++) {ix[$i]=i}} \
                NR>1 && $ix["info"]>0.7 {print $ix["rs_id"]}' >> $dirOut/${cohort}/${cohort}_info70.snps 
        else 
            awk 'NR==1{for(i=1;i<=NF;i++) {ix[$i]=i}} \
                NR>1 && $ix["info"]>0.7 {print $ix["rs_id"]}' $infoChr >> $dirOut/${cohort}/${cohort}_info70.snps
        fi
    done
}


prune(){
    cohort=$1
    eval gf='$'gf_${cohort}
    eval sampleKeep='$'sample_${cohort}

        for chr in {1..22}; do
        gfChr=$(echo "$gf" | sed "s/\[1:22\]/$chr/g")
        plink2 --pfile $gfChr \
            --exclude $dirCode/exclusion_regions_hg19.txt \
            --maf 0.01 \
            --rm-dup force-first \
            --extract $dirOut/${cohort}/${cohort}_info70.snps \
            --keep $sampleKeep \
            --indep-pairwise 1000 5 0.1 \
            --out $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr} &
    done
    wait

    rm -f  $dirOut/$cohort/mergelist
    for chr in {1..22}; do
        gfChr=$(echo "$gf" | sed "s/\[1:22\]/$chr/g")
        plink2 --pfile $gfChr \
            --extract $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr}.prune.in \
            --keep $sampleKeep \
            --make-bed \
            --out $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr}_pruned &

        echo $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr}_pruned >> $dirOut/$cohort/mergelist
    done
    wait

    plink1.9 --merge-list $dirOut/$cohort/mergelist --make-bed --out $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned
    
    if [[ -f $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned-merge.missnp ]]; then
        for chr in {1..22}; do
            gfChr=$(echo "$gf" | sed "s/\[1:22\]/$chr/g")
            plink2 --bfile $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr}_pruned \
                --exclude $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned-merge.missnp \
                --make-bed \
                --out $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr}_pruned.tmp & 
        done
        wait

        sed -i 's/$/\.tmp/g' $dirOut/$cohort/mergelist

        plink1.9 --merge-list $dirOut/$cohort/mergelist --make-bed --out $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned

        rm $dirOut/$cohort/*tmp*
    fi 
    rm $dirOut/$cohort/${cohort}_*chr*_pruned* $dirOut/$cohort/*prune.* $dirOut/$cohort/mergelist
}


subset_unrelated(){
    cohort=$1
    
    plink1.9 --bfile $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned \
        --rel-cutoff 0.05 \
        --out $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned_rel05

    sed -i 's/\t/ /g' $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned_rel05.rel.id

    awk 'NR==FNR{a[$2]=$2;next} \
        ($2 in a){print $1,$2,"unrelated";next}{print $1,$2,"other"}' OFS="\t"  $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned_rel05.rel.id $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned.fam > $dirOut/$cohort/${cohort}_EUR.clusters
}

PCs(){
    cohort=$1

    plink1.9 --bfile $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned \
            --within $dirOut/$cohort/${cohort}_EUR.clusters \
            --pca 20 \
            --pca-cluster-names unrelated \
            --out $dirOut/$cohort/${cohort}_PCs

    rm $dirOut/$cohort/${cohort}_EUR.clusters \
        $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned.bim \
        $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned.bed \
        $dirOut/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned.fam
}

main(){
    for cohort in WLS HRS2; do
        filterInfo $cohort
        prune $cohort
        subset_unrelated $cohort
        PCs $cohort
    done

}

main