#!/bin/bash

source $mainDir/code/paths

for cohort in AH Dunedin EGCUT ELSA ERisk HRS2 HRS3 MCTFR Texas STRpsych STRtwge STRyatssstage WLS; do
    declare gf_${cohort}_plink2="$mainDir/derived_data/7_Genotypes/${cohort}/plink2/${cohort}_chr[1:22]"
    declare EURsample_${cohort}="$mainDir/derived_data/7_Genotypes/${cohort}/sampleQC/${cohort}_EUR_FID_IID.txt"
done


#--------------------------------------------------------------------------------------------------------#

filterInfo(){
    cohort=$1

    eval info='$'info_orig_${cohort}

    rm -f $mainDir/derived_data/8_PCs/${cohort}/${cohort}_info70.snps
    for chr in {1..22}; do
        infoChr=$(echo "$info" | sed "s/\[1:22\]/$chr/g")
        if [[ $infoChr == *.gz ]]; then 
            zcat $infoChr | awk 'NR==1{for(i=1;i<=NF;i++) {ix[$i]=i}} \
                    NR>1 && $ix["info"]>0.7 {print $ix["rs_id"]}' >> $mainDir/derived_data/8_PCs/${cohort}/${cohort}_info70.snps 
        else
            awk 'NR==1{for(i=1;i<=NF;i++) {ix[$i]=i}} \
                NR>1 && $ix["info"]>0.7 {print $ix["rs_id"]}' $infoChr >> $mainDir/derived_data/8_PCs/${cohort}/${cohort}_info70.snps                  
        fi
    done
}


prune(){
    cohort=$1
    eval gf='$'gf_${cohort}_plink2
    eval sampleKeep='$'EURsample_${cohort}

    if [[ $cohort == HRS2 || $cohort == WLS ]]; then
        for chr in {1..22}; do
            gfChr=$(echo "$mainDir/derived_data/7_Genotypes/${cohort}/plink2/${cohort}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
            plink2 --pfile $gfChr \
                --exclude $mainDir/code/8_PCs/exclusion_regions_hg19.txt \
                --maf 0.01 \
                --rm-dup force-first \
                --extract $mainDir/derived_data/8_PCs/${cohort}/${cohort}_info70.snps \
                --keep $sampleKeep \
                --indep-pairwise 1000 5 0.1 \
                --out $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr} &
        done
    else
        if [[ $cohort == "EGCUT" ]]; then
            info=AR2
        else
            info=R2
        fi

        for chr in {1..22}; do
            gfChr=$(echo "$mainDir/derived_data/7_Genotypes/${cohort}/plink2/${cohort}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
            plink2 --pfile $gfChr \
                --exclude $mainDir/code/8_PCs/exclusion_regions_hg19.txt \
                --maf 0.01 \
                --rm-dup force-first \
                --extract-if-info $info '>'= 0.7 \
                --keep $sampleKeep \
                --indep-pairwise 1000 5 0.1 \
                --out $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr} &
        done
    fi
    wait

    rm -f  $mainDir/derived_data/8_PCs/$cohort/mergelist
    for chr in {1..22}; do
        gfChr=$(echo "$mainDir/derived_data/7_Genotypes/${cohort}/plink2/${cohort}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
        plink2 --pfile $gfChr \
            --extract $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr}.prune.in \
            --keep $sampleKeep \
            --make-bed \
            --out $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr}_pruned &

        echo $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr}_pruned >> $mainDir/derived_data/8_PCs/$cohort/mergelist
    done
    wait

    plink1.9 --merge-list $mainDir/derived_data/8_PCs/$cohort/mergelist --make-bed --out $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned
    
    if [[ -f $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned-merge.missnp ]]; then
        for chr in {1..22}; do
            gfChr=$(echo "$mainDir/derived_data/7_Genotypes/${cohort}/plink2/${cohort}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
            plink2 --bfile $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr}_pruned \
                --exclude $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned-merge.missnp \
                --make-bed \
                --out $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr}_pruned.tmp & 
        done
        wait

        sed -i 's/$/\.tmp/g' $mainDir/derived_data/8_PCs/$cohort/mergelist

        plink1.9 --merge-list $mainDir/derived_data/8_PCs/$cohort/mergelist --make-bed --out $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned

        rm $mainDir/derived_data/8_PCs/$cohort/*tmp*
    fi 
    rm $mainDir/derived_data/8_PCs/$cohort/${cohort}_*chr*_pruned* $mainDir/derived_data/8_PCs/$cohort/*prune.* $mainDir/derived_data/8_PCs/$cohort/mergelist
}


subset_unrelated(){
    cohort=$1
    
    plink1.9 --bfile $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned \
        --rel-cutoff 0.05 \
        --out $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned_rel05

    sed -i 's/\t/ /g' $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned_rel05.rel.id

    awk 'NR==FNR{a[$2]=$2;next} \
        ($2 in a){print $1,$2,"unrelated";next}{print $1,$2,"other"}' OFS="\t"  $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned_rel05.rel.id $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned.fam > $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR.clusters
}

PCs(){
    cohort=$1

    plink1.9 --bfile $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned \
            --within $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR.clusters \
            --pca 20 \
            --pca-cluster-names unrelated \
            --out $mainDir/derived_data/8_PCs/$cohort/${cohort}_PCs

    rm $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR.clusters \
        $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned.bim \
        $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned.bed \
        $mainDir/derived_data/8_PCs/$cohort/${cohort}_EUR_maf01_info70_highLDexcluded_pruned.fam
}

main(){
    for cohort in HRS2 WLS; do
        mkdir -p $mainDir/derived_data/8_PCs/$cohort
        filterInfo $cohort
        prune $cohort
        subset_unrelated $cohort
        PCs $cohort
    done

    for cohort in HRS3 EGCUT AH MCTFR Texas STRyatssstage STRpsych STRtwge ELSA Dunedin ERisk; do
        mkdir -p $mainDir/derived_data/8_PCs/$cohort
        prune $cohort
        subset_unrelated $cohort
        PCs $cohort
    done
}

main