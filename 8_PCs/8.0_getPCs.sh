#!/bin/bash

source $PGI_Repo/code/paths

#--------------------------------------------------------------------------------------------------------#

filterInfo(){
    cohort=$1

    eval info='$'info_orig_${cohort}
    eval gf_dir='$'gf_dir_${cohort}
    eval pc_dir='$'pc_dir_${cohort}

    rm -f $pc_dir/${cohort}_info70.snps
    for chr in {1..22}; do
        infoChr=$(echo "$info" | sed "s/\[1:22\]/$chr/g")
        if [[ $infoChr == *.gz ]]; then 
            zcat $infoChr | awk 'NR==1{for(i=1;i<=NF;i++) {ix[$i]=i}} \
                    NR>1 && $ix["info"]>0.7 {print $ix["rs_id"]}' >> $pc_dir/${cohort}_info70.snps 
        else
            awk 'NR==1{for(i=1;i<=NF;i++) {ix[$i]=i}} \
                NR>1 && $ix["info"]>0.7 {print $ix["rs_id"]}' $infoChr >> $pc_dir/${cohort}_info70.snps                  
        fi
    done
}


prune(){
    cohort=$1
    eval gf_dir='$'gf_dir_${cohort}
    eval pc_dir='$'pc_dir_${cohort}
    eval sampleKeep='$'EURsample_${cohort}

    if [[ $cohort == HRS2 || $cohort == WLS ]]; then
        for chr in {1..22}; do
            gfChr=$(echo "${gf_dir}/plink2/${cohort}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
            plink2 --pfile ${gfChr} \
                --exclude ${PGI_Repo}/code/8_PCs/exclusion_regions_hg19.txt \
                --maf 0.01 \
                --rm-dup force-first \
                --extract ${pc_dir}/${cohort}_info70.snps \
                --keep ${gf_dir}/sampleQC/${cohort}_EUR_FID_IID.txt \
                --indep-pairwise 1000 5 0.1 \
                --out ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr} &
        done
    else
        if [[ $cohort == "EGCUT" ]]; then
            info=AR2
        else
            info=R2
        fi

        for chr in {1..22}; do
            gfChr=$(echo "${gf_dir}/plink2/${cohort}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
            plink2 --pfile $gfChr \
                --exclude $PGI_Repo/code/8_PCs/exclusion_regions_hg19.txt \
                --maf 0.01 \
                --rm-dup force-first \
                --extract-if-info $info '>'= 0.7 \
                --keep ${gf_dir}/sampleQC/${cohort}_EUR_FID_IID.txt \
                --indep-pairwise 1000 5 0.1 \
                --out ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr} &
        done
    fi
    wait

    rm -f  ${pc_dir}/mergelist

    for chr in {1..22}; do
        gfChr=$(echo "${gf_dir}/plink2/${cohort}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
        plink2 --pfile $gfChr \
            --extract ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr}.prune.in \
            --keep ${gf_dir}/sampleQC/${cohort}_EUR_FID_IID.txt \
            --make-bed \
            --out ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr}_pruned &

        echo ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr}_pruned >> ${pc_dir}/mergelist
    done
    wait

    plink1.9 --merge-list ${pc_dir}/mergelist --make-bed --out ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_pruned
    
    if [[ -f ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_pruned-merge.missnp ]]; then
        for chr in {1..22}; do
            gfChr=$(echo "${gf_dir}/plink2/${cohort}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
            plink2 --bfile ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr}_pruned \
                --exclude ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_pruned-merge.missnp \
                --make-bed \
                --out ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_chr${chr}_pruned.tmp & 
        done
        wait

        sed -i 's/$/\.tmp/g' ${pc_dir}/mergelist

        plink1.9 --merge-list ${pc_dir}/mergelist --make-bed --out ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_pruned

        rm ${pc_dir}/*tmp*
    fi 
    rm ${pc_dir}/${cohort}_*chr*_pruned* ${pc_dir}/*prune.* ${pc_dir}/mergelist
}


subset_unrelated(){
    cohort=$1
    
    plink1.9 --bfile ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_pruned \
        --rel-cutoff 0.05 \
        --out ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_pruned_rel05

    sed -i 's/\t/ /g' ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_pruned_rel05.rel.id

    awk 'NR==FNR{a[$2]=$2;next} \
        ($2 in a){print $1,$2,"unrelated";next}{print $1,$2,"other"}' OFS="\t"  ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_pruned_rel05.rel.id ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_pruned.fam > ${pc_dir}/${cohort}_EUR.clusters
}

PCs(){
    cohort=$1

    plink1.9 --bfile ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_pruned \
            --within ${pc_dir}/${cohort}_EUR.clusters \
            --pca 20 \
            --pca-cluster-names unrelated \
            --out ${pc_dir}/${cohort}_PCs

    rm ${pc_dir}/${cohort}_EUR.clusters \
        ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_pruned.bim \
        ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_pruned.bed \
        ${pc_dir}/${cohort}_EUR_maf01_info70_highLDexcluded_pruned.fam
}

main(){
    for cohort in HRS2 WLS; do
        mkdir -p ${pc_dir}
        filterInfo $cohort
        prune $cohort
        subset_unrelated $cohort
        PCs $cohort
    done

    for cohort in HRS3 EGCUT AH MCTFR Texas STRyatssstage STRpsych STRtwge ELSA Dunedin ERisk; do
        mkdir -p ${pc_dir}
        prune $cohort
        subset_unrelated $cohort
        PCs $cohort
    done
}

main