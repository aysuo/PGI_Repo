#!/bin/bash

filterInfo(){
    cohort=$1
    cutoff=$2 
    infoCol=$3
    out=$4

    eval info='$'info_orig_${cohort}

    rm -f ${out}
    for chr in {1..22}; do
        infoChr=$(echo "$info" | sed "s/\[1:22\]/$chr/g")
        if [[ $infoChr == *.gz ]]; then 
            zcat $infoChr | awk -v minInfo=$cutoff -v infoCol=$infoCol '!($1~"#"){k++} k==1{for(i=1;i<=NF;i++) {ix[$i]=i}} \
                    k>1 && $ix[infoCol]>minInfo {print $ix["rs_id"]}' >> ${out}
        else
            awk -v minInfo=$cutoff -v infoCol=$infoCol '!($1~"#"){k++} k==1{for(i=1;i<=NF;i++) {ix[$i]=i}} \
                k>1 && $ix[infoCol]>minInfo {print $ix["rs_id"]}' $infoChr >> ${out}               
        fi
    done
}


prune(){
    cohort=$1
    excludeSNPs=$2
    pre_filterInfo=$3
    infoCol=$4
    infoCutoff=$5
    mafCutoff=$6
    pruneWindow=$7
    pruneShift=$8
    pruneR2=$9
    sampleKeep=${10}
    out=${11}


    mkdir -p logs

    if [[ $excludeSNPs == "NA" ]]
        then
            touch excludeSNPs.txt
            excludeSNPs=excludeSNPs.txt
    fi

    if [[ $pre_filterInfo == "yes" ]]; then
        for chr in {1..22}; do
            gfChr=$(echo "${gf_dir}/plink2/${cohort}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
            plink2 --pfile ${gfChr} \
                --exclude bed1 ${excludeSNPs} \
                --maf $mafCutoff \
                --rm-dup force-first \
                --extract ${cohort}_infofiltered.snps \
                --keep ${sampleKeep} \
                --indep-pairwise ${pruneWindow} ${pruneShift} ${pruneR2}  \
                --out ${out}_chr${chr} 
        done
    else
        for chr in {1..22}; do
            gfChr=$(echo "${gf_dir}/plink2/${cohort}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
            plink2 --pfile $gfChr \
                --exclude bed1 ${excludeSNPs} \
                --maf $mafCutoff \
                --rm-dup force-first \
                --extract-if-info $infoCol '>'= $infoCutoff \
                --keep ${sampleKeep} \
                --indep-pairwise ${pruneWindow} ${pruneShift} ${pruneR2} \
                --out ${out}_chr${chr} 
        done
    fi
    wait

    rm -f mergelist ${cohort}_infofiltered.snps
    echo $sampleKeep
    for chr in {1..22}; do
        gfChr=$(echo "${gf_dir}/plink2/${cohort}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
        plink2 --pfile $gfChr \
            --extract ${out}_chr${chr}.prune.in \
            --keep ${sampleKeep} \
            --make-bed \
            --out ${out}_chr${chr}_pruned 

        echo ${out}_chr${chr}_pruned >> mergelist
    done
    wait

    plink1.9 --merge-list mergelist --make-bed --out ${out}_pruned
    
    if [[ -f ${out}_pruned-merge.missnp ]]; then
        for chr in {1..22}; do
            gfChr=$(echo "${gf_dir}/plink2/${cohort}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
            plink2 --bfile ${out}_chr${chr}_pruned \
                --exclude ${out}_pruned-merge.missnp \
                --make-bed \
                --out ${out}_chr${chr}_pruned.tmp 
        done
        wait

        sed -i 's/$/\.tmp/g' mergelist

        plink1.9 --merge-list mergelist --make-bed --out ${out}_pruned

        rm *tmp*
    fi 
    rm -f ${cohort}_*chr*_pruned* *prune.* mergelist excludeSNPs.txt
    mv *.log logs/
    mv *.missnp logs/
}


subset_unrelated(){
    gf=$1
    relCutoff=$2

    plink1.9 --bfile ${gf} \
        --rel-cutoff ${relCutoff} \
        --out ${gf}_rel05

    sed -i 's/\t/ /g' ${gf}_rel05.rel.id

    awk 'NR==FNR{a[$2]=$2;next} \
        ($2 in a){print $1,$2,"unrelated";next}{print $1,$2,"other"}' OFS="\t"  ${gf}_rel05.rel.id ${gf}.fam > ${gf}.clusters
}

PCs(){
    gf=$1
    out=$2

    if [[ $(grep -c "other" ${gf}.clusters) -gt 100 ]]
    then
        plink1.9 --bfile ${gf} \
                --within ${gf}.clusters \
                --keep-cluster-names other \
                --mac 1 \
                --make-just-bim \
                --out ${gf}_other_mac1
        extractFlag="--extract ${gf}_other_mac1.bim"
    else
        extractFlag=""
    fi

    cut -f2 ${gf}_other_mac1.bim | sort -u > ${gf}_other_mac1.snps

    plink1.9 --bfile ${gf} \
            --within ${gf}.clusters \
            --pca 20 \
            --pca-cluster-names unrelated \
            --out ${out} \
            $extractFlag

    rm ${gf}.clusters  \
        ${gf}*.bim \
        ${gf}.bed \
        ${gf}.fam \
        ${gf}*.nosex \
        ${gf}*.snps

    mv ${gf}_* logs/
}
