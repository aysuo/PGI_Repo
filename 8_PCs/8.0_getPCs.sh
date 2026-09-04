#!/bin/bash

source $PGI_Repo/code/7_Genotypes/7.6.0_formatConversion.sh

getInfoFilteredSNPs(){
    coh=$1
    cutoff=$2 
    infoCol=$3
    out=$4

    eval info='$'info_orig_${coh}

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
    coh=$1
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
            touch ${pc_dir}/excludeSNPs.txt
            excludeSNPs=${pc_dir}/excludeSNPs.txt
    fi

    if [[ $coh == "1000G" ]]; then
        cohF="1000Gph3"
    else 
        cohF=$coh
    fi
    eval gf_dir='$'gf_dir_${coh}

    if [[ $pre_filterInfo == "yes" ]]; then
        echo "Filtering SNPs based on info score from info files"
        for chr in {1..22}; do
            gfChr=$(echo "${gf_dir}/plink2/${cohF}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
            plink2 --pfile ${gfChr} \
                --exclude bed1 ${excludeSNPs} \
                --maf $mafCutoff \
                --rm-dup force-first \
                --set-missing-var-ids @:# \
                --extract ${pc_dir}/${coh}_infofiltered.snps \
                --keep ${sampleKeep} \
                --indep-pairwise ${pruneWindow} ${pruneShift} ${pruneR2}  \
                --out ${out}_chr${chr} 
        done
        rm -f ${pc_dir}/${coh}_infofiltered.snps
    else
        if [[  $infoCol == "NA" || $infoCutoff == "NA" ]]; then
            echo "Skipping info filtering because info column or cutoff not provided"
            for chr in {1..22}; do
                gfChr=$(echo "${gf_dir}/plink2/${cohF}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
                plink2 --pfile $gfChr \
                    --exclude bed1 ${excludeSNPs} \
                    --maf $mafCutoff \
                    --set-missing-var-ids @:# \
                    --rm-dup force-first \
                    --keep ${sampleKeep} \
                    --indep-pairwise ${pruneWindow} ${pruneShift} ${pruneR2} \
                    --out ${out}_chr${chr} 
            done
        else
            echo "Filtering SNPs based on info score from plink2 files"
            for chr in {1..22}; do
                gfChr=$(echo "${gf_dir}/plink2/${cohF}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
                plink2 --pfile $gfChr \
                    --exclude bed1 ${excludeSNPs} \
                    --maf $mafCutoff \
                    --set-missing-var-ids @:# \
                    --rm-dup force-first \
                    --extract-if-info $infoCol '>'= $infoCutoff \
                    --keep ${sampleKeep} \
                    --indep-pairwise ${pruneWindow} ${pruneShift} ${pruneR2} \
                    --out ${out}_chr${chr} 
            done
        fi
    fi
    wait

    for chr in {1..22}; do
        gfChr=$(echo "${gf_dir}/plink2/${cohF}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
        plink2 --pfile $gfChr \
            --extract ${out}_chr${chr}.prune.in \
            --max-alleles 2 \
            --keep ${sampleKeep} \
            --make-bed \
            --out ${out}_chr${chr}_pruned 
    done
    wait

    mergePlink ${out}_chr[1:22]_pruned ${out}_pruned
    rm -f ${out}_chr*_pruned* ${out}_*prune.* ${pc_dir}/excludeSNPs.txt
    mv ${pc_dir}/*.log ${pc_dir}/logs/
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
        cut -f2 ${gf}_other_mac1.bim | sort -u > ${gf}_other_mac1.snps
    else
        extractFlag=""
    fi



    plink1.9 --bfile ${gf} \
            --within ${gf}.clusters \
            --pca 20 \
            --pca-cluster-names unrelated \
            --out ${out} \
            $extractFlag

    rm  ${gf}*.bim \
        ${gf}.bed \
        ${gf}.fam \
        ${gf}*.nosex \
        ${gf}*.snps \
        ${gf}.clusters \
        ${out}.nosex

    mv ${gf}_* ${pc_dir}/logs/
}

PC_project(){   
    gfRef=$1      # Reference genotype data (plink format) to estimate PC weights
    pcRefOut=$2   # PC output path for ref data
    coh=$3     # Cohort to project onto reference panel
    sampleKeep=$4   # List of individuals to keep in the coh (individuals of a specific ancestry)  
    out=$5          # Output file path

    cut -f2 $gfRef.bim  > ${gfRef}.snps
    for chr in {1..22}; do
        gfChr=$(echo "${gf_dir}/plink2/${coh}_chr[1:22]" | sed "s/\[1:22\]/$chr/g")
        plink2 --pfile $gfChr \
            --keep ${sampleKeep} \
            --mac 1 \
            --extract ${gfRef}.snps \
            --make-bed \
            --out ${pc_dir}/${coh}_${ancestry}_chr${chr} 
    done
    mergePlink ${pc_dir}/${coh}_${ancestry}_chr[1:22] ${pc_dir}/${coh}_${ancestry}
    rm ${pc_dir}/${coh}_${ancestry}_chr*

    echo "Estimating PC weights.."
    plink2 --bfile ${gfRef} \
        --rm-dup force-first \
        --freq counts \
        --pca 20 allele-wts vcols=chrom,ref,alt \
        --out ${pcRefOut}
    
    awk -F"\t" 'NR==1{print;next} $6>0{print}' ${pcRefOut}.acount > ${pcRefOut}.acount.tmp
    mv ${pcRefOut}.acount.tmp ${pcRefOut}.acount

    
    echo "Projecting coh data onto estimated weights.."
    plink2 --bfile ${pc_dir}/${coh}_${ancestry} \
        --read-freq ${pcRefOut}.acount \
        --score ${pcRefOut}.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
        --score-col-nums 6-25 \
        --out ${out} 

    rm ${pc_dir}/${coh}_${ancestry}.bed \
        ${pc_dir}/${coh}_${ancestry}.bim \
        ${pc_dir}/${coh}_${ancestry}.fam \
        ${pc_dir}/${coh}_${ancestry}.nosex

    mv ${pc_dir}/*.log ${pc_dir}/logs/
}