#!/bin/bash
dirCode="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/code/7_Genotypes"

pop1000G="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/ref_data/1000G_ph3/igsr_samples.tsv"

dirOut="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes"
gf_HRS2="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes/HRS2/plink/HM3/HRS2_HM3"
gf_WLS="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes/WLS/plink/HM3/WLS_HM3"
gf_1000G="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes/1000G/plink/HM3/1000Gph3_HM3"
#---------------------------------------------------------------------------------------#


PCA(){
    cohort=$1
    sampleKeep=$2
    eval gf='$'gf_${cohort}

    mkdir -p $dirOut/$cohort/sampleQC

    plink1.9 --bfile $gf_1000G \
        --bmerge $gf \
        --geno 0.99 \
        --make-bed \
        --out $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3

    if [[ -f $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3-merge.missnp ]]; then
        plink2 --bfile $gf_1000G \
            --exclude $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3-merge.missnp \
            --make-bed \
            --out $gf_1000G.tmp

        plink1.9 --bfile $gf_1000G.tmp \
            --bmerge $gf \
            --geno 0.99 \
            --make-bed \
            --out $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3
        
        rm $gf_1000G.tmp*
    fi

    awk -v cohort=$cohort '$1=="0"{print $1,$2,"1000G";next} \
        {print $1,$2,cohort}' $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3.fam > $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3.clusters

     # Estimate 4PCs in 1000G
    if [[ $sampleKeep == "NA" ]]; then
        plink1.9 --bfile $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3 \
            --maf 0.01 \
            --within $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3.clusters \
            --pca 4 \
            --pca-cluster-names 1000G \
            --out $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3_PCs
        sed -i 's/ /\t/g' $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec
    else
        plink1.9 --bfile $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3 \
            --keep $sampleKeep \
            --maf 0.01 \
            --within $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3.clusters \
            --pca 4 \
            --pca-cluster-names 1000G \
            --out $dirOut/$cohort/sampleQC/${cohort}_EUR_1kG_HM3_PCs
        sed -i 's/ /\t/g' $dirOut/$cohort/sampleQC/${cohort}_EUR_1kG_HM3_PCs.eigenvec
    fi

    rm $dirOut/$cohort/sampleQC/*bed \
        $dirOut/$cohort/sampleQC/*bim \
        $dirOut/$cohort/sampleQC/*fam \
        $dirOut/$cohort/sampleQC/*nosex \
        $dirOut/$cohort/sampleQC/*-merge.missnp \
        $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3.clusters
}

plotPCs(){
    cohort=$1
    PCs=$2
    out=$3

    # Annotate PCs with ancestry info
    awk -F"\t" -v cohort=$cohort 'BEGIN{OFS="\t";print "FID","IID","PC1","PC2","PC3","PC4","POP","SUPERPOP"} \
        NR==FNR{a[$1]=$4OFS$6;next} \
        $2 in a{print $0,a[$2];next} \
        !($1~"#"){print $0,cohort,cohort}' OFS="\t" $pop1000G $PCs > ${PCs}.annotated

    # Plot
    Rscript $dirCode/7.7.0_plotPCs.R "${PCs}.annotated" $out
}

extractEUR(){
    cohort=$1

    low1=$2
    up1=$3
    low2=$4
    up2=$5
    low3=$6
    up3=$7
    low4=$8
    up4=$9
        
    #low1=$(awk -F"\t" 'BEGIN{min1=9999} $8=="EUR" && NR>1 {if ($3<min1) min1=$3}END{print min1}' $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated)
    #up1=$(awk -F"\t" 'BEGIN{max1=-9999} $8=="EUR" && NR>1 {if ($3>max1) max1=$3}END{print max1}' $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated)
    #low2=$(awk -F"\t" 'BEGIN{min2=9999} $8=="EUR" && NR>1 {if ($4<min2) min2=$4}END{print min2}' $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated)
    #up2=$(awk -F"\t" 'BEGIN{max2=-9999} $8=="EUR" && NR>1 {if ($4>max2) max2=$4}END{print max2}' $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated)
    #low3=$(awk -F"\t" 'BEGIN{min3=9999} $8=="EUR" && NR>1 {if ($5<min3) min3=$5}END{print min3}' $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated)
    #up3=$(awk -F"\t" 'BEGIN{max3=-9999} $8=="EUR" && NR>1 {if ($5>max3) max3=$5}END{print max3}' $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated)
    #low4=$(awk -F"\t" 'BEGIN{min4=9999} $8=="EUR" && NR>1 {if ($6<min4) min4=$6}END{print min4}' $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated)
    #up4=$(awk -F"\t" 'BEGIN{max4=-9999} $8=="EUR" && NR>1 {if ($6>max4) max4=$6}END{print max4}' $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated)
    
    awk -F"\t" -v l1=$low1 -v u1=$up1 -v l2=$low2 -v u2=$up2 -v l3=$low3 -v u3=$up3 -v l4=$low4 -v u4=$up4  \
        '$1=="0" || ($3>l1 && $3<u1 && $4>l2 && $4<u2 && $5>l3 && $5<u3 && $6>l4 && $6<u4) {print $1,$2,$3,$4,$5,$6}' OFS="\t" $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated > $dirOut/$cohort/sampleQC/${cohort}_EUR_1kG_HM3_PCs.eigenvec

    awk -F"\t" '$1!=0{print $1,$2}' OFS="\t" $dirOut/$cohort/sampleQC/${cohort}_EUR_1kG_HM3_PCs.eigenvec > $dirOut/$cohort/sampleQC/${cohort}_EUR_FID_IID.txt

}


main(){
    for cohort in WLS HRS2; do
        PCA $cohort "NA"
        plotPCs $cohort $dirOut/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec $dirOut/$cohort/sampleQC/${cohort}_PCA.pdf

        extractEUR WLS -0.01 0 -0.005 0.005 -0.01 0.01 -0.01 0.01 
        extractEUR HRS2 -0.01 0 0 0.01 -0.008 0.01 -0.005 0.01

        plotPCs $cohort $dirOut/$cohort/sampleQC/${cohort}_EUR_1kG_HM3_PCs.eigenvec $dirOut/$cohort/sampleQC/${cohort}_EUR_PCA.pdf
    done
}

main

