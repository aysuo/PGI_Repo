#!/bin/bash

for cohort in AH Dunedin EGCUT ERisk ELSA HRS2 HRS3 MCTFR STRpsych STRtwge STRyatssstage Texas WLS; do
    declare gf_${cohort}_hm3="$PGI_Repo/derived_data/7_Genotypes/$cohort/plink/HM3/${cohort}_HM3"
done

gf_1000G_hm3_rs="$PGI_Repo/derived_data/7_Genotypes/1000G/plink/HM3/1000Gph3_HM3"
gf_1000G_hm3_chrpos="$PGI_Repo/derived_data/7_Genotypes/1000G/plink/HM3/1000Gph3_HM3_chrpos"
pop1000G="$PGI_Repo/original_data/ref_data/1000G_ph3/igsr_samples.tsv"

#---------------------------------------------------------------------------------------#


PCA(){
    cohort=$1
    sampleKeep=$2
    snpid=$3

    eval gf='$'gf_${cohort}_hm3
    eval gf_1000G='$'gf_1000G_hm3_${snpid}

    mkdir -p $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC

    plink1.9 --bfile $gf_1000G \
        --bmerge $gf \
        --geno 0.99 \
        --make-bed \
        --out $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3

    if [[ -f $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3-merge.missnp ]]; then
        plink2 --bfile $gf_1000G \
            --exclude $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3-merge.missnp \
            --make-bed \
            --out $gf_1000G.tmp

        plink1.9 --bfile $gf_1000G.tmp \
            --bmerge $gf \
            --geno 0.99 \
            --make-bed \
            --out $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3
        
        rm $gf_1000G.tmp*
    fi

    awk -v cohort=$cohort '$1=="0"{print $1,$2,"1000G";next} \
        {print $1,$2,cohort}' $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3.fam > $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3.clusters

    if [[ $sampleKeep == "NA" ]]; then
        plink1.9 --bfile $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3 \
            --maf 0.01 \
            --within $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3.clusters \
            --pca 4 \
            --pca-cluster-names 1000G \
            --out $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3_PCs
        sed -i 's/ /\t/g' $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec
    else
        plink1.9 --bfile $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3 \
            --keep $sampleKeep \
            --maf 0.01 \
            --within $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3.clusters \
            --pca 4 \
            --pca-cluster-names 1000G \
            --out $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_EUR_1kG_HM3_PCs
        sed -i 's/ /\t/g' $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_EUR_1kG_HM3_PCs.eigenvec
    fi

    rm $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/*bed \
        $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/*bim \
        $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/*fam \
        $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/*nosex \
        $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/*-merge.missnp \
        $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3.clusters
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
    Rscript $PGI_Repo/code/7_Genotypes/7.7.0_plotPCs.R "${PCs}.annotated" $out
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
        
    # Code to get min/max of 1kG EUR PCs if we want to automate the filter instead of visual inspection    
    #low1=$(awk -F"\t" 'BEGIN{min1=9999} $8=="EUR" && NR>1 {if ($3<min1) min1=$3}END{print min1}' $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated)
    #up1=$(awk -F"\t" 'BEGIN{max1=-9999} $8=="EUR" && NR>1 {if ($3>max1) max1=$3}END{print max1}' $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated)
    #low2=$(awk -F"\t" 'BEGIN{min2=9999} $8=="EUR" && NR>1 {if ($4<min2) min2=$4}END{print min2}' $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated)
    #up2=$(awk -F"\t" 'BEGIN{max2=-9999} $8=="EUR" && NR>1 {if ($4>max2) max2=$4}END{print max2}' $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated)
    #low3=$(awk -F"\t" 'BEGIN{min3=9999} $8=="EUR" && NR>1 {if ($5<min3) min3=$5}END{print min3}' $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated)
    #up3=$(awk -F"\t" 'BEGIN{max3=-9999} $8=="EUR" && NR>1 {if ($5>max3) max3=$5}END{print max3}' $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated)
    #low4=$(awk -F"\t" 'BEGIN{min4=9999} $8=="EUR" && NR>1 {if ($6<min4) min4=$6}END{print min4}' $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated)
    #up4=$(awk -F"\t" 'BEGIN{max4=-9999} $8=="EUR" && NR>1 {if ($6>max4) max4=$6}END{print max4}' $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated)
    
    awk -F"\t" -v l1=$low1 -v u1=$up1 -v l2=$low2 -v u2=$up2 -v l3=$low3 -v u3=$up3 -v l4=$low4 -v u4=$up4  \
        '$1=="0" || ($3>l1 && $3<u1 && $4>l2 && $4<u2 && $5>l3 && $5<u3 && $6>l4 && $6<u4) {print $1,$2,$3,$4,$5,$6}' OFS="\t" $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated > $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_EUR_1kG_HM3_PCs.eigenvec

    awk -F"\t" '$1!=0{print $1,$2}' OFS="\t" $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_EUR_1kG_HM3_PCs.eigenvec > $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_EUR_FID_IID.txt

}


main(){
    for cohort in AH Dunedin EGCUT ELSA ERisk HRS2 HRS3 MCTFR Texas STRpsych STRtwge STRyatssstage WLS; do
        if [[ $cohort == "HRS2" || $cohort == "WLS" ]]; then
            snpidtype="rs"
        else
            snpidtype="chrpos"
        fi 
        PCA $cohort "NA" $snpidtype
        plotPCs $cohort $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_PCA.pdf
    done

    extractEUR WLS -0.01 0 -0.005 0.005 -0.01 0.01 -0.01 0.01 
    extractEUR HRS2 -0.01 -0.0025 0 0.01 -0.006 0.005 -0.01 0
    extractEUR EGCUT -0.005 1 -1 0.005 -0.013 0.005 -0.01 0.01
    extractEUR ELSA -1 1 -0.005 1 -0.01 1 -1 1
    extractEUR MCTFR -1 0 0 1 -0.005 0.005 -0.005 0.01
    extractEUR STRtwge -1 1 -0.003 1 -1 1 -1 1
    extractEUR STRyatssstage -1 0.0025 -0.005 1 -0.01 1 -1 0.005
    extractEUR STRpsych  
    extractEUR Texas -1 1 0.015 1 0 0.015 -1 -0.005
    extractEUR Dunedin -1 1 0.01 1 -1 1 -1 1
    extractEUR ERisk 0 1 -1 -0.005 -0.005 0.01 -1 0
    extractEUR STRpsych -1 0.003 -1 0.005 -0.008 1 -0.005 0.005
    extractEUR AH -1 1 -1 1 -1 1 -1 1
    extractEUR HRS3 -1 1 -1 1 -1 1 -1 1
        
    for cohort in AH Dunedin EGCUT ELSA ERisk HRS2 HRS3 MCTFR Texas STRpsych STRtwge STRyatssstage WLS; do
        plotPCs $cohort $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_EUR_1kG_HM3_PCs.eigenvec $PGI_Repo/derived_data/7_Genotypes/$cohort/sampleQC/${cohort}_EUR_PCA.pdf
    done
}

main

