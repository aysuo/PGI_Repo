#!/bin/bash

source paths

#---------------------------------------------------------------------------------------#

PCA(){
    cohort=$1
    sampleKeep=$2
    snpid=$3

    eval gf_dir='$'gf_dir_${cohort}
    gf="$gf_dir/plink/HM3/${cohort}_HM3"
    
    if [[ $snpid == "rs" ]]
        then
            gf_1000G="$gf_out_1000G/plink/HM3/1000Gph3_HM3"
        else
            gf_1000G="$gf_out_1000G/plink/HM3/1000Gph3_HM3_chrpos"
    fi

    mkdir -p $gf_dir/sampleQC

    plink1.9 --bfile $gf_1000G \
        --bmerge $gf \
        --geno 0.99 \
        --make-bed \
        --out $gf_dir/sampleQC/${cohort}_1kG_HM3

    if [[ -f $gf_dir/sampleQC/${cohort}_1kG_HM3-merge.missnp ]]; then
        plink2 --bfile $gf_1000G \
            --exclude $gf_dir/sampleQC/${cohort}_1kG_HM3-merge.missnp \
            --make-bed \
            --out $gf_dir/sampleQC/1000G.tmp

        plink1.9 --bfile $gf_dir/sampleQC/1000G.tmp \
            --bmerge $gf \
            --geno 0.99 \
            --make-bed \
            --out $gf_dir/sampleQC/${cohort}_1kG_HM3
        
        rm $gf_dir/sampleQC/1000G.tmp*
    fi

    awk -v cohort=$cohort '$1=="0"{print $1,$2,"1000G";next} \
        {print $1,$2,cohort}' $gf_dir/sampleQC/${cohort}_1kG_HM3.fam > $gf_dir/sampleQC/${cohort}_1kG_HM3.clusters

    if [[ $sampleKeep == "NA" ]]; then
        plink1.9 --bfile $gf_dir/sampleQC/${cohort}_1kG_HM3 \
            --maf 0.01 \
            --within $gf_dir/sampleQC/${cohort}_1kG_HM3.clusters \
            --pca 4 \
            --pca-cluster-names 1000G \
            --out $gf_dir/sampleQC/${cohort}_1kG_HM3_PCs
        sed -i 's/ /\t/g' $gf_dir/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec
    else
        plink1.9 --bfile $gf_dir/sampleQC/${cohort}_1kG_HM3 \
            --keep $sampleKeep \
            --maf 0.01 \
            --within $gf_dir/sampleQC/${cohort}_1kG_HM3.clusters \
            --pca 4 \
            --pca-cluster-names 1000G \
            --out $gf_dir/sampleQC/${cohort}_EUR_1kG_HM3_PCs
        sed -i 's/ /\t/g' $gf_dir/sampleQC/${cohort}_EUR_1kG_HM3_PCs.eigenvec
    fi

    rm $gf_dir/sampleQC/*bed \
        $gf_dir/sampleQC/*bim \
        $gf_dir/sampleQC/*fam \
        $gf_dir/sampleQC/*nosex \
        $gf_dir/sampleQC/*-merge.missnp \
        $gf_dir/sampleQC/${cohort}_1kG_HM3.clusters
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
    margin=$2

    means=($(awk -F"\t" '$8=="EUR" && NR>1 {count++; sum1=sum1+$3; sum2=sum2+$4; sum3=sum3+$5; sum4=sum4+$6}
            END{print sum1/count, sum2/count, sum3/count, sum4/count}' $gf_dir/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated))
      
    limits=($(awk -F"\t" -v k=$margin -v mean1=${means[0]} -v mean2=${means[1]} -v mean3=${means[2]} -v mean4=${means[3]} '
            $8=="EUR" && NR>1 {count++; sum1=sum1+($3-mean1)^2;  sum2=sum2+($4-mean2)^2; sum3=sum3+($5-mean3)^2; sum4=sum4+($6-mean4)^2} 
            END{print mean1-k*sqrt(sum1/count), mean1+k*sqrt(sum1/count), mean2-k*sqrt(sum2/count), mean2+k*sqrt(sum2/count), mean3-k*sqrt(sum3/count), mean3+k*sqrt(sum3/count), mean4-k*sqrt(sum4/count), mean4+k*sqrt(sum4/count)}' $gf_dir/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated))
   
    awk -F"\t" -v l1=${limits[0]} -v u1=${limits[1]} -v l2=${limits[2]} -v u2=${limits[3]} -v l3=${limits[4]} -v u3=${limits[5]} -v l4=${limits[6]} -v u4=${limits[7]}  \
        '$1=="0" || ($3>l1 && $3<u1 && $4>l2 && $4<u2 && $5>l3 && $5<u3 && $6>l4 && $6<u4) {print $1,$2,$3,$4,$5,$6}' OFS="\t" $gf_dir/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated > $gf_dir/sampleQC/${cohort}_EUR_1kG_HM3_PCs.eigenvec
    
    awk -F"\t" '$1!=0{print $1,$2}' OFS="\t" $gf_dir/sampleQC/${cohort}_EUR_1kG_HM3_PCs.eigenvec > $gf_dir/sampleQC/${cohort}_EUR_FID_IID.txt
}


main(){
    for cohort in AH Dunedin EGCUT ELSA ERisk HRS2 HRS3 MCS MCTFR NCDS Texas STRpsych STRtwge STRyatssstage WLS; do
        if [[ $cohort == "HRS2" || $cohort == "WLS" ]]; then
            snpidtype="rs"
        else
            snpidtype="chrpos"
        fi 
        eval gf_dir='$'gf_dir_${cohort}
        PCA $cohort "NA" $snpidtype
        plotPCs $cohort $gf_dir/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec $gf_dir/sampleQC/${cohort}_PCA.pdf
        extractEUR $cohort 5
        plotPCs $cohort $gf_dir/sampleQC/${cohort}_EUR_1kG_HM3_PCs.eigenvec $gf_dir/sampleQC/${cohort}_EUR_PCA.pdf
    done
}

main

