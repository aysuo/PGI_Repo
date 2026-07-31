#!/bin/bash


PCA(){
    cohort=$1       # Name of cohort    
    sampleKeep=$2   # List of individuals to restricti the sample
    snpid=$3        # Type of snpid: options are rsID, ChrPosIDhg19 (Chr:BP build37) or ChrPosIDhg38 (Chr:BP build37)
    out=$4          # Output file path

    # Path to cohort bed/bim/fam files with only HapMap3 SNPs
    gf="$gf_dir/plink/HM3/${cohort}_HM3"
    
    # Path to 1000 Genomes Phase 3 bed/bim/fam files with only HapMap3 SNPs, depending on type of snpid
    if [[ $snpid == "rsID" ]]
        then
            gf_1000G="$gf_dir_1000G/plink/HM3/1000Gph3_HM3"
    elif [[ $snpid == "ChrPosIDhg19" ]]
        then
            gf_1000G="$gf_dir_1000G/plink/HM3/1000Gph3_HM3_chrpos"
    elif [[ $snpid == "ChrPosIDhg38" ]]
        then
            gf_1000G="$gf_dir_1000G/plink/HM3/1000Gph3_HM3_chrpos_hg38"
    fi
    
    mkdir -p $gf_dir/sampleQC

    # Merge cohort HM3 plink files with 1000 Genomes, drop SNPs with missingness > 1% so that only overlapping SNPs remain. 
    plink1.9 --bfile $gf_1000G \
        --bmerge $gf \
        --geno 0.001 \
        --make-bed \
        --out $gf_dir/sampleQC/${cohort}_1kG_HM3

    if [[ -f $gf_dir/sampleQC/${cohort}_1kG_HM3-merge.missnp ]]; then
        plink2 --bfile $gf_1000G \
            --exclude $gf_dir/sampleQC/${cohort}_1kG_HM3-merge.missnp \
            --make-bed \
            --out $gf_dir/sampleQC/1000G.tmp

        plink1.9 --bfile $gf_dir/sampleQC/1000G.tmp \
            --bmerge $gf \
            --geno 0.001 \
            --make-bed \
            --out $gf_dir/sampleQC/${cohort}_1kG_HM3
        
        rm $gf_dir/sampleQC/1000G.tmp*
    fi

    # Using the fam file from merged dataset, create file with FID,IID and cluster name (i.e. 1000G or cohort name) 
    awk -v cohort=$cohort '$1=="0"{print $1,$2,"1000G";next} \
        {print $1,$2,cohort}' $gf_dir/sampleQC/${cohort}_1kG_HM3.fam > $gf_dir/sampleQC/${cohort}_1kG_HM3.clusters

    # Create 4 PCs using the 1000G cluster, project cohort individuals on those loadings
    if [[ $sampleKeep == "NA" ]]; then
        plink1.9 --bfile $gf_dir/sampleQC/${cohort}_1kG_HM3 \
            --maf 0.01 \
            --within $gf_dir/sampleQC/${cohort}_1kG_HM3.clusters \
            --pca 4 \
            --pca-cluster-names 1000G \
            --out $gf_dir/sampleQC/${cohort}_1kG_HM3_PCs
        sed -i 's/ /\t/g' $gf_dir/sampleQC/${out}.eigenvec
    else
        plink1.9 --bfile $gf_dir/sampleQC/${cohort}_1kG_HM3 \
            --keep $sampleKeep \
            --maf 0.01 \
            --within $gf_dir/sampleQC/${cohort}_1kG_HM3.clusters \
            --pca 4 \
            --pca-cluster-names 1000G \
            --out $gf_dir/sampleQC/${out}
        sed -i 's/ /\t/g' $gf_dir/sampleQC/${out}.eigenvec
    fi

    rm $gf_dir/sampleQC/*-merge.missnp \
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
    Rscript $PGI_Repo/code/7_Genotypes/7.7.1_plotPCs.R "${PCs}.annotated" $out
}

extractAncestry(){
    cohort=$1
    margin=$2 # How many SD's away from the mean of an ancestry still counts as that ancestry? 
    ancestry=$3 # EUR, AFR, AMR, etc.

    # Obtain mean of each PC within 1000G ancestry
    means=($(awk -F"\t" -v ancestry=$ancestry '$8==ancestry && NR>1 {count++; sum1=sum1+$3; sum2=sum2+$4; sum3=sum3+$5; sum4=sum4+$6}
            END{print sum1/count, sum2/count, sum3/count, sum4/count}' $gf_dir/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated))

    # Define limits for ancestry based on margin (mean +- $margin * SD)
    limits=($(awk -F"\t" -v ancestry=$ancestry -v k=$margin -v mean1=${means[0]} -v mean2=${means[1]} -v mean3=${means[2]} -v mean4=${means[3]} '
            $8==ancestry && NR>1 {count++; sum1=sum1+($3-mean1)^2;  sum2=sum2+($4-mean2)^2; sum3=sum3+($5-mean3)^2; sum4=sum4+($6-mean4)^2} 
            END{print mean1-k*sqrt(sum1/(count-1)), mean1+k*sqrt(sum1/(count-1)), mean2-k*sqrt(sum2/(count-1)), mean2+k*sqrt(sum2/(count-1)), mean3-k*sqrt(sum3/(count-1)), mean3+k*sqrt(sum3/(count-1)), mean4-k*sqrt(sum4/(count-1)), mean4+k*sqrt(sum4/(count-1))}' $gf_dir/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated))
    
    # Drop individuals that are outside the limits from PC file
    awk -F"\t" -v l1=${limits[0]} -v u1=${limits[1]} -v l2=${limits[2]} -v u2=${limits[3]} -v l3=${limits[4]} -v u3=${limits[5]} -v l4=${limits[6]} -v u4=${limits[7]}  \
        '$1=="0" || ($3>l1 && $3<u1 && $4>l2 && $4<u2 && $5>l3 && $5<u3 && $6>l4 && $6<u4) {print $1,$2,$3,$4,$5,$6}' OFS="\t" $gf_dir/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated > $gf_dir/sampleQC/${cohort}_${ancestry}_1kG_HM3_PCs.eigenvec


    # Drop 1000G samples to get a list (FID IID) of individuals belonging to the requested ancestry in the cohort
    awk -F"\t" '$1!=0{print $1,$2}' OFS="\t" $gf_dir/sampleQC/${cohort}_${ancestry}_1kG_HM3_PCs.eigenvec > $gf_dir/sampleQC/${cohort}_${ancestry}_FID_IID.txt
}


Mahalanobis(){
    cohort=$1
    delta_EUR=$2
    delta_AFR=$3
    delta_EAS=$4
    delta_AMR=$5
    delta_SAS=$6

    if [[ $cohort == "PSID" ]]
    then
        output_dir=$gf_dir/sampleQC/
    else
        output_dir=$gf_dir/sampleQC/mahalanobis/
        mkdir -p $gf_dir/sampleQC/mahalanobis
    fi

    if [[ $cohort == "AH1kG" ]]
    then 
        EURlist=$gf_dir_AH/sampleQC/AH_EUR_FID_IID.txt
    elif [[ $cohort == "WLS1kG" ]]
    then
        EURlist=$gf_dir_WLS/sampleQC/WLS_EUR_FID_IID.txt
    else 
        EURlist=$gf_dir/sampleQC/${cohort}_EUR_FID_IID.txt
    fi

    echo "Extracting individuals that do not belong to EUR ancestry based on standard PCA results."
    awk 'NR==FNR{a[$2];next} FNR==1 {print;next} !($2 in a) {print $0}' OFS="\t" $EURlist $gf_dir/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec.annotated > $gf_dir/sampleQC/mahalanobis/${cohort}_nonEUR_1kG_HM3_PCs.eigenvec.annotated

    echo "Calculating Mahalanobis distance for each individual in cohort $cohort to the 1000G EUR cluster based on the 4 PCs."
    Rscript $PGI_Repo/code/7_Genotypes/7.7.2_Mahalanobis.R $gf_dir/sampleQC/mahalanobis/${cohort}_nonEUR_1kG_HM3_PCs.eigenvec.annotated ${output_dir}/${cohort}_Mahalanobis.txt $delta_EUR $delta_AFR $delta_EAS $delta_AMR $delta_SAS
    Rscript $PGI_Repo/code/7_Genotypes/7.7.3_plot_Mahalanobis.R ${output_dir}/${cohort}_Mahalanobis.txt ${output_dir}/${cohort}_Mahalanobis $delta_EUR $delta_AFR $delta_EAS $delta_AMR $delta_SAS

}

extractAncestryMahalanobis(){
    cohort=$1
    ancestry=$2
    delta_EUR=$3
    delta_AFR=$4
    delta_EAS=$5
    delta_AMR=$6
    delta_SAS=$7
    admixture_filter=$8 # Whether to further filter individuals based on admixture results after Mahalanobis distance filtering (NA or number)


    if ! [[ -f $gf_dir/sampleQC/mahalanobis/${cohort}_Mahalanobis.txt ]]
    then
        echo "Mahalanobis distance file not found for cohort $cohort. Running Mahalanobis distance calculation first."
        Mahalanobis $cohort $delta_EUR $delta_AFR $delta_EAS $delta_AMR $delta_SAS
    fi
    
    # Filter for the requested ancestry
    awk -F"\t" -v anc=$ancestry '$9 == anc {print $1,$2}' OFS="\t" $gf_dir/sampleQC/mahalanobis/${cohort}_Mahalanobis.txt > $gf_dir/sampleQC/${cohort}_${ancestry}_FID_IID.txt

    if [[ admixture_filter != "NA" ]]
    then
        echo "Further filtering individuals based on admixture results with threshold $admixture_filter."

        case $ancestry in
            "EAS")
                anc_comp=2
                ;;
            "AMR")
                anc_comp=3
                ;;
            "SAS")
                anc_comp=4
                ;;
            "AFR")
                anc_comp=5
                ;;
        esac          
        anc_comp=$(($anc_comp+3))

        awk -F"\t" -v col=$anc_comp -v cutoff=$admixture_filter 'NR==FNR{a[$2]=$col;next}(a[$2] > cutoff){print $1,$2}' OFS="\t" $gf_dir/sampleQC/admixture/${cohort}_1kG_ancestry_proportions.txt $gf_dir/sampleQC/${cohort}_${ancestry}_FID_IID.txt > $gf_dir/sampleQC/${cohort}_${ancestry}_FID_IID.tmp
        
        # Count individuals that failed the criterion
        total_input=$(wc -l < $gf_dir/sampleQC/${cohort}_${ancestry}_FID_IID.txt)
        passed=$(wc -l < $gf_dir/sampleQC/${cohort}_${ancestry}_FID_IID.tmp)
        failed=$((total_input - passed))
        echo "Number of individuals failing criterion ($ancestry ancestry component > $cutoff): $failed"
        
        mv $gf_dir/sampleQC/${cohort}_${ancestry}_FID_IID.tmp $gf_dir/sampleQC/${cohort}_${ancestry}_FID_IID.txt
    fi

    # Create PC file for the ancestry
    awk -F"\t" -v anc=$ancestry 'NR==FNR{a[$2];next} $1=="0" || $2 in a {print $1,$2,$3,$4,$5,$6}' OFS="\t" $gf_dir/sampleQC/${cohort}_${ancestry}_FID_IID.txt $gf_dir/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec > $gf_dir/sampleQC/${cohort}_${ancestry}_1kG_HM3_PCs.eigenvec
}  

admixture(){
    cohort=$1
    K=$2

    mkdir -p $gf_dir/sampleQC/admixture
    cd $gf_dir/sampleQC/admixture

    # Run ADMIXTURE on the merged cohort + 1000G dataset with only HapMap3 SNPs, using 1000G as reference panel to estimate ancestry proportions in the cohort
    echo "Will run ADMIXTURE for cohort $cohort with K=$K"

    echo "First pruning SNPs for LD."
    plink1.9 --bfile $gf_dir/sampleQC/${cohort}_1kG_HM3 \
        --indep-pairwise 50 10 0.1 \
        --out $gf_dir/sampleQC/admixture/${cohort}_1kG_HM3_prune

    plink1.9 --bfile $gf_dir/sampleQC/${cohort}_1kG_HM3 \
        --extract $gf_dir/sampleQC/admixture/${cohort}_1kG_HM3_prune.prune.in \
        --make-bed \
        --out $gf_dir/sampleQC/admixture/${cohort}_1kG_HM3_pruned

    # Annotate PCs with ancestry info
    echo "Annotating 1kG samples with ancestry info."
    sed 's/ /\t/g' $gf_dir/sampleQC/admixture/${cohort}_1kG_HM3_pruned.fam > $gf_dir/sampleQC/admixture/${cohort}_1kG_HM3_pruned.fam.tab
    awk -F"\t" 'NR==FNR{a[$1]=$6;next} \
        $2 in a{print a[$2];next} \
        !($1~"#"){print "?"}' $pop1000G $gf_dir/sampleQC/admixture/${cohort}_1kG_HM3_pruned.fam.tab > $gf_dir/sampleQC/admixture/${cohort}_1kG_HM3_pruned.pop


    echo "Running ADMIXTURE with K=$K. This may take a while."
    $admixture --supervised -j10 --cv $gf_dir/sampleQC/admixture/${cohort}_1kG_HM3_pruned.bed $K | tee $gf_dir/sampleQC/admixture/${cohort}_admixture_K5.log
    
    echo "Processing ADMIXTURE output to get ancestry proportions for each individual in the cohort."
    paste $gf_dir/sampleQC/admixture/${cohort}_1kG_HM3_pruned.fam.tab \
        $gf_dir/sampleQC/admixture/${cohort}_1kG_HM3_pruned.pop \
        $gf_dir/sampleQC/admixture/${cohort}_1kG_HM3_pruned.5.Q \
        $gf_dir/sampleQC/admixture/${cohort}_1kG_HM3_pruned.5.P | sed 's/ /\t/g' | cut --complement -f3-6 > $gf_dir/sampleQC/admixture/tmp
        
    echo -e "FID\tIID\tPOP\tK1\tK2\tK3\tK4\tK5\tP1\tP2\tP3\tP4\tP5" | cat - $gf_dir/sampleQC/admixture/tmp > $gf_dir/sampleQC/admixture/${cohort}_1kG_ancestry_proportions.txt
    rm $gf_dir/sampleQC/admixture/tmp \
        $gf_dir/sampleQC/*bed \
        $gf_dir/sampleQC/*bim \
        $gf_dir/sampleQC/*fam \
        $gf_dir/sampleQC/*nosex \
        $gf_dir/sampleQC/${cohort}_1kG_HM3_prune*
        
    
    echo "Done with ADMIXTURE for cohort $cohort."
}

extractAncestryAdmixture(){
    cohort=$1
    ancestry=$2
    k=$3
    cutoff=$4

    if ! [[ -f $gf_dir/sampleQC/admixture/${cohort}_1kG_ancestry_proportions.txt ]]
    then 
        admixture $cohort $k
    fi
    echo ""
    echo "Extracting individuals with ${ancestry} ancestry from cohort $cohort based on ADMIXTURE results."
    
    echo -e "Ancestry\tMean_K1\tMean_K2\tMean_K3\tMean_K4\tMean_K5" > $gf_dir/sampleQC/admixture/Ancestry_mean_proportions_1kG.txt

    for i in EUR AFR EAS SAS AMR
    do
        echo "Obtaining mean of each ancestry component for 1kG $i ancestry."
        means=($(awk -F"\t" -v ancestry=$i '$3==ancestry{sum1=sum1+$4; sum2=sum2+$5; sum3=sum3+$6; sum4=sum4+$7 ; sum5=sum5+$8; k++} 
            END {print sum1/k, sum2/k, sum3/k, sum4/k, sum5/k}' OFS="\t" $gf_dir/sampleQC/admixture/${cohort}_1kG_ancestry_proportions.txt)) 
        echo -e "$i\t${means[0]}\t${means[1]}\t${means[2]}\t${means[3]}\t${means[4]}" >> $gf_dir/sampleQC/admixture/Ancestry_mean_proportions_1kG.txt
    done

    cat $gf_dir/sampleQC/admixture/Ancestry_mean_proportions_1kG.txt
    
    echo "Determining which ancestry component corresponds to ${ancestry} ancestry in the ADMIXTURE output."
    Kcol=$(awk -F"\t" -v ancestry=$ancestry '$1==ancestry{
        max=$2; col=2
        for(i=3;i<=NF;i++) if($i+0>max+0){max=$i; col=i}
        print col; exit
        }' $gf_dir/sampleQC/admixture/Ancestry_mean_proportions_1kG.txt)
    
    echo "$ancestry ancestry corresponds to component $(($Kcol-1)) in the ADMIXTURE output."
    
    Kcol=$(($Kcol+2)) # Column in the ancestry proportions file that corresponds to the ancestry of interest 


    if ! [[ -f $gf_dir/sampleQC/${cohort}_EUR_FID_IID.txt ]]; then
        echo "${cohort}_EUR_FID_IID.txt does not exist. Extracting EUR ancestry for cohort $cohort based on PCA results first to exclude from potential AMR ancestry."
        extractAncestry $cohort 5 EUR
    fi
    
    echo "Extracting non-EUR individuals with ${ancestry} ancestry component above $cutoff."
    awk -v K=$Kcol -v cutoff=$cutoff 'NR==FNR{a[$2];next} !($2 in a) && !($1==0) && $K>cutoff{print $1,$2}' OFS="\t" $gf_dir/sampleQC/${cohort}_EUR_FID_IID.txt $gf_dir/sampleQC/admixture/${cohort}_1kG_ancestry_proportions.txt > $gf_dir/sampleQC/admixture/${cohort}_${ancestry}_FID_IID.txt

    # Create PC file for individuals of the requested ancestry, to be plotted and checked manually to see if the ancestry assignment makes sense based on PCA results.
    awk -F"\t" 'NR==FNR{a[$2]=$2;next}($1==0 || $2 in a){print}' OFS="\t" $gf_dir/sampleQC/admixture/${cohort}_${ancestry}_FID_IID.txt $gf_dir/sampleQC/${cohort}_1kG_HM3_PCs.eigenvec > $gf_dir/sampleQC/admixture/${cohort}_${ancestry}_1kG_HM3_PCs.eigenvec
    
    echo "Done."
}
