#!/bin/bash

#---------------------------------------------------------------------------------------#

filterVCF(){
    cohort=$1
    colname_Rsq=$2
    cutoff_Rsq=$3
    cutoff_AvgCall=$4
    cutoff_hwe=$5
    cutoff_maf=$6
    snpidtype=$7
    phased=$8

    eval gfIn='$'gf_orig_${cohort} 
    eval info='$'info_orig_${cohort}

    echo "filterVCF() Step 1: Getting list of SNPs with INFO>$cutoff_Rsq and AvgCall>$cutoff_AvgCall ".
    start1=$(date +%s)
    for chr in {1..22}
    do
        infoChr=$(echo "$info" | sed "s/\[1:22\]/${chr}/g")
        if [[ "$infoChr" = *.gz ]]
            then
                zcat "$infoChr" | awk -F"\t" -v colRsq=$colname_Rsq -v cRsq=$cutoff_Rsq -v cAvgC=$cutoff_AvgCall 'NR==1 { for (i=1; i<=NF; i++) { ix[$i] = i }; next } $ix[colRsq]>cRsq && $ix["AvgCall"]>cAvgC {print $1}' > $gfDir/parental/tmp/${cohort}_chr$chr.avgcall99.rsq99.snps &
            else
                awk -F"\t" -v colRsq=$colname_Rsq -v cRsq=$cutoff_Rsq -v cAvgC=$cutoff_AvgCall 'NR==1 { for (i=1; i<=NF; i++) { ix[$i] = i }; next } $ix[colRsq]>cRsq && $ix["AvgCall"]>cAvgC {print $1}' "$infoChr" > $gfDir/parental/tmp/${cohort}_chr$chr.avgcall99.rsq99.snps &
        fi
    done
    wait

    end1=$(date +%s)
    echo "filterVCF() Step 1: Done. (Time: $(( ($end1 - $start1)/60 )) minutes)"
    
    echo "filterVCF() Step 2: Filtering vcf files for INFO>$cutoff_Rsq , AvgCall>$cutoff_AvgCall, HWE P-value > $cutoff_hwe , MAF> $cutoff_maf ; dropping INDELs, multi-allelic SNPs; dropping non-EUR ancestry."
    start2=$(date +%s)
    i=0

    if [[ $phased == 0 ]]; then
        phasedFlag=$(echo -e "--recode\n")
    else
        phasedFlag=$(echo -e "--recode\n--phased")
    fi

    for chr in {1..22}
    do
        gfChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
        $vcftools --gzvcf "$gfChr" \
            $phasedFlag \
            --keep $gfDir/sampleQC/${cohort}_EUR_FID_IID.txt \
            --snps $gfDir/parental/tmp/${cohort}_chr$chr.avgcall99.rsq99.snps \
            --stdout \
            --min-alleles 2 \
            --max-alleles 2 \
            --remove-indels \
            --maf $cutoff_maf \
            --hwe $cutoff_hwe > $gfDir/parental/tmp/${cohort}_chr${chr}.vcf &
        let i=i+1

        if [[ $i == 5 ]]
        then
            wait
            i=0
        fi

    done
    wait

    if [[ $snpidtype == "Chr:BP:A1:A2" ]]
    then
        echo "filterVCF() Step 2.1: SNPIDs are in Chr:BP:A1:A2 format. Converting to Chr:BP format."
        start21=$(date +%s)
        i=0
        for chr in {1..22}
        do
            awk -F"\t" '$1~"#"{print;next} 
                {gsub("chr","",$3) ; split($3,a,":"); $3=a[1]":"a[2]; print}' OFS="\t" $gfDir/parental/tmp/${cohort}_chr${chr}.vcf > $gfDir/parental/tmp/${cohort}_chr${chr}_ChrPosID.vcf &
    
            let i=i+1

            if [[ $i == 5 ]]
            then
                wait
                i=0
            fi
        done
        wait

        for chr in {1..22}
        do
            mv $gfDir/parental/tmp/${cohort}_chr${chr}_ChrPosID.vcf $gfDir/parental/tmp/${cohort}_chr${chr}.vcf 
        done
        
        end21=$(date +%s)
        echo "filterVCF() Step 2.1: Done. (Time: $(( ($end21 - $start21)/60 )) minutes)"
    fi

    for chr in {1..22}
    do
        gzip $gfDir/parental/tmp/${cohort}_chr$chr.vcf &
    done
    wait

    end2=$(date +%s) 
    echo "filterVCF() Step 2: Done. (Time: $(( ($end2 - $start2)/60 )) minutes)"
}

#            --set-all-var-ids @:# \

filterBGEN(){
    cohort=$1
    colname_Rsq=$2
    cutoff_Rsq=$3
    cutoff_AvgCall=$4
    cutoff_hwe=$5
    cutoff_maf=$6
    snpidtype=$7

    mkdir -p $gfDir/parental/input/bgen

    eval gfIn='$'gf_orig_${cohort} 
    eval info='$'info_orig_${cohort}
    eval sample='$'sample_orig_${cohort} 

    echo "filterBGEN() Step 1.0: Getting list of SNPs with INFO>$cutoff_Rsq , AvgCall>$cutoff_AvgCall , HWE P-value > $cutoff_hwe , MAF> $cutoff_maf, excluding INDELs and multi-allelic SNPs."
    start1=$(date +%s)
    for chr in {1..22}
    do
        infoChr=$(echo "$info" | sed "s/\[1:22\]/${chr}/g")
        if [[ $infoChr = *.gz ]]
            then
                zcat $infoChr | awk -v colRsq=$colname_Rsq -v cRsq=$cutoff_Rsq -v cAvgC=$cutoff_AvgCall -v cHWE=$cutoff_hwe -v cMAF=$cutoff_maf ' !($1~"#"){k++} k==1{ for (i=1; i<=NF; i++) { ix[$i] = i }; next } 
                    seen[$1]++ || $ix[colRsq]<=cRsq || $ix["AvgCall"]<=cAvgC || $ix["HW_exact_p_value"]<=cHWE || $ix["minor_allele_frequency"]<=cMAF {print $1}' > $gfDir/parental/tmp/${cohort}_chr$chr.avgcall.rsq.maf.hwe.snps &
            else
                awk -v colRsq=$colname_Rsq -v cRsq=$cutoff_Rsq -v cAvgC=$cutoff_AvgCall -v cHWE=$cutoff_hwe -v cMAF=$cutoff_maf '!($1~"#"){k++} k==1{ for (i=1; i<=NF; i++) { ix[$i] = i }; next } 
                    seen[$1]++ || $ix[colRsq]<=cRsq || $ix["AvgCall"]<=cAvgC || $ix["HW_exact_p_value"]<=cHWE || $ix["minor_allele_frequency"]<=cMAF {print $1}' $infoChr > $gfDir/parental/tmp/${cohort}_chr$chr.avgcall.rsq.maf.hwe.snps &
        fi
    done
    wait

    echo "filterBGEN() Step 1.1: Getting list of non-EUR samples to drop."
    awk 'NR==FNR{FS="\t";a[$2]=$2;next} {FS=" "} !($2 in a) && FNR>2{print $1,$2}' $gfDir/sampleQC/${cohort}_EUR_FID_IID.txt $sample > $gfDir/parental/tmp/${cohort}_nonEUR.txt

    end1=$(date +%s)
    echo "filterBGEN() Step 1: Done. (Time: $(( ($end1 - $start1)/60 )) minutes)"
    
    echo "filterBGEN() Step 2: Filtering BGEN files for INFO>$cutoff_Rsq , AvgCall>$cutoff_AvgCall, HWE P-value > $cutoff_hwe , MAF> $cutoff_maf ; dropping INDELs, multi-allelic SNPs; dropping non-EUR ancestry."
    start2=$(date +%s)
    i=0
    for chr in {1..22}
    do
        gfChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")

        if ! [[ -f $gfChr ]]
        then 
            gfChr=$(echo "$gfIn" | sed "s/\[1:22\]/0${chr}/g")
        fi

        $qctool -g $gfChr \
            -s $sample \
            -excl-rsids $gfDir/parental/tmp/${cohort}_chr$chr.avgcall.rsq.maf.hwe.snps \
            -excl-samples $gfDir/parental/tmp/${cohort}_nonEUR.txt \
            -og $gfDir/parental/input/bgen/${cohort}_chr${chr}.bgen \
            -os $gfDir/parental/input/bgen/${cohort}_chr${chr}.sample &

        let i=i+1

        if [[ $i == 5 ]]
        then
            wait
            i=0
        fi    
    done
    wait

    end2=$(date +%s) 
    echo "filterBGEN() Step 2: Done. (Time: $(( ($end2 - $start2)/60 )) minutes)"
}


kingRel(){
    cohort=$1

    mkdir -p $gfDir/parental/input/ibd
    
    echo "kingRel(): Obtaining relationships with KING"
    start=$(date +%s)

    $king -b $gfDir/parental/tmp/${cohort}_chr1_22.bed \
        --related \
        --prefix $gfDir/parental/input/ibd/${cohort}

    end=$(date +%s) 
    echo "kingRel(): Done. (Time: $(( ($end - $start)/60 )) minutes)"
}

ibd(){
    cohort=$1

    echo "ibd() Step 1: Pruning per chr plink files"
    start1=$(date +%s)

    for chr in {1..22}
    do
        plink2 --bfile $gfDir/parental/input/plink/${cohort}_chr${chr} \
            --indep-pairwise 1000kb 1 0.3 \
            --out $gfDir/parental/tmp/${cohort}_chr${chr}_pruned 
    done
    wait

    for chr in {1..22}
    do
        plink2 --bfile $gfDir/parental/input/plink/${cohort}_chr${chr} \
            --extract $gfDir/parental/tmp/${cohort}_chr${chr}_pruned.prune.in \
            --make-bed \
            --out $gfDir/parental/tmp/${cohort}_chr${chr}_pruned 
    done
    wait

    end1=$(date +%s) 
    echo "ibd() Step 1: Done. (Time: $(( ($end1 - $start1)/60 )) minutes)"

    echo "ibd() Step 2: Obtaining ibd segments"
    start2=$(date +%s)

    N_PO=$(grep "PO" -c $gfDir/parental/input/ibd/${cohort}.kin0)

    # Note: Using p_error 0.001 for WLS because all the PO pairs in the data are HapMap controls
    if [[ $N_PO -gt 0 && ${cohort} != "WLS" && ${cohort} != "PSID" ]]
    then
        eval agesex='$'{phen_dir_${cohort}}/${cohort}.agesex

        ibd.py --bed $gfDir/parental/tmp/${cohort}_chr@_pruned \
            --king $gfDir/parental/input/ibd/${cohort}.kin0 \
            --agesex $agesex \
            --out $gfDir/parental/input/ibd/${cohort}_chr@ \
            --threads 10 \
            --ld_out
    else
        ibd.py --bed $gfDir/parental/tmp/${cohort}_chr@_pruned \
            --king $gfDir/parental/input/ibd/${cohort}.kin0 \
            --out $gfDir/parental/input/ibd/${cohort}_chr@ \
            --threads 10 \
            --p_error 0.001 \
            --ld_out
    fi

    end2=$(date +%s) 
    echo "ibd() Step 2: Done. (Time: $(( ($end2 - $start2)/60 )) minutes)"
}

impute(){
    cohort=$1
    phased=$2

    mkdir -p $gfDir/parental/output
    eval agesex='$'{phen_dir_${cohort}}/${cohort}.agesex
    start=$(date +%s)

    if  [[ $phased == 1 ]]
    then
        echo "impute(): Imputing parental genotypes from phased data"
        impute.py \
            --bgen $gfDir/parental/input/bgen/${cohort}_chr@ \
            --king $gfDir/parental/input/ibd/${cohort}.kin0 \
            --ibd $gfDir/parental/input/ibd/${cohort}_chr@.ibd \
            --agesex $agesex \
            --chr_range 1-22 \
            --out $gfDir/parental/output/${cohort}_parental_chr@      
    elif [[ $phased == 0 ]]  
    then
        echo "impute(): Imputing parental genotypes from unphased data"
        impute.py \
            --bed $gfDir/parental/input/plink/${cohort}_chr@ \
            --king $gfDir/parental/input/ibd/${cohort}.kin0 \
            --ibd $gfDir/parental/input/ibd/${cohort}_chr@.ibd \
            --agesex $agesex \
            --chr_range 1-22 \
            --out $gfDir/parental/output/${cohort}_parental_chr@
    fi     

    end=$(date +%s) 
    echo "impute(): Done. (Time: $(( ($end - $start)/60 )) minutes)"
    
}
