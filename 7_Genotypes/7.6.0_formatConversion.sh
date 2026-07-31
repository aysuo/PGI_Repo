#!/bin/bash

source $PGI_Repo/code/paths

oxford2plink2(){
    gfIn=$1
    sampleIn=$2
    gfOut=$3

    for chr in {1..22}; do
        gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
        sampleInChr=$(echo "$sampleIn" | sed "s/\[1:22\]/${chr}/g")

        plink2 --gen $gfInChr ref-first \
            --sample $sampleInChr \
            --missing-code -9 \
            --oxford-single-chr $chr \
            --make-pgen \
            --set-all-var-ids @:# \
            --out ${gfOut}_chr$chr &
    done 
    wait
}


vcf2plink2(){
    gfIn=$1
    gfOut=$2
    snpidtype=$3

    i=0

    if [[ $snpidtype = ChrPosID* ]]; then
        for chr in {1..22} ; 
        do
            if ! [[ -f ${gfOut}_chr$chr.pgen ]]
            then
                gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
                plink2 --vcf "$gfInChr" \
                    --make-pgen \
                    --set-all-var-ids @:# \
                    --double-id \
                    --out ${gfOut}_chr$chr &
                let i=i+1

                if [[ $i == 5 ]]
                then
                    wait
                    i=0
                fi
            fi
        done

    else
        for chr in {1..22};
        do
            if ! [[ -f ${gfOut}_chr$chr.pgen ]]
            then
                gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
                plink2 --vcf $gfInChr \
                    --double-id \
                    --make-pgen \
                    --out ${gfOut}_chr$chr &

                let i=i+1

                if [[ $i == 5 ]]
                then
                    wait
                    i=0
                fi
            fi
        done
    fi
    wait
}

# Used in parental imputation
vcf2plink(){
    gfIn=$1
    gfOut=$2

    i=0
    
    echo "vcf2plink(): Converting vcf files to plink format."
    start1=$(date +%s)

    for chr in {1..22}
    do
        gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
        plink2 --vcf $gfInChr \
            --make-bed \
            --double-id \
            --rm-dup force-first \
            --out ${gfOut}_chr${chr} &
        
        let i=i+1
    
        if [[ $i == 5 ]]
        then
            wait
            i=0
        fi

    done
    wait

    end1=$(date +%s) 
    echo "vcf2plink(): Done. (Time: $(( ($end1 - $start1)/60 )) minutes)"
}


mergePlink(){
    gfIn=$1
    gfOut=$2

    echo "mergePlink(): Merging per chr plink files."
    start=$(date +%s) 

    rm -f ${gfOut}_mergelist
    for chr in {1..22}
    do
        gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
        echo $gfInChr >> ${gfOut}_mergelist
    done

    plink1.9 --merge-list ${gfOut}_mergelist \
        --make-bed \
        --out ${gfOut}

    if [[ -f ${gfOut}-merge.missnp ]]; then
        for chr in {1..22}; do
            gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
            plink2 --bfile $gfInChr \
                --exclude ${gfOut}-merge.missnp \
                --make-bed \
                --out $gfOut.chr$chr.tmp
        done
        #wait

        rm -f ${gfOut}_mergelist
        for chr in {1..22}; do
            echo $gfOut.chr$chr.tmp >> ${gfOut}_mergelist
        done

        plink1.9 --merge-list ${gfOut}_mergelist --make-bed --out ${gfOut}

        rm ${gfOut}.chr*.tmp*
    fi 
    rm ${gfOut}_chr*
    end=$(date +%s) 
    echo "mergePlink() Done. (Time: $(( ($end - $start)/60 )) minutes)"
}

# Used in parental imputation
vcf2bgen(){
    gfIn=$1
    gfOut=$2

    echo "vcf2bgen(): Converting vcf files to bgen format."
    start=$(date +%s)

    i=0
    for chr in {1..22}
    do
        gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")

        $qctool -g $gfInChr \
            -og ${gfOut}_chr${chr}.bgen \
            -os ${gfOut}_chr${chr}.sample &

        let i=i+1

        if [[ $i == 5 ]]
        then
            wait
            i=0
        fi    
    done
    wait

    end=$(date +%s) 
    echo "vcf2bgen(): Done. (Time: $(( ($end - $start)/60 )) minutes)"
}

bgen2plink2(){
    gfIn=$1
    gfOut=$2
    snpidtype=$3
    sampleIn=$4

    if [[ $snpidtype = ChrPosID* ]]; then
        for chr in {1..22}; 
        do
            if ! [[ -f ${gfOut}_chr$chr.pgen ]]
            then
                echo $gfIn
                if [[ $gfIn == *alspac* ]] && [[ $chr -lt 10 ]]
                then
                    gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/0${chr}/g")
                else
                    gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
                fi
                
                plink2 --bgen "$gfInChr" ref-unknown \
                    --sample $sampleIn \
                    --make-pgen \
                    --double-id \
                    --out ${gfOut}_chr$chr \
                    --set-all-var-ids @:# 
            fi
        done
    else
        for chr in {1..22};
        do
            if ! [[ -f ${gfOut}_chr$chr.pgen ]]
            then
                gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
                plink2 --bgen "$gfInChr" ref-unknown \
                    --sample $sampleIn \
                    --double-id \
                    --make-pgen \
                    --out ${gfOut}_chr$chr &
            fi
        done
    fi

    wait
}

# Used in parental imputation
bgen2plink(){
    gfIn=$1
    gfOut=$2
    snpidtype=$3
    sampleIn=$4
    sampleKeep=$5

    echo "bgen2plink(): Converting bgen files to plink format."
    start=$(date +%s)

    if [[ -z $sampleKeep ]]; then
        awk 'NR>2{print $1,$2}' $sampleIn > ${gfOut}_sampleKeep
        sampleKeep=${gfOut}_sampleKeep
    fi

    i=0
    if [[ $snpidtype = ChrPosID* ]]; then
        for chr in {1..22}; 
        do
            if ! [[ -f ${gfOut}_chr$chr.bed ]]
            then
                if [[ $gfIn == *alspac* ]] && [[ $chr -lt 10 ]]
                then
                    gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/0${chr}/g")
                else
                    gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
                fi
                
                plink2 --bgen "$gfInChr" ref-unknown \
                    --sample $sampleIn \
                    --keep $sampleKeep \
                    --make-bed \
                    --out ${gfOut}_chr$chr \
                    --set-all-var-ids @:# &

                let i=i+1

                if [[ $i == 5 ]]
                then
                    wait
                    i=0
                fi
            fi
        done
    else
        for chr in {1..22};
        do
            if ! [[ -f ${gfOut}_chr$chr.bed ]]
            then
                gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
                gfInChr=$(echo "$gfInChr" | sed "s/{1:22}/${chr}/g")

                plink2 --bgen "$gfInChr" ref-unknown \
                    --sample $sampleIn \
                    --keep $sampleKeep \
                    --make-bed \
                    --out ${gfOut}_chr$chr &

                let i=i+1

                if [[ $i == 5 ]]
                then
                    wait
                    i=0
                fi
            fi
        done
    fi

    wait

    rm -f ${gfOut}_sampleKeep

    end=$(date +%s) 
    echo "bgen2plink(): Done. (Time: $(( ($end - $start)/60 )) minutes)"
}


subsetSNPs(){
    gfIn=$1
    SNPlist=$2
    gfOut=$3
    

    for chr in {1..22}; do
        gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")

        if [[ -f $gfInChr.pgen ]]
        then
            plink2 --pfile $gfInChr \
                --extract $SNPlist \
                --max-alleles 2 \
                --make-bed \
                --out ${gfOut}_chr$chr 
        else
            plink2 --bfile $gfInChr \
                --extract $SNPlist \
                --max-alleles 2 \
                --make-bed \
                --out ${gfOut}_chr$chr 
        fi
    done
    wait
}


rs2chrpos(){
    gfIn=$1
    gfOut=$2

    plink2 --bfile $gfIn \
        --set-all-var-ids @:# \
        --make-bed \
        --out $gfOut
}

liftOver(){
    gfIn=$1
    gfOut=$2
    chain=$3

    awk '{print "chr"$1,$4-1,$4,$2}' OFS="\t" $gfIn.bim > $gfOut.lift.bed

    $liftover $gfOut.lift.bed $chain $gfOut.lifted.bed $gfOut.unlifted.bed
    awk -F"\t" '{gsub(/chr/,"",$1); print $4,$1":"$3,$1,$3}' $gfOut.lifted.bed > $gfOut.update.map 
    sed -i '/alt/d' $gfOut.update.map
    cut -f2 -d" " $gfOut.update.map > $gfOut.keep.snps

    plink1.9 --bfile $gfIn \
        --update-name $gfOut.update.map 2 1 \
        --make-bed \
        --out $gfOut.tmp

    plink1.9 --bfile $gfOut.tmp \
        --extract $gfOut.keep.snps \
        --update-chr $gfOut.update.map 3 2 \
        --update-map $gfOut.update.map 4 2 \
        --make-bed \
        --out $gfOut
        
    rm $gfOut.tmp* $gfOut.update.map $gfOut.keep.snps $gfOut.lift.bed $gfOut.lifted.bed
}

kgp2rs(){
    gfIn=$1
    gfOut=$2

    awk -F"\t" '$5=="TRUE" && $6=="TRUE"{print $1,$2}' OFS="\t" $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv > $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.both
    awk -F"\t" '$5=="TRUE" && $6=="FALSE"{print $1,$2}' OFS="\t" $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv > $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.illumina
    awk -F"\t" '$5=="FALSE" && $6=="TRUE"{print $1,$2}' OFS="\t" $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv > $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.gcc

    awk -F"\t" 'NR==FNR{a[$1]=$1;print;next}!($1 in a){print}' OFS="\t" $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.both $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.illumina > $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.both.illumina
    awk -F"\t" 'NR==FNR{a[$1]=$1;print;next}!($1 in a){print}' OFS="\t" $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.both.illumina $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.gcc > $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.preferred
    awk -F"\t" '!seen[$1]++{print}' OFS="\t" $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.preferred > $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.preferred.nodups
    rm $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.both* $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.illumina $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.gcc $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.preferred

    for chr in {1..22}; do
        plink2 --pfile ${gfIn}_chr$chr \
            --update-name $PGI_Repo/original_data/ref_data/SNP_kgpID2rsID.csv.preferred.nodups 2 1 \
            --make-pgen \
            --out ${gfOut}_chr$chr
    done
}
