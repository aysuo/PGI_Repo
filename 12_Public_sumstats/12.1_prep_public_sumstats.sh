#!/bin/bash

source $PGI_Repo/code/paths

cd $PGI_Repo/derived_data/12_Public_sumstats
mkdir -p single multi public

gunzip -c $HRC_EasyQC_afref > afref.tmp
sed 's/\t/ /g' $HRC_rsid2chrpos_map > cptref.tmp

# Clumping function
clump(){
    pheno=$1
    ss=$2
    version=$3

    for chr in {1..22}; do
        plink1.9 --bfile ${HRC_LDgf_full}_chr$chr \
        --clump $ss \
        --clump-snp-field SNPID \
        --clump-p1 0.01 \
        --clump-p2 0.01 \
        --clump-field PVALUE \
        --clump-r2 0.1 \
        --clump-kb 1000000 \
        --out $version/$pheno/chr$chr &
    done
    wait
  
    awk '{print $3,$1,$4,$5}' OFS="\t" $version/$pheno/chr*.clumped > $version/$pheno/${pheno}_lead_SNPs_p1e-2.txt
    if [[ $version == "single" ]]; then
        awk -F"\t" 'NR==FNR{a[$1]=$4;next}
            FNR==1{$8="BETA" ; $11="N" ; $12="GWASeqN";print;next} \
            ($2 in a && a[$2]<=5e-8){$7=sprintf("%.5f",$7) ; $8=sprintf("%.5f",$8) ; $9=sprintf("%.5f",$9) ; $10=sprintf("%3.2e",$10) ; print}' OFS="\t" $version/$pheno/${pheno}_lead_SNPs_p1e-2.txt $ss | cut --complement -f1 > $version/$pheno/${pheno}_${version}_p5e-8_sumstats.txt
    else
        awk -F"\t" 'NR==FNR{a[$1]=$4;next}
            FNR==1{$8="BETA" ; $12="GWASeqN" ; print;next} \
            ($2 in a && a[$2]<=5e-8){$7=sprintf("%.5f",$7) ; $8=sprintf("%.5f",$8) ; $9=sprintf("%.5f",$9) ; $10=sprintf("%3.2e",$10) ; print}' OFS="\t" $version/$pheno/${pheno}_lead_SNPs_p1e-2.txt $ss | cut --complement -f1,11 > $version/$pheno/${pheno}_${version}_p5e-8_sumstats.txt
    fi
    cat $version/$pheno/chr*.log > $version/$pheno/clump.log
    rm $version/$pheno/chr* $version/$pheno/clump.log
}

############################################################################

### Single-trait sumstats
for SSversion in clump gwide LDpred; do
    rm -f $PGI_Repo/code/12_Public_sumstats/ss_${SSversion}_single
    
    while read pheno; do
        mkdir -p single/$pheno
        if [[ $SSversion == "LDpred" ]]; then
            phenoNoNum=$(echo $pheno | sed 's/[1-9]//g')
            if [[ -f  $PGI_Repo/derived_data/9_Scores/single/weights/AH_${phenoNoNum}-single_weights_LDpred_p1.0000e+00.txt ]]; then
                maxNsnp=0
                for cohort in AH Dunedin ELSA ERisk Texas WLS; do
                    Nsnp=$(wc -l $PGI_Repo/derived_data/9_Scores/single/weights/${cohort}_${phenoNoNum}-single_weights_LDpred_p1.0000e+00.txt | cut -d" " -f1)
                    if [[ $Nsnp > $maxNsnp ]]; then
                        maxNsnp=$Nsnp
                        maxNsnpCohort=$cohort
                    fi
                done
                ss="$PGI_Repo/derived_data/9_Scores/single/weights/${maxNsnpCohort}_${phenoNoNum}-single_weights_LDpred_p1.0000e+00.txt"
            else
                ss="NA"
            fi
        elif ! [[ $pheno = *_excl_23andMe ]]; then
            trait=$(awk -F"\t" -v pheno=$pheno '$1==pheno{print $2}' ../5_LDSC/singleMTAG/ER2_table_full.txt)
            ss="$PGI_Repo/derived_data/4_MTAG_single/$pheno/${pheno}_trait_${trait}_formatted.txt"
        else 
            trait=1
            ss="$PGI_Repo/derived_data/12_Public_sumstats/single/$pheno/${pheno}_trait_${trait}_formatted.txt"
        fi
       
        if ! [[ $SSversion == "LDpred" ]] && ! [[ -f $ss ]]; then
            ss=$(echo $ss | sed 's/trait_[1-9]/trait/g')
        fi

        # Write version & sumstats into a file to keep record 
        echo -e "$pheno\t$ss" >> $PGI_Repo/code/12_Public_sumstats/ss_${SSversion}_single

        # clumped
        if ! [[ -f single/$pheno/${pheno}_lead_SNPs_p1e-2.txt ]] && [[ $SSversion == "clump" ]]; then
            echo "Clumping single $pheno trait $trait .."
            clump $pheno $ss single > single/${pheno}.log
            echo "Clumping for $pheno trait $trait finished."
        
        # genome-wide, raw
        elif ! [[ -f single/$pheno/${pheno}_single_gwide_sumstats.txt.gz ]] && [[ $SSversion == "gwide" ]]; then
            echo "Copying genome-wide sumstats for single $pheno trait $trait .."
            mkdir -p single/$pheno
            if ! [[ $pheno = *_excl_23andMe ]]; then
                cp $ss single/$pheno/${pheno}_single_gwide_sumstats.txt
            else
                mv $ss single/$pheno/${pheno}_single_gwide_sumstats.txt
                rm single/${pheno}/*formatted*
            fi
            awk -F"\t" 'NR==FNR{ref[$1]=$2 ; alt[$1]=$3 ; raf[$1]=$4;next} \
                FNR==1{$7="EAF_HRC"; $8="BETA"; $11="N"; $12="GWASeqN"; print; next} \
                $5 == ref[$1] && $6 == alt[$1] {$7=sprintf("%.5f",raf[$1])} \
                $5 == alt[$1] && $6 == ref[$1] {$7=sprintf("%.5f",1-raf[$1])} \
                !( $5 == ref[$1] && $6 == alt[$1]) && !($5 == alt[$1] && $6 == ref[$1]) {$7=sprintf("%.5f",$7)} 
                {$8=sprintf("%.5f",$8) ; $9=sprintf("%.5f",$9) ; $10=sprintf("%3.2e",$10);print}' OFS="\t" afref.tmp single/$pheno/${pheno}_single_gwide_sumstats.txt | cut --complement -f1 > single/$pheno/tmp
            mv single/$pheno/tmp single/$pheno/${pheno}_single_gwide_sumstats.txt
            gzip single/$pheno/${pheno}_single_gwide_sumstats.txt
            echo "Genome-wide sumstats for single $pheno trait $trait copied."
        
        # genome-wide, LDpred
        elif ! [[ -f $PGI_Repo/derived_data/12_Public_sumstats/single/${pheno}/${pheno}_single_LDpred_weights.txt.gz ]] && [[ $SSversion == "LDpred" ]]; then
            if ! [[ $ss == "NA" ]]; then
                echo "Copying LDpred weights for single $pheno .."

                awk 'NR==FNR{a[$2]=$1; next} \
                    FNR==1{print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","LDPRED_BETA"; next} \
                    $3 in a && $3 ~ ":" {$3=a[$3]} \
                    {gsub(/chrom_/,"",$1); print $3,$1,$2,$4,$5,$7}' OFS="\t" cptref.tmp $ss > $PGI_Repo/derived_data/12_Public_sumstats/single/${pheno}/${pheno}_single_LDpred_weights.txt
                gzip $PGI_Repo/derived_data/12_Public_sumstats/single/${pheno}/${pheno}_single_LDpred_weights.txt

            else
                echo "There are no LDpred weights for single $pheno, skipping.."
            fi            
         fi
    done < $PGI_Repo/code/12_Public_sumstats/version_${SSversion}_single
done    


##############################################################################

### Multi-trait sumstats
while read pheno; do
    mkdir -p multi/$pheno
    ss=$(echo $PGI_Repo/derived_data/6_MTAG_multi/$pheno/${pheno}_trait_1_formatted.txt)

    if ! [[ -f multi/$pheno/${pheno}_lead_SNPs_p1e-2.txt ]]; then
        echo "Clumping $pheno - multi.."
        clump $pheno $ss multi > multi/${pheno}.log
        echo "Clumping for $pheno - multi finished."
    fi
done < $PGI_Repo/code/12_Public_sumstats/version_clump_multi

##############################################################################

### Public sumstats
cat $PGI_Repo/code/12_Public_sumstats/../9_Scores/version_public_* | sort | uniq > $PGI_Repo/code/12_Public_sumstats/version_LDpred_public 
rm -f $PGI_Repo/code/12_Public_sumstats/ss_LDpred_public

while read pheno; do
    mkdir -p public/$pheno
    maxNsnp=0
    for file in $PGI_Repo/derived_data/9_Scores/public/weights/*_${pheno}_weights_LDpred_p1.0000e+00.txt; do
        Nsnp=$(wc -l $file | cut -d" " -f1)
        if [[ $Nsnp > $maxNsnp ]]; then
            maxNsnp=$Nsnp
            maxNsnpFile=$file
        fi
    done
    ss=$file

    # Write version & sumstats into a file to keep record 
    echo -e "$pheno\t$ss" >> $PGI_Repo/code/12_Public_sumstats/ss_${SSversion}_public
        
    # genome-wide, LDpred
    if ! [[ -f $PGI_Repo/derived_data/12_Public_sumstats/public/${pheno}/${pheno}_public_LDpred_weights.txt ]]; then
        echo "Copying LDpred weights for public $pheno .."
        awk 'NR==FNR{a[$2]=$1; next} \
            FNR==1{print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","LDPRED_BETA"; next} \
            $3 in a && $3 ~ ":" {$3=a[$3]} \
            {gsub(/chrom_/,"",$1); print $3,$1,$2,$4,$5,$7}' OFS="\t" cptref.tmp $ss > $PGI_Repo/derived_data/12_Public_sumstats/public/${pheno}/${pheno}_public_LDpred_weights.txt
    fi
done < $PGI_Repo/code/12_Public_sumstats/version_LDpred_public


##############################################################################
