#!/bin/bash

gf_HRS2="/disk/genetics/dbgap/data/HRS_Jun19_2018/HRS/files/62567/PhenoGenotypeFiles/RootStudyConsentSet_phs000428.CIDR_Aging_Omni1.v2.p2.c1.NPR/GenotypeFiles/phg000515.v1.HRS_phase123_imputation.genotype-imputed-data.c1/imputed/HRS2_chr[1:22].gprobs.gz"
sample_HRS2="/disk/genetics/dbgap/data/HRS_Jun19_2018/HRS2.sample"

gf_WLS="/disk/genetics/WLS_DBGAP/derived/gen/WLS_chr[1:22].gen.gz"
sample_WLS="/disk/genetics/WLS_DBGAP/derived/gen/WLS_chr[1:22].sample"

gf_AH="/disk/genetics/dbgap/data/addhealth/unpacked/phg001099.v1.AddHealth.genotype-imputed-data.c1.GRU-IRB-PUB-GSO.set2/chr[1:22].dbGaP.dose.vcf.gz"
gf_Dunedin="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/genotype_data/Dunedin/imputed/vcf/chr[1:22].dose.vcf.gz"
gf_ERisk="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/genotype_data/ERisk/imputed/vcf/chr[1:22].dose.vcf.gz"
gf_STRpsych="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/genotype_data/STR/imputed/psych/decrypted/chr[1:22].dose.vcf.gz"
gf_STRyatssstage="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/genotype_data/STR/imputed/yatss_stage/decrypted/chr[1:22].dose.vcf.gz"
gf_STRtwge="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/genotype_data/STR/imputed/twge/decrypted/chr[1:22].dose.vcf.gz"
gf_ELSA="/disk/genetics/ELSA/data/imputed_data/chr[1:22].dose.vcf.gz"
gf_EGCUT="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/genotype_data/EGCUT/imputed/EGCUT_chr[1:22].vcf.gz"
gf_Texas="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/genotype_data/Texas/imputed/chr[1:22].dose.vcf.gz"
gf_MCTFR="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/genotype_data/MCTFR/imputed/decrypted/chr[1:22].dose.vcf.gz"
gf_HRS3="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/genotype_data/HRS3/imputed/HRS.EA_chr[1:22].vcf.gz"

HM3_ChrPosID="/disk/genetics/ukb/aokbay/reffiles/HM3_SNPs_ChrPosID.txt"
HM3_rsID="/disk/genetics/ukb/aokbay/reffiles/w_hm3.snplist"

gf1000G="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/ref_data/1000G_ph3/ALL.chr[1:22].phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

dirOut="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/7_Genotypes"

#---------------------------------------------------------------------------------#

oxford2plink2(){
    gfIn=$1
    sampleIn=$2
    gfOut=$3

    for chr in {1..22}; do
        gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
        sampleInChr=$(echo "$sampleIn" | sed "s/\[1:22\]/${chr}/g")
        plink2 --gen $gfInChr \
            --sample $sampleInChr \
            --oxford-single-chr $chr \
            --make-pgen \
            --out ${gfOut}_chr$chr &
    done 
    wait
}

vcf2plink2(){
    gfIn=$1
    gfOut=$2
    snpidtype=$3

    if [[ $snpidtype == "ChrPosID" ]]; then
        for chr in {1..22}; 
        do
            gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
            plink2 --vcf $gfInChr \
                --make-pgen \
                --set-all-var-ids @:# \
                --out ${gfOut}_chr$chr &
        done
    else
        for chr in {1..22};
        do
            gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
            plink2 --vcf $gfInChr \
                --make-pgen \
                --out ${gfOut}_chr$chr &
        done
    fi

    wait
}

subsetHM3(){
    gfIn=$1
    snpidtype=$2
    gfOut=$3
    
    eval hm3='$'HM3_${snpidtype}

    for chr in {1..22}; do
        gfInChr=$(echo "$gfIn" | sed "s/\[1:22\]/${chr}/g")
        plink2 --pfile $gfInChr \
            --extract $hm3 \
            --max-alleles 2 \
            --make-bed \
            --out ${gfOut}_chr$chr &
    done
    wait
    
    for chr in {1..22}; do
        echo ${gfOut}_chr$chr >> ${gfOut}_merge
    done

    plink1.9 --merge-list ${gfOut}_merge --make-bed --out ${gfOut}

    if [[ -f ${gfOut}-merge.missnp ]]; then
        for chr in {1..22}; do
            plink2 --bfile ${gfOut}_chr$chr --exclude ${gfOut}-merge.missnp --make-bed --out ${gfOut}_chr$chr.tmp &
        done
        wait

        rm ${gfOut}_merge
        for chr in {1..22}; do
            echo ${gfOut}_chr$chr.tmp >> ${gfOut}_merge
        done

        plink1.9 --merge-list ${gfOut}_merge --make-bed --out ${gfOut}
    fi

    rm ${gfOut}_merge ${gfOut}_chr*
}

rs2chrpos(){
    gfIn=$1
    gfOut=$2

    plink2 --bfile $gfIn \
        --set-all-var-ids @:# \
        --make-bed \
        --out $gfOut
}

#---------------------------------------------------------------------------------#

main () {
    for cohort in HRS2 WLS AH Dunedin ERisk EGCUT MCTFR STRtwge STRpsych STRyatssstage ELSA; do
        mkdir -p $dirOut/$cohort/plink2
        mkdir -p $dirOut/$cohort/plink/HM3
    done

    #vcf2plink2 "${gf1000G}" $dirOut/1000G/plink2/1000Gph3 rsID
    #subsetHM3  $dirOut/1000G/plink2/1000Gph3_chr[1:22] rsID $dirOut/1000G/plink/HM3/1000Gph3_HM3
    #rs2chrpos $dirOut/1000G/plink/HM3/1000Gph3_HM3 $dirOut/1000G/plink/HM3/1000Gph3_HM3_chrpos

    for cohort in HRS2 WLS; do
        eval gf='$'gf_${cohort} 
        eval sample='$'sample_${cohort}
        oxford2plink2 "$gf" $sample $dirOut/$cohort/plink2/$cohort
        subsetHM3 $dirOut/$cohort/plink2/${cohort}_chr[1:22] rsID $dirOut/$cohort/plink/HM3/${cohort}_HM3
    done

    for cohort in AH Dunedin ERisk EGCUT MCTFR STRtwge STRpsych STRyatssstage ELSA; do
        eval gf='$'gf_${cohort} 
        vcf2plink2 "${gf}" $dirOut/$cohort/plink2/$cohort ChrPosID
        subsetHM3 $dirOut/$cohort/plink2/${cohort}_chr[1:22] ChrPosID $dirOut/$cohort/plink/HM3/${cohort}_HM3
    done

}

main