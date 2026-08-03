#!/bin/bash

source $PGI_Repo/code/paths

#--------------------------------------------------------------------------------------------------------#

PGI(){
    cohort=$1
    method=$2
    gen=$3
    snps=$4
    ancestry=$5
    sBayesRsnps=$6

    eval gf_dir='$'gf_dir_${cohort}

    case $cohort in
	"UKB1" | "UKB2" | "UKB3")
		snpidtype="rs"
		;;
	"STRgsa" | "MCS" | "NCDS" | "BCS70" | "NSHD" | "MCTFR") 
		snpidtype="ChrPosIDhg38"
		;;
	*)
		snpidtype="ChrPosIDhg19"
		;;
    esac	

    case $cohort in
	"UKB1" | "UKB2" | "UKB3")
		valgf=$gf_plink2_UKB
		part=$(echo $cohort | sed 's/UKB//g')
		sample="$phen_dir_UKB/GWAS/partitions/UKB_part${part}_eid.txt"
		;;
	*)		
		valgf="${gf_dir}/plink2/${cohort}_chr[1:22]"
        sample="${gf_dir}/sampleQC/${cohort}_${ancestry}_FID_IID.txt"
		;;
    esac	

    if [[ $method == "LDpred" ]]
    then
        case $cohort in
        "UKB1" | "UKB2" | "UKB3")
            valbim="${gf_dir}/plink/HM3/UKB_HM3"
            ;;
	    *)		
            valbim="${gf_dir}/plink/HM3/${cohort}_HM3"
            ;;
        esac
    else
        valbim="NA"	
    fi

    if [[ $gen == "parent" ]]
    then
        case $cohort in 
        "ALSPAC" | "AH" | "UKB1" | "PSID" )
            phased=0
            ;;
	    *)
            phased=1
            ;;
        esac
    else
        phased="NA"
    fi

    if [[ "${cohorts_nogf[@]}" =~ $cohort  ]]
    then
        onlyweights=1
    else
        onlyweights=0
    fi

    rm -f $PGI_Repo/code/9_Scores/ss_single_${cohort}_${gen}_${ancestry}
    # Get list of sumstats: Pheno name on first column (e.g. SWB-Okbay), file path on second
    for pheno in $(cat $PGI_Repo/code/9_Scores/versions/version_single_${cohort}_${ancestry}); do
        if [[ $pheno != \#* ]]
        then

            if ls $PGI_Repo/derived_data/4_MTAG_single/$pheno/*_trait_formatted*  1> /dev/null 2>&1
            then 
                path="$PGI_Repo/derived_data/4_MTAG_single/$pheno/${pheno}_trait_formatted.txt"
            else
                path="$PGI_Repo/derived_data/4_MTAG_single/$pheno/${pheno}_trait_1_formatted.txt"
            fi

            if [[ $method == "LDpred" ]]; then
                pheno=$(echo "$pheno-single" | sed 's/[1-9]//g')
            elif [[ $method == "SBayesR" ]]; then
                pheno=$(echo "$pheno-single")
                path=$(echo $path | sed 's/formatted/formatted_SBayesR/g')
            fi
            echo $pheno $path >> $PGI_Repo/code/9_Scores/ss_single_${cohort}_${gen}_${ancestry}
        fi
    done
    sh $PGI_Repo/code/9_Scores/9.0_PGI.sh single $cohort $method $gen $snps $snpidtype $valgf $sample $ancestry $valbim $phased $onlyweights $sBayesRsnps
}


# Using --impute-n and --robust for CPD1, CPD4-6, EVERSMOKE1, EVERSMOKE4-6 
EURcohorts=(AH ALSPAC BCS70 Dunedin ELSA ERisk GS GSOEP HRS MCTFR MCS MIDUS NCDS NSHD PSID STRtwge STRpsych STRgsa Texas WLS UKB1 UKB2 UKB3 FinnGen EstBB TwinLife)
AFRcohorts=(AH1kG HRS MCS MIDUS PSID) #UKB
EAScohorts=(AH1kG HRS MCTFR MIDUS) #UKB
SAScohorts=(BCS70 MCS) #UKB
AMRcohorts=(AH1kG HRS MCTFR PSID) #UKB
cohorts_nogf=("UKB1" "UKB2" "UKB3" "UKB" "FinnGen" "EstBB" "TwinLife" "Dunedin")

for ancestry in EUR AFR EAS SAS AMR EUR
do
    eval "cohorts=(\${${ancestry}cohorts[@]})"
    for cohort in "${cohorts[@]}"
    do
        echo "==================================="
        echo COHORT: $cohort, ANCESTRY: $ancestry
        PGI $cohort SBayesR proband NA $ancestry 2.9m > $PGI_Repo/code/9_Scores/9.2_PGI_single_${cohort}_${ancestry}.log
        echo "==================================="
    done
done
wait

# HapMap3 SNPs only PGIs - EUR ancestry
for ancestry in EUR
do
    eval "cohorts=(\${${ancestry}cohorts[@]})"
    for cohort in "${cohorts[@]}"
    do
        echo "==================================="
        echo COHORT: $cohort, ANCESTRY: $ancestry, HapMap3 SNPs only
        PGI $cohort SBayesR proband NA $ancestry HM3 >> $PGI_Repo/code/9_Scores/9.2_PGI_single_${cohort}_${ancestry}_HM3.log
        echo "==================================="
    done
done
wait

# Parental PGIs - EUR ancestry
for cohort in WLS AH ALSPAC ERisk GS MCTFR MCS STRtwge STRgsa STRpsych Texas UKB1 
do 
    echo "==================================="
    echo COHORT: $cohort
    PGI $cohort SBayesR parent NA EUR 2.9m>> $PGI_Repo/code/9_Scores/9.2_PGI_single_${cohort}_parent.log
    echo "==================================="
done

# Restricted SNP PGIs (for direct effect - population association comparison)
sniparSNPs_WLS=$PGI_Repo/derived_data/13_Parental_imputation/QC/SNPlists/WLS/WLS.snps
PGI WLS SBayesR proband $sniparSNPs_WLS > $PGI_Repo/code/9_Scores/9.2_PGI_single_WLS_sniparSNPs.log

# Pre-RAP code
# sniparSNPs_UKB1=$PGI_Repo/derived_data/13_Parental_imputation/QC/SNPlists/UKB1/UKB1.snps
# PGI UKB3 SBayesR proband $sniparSNPs_UKB1 > $PGI_Repo/code/9_Scores/9.2_PGI_single_UKB3_UKB1sniparSNPs.log
