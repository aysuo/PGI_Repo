
source $PGI_Repo/code/paths
source $snipar_venv
source $PGI_Repo/code/13_Parental_imputation/13.1.0_parental_imputation.sh
source $PGI_Repo/code/7_Genotypes/7.6.0_formatConversion.sh

colRsq_MCS="Rsq"
rsq_MCS=0.99
avgCall_MCS=0.99
hwe_MCS=0.00001
snpidtype_MCS="Chr:BP:A1:A2"
phased_MCS=1

colRsq_WLS="R2"
rsq_WLS=0.99
avgCall_WLS=0
hwe_WLS=0.0001
snpidtype_WLS="ChrPosID"
phased_WLS=1

colRsq_ERisk="Rsq"
rsq_ERisk=0.99
avgCall_ERisk=0.99
hwe_ERisk=0.0001
snpidtype_ERisk="Chr:BP:A1:A2"
phased_ERisk=1

colRsq_Texas="Rsq"
rsq_Texas=0.99
avgCall_Texas=0.99
hwe_Texas=0.0001
snpidtype_Texas="ChrPosID"
phased_Texas=1

colRsq_ALSPAC="info"
rsq_ALSPAC=0.99
avgCall_ALSPAC=0
hwe_ALSPAC=0.00001
snpidtype_ALSPAC="ChrPosID"
phased_ALSPAC=0

colRsq_AH="R2"
rsq_AH=0.99
avgCall_AH=0
hwe_AH=0.0001
snpidtype_AH="ChrPosID"
phased_AH=0

colRsq_GS="Rsq"
rsq_GS=0.99
avgCall_GS=0.99
hwe_GS=0.00001
snpidtype_GS="Chr:BP:A1:A2"
phased_GS=1

colRsq_STRtwge="Rsq"
rsq_STRtwge=0.99
avgCall_STRtwge=0.99
hwe_STRtwge=0.0001
snpidtype_STRtwge="ChrPosID"
phased_STRtwge=1

colRsq_STRgsa="Rsq"
rsq_STRgsa=0.99
avgCall_STRgsa=0.99
hwe_STRgsa=0.00001
snpidtype_STRgsa="Chr:BP:A1:A2"
phased_STRgsa=1

colRsq_STRpsych="Rsq"
rsq_STRpsych=0.99
avgCall_STRpsych=0.99
hwe_STRpsych=0.00001
snpidtype_STRpsych="Chr:BP:A1:A2"
phased_STRpsych=1

colRsq_MCTFR="Rsq"
rsq_MCTFR=0.99
avgCall_MCTFR=0.99
hwe_MCTFR=0.0001
snpidtype_MCTFR="Chr:BP:A1:A2"
phased_MCTFR=1

colRsq_GSOEP="Rsq"
rsq_GSOEP=0.95
avgCall_GSOEP=0.95
hwe_GSOEP=0.0001
snpidtype_GSOEP="Chr:BP:A1:A2"
phased_GSOEP=1

# colRsq_PSID="R2"
# rsq_PSID=0.95
# avgCall_PSID=0.95
# hwe_PSID=0.0001
# snpidtype_PSID="Chr:BP:A1:A2"
# phased_PSID=0



for cohort in AH ALSPAC ERisk GS GSOEP MCTFR MCS STRtwge STRgsa STRpsych Texas WLS UKB1
do
    eval colRsq='$'{colRsq_${cohort}}
    eval rsq='$'{rsq_${cohort}}
    eval avgCall='$'{avgCall_${cohort}}
    eval hwe='$'{hwe_${cohort}}
    eval snpidtype='$'{snpidtype_${cohort}}
    eval phased='$'{phased_${cohort}}

    eval gfDir='$'{gf_dir_${cohort}}
    mkdir -p $gfDir/parental/tmp
    mkdir -p $gfDir/parental/input/bgen
    mkdir -p $gfDir/parental/input/plink

    if [[ $cohort == "ALSPAC" ]]
        then
            filterBGEN $cohort $colRsq $rsq $avgCall $hwe 0.01 $snpidtype
            bgen2plink "$gfDir/parental/input/bgen/${cohort}_chr[1:22].bgen" \
                $gfDir/parental/input/plink/${cohort} \
                $snpidtype \
                $gfDir/parental/input/bgen/${cohort}_chr1.sample
        else
            filterVCF $cohort $colRsq $rsq $avgCall $hwe 0.01 $snpidtype $phased
            wait
            if [[ $phased == 1 ]]
            then
                vcf2bgen "$gfDir/parental/tmp/${cohort}_chr[1:22].vcf.gz" "$gfDir/parental/input/bgen/${cohort}"
                wait
            fi        
            vcf2plink "$gfDir/parental/tmp/${cohort}_chr[1:22].vcf.gz" "$gfDir/parental/input/plink/${cohort}"
            wait
    fi
    mergePlink $gfDir/parental/input/plink/${cohort}_chr[1:22] $gfDir/parental/tmp/${cohort}_chr1_22
    wait
    kingRel $cohort
    wait
    ibd $cohort
    wait
    impute $cohort $phased
done
