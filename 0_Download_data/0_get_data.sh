#!/bin/bash

#------------------------------------------------------------------------------------------------------------#
# ------------------------------------------ For Repository PGI----------------------------------------------#
#------------------------------------------------------------------------------------------------------------#
cd $PGI_Repo/original_data/public

#------------------------------------#
# Barban et al - Fertility
wget http://sociogenome.com/material/GWASresults/AgeFirstBirth_Pooled.txt.gz
wget http://sociogenome.com/material/GWASresults/NumberChildrenEverBorn_Female.txt.gz
wget http://sociogenome.com/material/GWASresults/NumberChildrenEverBorn_Male.txt.gz
wget http://ssgac.org/documents/readme_reproductivebehavior.txt
#------------------------------------#

#------------------------------------#
# Day et al - Age at first menses
wget https://www.reprogen.org/Menarche_1KG_NatGen2017_WebsiteUpload.zip
#------------------------------------#

#------------------------------------#
# Demontis et al - ADHD
# Download from here: https://www.med.unc.edu/pgc/download-results/adhd/
#------------------------------------#

#------------------------------------#
# Doherty - Physical activity
# Download from here: https://ora.ox.ac.uk/objects/uuid:ff479f44-bf35-48b9-9e67-e690a2937b22
#------------------------------------#

#------------------------------------#
# Duncan - Anorexia Nervosa (excl UKB)
wget https://figshare.com/ndownloader/files/28169259 -O ANOREX-Duncan.txt.gz
# ReadMe
wget https://figshare.com/ndownloader/files/28169256 -O ANOREX-Duncan_ReadMe.pdf
#------------------------------------#

#------------------------------------#
# Evangelou - Blood pressure

# ICBP+UKB
# Systolic
# Data: 
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006624/Evangelou_30224653_SBP.txt.gz
# ReadMe: 
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006624/README_Evangelou_30224653.txt

# Diastolic
# Data: 
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006630/Evangelou_30224653_DBP.txt.gz
# ReadMe: 
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006630/README_Evangelou_30224653.txt

# Pulse pressure
# Data: 
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006629/Evangelou_30224653_PP.txt.gz
# ReadMe: 
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006629/README_Evangelou_30224653.txt

# ICBP-only sumstats obtained by email
# ICBP_DBP_02082017.txt.gz, ICBP_PP_02082017.txt.gz, ICBP_SBP_02082017.txt.gz
#------------------------------------#

#------------------------------------#
# Ferreira et al - Asthma/Eczema/Rhinitis
wget https://genepi.qimr.edu.au/staff/manuelF/gwas_results/SHARE-without23andMe.LDSCORE-GC.SE-META.v0.gz
#------------------------------------#

#------------------------------------#
# Gormley-Migraine obtained by email 
#------------------------------------#

#------------------------------------#
# Grove - Autism
wget https://figshare.com/ndownloader/files/28169292 -O ASD-Grove.txt.gz
# ReadMe
wget https://figshare.com/ndownloader/files/28169289 -O ASD-Grove_ReadMe.pdf
#------------------------------------#

#------------------------------------#
# Hautakangas-Migraine obtained by email
#------------------------------------#

#------------------------------------#
# Hysi - Nearsightedness
wget ftp://twinr-ftp.kcl.ac.uk/Refractive_Error_MetaAnalysis_2020/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020.txt.gz
wget ftp://twinr-ftp.kcl.ac.uk/Refractive_Error_MetaAnalysis_2020/ReadMe.txt -O NEARSIGHTED-Hysi_ReadMe.txt
#------------------------------------#

#------------------------------------#
# Jansen - Insomnia
wget https://ctg.cncr.nl/documents/p1651/Insomnia_sumstats_Jansenetal.txt.gz
#------------------------------------#

#------------------------------------#
# Kunkle et al - Alzheimer's
wget https://www.niagads.org/system/tdf/public_docs/Kunkle_etal_2019_IGAP_summary_statistics_README_0.docx?file=1&type=field_collection_item&id=120&force=
wget https://www.niagads.org/system/tdf/public_docs/Kunkle_etal_Stage1_results.txt?file=1&type=field_collection_item&id=121&force=
wget https://www.niagads.org/system/tdf/public_docs/Kunkle_etal_Stage2_results.txt?file=1&type=field_collection_item&id=122&force=
#------------------------------------#

#------------------------------------#
# Liu et al 
# Age started smoking
wget https://conservancy.umn.edu/bitstream/handle/11299/201564/AgeofInitiation.txt.gz
# W/o UKB
wget https://conservancy.umn.edu/bitstream/handle/11299/201564/AgeOfInitiation.WithoutUKB.txt.gz

# Cigarettes per day
wget https://conservancy.umn.edu/bitstream/handle/11299/201564/CigarettesPerDay.txt.gz?sequence=31&isAllowed=y
# W/o UKB
wget https://conservancy.umn.edu/bitstream/handle/11299/201564/CigarettesPerDay.WithoutUKB.txt.gz

# Drinks per week
wget https://conservancy.umn.edu/bitstream/handle/11299/201564/DrinksPerWeek.txt.gz?sequence=32&isAllowed=y
# W/o UKB
wget https://conservancy.umn.edu/bitstream/handle/11299/201564/DrinksPerWeek.WithoutUKB.txt.gz

# Ever smoker
wget https://conservancy.umn.edu/bitstream/handle/11299/201564/SmokingInitiation.txt.gz?sequence=34&isAllowed=y
# W/o UKB
wget https://conservancy.umn.edu/bitstream/handle/11299/201564/SmokingInitiation.WithoutUKB.txt.gz

# Smoking cessation
wget https://conservancy.umn.edu/bitstream/handle/11299/201564/SmokingCessation.txt.gz
# W/o UKB
wget https://conservancy.umn.edu/bitstream/handle/11299/201564/SmokingCessation.WithoutUKB.txt.gz
#------------------------------------#

#------------------------------------#
# Liu et al. - Inflammatory bowel disease
wget ftp://ftp.sanger.ac.uk/pub/consortia/ibdgenetics/iibdgc-trans-ancestry-filtered-summary-stats.tgz
#------------------------------------#

#------------------------------------#
# Locke et al - BMI
# https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files
wget https://portals.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz
#------------------------------------#

#------------------------------------#
# Mahajan et al. - T2D (unadjusted for BMI)
# Download from http://www.diagram-consortium.org/downloads.html
#------------------------------------#

#------------------------------------#
# Michailidou - Breast cancer
# Obtained by email (database was down)
# oncoarray_bcac_public_release_oct17.txt.gz
#------------------------------------#

#------------------------------------#
# Mills et al 
# Age First Birth
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90000001-GCST90001000/GCST90000050/GCST90000050_buildGRCh37.tsv.gz -O AFB-Mills.txt.gz

# Age First Sex
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90000001-GCST90001000/GCST90000047/GCST90000047_buildGRCh37.tsv.gz
#------------------------------------#

#------------------------------------#
# Mullins - Bipolar
wget https://figshare.com/ndownloader/files/26603681 -O BIPOLAR-Mullins.txt.gz
#------------------------------------#

#------------------------------------#
# Nelson et al - CAD
wget http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz
#------------------------------------#

#------------------------------------#
# Nikpey et al. CAD
wget http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/cad.additive.Oct2015.pub.zip
#------------------------------------#

#------------------------------------#
# Pasman et al - Cannabis (lifetime use)
# Download from here: https://www.ru.nl/bsi/research/group-pages/substance-use-addiction-food-saf/vm-saf/genetics/international-cannabis-consortium-icc/
wget https://www.ru.nl/publish/pages/898181/cannabis_readme_1.docx
#------------------------------------#

#------------------------------------#
# Perry et al - Age at first menses (excludes UKB)
wget https://www.reprogen.org/Menarche_Nature2014_GWASMetaResults_17122014.zip
#------------------------------------#

#------------------------------------#
# PGC3 - Schizophrenia 
wget https://figshare.com/ndownloader/files/28169757 -O SCZ-PGC3.txt.gz
#------------------------------------#

#------------------------------------#
# Schumacher et al. - Prostate cancer
wget http://practical.icr.ac.uk/blog/wp-content/uploads/uploadedfiles/oncoarray/MetaSummaryData/meta_v3_onco_euro_overall_ChrAll_1_release.zip
#------------------------------------#

#------------------------------------#
# Scott et al. - T2D (unadjusted for BMI)
# Download from http://www.diagram-consortium.org/downloads.html
#------------------------------------#

#------------------------------------#
# Sinnott-Armstrong et al. - HDL cholesterol, LDL cholesterol, total cholesterol, tryglycerids  
# HDL cholesterol
wget https://nih.figshare.com/ndownloader/files/22770203 -O CHOHDL-SinnottArmstrong.txt.gz
# LDL cholesterol
wget https://nih.figshare.com/ndownloader/files/22770218 -O CHOLDL-SinnottArmstrong.txt.gz
# Total cholesterol
wget https://nih.figshare.com/ndownloader/files/22770140 -O CHOTOT-SinnottArmstrong.txt.gz
# Tryglycerids
wget https://nih.figshare.com/ndownloader/files/22770296 -O TRYGL-SinnottArmstrong.txt.gz
#------------------------------------#

#------------------------------------#
# Stahl - Biploar (excl UKB)
wget https://figshare.com/ndownloader/files/28169307 -O BIPOLAR-Stahl.txt.gz
#------------------------------------#

#------------------------------------#
# Stringer et al - Cannabis (lifetime use, excludes UKB)
# Download from here: https://www.ru.nl/bsi/research/group-pages/substance-use-addiction-food-saf/vm-saf/genetics/international-cannabis-consortium-icc/
wget https://www.ru.nl/publish/pages/898181/paper2_readme.docx
#------------------------------------#

#------------------------------------#
# TAG - Cigarettes per day
# Download from here: https://www.med.unc.edu/pgc/download-results/tag/
#------------------------------------#

#------------------------------------#
# TAG - Ever smoker
# Download from here: https://www.med.unc.edu/pgc/download-results/tag/
#------------------------------------#

#------------------------------------#
# Tsuo - Asthma
# EUR: 
wget https://gbmi-sumstats.s3.amazonaws.com/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz  -O Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz
# ALL: 
wget https://gbmi-sumstats.s3.amazonaws.com/Asthma_Bothsex_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz  -O Asthma_Bothsex_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz
#------------------------------------#

#------------------------------------#
# Tsuo - COPD
# EUR: 
wget https://gbmi-sumstats.s3.amazonaws.com/COPD_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz  -O COPD_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz
# ALL: 
wget https://gbmi-sumstats.s3.amazonaws.com/COPD_Bothsex_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz  -O COPD_Bothsex_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz
#------------------------------------#

#------------------------------------#
# Van den Berg et al - Extraversion
wget https://www.dropbox.com/s/bk2jn41vrfl3zna/GPC-2.EXTRAVERSION.zip?dl=0
#------------------------------------#

#------------------------------------#
# Walters - Alcohol dependence
# ReadMe
wget https://figshare.com/ndownloader/files/28169775 -O AUDIT-Walters_ReadMe.txt
# Data
wget https://figshare.com/ndownloader/files/28169778 -O AUDIT-Walters.tar.gz
#------------------------------------#

#------------------------------------#
# Watson - Anorexia Nervosa
wget https://figshare.com/ndownloader/files/28169271 -O ANOREX-Watson.txt.gz
# ReadMe
wget https://figshare.com/ndownloader/files/28169268 -O ANOREX-Watson_ReadMe.pdf
#------------------------------------#

#------------------------------------#
# Wightman - Alzheimer
wget https://ctg.cncr.nl/documents/p1651/PGCALZ2sumstatsExcluding23andMe.txt.gz
#------------------------------------#

#------------------------------------#
# Wood et al - Height
# https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files
wget https://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz
#------------------------------------#

#------------------------------------#
# Zhang - Breast cancer
# Obtained by email (database was down)
# icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt.zip
# icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt.zip
# data_dictionary_summary_stats_bc_risk_2020.xlsx
#------------------------------------#

#------------------------------------#
# In-house results 
# sh get_data_inhouse.sh
#------------------------------------#


#------------------------------------------------------------------------------------------------------------#
# -------------------------------------------- For Public PGI------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#
cd $PGI_Repo/original_data/public_scores
#------------------------------------#
# ADHD
# See above
#------------------------------------#

#------------------------------------#
# Asthma/Eczema/Rhinitis
# See above
#------------------------------------#

#------------------------------------#
# Barban - Age first birth 
# See above
#------------------------------------#

#------------------------------------#
# Barban et al - NEBwomen
# See above
#------------------------------------#

#------------------------------------#
# Day et al - Age first menses
# See above
#------------------------------------#

#------------------------------------#
# De Moor et al - Openness
wget https://www.dropbox.com/s/df3bifmk238ks2y/GPC-1.BigFiveNEO.zip?dl=0 -O GPC-1.BigFiveNEO.zip
#------------------------------------#

#------------------------------------#
# Doherty et al - Physical activity
# Download from here: https://ora.ox.ac.uk/objects/uuid:ff479f44-bf35-48b9-9e67-e690a2937b22
#------------------------------------#

#------------------------------------#
# Jones et al - Morning person
wget https://personal.broadinstitute.org/mvon/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.gz
wget https://s3.amazonaws.com/broad-portal-resources/sleep/chronotype_raw_README.txt
#------------------------------------#

#------------------------------------#
# Harris et al - Self-rated health
wget https://datashare.is.ed.ac.uk/bitstream/handle/10283/3342/Harris2016_UKB_self_rated_health_summary_results_10112016.txt?sequence=2&isAllowed=y
wget https://datashare.is.ed.ac.uk/bitstream/handle/10283/3342/Harris2016_README.txt?sequence=1&isAllowed=y
#------------------------------------#

#------------------------------------#
# Howard et al - Depression
wget https://datashare.is.ed.ac.uk/bitstream/handle/10283/3203/PGC_UKB_depression_genome-wide.txt?sequence=3&isAllowed=y
wget https://datashare.is.ed.ac.uk/bitstream/handle/10283/3203/ReadMe.txt?sequence=4&isAllowed=y
#------------------------------------#

#------------------------------------#
# Lee et al - EA
wget https://www.dropbox.com/s/ho58e9jmytmpaf8/GWAS_EA_excl23andMe.txt?dl=0
#------------------------------------#

#------------------------------------#
# Lee et al - Intelligence
wget https://www.dropbox.com/s/ibjoh0g5e3sdd8t/GWAS_CP_all.txt?dl=0
#------------------------------------#

#------------------------------------#
# Linner et al - Drinks per week
wget https://www.dropbox.com/s/7hjxdhlxlwa482n/DRINKS_PER_WEEK_GWAS.txt?dl=0
#------------------------------------#

#------------------------------------#
# Linner et al - Ever-smoker
wget https://www.dropbox.com/s/o7wgwhnhjgt3eyn/EVER_SMOKER_GWAS_MA_UKB%2BTAG.txt?dl=0
#------------------------------------#

#------------------------------------#
# Linner et al - Risk
wget https://www.dropbox.com/s/il1d7vabk5283dm/RISK_GWAS_MA_UKB%2Breplication.txt?dl=0
#------------------------------------#

#------------------------------------#
# Liu et al - Cigarettes per day 
# See above
#------------------------------------#

#------------------------------------#
# Nagel et al - Neuroticism
wget https://ctg.cncr.nl/documents/p1651/sumstats_neuro_sum_ctg_format.txt.gz
wget https://ctg.cncr.nl/documents/p1651/readme_neuro_items_ctg
#------------------------------------#

#------------------------------------#
# Neale Lab

# Pheno description files
wget https://www.dropbox.com/s/d4mlq9ly93yhjyt/phenotypes.both_sexes.tsv.bgz -O phenotypes.both_sexes.tsv.gz
wget https://www.dropbox.com/s/r7idtoulxfpyjss/phenotypes.female.tsv.bgz -O phenotypes.female.tsv.gz
wget https://www.dropbox.com/s/mywldevz4nsla2r/phenotypes.male.tsv.bgz -O phenotypes.male.tsv.gz

# Variant annotation file
wget https://www.dropbox.com/s/puxks683vb0omeg/variants.tsv.bgz?dl=0 -O variants.tsv.bgz

# Asthma (multiple versions were available, chose the GWAS with largest Neff)
# N_controls=319207, N_cases=41934 --> N_eff=148259.282
wget https://www.dropbox.com/s/kp9bollwekaco0s/20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz

# Cannabis
wget https://www.dropbox.com/s/m6frhnk184rd7uz/20453.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 20453.gwas.imputed_v3.both_sexes.tsv.bgz &

# COPD
wget https://www.dropbox.com/s/hxx4huusrq7nwqq/22130.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 22130.gwas.imputed_v3.both_sexes.tsv.bgz

# Family satisfaction
wget https://www.dropbox.com/s/tsn2m8uzgwnkkr6/4559.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 4559.gwas.imputed_v3.both_sexes.tsv.bgz

# Financial satisfaction
wget https://www.dropbox.com/s/pezhe9cgcfehv9a/4581.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 4581.gwas.imputed_v3.both_sexes.tsv.bgz

# Friend satisfaction
wget https://www.dropbox.com/s/274tmv6f3q3z4q1/4570.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 4570.gwas.imputed_v3.both_sexes.tsv.bgz

# Hayfever (multiple versions were available, chose the GWAS with largest Neff)
# Hayfever/eczema (N_controls=277120, N_cases=83407 --> 256444.0149)
wget https://www.dropbox.com/s/strqg2hhhgl09rj/6152_9.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 6152_9.gwas.imputed_v3.both_sexes.tsv.bgz

# Loneliness
wget https://www.dropbox.com/s/nf4jl3mdppu1ng8/2020.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 2020.gwas.imputed_v3.both_sexes.tsv.bgz

# Migraine
wget https://www.dropbox.com/s/31xvr94v67qbt7e/G43.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O G43.gwas.imputed_v3.both_sexes.tsv.bgz

# Nearsightedness
wget https://www.dropbox.com/s/gmh9g65q5iw93ee/6147_1.gwas.imputed_v3.both_sexes.v2.tsv.bgz?dl=0 -O 6147_1.gwas.imputed_v3.both_sexes.v2.tsv.bgz 

# NEBmen
wget https://www.dropbox.com/s/uhj063tukv01d8m/2405.gwas.imputed_v3.male.tsv.bgz?dl=0 -O 2405.gwas.imputed_v3.male.tsv.bgz

# NEB women
wget https://www.dropbox.com/s/jdns4h91nip8qvl/2734.gwas.imputed_v3.female.tsv.bgz?dl=0 -O 2734.gwas.imputed_v3.female.tsv.bgz

# Religious attendance
wget https://www.dropbox.com/s/0ruqyyrqvs1osbp/6160_3.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 6160_3.gwas.imputed_v3.both_sexes.tsv.bgz

# Self-rated health
wget https://www.dropbox.com/s/aawh07hlhldbckc/2178.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 2178.gwas.imputed_v3.both_sexes.tsv.bgz

# SWB
wget https://www.dropbox.com/s/0u43m2w695pxk3x/4526.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 4526.gwas.imputed_v3.both_sexes.tsv.bgz

# Work satisfaction
wget https://www.dropbox.com/s/1o3l8jdcsugwa0e/4537.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 4537.gwas.imputed_v3.both_sexes.tsv.bgz
#------------------------------------#

#------------------------------------#
# Okbay et al - EA
wget http://ssgac.org/documents/EduYears_Main.txt.gz
#------------------------------------#

#------------------------------------#
# Okbay et al - SWB
wget http://ssgac.org/documents/SWB_Full.txt.gz
#------------------------------------#

#------------------------------------#
# Perry et al - Age first menses
# See above
#------------------------------------#

#------------------------------------#
# Rietveld et al - EA
wget http://ssgac.org/documents/MA_EA_1st_stage.txt.gz
#------------------------------------#

#------------------------------------#
# Savage et al - Intelligence
wget https://ctg.cncr.nl/documents/p1651/SavageJansen_IntMeta_sumstats.zip
#------------------------------------#

#------------------------------------#
# Van den Berg et al - Extraversion
# See above
#------------------------------------#

#------------------------------------#
# Yengo et al - BMI
wget https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz
#------------------------------------#

#------------------------------------#
# Yengo et al - Height
wget https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz
#------------------------------------#





#------------------------------------------------------------------------------------------------------------#
# -------------------------------------------- Reference data -----------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

#------------------------------------#
# 1000 Genomes
cd $PGI_Repo/original_data/ref_data/1000G_ph3
for chr in {1..22}; do
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz &
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi &
done
#------------------------------------#

#------------------------------------#
# HRC ref data for pre-imputation qc
cd $PGI_Repo/original_data/ref_data/HRC_imputation_qc
wget https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.10.zip
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
gunzip *.gz
unzip HRC-1000G-check-bim-v4.2.10.zip
rm HRC-1000G-check-bim-v4.2.10.zip
#------------------------------------#