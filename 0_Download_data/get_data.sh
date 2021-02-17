# Download the public data in this folder
# Aysu Okbay 21/02/2020

cd /disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/original_data/public

#------------------------------------#
# Wood et al height GWAS
# https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files
#------------------------------------#
wget https://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz

#------------------------------------#
# Locke BMI data
# https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files
#------------------------------------#
wget https://portals.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz

#---------------------
# Van den Berg et al. - Extraversion
#---------------------
wget https://www.dropbox.com/s/bk2jn41vrfl3zna/GPC-2.EXTRAVERSION.zip?dl=0

#------------------------------------#
# Demontis et al. - ADHD
# Download from here: https://www.med.unc.edu/pgc/download-results/adhd/
#------------------------------------#

#------------------------------------#
# Ferreira et al.- Asthma/Eczema/Rhinitis
#------------------------------------#
wget https://genepi.qimr.edu.au/staff/manuelF/gwas_results/SHARE-without23andMe.LDSCORE-GC.SE-META.v0.gz

#------------------------------------#
# Day et al - Age at menarche
#------------------------------------#
wget https://www.reprogen.org/Menarche_1KG_NatGen2017_WebsiteUpload.zip

#------------------------------------#
# Perry et al - Age at menarche (excludes UKB)
wget https://www.reprogen.org/Menarche_Nature2014_GWASMetaResults_17122014.zip
#------------------------------------#

#------------------------------------#
# TAG - Cigarettes per day
# Download from here: https://www.med.unc.edu/pgc/download-results/tag/
#------------------------------------#

#------------------------------------#
# TAG - Eversmoke
# Download from here: https://www.med.unc.edu/pgc/download-results/tag/
#------------------------------------#

#------------------------------------#
# Pasman et al - Cannabis (lifetime use)
# Download from here: https://www.ru.nl/bsi/research/group-pages/substance-use-addiction-food-saf/vm-saf/genetics/international-cannabis-consortium-icc/
# wget https://www.ru.nl/publish/pages/898181/cannabis_readme_1.docx
#------------------------------------#

#------------------------------------#
# Stringer et al - Cannabis (lifetime use, exlcudes UKB)
# Download from here: https://www.ru.nl/bsi/research/group-pages/substance-use-addiction-food-saf/vm-saf/genetics/international-cannabis-consortium-icc/
# https://www.ru.nl/publish/pages/898181/paper2_readme.docx
#------------------------------------#

#------------------------------------#
# Demontis et al - Cannabis use disorder
# Download from here: https://ipsych.dk/en/research/downloads/data-download-agreement-ipsych-secondary-phenotypes-cannabis-2019/
wget https://ipsych.dk/fileadmin/ipsych.dk/Downloads/README_CUD_GWAS.pdf
#------------------------------------#

#------------------------------------#
# AUDIT
# Sanchez-Roige - GWAS summary statistics for the UKB GWAS of AUDIT scores will be available on request.
#------------------------------------#

#------------------------------------#
# ALZHEIMER'S IGAP

# OLD
# Download from here: http://web.pasteur-lille.fr/en/recherche/u744/igap/igap_download.php

# Kunkle 
wget https://www.niagads.org/system/tdf/public_docs/Kunkle_etal_2019_IGAP_summary_statistics_README_0.docx?file=1&type=field_collection_item&id=120&force=
wget https://www.niagads.org/system/tdf/public_docs/Kunkle_etal_Stage1_results.txt?file=1&type=field_collection_item&id=121&force=
wget https://www.niagads.org/system/tdf/public_docs/Kunkle_etal_Stage2_results.txt?file=1&type=field_collection_item&id=122&force=

#------------------------------------#

#------------------------------------#
# Barban et al fertility results 
#------------------------------------#
# wget http://sociogenome.com/material/GWASresults/AgeFirstBirth_Pooled.txt.gz
wget http://sociogenome.com/material/GWASresults/NumberChildrenEverBorn_Pooled.txt.gz

wget http://sociogenome.com/material/GWASresults/AgeFirstBirth_Female.txt.gz
wget http://sociogenome.com/material/GWASresults/AgeFirstBirth_Male.txt.gz

wget http://sociogenome.com/material/GWASresults/NumberChildrenEverBorn_Female.txt.gz
wget http://sociogenome.com/material/GWASresults/NumberChildrenEverBorn_Male.txt.gz

wget http://ssgac.org/documents/readme_reproductivebehavior.txt


#------------------------------------#
# Physical activity
# Download from here: https://ora.ox.ac.uk/objects/uuid:ff479f44-bf35-48b9-9e67-e690a2937b22
#------------------------------------#

#------------------------------------------------------------------------------------------------------------#
# ------------------------------------------ In-house results -----------------------------------------------#

#------------------------------------#
# Alzheimer's 
# Visscher meta analysis of IGAP + proxycase (use this - Riccardo Marioni used approx N of 450K for this)
cp /var/genetics/ukb/edkong/E3_ProxyPhenotype/INPUT/AD.meta.assoc_cleaned .
#------------------------------------#

#------------------------------------#
# EA
# Aysu ran EA3 meta-analysis excluding HRS, TEDS, WLS, AH, MCTFR, ELSA, STR, EGCUT
cp /disk/genetics4/PGS/Aysu/META/EA/EA3_PGSrepo.meta .
#------------------------------------#

#------------------------------------#
# EA excluding UKB
# Aysu ran EA3 meta-analysis excluding UKB
cp /disk/genetics/PGS/Aysu/META/EA/EA3_excl_UKB/EA3_excl_UKB.meta .
#------------------------------------#

#------------------------------------#
# SWB
# Aysu ran SWB meta-analysis excluding 23andMe, EGCUT, ELSA, HRS, MCTFR, STR, TEDS, UKB
cp /disk/genetics/PGS/Aysu/META/SWB/SWB_excl_PGSrepo_ldscGC.meta.gz .
#------------------------------------#

#------------------------------------#
# DEP
# Aysu ran meta-analysis of UKB and GERA (not done at the within-trait MTAG stage because too many SNPs get lost when taking overlapping )
# cp /disk/genetics/PGS/Aysu/META/DEP/DEP_UKB_GERA.meta.gz .

cp /disk/genetics/ukb/aokbay/SWB/Depression/Kaiser/KP_DEPR_BETA_EAF.txt.gz .
#------------------------------------#

#------------------------------------#
# Intelligence
# COGENT results from EA3
cp /disk/genetics3/EA3/data/COHORT_LEVEL/NEW/ORIGINAL/cogent.hrc.meta.chr.bp.rsid.assoc.full.gz .
# UKB IQ results
cp /disk/genetics3/EA3/data/COHORT_LEVEL/NEW/ORIGINAL/OUTPUT_BOLT_LMM_UKB_v2_IQ_2017_07_31_STRICTER/IQ_BOLT_LMM_UKB_v2_BGEN.tar.gz .
#------------------------------------#

#------------------------------------#
# Neuroticism - de Moor
cp /disk/genetics/ukb/aokbay/SWB/Neuro/GPC/META_NEUROTICISM_ALLIwv_iwv_20150402_1.dat .
#------------------------------------#

#------------------------------------#
# Risk
# Richard ran meta-analysis excluding STR
cp /disk/genetics4/PGS/PGS_Repo/data/other_sumstats/RISK_GWAS_MA_Nweighted_ID15_2017_08_06.tbl.gz .
# We will need results excluding UKB!!!
#------------------------------------#

#------------------------------------#
# Ever-smoker
# Got from the risk paper (Ed had it?)
# cp /var/genetics/23andmess/RISK_SUMSTATS/RISK_META/SMOKE_EVER_MA_Nweighted_2017_08_28.tbl.gz .
# DO TAG + UKB instead!!
#------------------------------------#




#------------------------------------------------------------------------------------------------------------#
# ------------------------------------------ For public scores ----------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

cd /disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/original_data/public_scores

# Neale Lab pheno description files
wget https://www.dropbox.com/s/d4mlq9ly93yhjyt/phenotypes.both_sexes.tsv.bgz -O phenotypes.both_sexes.tsv.gz
wget https://www.dropbox.com/s/r7idtoulxfpyjss/phenotypes.female.tsv.bgz -O phenotypes.female.tsv.gz
wget https://www.dropbox.com/s/mywldevz4nsla2r/phenotypes.male.tsv.bgz -O phenotypes.male.tsv.gz
# Neale Lab variant annotation file
wget https://www.dropbox.com/s/puxks683vb0omeg/variants.tsv.bgz?dl=0 -O variants.tsv.bgz

# BMI
wget https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz &

# HEIGHT
wget https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz &

# Extraversion
# See above

# Morning person
wget https://personal.broadinstitute.org/mvon/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.gz &
wget https://s3.amazonaws.com/broad-portal-resources/sleep/chronotype_raw_README.txt &

# Risk
wget https://www.dropbox.com/s/il1d7vabk5283dm/RISK_GWAS_MA_UKB%2Breplication.txt?dl=0 &

# Neuroticism
wget https://ctg.cncr.nl/documents/p1651/sumstats_neuro_sum_ctg_format.txt.gz
wget https://ctg.cncr.nl/documents/p1651/readme_neuro_items_ctg

# Religiosity
wget https://www.dropbox.com/s/0ruqyyrqvs1osbp/6160_3.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 6160_3.gwas.imputed_v3.both_sexes.tsv.bgz

# Openness
wget https://www.dropbox.com/s/df3bifmk238ks2y/GPC-1.BigFiveNEO.zip?dl=0

# EA
wget https://www.dropbox.com/s/ho58e9jmytmpaf8/GWAS_EA_excl23andMe.txt?dl=0
wget http://ssgac.org/documents/EduYears_Main.txt.gz
wget http://ssgac.org/documents/MA_EA_1st_stage.txt.gz

# Friend satisfaction
wget https://www.dropbox.com/s/274tmv6f3q3z4q1/4570.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 4570.gwas.imputed_v3.both_sexes.tsv.bgz

# Family satisfaction
wget https://www.dropbox.com/s/tsn2m8uzgwnkkr6/4559.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 4559.gwas.imputed_v3.both_sexes.tsv.bgz

# Work satisfaction
wget https://www.dropbox.com/s/1o3l8jdcsugwa0e/4537.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 4537.gwas.imputed_v3.both_sexes.tsv.bgz

# Financial satisfaction
wget https://www.dropbox.com/s/pezhe9cgcfehv9a/4581.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 4581.gwas.imputed_v3.both_sexes.tsv.bgz

# Loneliness
wget https://www.dropbox.com/s/nf4jl3mdppu1ng8/2020.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 2020.gwas.imputed_v3.both_sexes.tsv.bgz

# COPD
wget https://www.dropbox.com/s/hxx4huusrq7nwqq/22130.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 22130.gwas.imputed_v3.both_sexes.tsv.bgz


#---------------------------------------------------------#
# Intelligence
# Savage
wget https://ctg.cncr.nl/documents/p1651/SavageJansen_IntMeta_sumstats.zip
# Lee
wget https://www.dropbox.com/s/ibjoh0g5e3sdd8t/GWAS_CP_all.txt?dl=0
#---------------------------------------------------------#

# SWB
wget http://ssgac.org/documents/SWB_Full.txt.gz
wget https://www.dropbox.com/s/0u43m2w695pxk3x/4526.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 4526.gwas.imputed_v3.both_sexes.tsv.bgz

# Depression
wget https://datashare.is.ed.ac.uk/bitstream/handle/10283/3203/PGC_UKB_depression_genome-wide.txt?sequence=3&isAllowed=y
wget https://datashare.is.ed.ac.uk/bitstream/handle/10283/3203/ReadMe.txt?sequence=4&isAllowed=y

# ADHD
# See above

#---------------------------------------------------------#
# Self-rated health
# Harris
wget https://datashare.is.ed.ac.uk/bitstream/handle/10283/3342/Harris2016_UKB_self_rated_health_summary_results_10112016.txt?sequence=2&isAllowed=y
wget https://datashare.is.ed.ac.uk/bitstream/handle/10283/3342/Harris2016_README.txt?sequence=1&isAllowed=y
# Neale
wget https://www.dropbox.com/s/aawh07hlhldbckc/2178.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 2178.gwas.imputed_v3.both_sexes.tsv.bgz
#---------------------------------------------------------#

#---------------------------------------------------------#
# Physical activity
# Download from here: https://ora.ox.ac.uk/objects/uuid:ff479f44-bf35-48b9-9e67-e690a2937b22
#---------------------------------------------------------#

# Migraine
wget https://www.dropbox.com/s/31xvr94v67qbt7e/G43.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O G43.gwas.imputed_v3.both_sexes.tsv.bgz

# Ever-smoker
wget https://conservancy.umn.edu/bitstream/handle/11299/201564/SmokingInitiation.txt.gz?sequence=34&isAllowed=y
wget https://www.dropbox.com/s/o7wgwhnhjgt3eyn/EVER_SMOKER_GWAS_MA_UKB%2BTAG.txt?dl=0

# ---------------------------------------------------------- #
# Asthma
# Choose the one with largest effective sample size 

# N_controls=359201, N_cases=1993 -->  N_eff=7928.012016
# wget https://www.dropbox.com/s/9fo51y4uaho1i7y/J10_ASTHMA_MAIN.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O J10_ASTHMA_MAIN.gwas.imputed_v3.both_sexes.tsv.bgz
# wget https://www.dropbox.com/s/zc6mrd6ztbnnzf7/J10_ASTHMA.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O J10_ASTHMA.gwas.imputed_v3.both_sexes.tsv.bgz

# N_controls=359501, N_cases=1693 --> N_eff=6740.258066
# wget https://www.dropbox.com/s/ft8wkalnzhtwjrg/J45.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O J45.gwas.imputed_v3.both_sexes.tsv.bgz

# N_controls=319207, N_cases=41934 --> N_eff=148259.282
wget https://www.dropbox.com/s/kp9bollwekaco0s/20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz

# N_controls=80070, N_cases=11717 --> N_eff=40885.10094
# wget https://www.dropbox.com/s/26tsq3xfzi0bqcu/22127.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 22127.gwas.imputed_v3.both_sexes.tsv.bgz

# N_cases=318894, N_controls=41633 --> N_eff=147301.1886
# wget https://www.dropbox.com/s/4c03qds07w2gpii/6152_8.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 6152_8.gwas.imputed_v3.both_sexes.tsv.bgz
# ---------------------------------------------------------- #


# ---------------------------------------------------------- #
# Hayfever
# Choose the one with largest effective sample size 
# Do both Hayfever/Eczema (6152) and largest hayfever only (20002)

# N_controls=277120, N_cases=83407 --> 256444.0149
wget https://www.dropbox.com/s/strqg2hhhgl09rj/6152_9.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 6152_9.gwas.imputed_v3.both_sexes.tsv.bgz

# N_controls=340474, N_cases=20667 --> 77937.16203
wget https://www.dropbox.com/s/9zrzjw9fp8seo4m/20002_1387.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 20002_1387.gwas.imputed_v3.both_sexes.tsv.bgz &

# N_controls=70883 , N_cases=20904 --> 64572.9017
# wget https://www.dropbox.com/s/d9abufm8pajk79d/22126.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 22126.gwas.imputed_v3.both_sexes.tsv.bgz
# ---------------------------------------------------------- #

# Asthma/Eczema/Rhinitis
# See above

# Cannabis
# Ask Joel where we got Pasman incl UKB
wget https://www.dropbox.com/s/m6frhnk184rd7uz/20453.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 20453.gwas.imputed_v3.both_sexes.tsv.bgz &

# Nearsightedness
wget https://www.dropbox.com/s/lawm9hw7oui862a/6147_1.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 6147_1.gwas.imputed_v3.both_sexes.tsv.bgz &
wget https://www.dropbox.com/s/gmh9g65q5iw93ee/6147_1.gwas.imputed_v3.both_sexes.v2.tsv.bgz?dl=0 -O 6147_1.gwas.imputed_v3.both_sexes.v2.tsv.bgz 

# Cigarettes per day 
wget https://conservancy.umn.edu/bitstream/handle/11299/201564/CigarettesPerDay.txt.gz?sequence=31&isAllowed=y

# Drinks per week
wget https://conservancy.umn.edu/bitstream/handle/11299/201564/DrinksPerWeek.txt.gz?sequence=32&isAllowed=y
wget https://www.dropbox.com/s/7hjxdhlxlwa482n/DRINKS_PER_WEEK_GWAS.txt?dl=0

# NEB women
# See above for Barban
# UKB:
wget https://www.dropbox.com/s/jdns4h91nip8qvl/2734.gwas.imputed_v3.female.tsv.bgz?dl=0 -O 2734.gwas.imputed_v3.female.tsv.bgz

# NEBmen
wget https://www.dropbox.com/s/uhj063tukv01d8m/2405.gwas.imputed_v3.male.tsv.bgz?dl=0 -O 2405.gwas.imputed_v3.male.tsv.bgz


# AFB 
# See above

# Menarche
# See above


#------------------------------------------------------------------------------------------------------------#
# -------------------------------------------- Reference data -----------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

# Get 1000 Genomes data
cd /disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/ref_data/1000G_ph3
for chr in {1..22}; do
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz &
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi &
done


# Get HRC imputation qc ref data
cd /disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/ref_data/HRC_imputation_qc
wget https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.10.zip
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
gunzip *.gz
unzip HRC-1000G-check-bim-v4.2.10.zip
rm HRC-1000G-check-bim-v4.2.10.zip