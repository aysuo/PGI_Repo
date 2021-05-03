#!/bin/bash

source dirs0
cd $dir0_PublicData

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
# Depression - GERA
cp /disk/genetics/ukb/aokbay/SWB/Depression/Kaiser/KP_DEPR_BETA_EAF.txt.gz .
#------------------------------------#

#------------------------------------#
# Intelligence
# COGENT (from EA3)
cp /disk/genetics3/EA3/data/COHORT_LEVEL/NEW/ORIGINAL/cogent.hrc.meta.chr.bp.rsid.assoc.full.gz .
# UKB
cp /disk/genetics3/EA3/data/COHORT_LEVEL/NEW/ORIGINAL/OUTPUT_BOLT_LMM_UKB_v2_IQ_2017_07_31_STRICTER/IQ_BOLT_LMM_UKB_v2_BGEN.tar.gz .
#------------------------------------#

#------------------------------------#
# Neuroticism - de Moor
cp /disk/genetics/ukb/aokbay/SWB/Neuro/GPC/META_NEUROTICISM_ALLIwv_iwv_20150402_1.dat .
#------------------------------------#

#------------------------------------#
# Risk - Linner excluding STR
cp /disk/genetics4/PGS/PGS_Repo/data/other_sumstats/RISK_GWAS_MA_Nweighted_ID15_2017_08_06.tbl.gz .
#------------------------------------#