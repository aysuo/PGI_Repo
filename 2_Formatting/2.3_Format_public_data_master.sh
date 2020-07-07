#!/bin/bash

dirIn="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/public"
dirOut="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public"
dirCode="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code/2_Formatting"
export R_LIBS=/homes/nber/aokbay/R/x86_64-redhat-linux-gnu-library/3.1/:$R_LIBS
cd $dirOut


#-----------------------------------#
#-------- PRE-PROCESSING -----------#
#-----------------------------------#
mkdir tmp

#-----------------------------------#
# Alzheimer's: IGAP + UKB
# What needs to be done: There's an extra first column containing line numbers, needs to be removed because there's no header for that column
awk -F"\t" 'NR==1{print}NR>1{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' OFS="\t" $dirIn/Alzheimer-AD.meta.assoc_cleaned > tmp/ALZ-UKB_IGAP.txt
#-----------------------------------#

#-----------------------------------#
# Lifetime cannabis use (Pasman)
# What needs to be done: Per-chromosome files, no header (get from readme)
# Header: SNP (rs-number, or CHR:BP position on build GRCh37 if no rs-number was available), Allele1 (effect allele), Allele2 (reference allele), MAF (minor allele frequency), Effect (beta regression coefficient from the meta-analysis), StdErr (standard error), P (p-value), Direction (direction of the effect per sample, order: ICC, [23andMe,] UKB), Chr (position on chromosome), Bp (position in basepairs), and N (sample size, depending on in how many samples this SNP was present).
awk 'BEGIN{OFS="\t";print "SNPID","EFFECT_ALLELE","OTHER_ALLELE","MAF","BETA","SE","P","CHR","BP","N"} \
{print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11}' $dirIn/Pasman_cannabis/cannabis_icc_ukb_chr*.txt > tmp/Cannabis-Pasman.txt
#-----------------------------------#

#-----------------------------------#
# Intelligence
# What needs to be done: Tar file, per chromosome
tar -xzf $dirIn/Intelligence-IQ_BOLT_LMM_UKB_v2_BGEN.tar.gz
awk -F"\t" 'NR==1{print}NR>1 && $1!="SNP"{print}' OFS="\t" IQ_BOLT_LMM_UKB_v2_BGEN_Chr* > tmp/IQ-UKB.txt
rm IQ_BOLT_LMM_UKB_v2_BGEN_Chr*
#-----------------------------------#

#-----------------------------------#
# ADHD, AFB (men, women, pooled), Asthma/Eczema/Rhinitis, BMI, Cannabis (disorder), CPD, Eversmoke, DEP, Height,
# NEB (men, women, pooled), Risk, SWB
# What needs to be done: Decompress, rename

gunzip -c $dirIn/ADHD-adhd_eur_jun2017.gz > tmp/ADHD-Demontis.txt &
gunzip -c $dirIn/AFBmen-AgeFirstBirth_Male.txt.gz > tmp/AFBmen-Barban.txt &
gunzip -c $dirIn/AFBpooled-AgeFirstBirth_Pooled.txt.gz > tmp/AFBpooled-Barban.txt &
gunzip -c $dirIn/AFBwomen-AgeFirstBirth_Female.txt.gz > tmp/AFBwomen-Barban.txt &
gunzip -c $dirIn/AstEczRhi-SHARE-without23andMe.LDSCORE-GC.SE-META.v0.gz > tmp/AstEczRhi-Ferreira.txt &
gunzip -c $dirIn/BMI-SNP_gwas_mc_merge_nogc.tbl.uniq.gz > tmp/BMI-Locke.txt &
gunzip -c $dirIn/CannabisDisorder-CUD_GWAS_iPSYCH_June2019.gz > tmp/Cannabis-Demontis.txt &
gunzip -c $dirIn/CPD-tag.cpd.tbl.gz > tmp/CPD-TAG.txt &
gunzip -c $dirIn/EVERSMOKE-tag.evrsmk.tbl.gz > tmp/Eversmoke-TAG.txt &
#gunzip -c $dirIn/DEP_UKB_GERA.meta.gz > tmp/DEP-UKBGERA.txt &
gunzip -c $dirIn/KP_DEPR_BETA_EAF.txt.gz > tmp/DEP-GERA.txt &
gunzip -c $dirIn/daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.gz > tmp/DEP-PGC.txt &
gunzip -c $dirIn/Doherty-2018-NatureComms-overall-activity-sexBMI.csv.gz > tmp/ActivityAdj-Doherty.txt &
gunzip -c $dirIn/Doherty-2018-NatureComms-overall-activity.csv.gz > tmp/Activity-Doherty.txt &
gunzip -c $dirIn/HEIGHT-GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz > tmp/Height-Wood.txt &
gunzip -c $dirIn/Intelligence-cogent.hrc.meta.chr.bp.rsid.assoc.full.gz > tmp/IQ-COGENT.txt &
gunzip -c $dirIn/NEBmen-NumberChildrenEverBorn_Male.txt.gz > tmp/NEBmen-Barban.txt &
gunzip -c $dirIn/NEBpooled-NumberChildrenEverBorn_Pooled.txt.gz > tmp/NEBpooled-Barban.txt &
gunzip -c $dirIn/NEBwomen-NumberChildrenEverBorn_Female.txt.gz > tmp/NEBwomen-Barban.txt &
gunzip -c $dirIn/RISK_GWAS_MA_Nweighted_ID15_2017_08_06.tbl.gz > tmp/Risk-Linner.txt &
gunzip -c $dirIn/SWB_excl_PGSrepo_ldscGC.meta.gz > tmp/SWB-Okbay.txt &

wait
#-----------------------------------#


#-----------------------------------#
# AgeFirstMenses-Day, AgeFirstMenses-Perry, Alzheimer-IGAP, EA3 excl UKB, EA3 excl PGSrepo, Extraversion
# What needs to be done: Rename, copy into tmp

cp $dirIn/AGEFIRSTMENSES-Menarche_1KG_NatGen2017_WebsiteUpload.txt tmp/AgeFirstMenses-Day.txt &
cp $dirIn/AGEFIRSTMENSES-Menarche_Nature2014_GWASMetaResults_17122014.txt tmp/AgeFirstMenses-Perry.txt &
#cp $dirIn/IGAP_stage_1.txt tmp/ALZ-IGAP.txt &
cp $dirIn/Kunkle_etal_Stage1_results.txt?file=1 tmp/ALZ-Kunkle.txt &
cp $dirIn/EA3_excl_UKB.meta tmp/EA-LeeExclUKB.txt &
cp $dirIn/EA3_PGSrepo.meta tmp/EA-LeeExclPGSrepo.txt &
cp $dirIn/GPC-2.EXTRAVERSION.full.txt tmp/Extraversion-GPC2.txt &
cp $dirIn/META_NEUROTICISM_ALLIwv_iwv_20150402_1.dat tmp/Neuro-GPC2.txt &
wait
#-----------------------------------#

#-----------------------------------#
# ALZ-Kunkle, Cannabis-Stringer, AstEczRhi-Ferreira, Cannabis-Demontis, Cannabis-Pasman, IQ-COGENT
# What needs to be done: Convert to tab-delimited
sed -i 's/ /\t/g' tmp/ALZ-Kunkle.txt $dirIn/CCI_discovery_MA_13_samples_07-18-2015_rs.txt.txt tmp/AstEczRhi-Ferreira.txt tmp/Cannabis-Demontis.txt tmp/Cannabis-Pasman.txt tmp/IQ-COGENT.txt
sed -i 's/,/\t/g' tmp/Activity-Doherty.txt tmp/ActivityAdj-Doherty.txt
#-----------------------------------#

#--------------------------------------------------------------------------------------#
########################################################################################
#--------------------------------------------------------------------------------------#

## FORMAT EACH FILE, IGNORE SNPID FORMAT OR EAF ISSUES 

# ADHD: 20,183 cases, 35,191 controls
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=20183+35191;print $2,$1,$3,$4,$5,"NA",log($7),$8,$9,$6,N,1,1,1,"A"}' tmp/ADHD-Demontis.txt > tmp/tmp_ADHD-Demontis.txt &

#---------------------------------------#

# AFB
# women (N=189,656)-9,370=180286
# all (N=251,151)-9,370=241781
# men (N=48,408)-0

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1 && $6>0 && $6<1{N=180286;print $1,$2,$3,$4,$5,$6,$7/sqrt(2*N*$6*(1-$6)),1/sqrt(2*N*$6*(1-$6)),$8,1,N,1,1,1,"A"}'  tmp/AFBwomen-Barban.txt > AFBwomen-Barban.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1 && $6>0 && $6<1{N=241781;print $1,$2,$3,$4,$5,$6,$7/sqrt(2*N*$6*(1-$6)),1/sqrt(2*N*$6*(1-$6)),$8,1,N,1,1,1,"A"}' tmp/AFBpooled-Barban.txt > AFBpooled-Barban.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1 && $6>0 && $6<1{N=48408;print $1,$2,$3,$4,$5,$6,$7/sqrt(2*N*$6*(1-$6)),1/sqrt(2*N*$6*(1-$6)),$8,1,N,1,1,1,"A"}' tmp/AFBmen-Barban.txt > AFBmen-Barban.txt &


#---------------------------------------#

# AGE FIRST MENSES
# 329,345 = ReproGen consortium (N = 179,117) + 23andMe (N = 76,831) + UK Biobank (N = 73,397) studies 
# Npublic = 329,345 - 76,831

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=252514;split($1,a,":");print a[1]":"a[2],a[1],a[2],toupper($2),toupper($3),"NA",$4,"NA",$5,1,N,1,1,1,"A"}' tmp/AgeFirstMenses-Day.txt | sed 's/chr//g' >  tmp/tmp_AgeFirstMenses-Day.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=132989;print $1,"NA","NA",toupper($2),toupper($3),$4,$5,"NA",$6,1,N,1,1,1,"A"}' tmp/AgeFirstMenses-Perry.txt > tmp/tmp_AgeFirstMenses-Perry.txt &

#---------------------------------------#

# ALZHEIMERS
# Riccardo assumed N=450k
# N_UKB= 314,278 (27,696 maternal cases, 14,338 paternal cases), IGAP (n = 74,046 with 25,580 cases)
# In total, 67,614‬ cases + 320710 controls --> N=388,324‬, N_eff=223,364
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=314278+74046;print $2,$1,$3,$4,$5,$11,$8,$9,$10,1,N,1,1,1,"A"}' tmp/ALZ-UKB_IGAP.txt > ALZ-UKB_IGAP.txt &

# IGAP-Lambert
# (Build 37, Assembly Hg19)
#awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
#NR>1{N=25580+48466;print $3,$1,$2,$4,$5,"NA",$6,$7,$8,1,N,1,1,1,"A"}' tmp/ALZ-IGAP.txt > tmp/tmp_ALZ-IGAP.txt &

# IGAP-Kunkle - 21,982 cases + 41,944 controls 
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=21892+41944;print $3,$1,$2,$4,$5,"NA",$6,$7,$8,1,N,1,1,1,"A"}' tmp/ALZ-Kunkle.txt > tmp/tmp_ALZ-Kunkle.txt &


#---------------------------------------#

# Asthma/Eczema/Rhinitis								
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1 && toupper($4)==$11{N=$NF;print $10,$2,$3,toupper($4),toupper($5),$12,$6,$7,$8,1,N,1,1,1,"A"} \
NR>1 && toupper($5)==$11{N=$NF;print $10,$2,$3,toupper($4),toupper($5),1-$12,$6,$7,$8,1,N,1,1,1,"A"}' tmp/AstEczRhi-Ferreira.txt > AstEczRhi-Ferreira.txt &

#---------------------------------------#

# BMI
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=$NF;print $1,"NA","NA",$2,$3,$4,$5,$6,$7,1,N,1,1,1,"A"}' tmp/BMI-Locke.txt > BMI-Locke.txt &

#---------------------------------------#

# Cannabis use disorder
# base pair location (hg19)
# N=2,387 CUD cases and 48,985 controls
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=2387+48985;print $2,$1,$3,$4,$5,"NA",log($7),$8,$9,$6,N,1,1,1,"A"}' tmp/Cannabis-Demontis.txt > tmp/tmp_Cannabis-Demontis.txt &

# Cannabis use (lifetime) Pasman 
# ICC - N = 35,297, 42.8% cases --> 15107 cases, 20190 controls
# 23andMe - N=22,683, 43.2% cases
# UK-Biobank - N=126,785, 22.3% cases --> 28273 cases, 98512 controls
# ICC + UKB --> 43380 cases , 118702 controls
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=35297+126785;print $1,$8,$9,toupper($2),toupper($3),$4,$5,$6,$7,1,N,1,1,1,"A"}' tmp/Cannabis-Pasman.txt > Cannabis-Pasman.txt &

# Cannabis use (lifetime) Stringer
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>2{N=35297;print $2,$3,$4,toupper($5),toupper($6),$7,$11,$12,$13,1,N,1,1,1,"A"}' $dirIn/CCI_discovery_MA_13_samples_07-18-2015_rs.txt.txt > Cannabis-Stringer.txt &


#---------------------------------------#

# TAG
# chr=hg18, FRQ_A=FRQ_U (only had "freq1"), INFO=1 for all (unavailable), OR is not an OR !!! Is the linear regression beta for the 
# continuous variables and logistic regression beta for the discrete variables. 
# Sample sizes across the three consortia were n = 143,023 for smoking initiation, n = 73,853 for CPD

# CPD - N=38,181
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=38181;print $2,"NA","NA",$4,$5,$6,$9,$10,$11,$8,N,1,1,1,"A"}' tmp/CPD-TAG.txt > tmp/tmp_CPD-TAG.txt &

# EVERSMOKE - N=74,035, 41969 cases 32066 controls‬ 
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=41969+32066;print $2,"NA","NA",$4,$5,$6,$9,$10,$11,$8,N,1,1,1,"A"}' tmp/Eversmoke-TAG.txt > tmp/tmp_Eversmoke-TAG.txt &

###################

# GSCAN
# See formatting script for public scores
cp /disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public_scores/Eversmoke-GSCAN.txt.gz . &
cp /disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public_scores/DPW-GSCAN.txt.gz . &
cp /disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public_scores/CPD-GSCAN.txt.gz . &

#---------------------------------------#

# DEP
#awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
#NR>1 && $4>0 && $4<1 {N=$8;split($1,a,":");print $1,a[1],a[2],toupper($2),toupper($3),$4,$9/sqrt(2*N*$4*(1-$4)),1/sqrt(2*N*$4*(1-$4)),$10,1,N,1,1,1,"A"}' tmp/DEP-UKBGERA.txt > DEP-UKBGERA.txt &

# Ntot= 56,368, Ncases = 7,231, Ncontrols = 49,137
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
{N=56368; split($2, a, ":", seps)}{if ($2~/^[1-9]/ ||  $2~/^X:/) $2=a[1]":"a[2] ; else $2=a[1]} NR>1{print $2,$1,$3,$4,$5,$14,$13,$11,$12,$8,N,1,1,1,"A"}' tmp/DEP-GERA.txt > DEP-GERA.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=$17+$18;EAF=(($6*45396)+($7*97250))/(45396+97250); print $2,$1,$3,$4,$5,EAF,log($9),$10,$11,$8,N,1,1,1,"A"}' tmp/DEP-PGC.txt > DEP-PGC.txt &

#---------------------------------------#

# EA
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=$10;print $1,$2,$3,$4,$5,$6,$16,$18,$12,1,$10,1,1,1,"A"}' tmp/EA-LeeExclPGSrepo.txt > EA-LeeExclPGSrepo.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=$10;print $1,$2,$3,$4,$5,$6,$16,$18,$12,1,$10,1,1,1,"A"}' tmp/EA-LeeExclUKB.txt > EA-LeeExclUKB.txt &

#---------------------------------------#

# GPC 
# Extraversion, N=63,030
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=63030;print $1,$2,$3,toupper($4),toupper($5),"NA",$6,$7,$8,1,N,1,1,1,"A"}' tmp/Extraversion-GPC2.txt > tmp/tmp_Extraversion-GPC2.txt &

#---------------------------------------#

# HEIGHT 
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=$NF;print $1,"NA","NA",$2,$3,$4,$5,$6,$7,1,N,1,1,1,"A"}' tmp/Height-Wood.txt > Height-Wood.txt &

#---------------------------------------#

# INTELLIGENCE
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=$13;print $2,$3,$4,$6,$7,$8,$10,$11,$12,$9,N,1,1,1,"A"}' tmp/IQ-COGENT.txt > IQ-COGENT.txt

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=222543;print $1,$2,$3,$5,$6,$7,$9,$10,$11,$8,N,1,1,1,"A"}' tmp/IQ-UKB.txt > IQ-UKB.txt

#---------------------------------------#

# NEB
# NEB women (N=225,230)-13,635= 211,595
# NEB men (N=103,909)-10,974= 92,935
# NEB all (343,072)-13,635= 329,437
    
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1 && $6>0 && $6<1{N=92935;print $1,$2,$3,$4,$5,$6,$7/sqrt(2*N*$6*(1-$6)),1/sqrt(2*N*$6*(1-$6)),$8,1,N,1,1,1,"A"}' tmp/NEBmen-Barban.txt > NEBmen-Barban.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1 && $6>0 && $6<1{N=329437;print $1,$2,$3,$4,$5,$6,$7/sqrt(2*N*$6*(1-$6)),1/sqrt(2*N*$6*(1-$6)),$8,1,N,1,1,1,"A"}' tmp/NEBpooled-Barban.txt > NEBpooled-Barban.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1 && $6>0 && $6<1{N=211595;print $1,$2,$3,$4,$5,$6,$7/sqrt(2*N*$6*(1-$6)),1/sqrt(2*N*$6*(1-$6)),$8,1,N,1,1,1,"A"}' tmp/NEBwomen-Barban.txt > NEBwomen-Barban.txt &

#---------------------------------------#

# Neuroticism
# Setting EAF to NA because Marleen de Moor had said it's not reliable (mix of MAF and EAF)
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{print $3,$1,$2,$4,$5,"NA",$10,$11,$12,1,$14,1,1,1,"A"}' OFS="\t" tmp/Neuro-GPC2.txt > tmp/tmp_Neuro-GPC2.txt &

#---------------------------------------#

# Physical activity
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=91105;print $1,$2,$3,toupper($5),toupper($6),$7,$11,$12,$16,$8,N,1,1,1,"A"}' tmp/ActivityAdj-Doherty.txt > ActivityAdj-Doherty.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=91105;print $1,$2,$3,toupper($5),toupper($6),$7,$11,$12,$16,$8,N,1,1,1,"A"}' tmp/Activity-Doherty.txt > Activity-Doherty.txt &

#---------------------------------------#

# RISK 
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1 && $6>0 && $6<1{N=$10;split($3,a,":");print $2,a[1],a[2],toupper($4),toupper($5),$6,$11/sqrt(2*N*$6*(1-$6)),1/sqrt(2*N*$6*(1-$6)),$12,1,N,1,1,1,"A"}' tmp/Risk-Linner.txt > Risk-Linner.txt &

#---------------------------------------#

# SWB
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1 && $6>0 && $6<1{N=$10;print $1,$2,$3,toupper($4),toupper($5),$6,$11/sqrt(2*N*$6*(1-$6)),1/sqrt(2*N*$6*(1-$6)),$12,1,N,1,1,1,"A"}' tmp/SWB-Okbay.txt > SWB-Okbay.txt &

wait

#--------------------------------------------------------------------------------------#
########################################################################################
#--------------------------------------------------------------------------------------#

## HANDLE EAF ISSUES
## Files with rsID and missing EAF

for file in ALZ-Kunkle.txt ADHD-Demontis.txt ALZ-IGAP.txt Cannabis-Demontis.txt Extraversion-GPC2.txt; do
	awk -F"\t" 'NR==FNR{ref[$1]=$3;raf[$1]=$5;next} \
	FNR==1{print;next} \
	FNR>1 && $1 in ref && ref[$1]==$4{$6=raf[$1];print;next}
	FNR>1 && $1 in ref && ref[$1]==$5{$6=1-raf[$1];print}' OFS="\t" /disk/genetics/ukb/aokbay/reffiles/HRC/HRC_r1-1.GRCh37.wgs.mac5.maf001.rsID_cptid_alleles_raf tmp/tmp_$file > $file &
done
wait

## Files with ChrPosID and missing EAF (also adding rsID)
for file in Neuro-GPC2.txt AgeFirstMenses-Day.txt; do
	awk -F"\t" 'NR==FNR{rs[$2]=$1;ref[$2]=$3;raf[$2]=$5;next} \
	FNR==1{print;next} \
	FNR>1 && $1 in ref && ref[$1]==$4{$6=raf[$1];$1=rs[$1];print;next}
	FNR>1 && $1 in ref && ref[$1]==$5{$6=1-raf[$1];$1=rs[$1];print}' OFS="\t" /disk/genetics/ukb/aokbay/reffiles/HRC/HRC_r1-1.GRCh37.wgs.mac5.maf001.rsID_cptid_alleles_raf tmp/tmp_$file > $file &
done 
wait
mv AgeFirstMenses-Day.txt tmp/tmp_AgeFirstMenses-Day.txt

## Files with no rsID
#for file in DEP-GERA.txt; do
#	awk -F"\t" 'NR==FNR{rs[$2]=$1;ref[$2]=$3;alt[$2]=$4;next} \
#	FNR==1 || !($1 in ref) {print;next} \
#	FNR>1 && $1 in ref && ((ref[$1]==$4 && alt[$1]==$5) || (ref[$1]==$5 && alt[$1]==$4)) {$1=rs[$1];print;next}' OFS="\t" /disk/genetics/ukb/aokbay/reffiles/HRC/HRC_r1-1.GRCh37.wgs.mac5.maf001.rsID_cptid_alleles_raf tmp/tmp_$file > $file &
#done 
#wait

## Files with no Chr BP
for file in Eversmoke-TAG.txt CPD-TAG.txt; do
	awk -F"\t" 'NR==FNR{ChrPosID[$1]=$2;next} \
		FNR==1{print;next} FNR>1 && $1 in ChrPosID {split(ChrPosID[$1],a,":");$2=a[1];$3=a[2];print}' OFS="\t" /disk/genetics/ukb/aokbay/reffiles/HRC/HRC_r1-1.GRCh37.wgs.mac5.maf001.rsID_cptid_alleles_raf tmp/tmp_$file > $file &
done


## Files with no SE
for file in AgeFirstMenses-Day.txt AgeFirstMenses-Perry.txt; do
	Rscript $dirCode/1.1_Add_SE.R $dirOut/tmp/tmp_$file $dirOut/$file &
done
wait

gzip *
