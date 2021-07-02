#!/bin/bash

source $mainDir/code/paths
export R_LIBS=$Rlib:$R_LIBS

cd $mainDir/derived_data/2_Formatted/public
mkdir tmp

#------------------------------------------------------------------------------------------------------------#
# -------------------------------------------- PRE-PROCESSING -----------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

# Get openness from GPC1, remove rest 
mv $mainDir/original_data/public/GPC-1.BigFiveNEO.zip?dl=0 $mainDir/original_data/public/GPC-1.BigFiveNEO.zip
unzip $mainDir/original_data/public/GPC-1.BigFiveNEO.zip
rm GPC-1.NEO-CONSCIENTIOUSNESS.full.txt  GPC-1.NEO-NEUROTICISM.full.txt GPC-1.NEO-AGREEABLENESS.full.txt  GPC-1.NEO-EXTRAVERSION.full.txt

#-----------------------------------#
## Decompress and rename
gunzip -c $mainDir/original_data/public/ADHD-adhd_eur_jun2017.gz > tmp/ADHD-Demontis.txt &
gunzip -c $mainDir/original_data/public/AFBpooled-AgeFirstBirth_Pooled.txt.gz > tmp/AFB-Barban.txt &
gunzip -c $mainDir/original_data/public/AstEczRhi-SHARE-without23andMe.LDSCORE-GC.SE-META.v0.gz > tmp/ASTECZRHI-Ferreira.txt &
gunzip -c $mainDir/original_data/public/BMI-SNP_gwas_mc_merge_nogc.tbl.uniq.gz > tmp/BMI-Locke.txt &
gunzip -c $mainDir/original_data/public/CPD-tag.cpd.tbl.gz > tmp/CPD-Furberg.txt &
gunzip -c $mainDir/original_data/public/CigarettesPerDay.txt.gz?sequence=31 > tmp/CPD.txt.gz &
gunzip -c $mainDir/original_data/public/EVERSMOKE-tag.evrsmk.tbl.gz > tmp/EVERSMOKE-Furberg.txt &
gunzip -c $mainDir/original_data/public/SmokingInitiation.txt.gz?sequence=34 > tmp/EVERSMOKE-Liu.txt.gz
gunzip -c $mainDir/original_data/public/KP_DEPR_BETA_EAF.txt.gz > tmp/DEP-GERA.txt &
gunzip -c $mainDir/original_data/public/daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.gz > tmp/DEP-PGC.txt &
gunzip -c $mainDir/original_data/public/DrinksPerWeek.txt.gz?sequence=32 > tmp/DPW-Liu.txt.gz
gunzip -c $mainDir/original_data/public/Doherty-2018-NatureComms-overall-activity.csv.gz > tmp/ACTIVITY-Doherty.txt &
gunzip -c $mainDir/original_data/public/HEIGHT-GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz > tmp/HEIGHT-Wood.txt &
gunzip -c $mainDir/original_data/public/Intelligence-cogent.hrc.meta.chr.bp.rsid.assoc.full.gz > tmp/CP-COGENT.txt &
gunzip -c $mainDir/original_data/public/NEBmen-NumberChildrenEverBorn_Male.txt.gz > tmp/NEBmen-Barban.txt &
gunzip -c $mainDir/original_data/public/NEBwomen-NumberChildrenEverBorn_Female.txt.gz > tmp/NEBwomen-Barban.txt &
gunzip -c $mainDir/original_data/public/RISK_GWAS_MA_Nweighted_ID15_2017_08_06.tbl.gz > tmp/RISK-Linner.txt &
gunzip -c $mainDir/original_data/public/SWB_excl_PGSrepo_ldscGC.meta.gz > tmp/SWB-Okbay.txt &
wait
#-----------------------------------#

#-----------------------------------#
# Copy into tmp and rename
cp $mainDir/original_data/public/AGEFIRSTMENSES-Menarche_1KG_NatGen2017_WebsiteUpload.txt tmp/MENARCHE-Day.txt &
cp $mainDir/original_data/public/AGEFIRSTMENSES-Menarche_Nature2014_GWASMetaResults_17122014.txt tmp/MENARCHE-Perry.txt &
cp $mainDir/original_data/public/Kunkle_etal_Stage1_results.txt?file=1 tmp/ALZ-Kunkle.txt &
cp $mainDir/original_data/public/EA3_excl_UKB.meta tmp/EA-LeeExclUKB.txt &
cp $mainDir/original_data/public/EA3_PGSrepo.meta tmp/EA-LeeExclPGIrepo.txt &
cp $mainDir/original_data/public/GPC-2.EXTRAVERSION.full.txt tmp/EXTRA-vandenBerg.txt &
cp $mainDir/original_data/public/META_NEUROTICISM_ALLIwv_iwv_20150402_1.dat tmp/NEURO-deMoor.txt &
cp $mainDir/original_data/public/GPC-1.NEO-OPENNESS.full.txt tmp/OPEN-deMoor.txt &
cp $mainDir/original_data/public/CCI_discovery_MA_13_samples_07-18-2015_rs.txt.txt tmp/CANNABIS-Stringer.txt &
wait
#-----------------------------------#

#-----------------------------------#
# CANNABIS-Pasman
# Per-chromosome files, no header (get from readme)
# Header: SNP (rs-number, or CHR:BP position on build GRCh37 if no rs-number was available), Allele1 (effect allele), Allele2 (reference allele), MAF (minor allele frequency), Effect (beta regression coefficient from the meta-analysis), StdErr (standard error), P (p-value), Direction (direction of the effect per sample, order: ICC, [23andMe,] UKB), Chr (position on chromosome), Bp (position in basepairs), and N (sample size, depending on in how many samples this SNP was present).
awk 'BEGIN{OFS="\t";print "SNPID","EFFECT_ALLELE","OTHER_ALLELE","MAF","BETA","SE","P","CHR","BP","N"} \
{print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11}' $mainDir/original_data/public/Pasman_cannabis/cannabis_icc_ukb_chr*.txt > tmp/CANNABIS-Pasman.txt
#-----------------------------------#

#-----------------------------------#
# CP-UKB
# Tar file, per chromosome
tar -xzf $mainDir/original_data/public/Intelligence-IQ_BOLT_LMM_UKB_v2_BGEN.tar.gz
awk -F"\t" 'NR==1{print}NR>1 && $1!="SNP"{print}' OFS="\t" IQ_BOLT_LMM_UKB_v2_BGEN_Chr* > tmp/CP-UKB.txt
rm IQ_BOLT_LMM_UKB_v2_BGEN_Chr*
#-----------------------------------#

#-----------------------------------#
# ALZ-Kunkle, CANNABIS-Stringer, ASTECZRHI-Ferreira, CANNABIS-Pasman, CP-COGENT
# Convert to tab-delimited
sed -i 's/ /\t/g' tmp/ALZ-Kunkle.txt tmp/CANNABIS-Stringer.txt tmp/ASTECZRHI-Ferreira.txt tmp/CANNABIS-Pasman.txt tmp/CP-COGENT.txt
sed -i 's/,/\t/g' tmp/ACTIVITY-Doherty.txt
#-----------------------------------#




#------------------------------------------------------------------------------------------------------------#
# ---------------------- FORMAT EACH FILE, IGNORE SNPID FORMAT OR EAF ISSUES --------------------------------#
#------------------------------------------------------------------------------------------------------------#


# ACTIVITY
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=91105;print $1,$2,$3,toupper($5),toupper($6),$7,$11,$12,$16,$8,N,1,1,1,"A"}' tmp/ACTIVITY-Doherty.txt > ACTIVITY-Doherty.txt &

#---------------------------------------#

# ADHD: 20,183 cases, 35,191 controls
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=20183+35191;print $2,$1,$3,$4,$5,"NA",log($7),$8,$9,$6,N,1,1,1,"A"}' tmp/ADHD-Demontis.txt > tmp/tmp_ADHD-Demontis.txt &

#---------------------------------------#

# AFB
# all (N=251,151)-9,370=241781
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1 && $6>0 && $6<1{N=241781;print $1,$2,$3,$4,$5,$6,$7/sqrt(2*N*$6*(1-$6)),1/sqrt(2*N*$6*(1-$6)),$8,1,N,1,1,1,"A"}' tmp/AFB-Barban.txt > AFB-Barban.txt &

#---------------------------------------#

# AGE FIRST MENSES
# 329,345 = ReproGen consortium (N = 179,117) + 23andMe (N = 76,831) + UK Biobank (N = 73,397) studies 
# Npublic = 329,345 - 76,831
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=252514;split($1,a,":");print a[1]":"a[2],a[1],a[2],toupper($2),toupper($3),"NA",$4,"NA",$5,1,N,1,1,1,"A"}' tmp/MENARCHE-Day.txt | sed 's/chr//g' >  tmp/tmp_MENARCHE-Day.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=132989;print $1,"NA","NA",toupper($2),toupper($3),$4,$5,"NA",$6,1,N,1,1,1,"A"}' tmp/MENARCHE-Perry.txt > tmp/tmp_MENARCHE-Perry.txt &

#---------------------------------------#

# ALZHEIMERS
# IGAP-Kunkle - 21,982 cases + 41,944 controls 
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=21892+41944;print $3,$1,$2,$4,$5,"NA",$6,$7,$8,1,N,1,1,1,"A"}' tmp/ALZ-Kunkle.txt > tmp/tmp_ALZ-Kunkle.txt &


#---------------------------------------#

# Asthma/Eczema/Rhinitis								
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1 && toupper($4)==$11{N=$NF;print $10,$2,$3,toupper($4),toupper($5),$12,$6,$7,$8,1,N,1,1,1,"A"} \
NR>1 && toupper($5)==$11{N=$NF;print $10,$2,$3,toupper($4),toupper($5),1-$12,$6,$7,$8,1,N,1,1,1,"A"}' tmp/ASTECZRHI-Ferreira.txt > ASTECZRHI-Ferreira.txt &

#---------------------------------------#

# BMI
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=$NF;print $1,"NA","NA",$2,$3,$4,$5,$6,$7,1,N,1,1,1,"A"}' tmp/BMI-Locke.txt > tmp/tmp_BMI-Locke.txt &

#---------------------------------------#

# CANNABIS - Pasman
# ICC - N = 35,297, 42.8% cases --> 15107 cases, 20190 controls
# 23andMe - N=22,683, 43.2% cases
# UK-Biobank - N=126,785, 22.3% cases --> 28273 cases, 98512 controls
# ICC + UKB --> 43380 cases , 118702 controls
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=35297+126785;print $1,$8,$9,toupper($2),toupper($3),$4,$5,$6,$7,1,N,1,1,1,"A"}' tmp/CANNABIS-Pasman.txt > CANNABIS-Pasman.txt &

# CANNABIS - Stringer
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>2{N=35297;print $2,$3,$4,toupper($5),toupper($6),$7,$11,$12,$13,1,N,1,1,1,"A"}' tmp/CANNABIS-Stringer.txt > CANNABIS-Stringer.txt &

#---------------------------------------#

# CP
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=$13;print $2,$3,$4,$6,$7,$8,$10,$11,$12,$9,N,1,1,1,"A"}' tmp/CP-COGENT.txt > CP-COGENT.txt

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=222543;print $1,$2,$3,$5,$6,$7,$9,$10,$11,$8,N,1,1,1,"A"}' tmp/CP-UKB.txt > CP-UKB.txt

#---------------------------------------#

# CPD
# Furberg 
# chr=hg18, FRQ_A=FRQ_U (only had "freq1"), INFO=1 for all (unavailable), OR is not an OR !!! It's the linear regression beta for the continuous variables and logistic regression beta for the discrete variables. 
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=38181;print $2,"NA","NA",$4,$5,$6,$9,$10,$11,$8,N,1,1,1,"A"}' tmp/CPD-Furberg.txt > tmp/tmp_CPD-Furberg.txt &

# Liu
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{print $3,$1,$2,$5,$4,$6,$9,$10,$8,1,$11,1,1,1,"A"}' OFS="\t" tmp/CPD-GSCAN.txt > CPD-GSCAN.txt &

#---------------------------------------#

# DEP
# Ntot= 56,368, Ncases = 7,231, Ncontrols = 49,137
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
{N=56368; split($2, a, ":", seps)}{if ($2~/^[1-9]/ ||  $2~/^X:/) $2=a[1]":"a[2] ; else $2=a[1]} NR>1{print $2,$1,$3,$4,$5,$14,$13,$11,$12,$8,N,1,1,1,"A"}' tmp/DEP-GERA.txt > DEP-GERA.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=$17+$18;EAF=(($6*45396)+($7*97250))/(45396+97250); print $2,$1,$3,$4,$5,EAF,log($9),$10,$11,$8,N,1,1,1,"A"}' tmp/DEP-PGC.txt > DEP-PGC.txt &

#---------------------------------------#

# DPW
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{print $3,$1,$2,$5,$4,$6,$9,$10,$8,1,$11,1,1,1,"A"}' OFS="\t" tmp/DPW-GSCAN.txt > DPW-GSCAN.txt &

#---------------------------------------#

# EA
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=$10;print $1,$2,$3,$4,$5,$6,$16,$18,$12,1,$10,1,1,1,"A"}' tmp/EA-LeeExclPGSrepo.txt > EA-LeeExclPGSrepo.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=$10;print $1,$2,$3,$4,$5,$6,$16,$18,$12,1,$10,1,1,1,"A"}' tmp/EA-LeeExclUKB.txt > EA-LeeExclUKB.txt &

#---------------------------------------#

# EXTRA
# N=63,030
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=63030;print $1,$2,$3,toupper($4),toupper($5),"NA",$6,$7,$8,1,N,1,1,1,"A"}' tmp/EXTRA-vandenBerg.txt > tmp/tmp_EXTRA-vandenBerg.txt &

#---------------------------------------#

# EVERSMOKE
# Furberg 
# N=74,035, 41969 cases 32066 controls‬ 
# chr=hg18, FRQ_A=FRQ_U (only had "freq1"), INFO=1 for all (unavailable), OR is not an OR !!! It's the linear regression beta for the continuous variables and logistic regression beta for the discrete variables. 
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=41969+32066;print $2,"NA","NA",$4,$5,$6,$9,$10,$11,$8,N,1,1,1,"A"}' tmp/EVERSMOKE-Furberg.txt > tmp/tmp_EVERSMOKE-Furberg.txt &

# Liu
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{print $3,$1,$2,$5,$4,$6,$9,$10,$8,1,$11,1,1,1,"A"}' OFS="\t" tmp/EVERSMOKE-GSCAN.txt > EVERSMOKE-GSCAN.txt &

#---------------------------------------#

# HEIGHT 
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=$NF;print $1,"NA","NA",$2,$3,$4,$5,$6,$7,1,N,1,1,1,"A"}' tmp/HEIGHT-Wood.txt > tmp/tmp_HEIGHT-Wood.txt &

#---------------------------------------#

# NEB
# NEB women (N=225,230)-13,635= 211,595
# NEB men (N=103,909)-10,974= 92,935
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1 && $6>0 && $6<1{N=92935;print $1,$2,$3,$4,$5,$6,$7/sqrt(2*N*$6*(1-$6)),1/sqrt(2*N*$6*(1-$6)),$8,1,N,1,1,1,"A"}' tmp/NEBmen-Barban.txt > NEBmen-Barban.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1 && $6>0 && $6<1{N=211595;print $1,$2,$3,$4,$5,$6,$7/sqrt(2*N*$6*(1-$6)),1/sqrt(2*N*$6*(1-$6)),$8,1,N,1,1,1,"A"}' tmp/NEBwomen-Barban.txt > NEBwomen-Barban.txt &

#---------------------------------------#

# NEURO
# Setting EAF to NA because Marleen de Moor had said it's not reliable (mix of MAF and EAF)
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{print $3,$1,$2,$4,$5,"NA",$10,$11,$12,1,$14,1,1,1,"A"}' OFS="\t" tmp/NEURO-deMoor.txt > tmp/tmp_NEURO-deMoor.txt &

#---------------------------------------#

# OPEN
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	{N=17375;print $1,"NA","NA",toupper($4),toupper($5),"NA",$6,$7,$8,$9,N,1,1,1,"A"}' OFS="\t" tmp/OPEN-deMoor.txt > tmp/tmp_OPEN-deMoor.txt &

#---------------------------------------#

# RISK 
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1 && $6>0 && $6<1{N=$10;split($3,a,":");print $2,a[1],a[2],toupper($4),toupper($5),$6,$11/sqrt(2*N*$6*(1-$6)),1/sqrt(2*N*$6*(1-$6)),$12,1,N,1,1,1,"A"}' tmp/RISK-Linner.txt > RISK-Linner.txt &

#---------------------------------------#

# SWB
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1 && $6>0 && $6<1{N=$10;print $1,$2,$3,toupper($4),toupper($5),$6,$11/sqrt(2*N*$6*(1-$6)),1/sqrt(2*N*$6*(1-$6)),$12,1,N,1,1,1,"A"}' tmp/SWB-Okbay.txt > SWB-Okbay.txt &

wait


#------------------------------------------------------------------------------------------------------------#
# ----------------------------------- HANDLE SNPID FORMAT, EAF, SE ISSUES -----------------------------------#
#------------------------------------------------------------------------------------------------------------#

## EAF ISSUES

# Files with rsID and missing EAF
for file in ALZ-Kunkle.txt ADHD-Demontis.txt EXTRA-vandenBerg.txt OPEN-deMoor.txt ; do
	awk -F"\t" 'NR==FNR{ref[$1]=$3;raf[$1]=$5;next} \
	FNR==1{print;next} \
	FNR>1 && $1 in ref && ref[$1]==$4{$6=raf[$1];print;next}
	FNR>1 && $1 in ref && ref[$1]==$5{$6=1-raf[$1];print}' OFS="\t" $HRC_raf tmp/tmp_$file > $file &
done
wait
mv OPEN-GPC1.txt tmp/tmp_OPEN-GPC1.txt

## Files with ChrPosID and missing EAF (also adding rsID)
for file in NEURO-deMoor.txt MENARCHE-Day.txt; do
	awk -F"\t" 'NR==FNR{rs[$2]=$1;ref[$2]=$3;raf[$2]=$5;next} \
	FNR==1{print;next} \
	FNR>1 && $1 in ref && ref[$1]==$4{$6=raf[$1];$1=rs[$1];print;next}
	FNR>1 && $1 in ref && ref[$1]==$5{$6=1-raf[$1];$1=rs[$1];print}' OFS="\t" $HRC_raf tmp/tmp_$file > $file &
done 
wait
mv MENARCHE-Day.txt tmp/tmp_MENARCHE-Day.txt

#---------------------------------------#

## Missing CHR/BP
for file in BMI-Locke.txt MENARCHE-Perry.txt HEIGHT-Wood.txt EVERSMOKE-Furberg.txt CPD-Furberg.txt OPEN-deMoor.txt; do
	awk -F"\t" 'NR==FNR{ChrPosID[$1]=$2;next} \
		FNR==1{print;next} FNR>1 && $1 in ChrPosID {split(ChrPosID[$1],a,":");$2=a[1];$3=a[2];print}' OFS="\t" $HRC_raf tmp/tmp_$file > $file &
done
mv MENARCHE-Perry.txt tmp/tmp_MENARCHE-Perry.txt

#---------------------------------------#

## Missing SE
for file in MENARCHE-Day.txt MENARCHE-Perry.txt; do
	Rscript $mainDir/code/2_Formatting/2.3.1_Add_SE.R $mainDir/derived_data/2_Formatted/public/tmp/tmp_$file $mainDir/derived_data/2_Formatted/public/$file &
done
wait

rm tmp/*
gzip *
