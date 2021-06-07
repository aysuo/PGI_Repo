#!/bin/bash

source paths2
export R_LIBS=$Rlib:$R_LIBS

cd $mainDir/derived_data/2_Formatted/public_scores
mkdir -p tmp

#------------------------------------------------------------------------------------------------------------#
# -------------------------------------------- UNZIP & RENAME -----------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

# Rename *bgz to *gz
rename bgz gz $mainDir/original_data/public_scores/*bgz

# Get rid of weird extensions
mv $mainDir/original_data/public_scores/GWAS_CP_all.txt?dl=0 $mainDir/original_data/public_scores/GWAS_CP_all.txt
mv $mainDir/original_data/public_scores/DRINKS_PER_WEEK_GWAS.txt\?dl\=0 $mainDir/original_data/public_scores/DRINKS_PER_WEEK_GWAS.txt
mv $mainDir/original_data/public_scores/EVER_SMOKER_GWAS_MA_UKB+TAG.txt?dl=0 $mainDir/original_data/public_scores/EVER_SMOKER_GWAS_MA_UKB+TAG.txt

# Rename Neale Lab GWAS
cp $mainDir/original_data/public_scores/6160_3.gwas.imputed_v3.both_sexes.tsv tmp/RELIGATT-Neale.txt
cp $mainDir/original_data/public_scores/20002_1111.gwas.imputed_v3.both_sexes.tsv tmp/ASTHMA-Neale.txt
cp $mainDir/original_data/public_scores/6152_9.gwas.imputed_v3.both_sexes.tsv tmp/HAYFEVER-Neale.txt
cp $mainDir/original_data/public_scores/6147_1.gwas.imputed_v3.both_sexes.v2.tsv tmp/NEARSIGHTED-Neale.txt
cp $mainDir/original_data/public_scores/G43.gwas.imputed_v3.both_sexes.tsv tmp/MIGRAINE-Neale.txt
cp $mainDir/original_data/public_scores/20453.gwas.imputed_v3.both_sexes.tsv tmp/CANNABIS-Neale.txt
cp $mainDir/original_data/public_scores/2734.gwas.imputed_v3.female.tsv tmp/NEBwomen-Neale.txt
cp $mainDir/original_data/public_scores/2178.gwas.imputed_v3.both_sexes.tsv tmp/SELFHEALTH-Neale.txt
cp $mainDir/original_data/public_scores/4526.gwas.imputed_v3.both_sexes.tsv tmp/SWB-Neale.txt
cp $mainDir/original_data/public_scores/4559.gwas.imputed_v3.both_sexes.tsv tmp/FAMSAT-Neale.txt
cp $mainDir/original_data/public_scores/4570.gwas.imputed_v3.both_sexes.tsv tmp/FRIENDSAT-Neale.txt
cp $mainDir/original_data/public_scores/2405.gwas.imputed_v3.male.tsv tmp/NEBmen-Neale.txt
cp $mainDir/original_data/public_scores/4537.gwas.imputed_v3.both_sexes.tsv tmp/WORKSAT-Neale.txt
cp $mainDir/original_data/public_scores/4581.gwas.imputed_v3.both_sexes.tsv tmp/FINSAT-Neale.txt
cp $mainDir/original_data/public_scores/2020.gwas.imputed_v3.both_sexes.tsv tmp/LONELY-Neale.txt
cp $mainDir/original_data/public_scores/22130.gwas.imputed_v3.both_sexes.tsv tmp/COPD-Neale.txt

# Unzip everything
gunzip $mainDir/original_data/public_scores/*.gz

for file in *.zip; do
	unzip $file
done

# Clean-up
rm *.zip



#------------------------------------------------------------------------------------------------------------#
# ---------------------- FORMAT EACH FILE, IGNORE SNPID FORMAT OR EAF ISSUES --------------------------------#
#------------------------------------------------------------------------------------------------------------#

# Neale Lab - ASTHMA CANNABIS COPD FAMSAT FINSAT FRIENDSAT HAYFEVER LONELY MIGRAINE NEARSIGHTED NEBmen NEBwomen RELIGATT SELFHEALTH SWB WORKSAT   
for file in NEBmen LONELY COPD NEBwomen RELIGATT ASTHMA HAYFEVER NEARSIGHTED MIGRAINE; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR==FNR{chr[$1]=$2; bp[$1]=$3; ea[$1]=$5; oa[$1]=$4; rs[$1]=$6; info[$1]=$10; cr[$1]=$11; eaf[$1]=$13; hwe[$1]=$16;next} \
		(FNR>1){print rs[$1],chr[$1],bp[$1],ea[$1],oa[$1],eaf[$1],$9,$10,$12,info[$1],$6,1,cr[$1],hwe[$1],"A" }' OFS="\t" $mainDir/original_data/public_scores/variants.tsv tmp/${file}-Neale.txt > ${file}-Neale.txt &
done

for file in WORKSAT FINSAT FRIENDSAT FAMSAT SWB SELFHEALTH CANNABIS; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR==FNR{chr[$1]=$2; bp[$1]=$3; ea[$1]=$5; oa[$1]=$4; rs[$1]=$6; info[$1]=$10; cr[$1]=$11; eaf[$1]=$13; hwe[$1]=$16;next} \
		(FNR>1){print rs[$1],chr[$1],bp[$1],ea[$1],oa[$1],eaf[$1],$8,$9,$11,info[$1],$5,1,cr[$1],hwe[$1],"A" }' OFS="\t" $mainDir/original_data/public_scores/variants.tsv tmp/${file}-Neale.txt > ${file}-Neale.txt &
done

#-------------------------------------#

# ADHD
cp $mainDir/derived_data/2_Formatted/public/ADHD-Demontis.txt.gz . &

#-------------------------------------#

# AFB 
cp $mainDir/derived_data/2_Formatted/public/AFB-Barban.txt.gz . &

#-------------------------------------#

# ASTECZRHI
cp $mainDir/derived_data/2_Formatted/public/ASTECZRHI-Ferreira.txt.gz . &

#-------------------------------------#

# BMI 
# Yengo
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{print $3,$1,$2,$4,$5,$6,$7,$8,$9,"1",$10,1,1,1,"A"}' OFS="\t" $mainDir/original_data/public_scores/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt > BMI-Yengo.txt &

# Locke
cp $mainDir/derived_data/2_Formatted/public/BMI-Locke.txt.gz .

#-------------------------------------#

# CP
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{print $1,$3,$4,toupper($5),toupper($6),$7,$9,$10,$11,$13,$12,1,1,1,"A"}' OFS="\t" $mainDir/original_data/public_scores/SavageJansen_2018_intelligence_metaanalysis.txt > CP-Savage.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=257841; print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,N,1,1,1,"A"}' OFS="\t" $mainDir/original_data/public_scores/GWAS_CP_all.txt > CP-Lee.txt &

#-------------------------------------#

# CPD
cp $mainDir/derived_data/2_Formatted/public/CPD-Liu.txt.gz .

#-------------------------------------#

# EA
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=766345; print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,N,1,1,1,"A"}' OFS="\t" $mainDir/original_data/public_scores/GWAS_EA_excl23andMe.txt?dl=0 > EA-Lee.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=328917; print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,N,1,1,1,"A"}' OFS="\t" $mainDir/original_data/public_scores/EduYears_Main.txt > EA-Okbay.txt &

awk 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=126559; print $1,"NA","NA",toupper($2),toupper($3),"NA",$4,"NA",$5,1,N,1,1,1,"A"}' OFS="\t" $mainDir/original_data/public_scores/MA_EA_1st_stage.txt > tmp/EA-Rietveld.txt &

#-------------------------------------#

# DEP
# N=500,199 (170,756 cases and 329,443 controls)
awk 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=500199; print $1,"NA","NA",toupper($2),toupper($3),$4,$5,$6,$7,1,N,1,1,1,"A"}' OFS="\t" $mainDir/original_data/public_scores/PGC_UKB_depression_genome-wide.txt?sequence=3 > tmp/DEP-Howard.txt &

#-------------------------------------#

# DPW
# Linner
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=414343; print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,N,1,1,1,"A"}' OFS="\t" $mainDir/original_data/public_scores/DRINKS_PER_WEEK_GWAS.txt > DPW-Linner.txt &

# Liu
cp $mainDir/derived_data/2_Formatted/public/DPW-Liu.txt.gz .

#-------------------------------------#

# EXTRA
cp $mainDir/derived_data/2_Formatted/public/EXTRA-vandenBerg.txt.gz .

#-------------------------------------#

# EVERSMOKE
# Linner
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=518633; print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,N,1,1,1,"A"}' OFS="\t" $mainDir/original_data/public_scores/EVER_SMOKER_GWAS_MA_UKB+TAG.txt > EVERSMOKE-Linner.txt &

# Liu
cp $mainDir/derived_data/2_Formatted/public/EVERSMOKE-Liu.txt.gz .

# Furberg
cp $mainDir/derived_data/2_Formatted/public/EVERSMOKE-Furberg.txt.gz .

#-------------------------------------#

# HEIGHT
# Yengo
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{print $3,$1,$2,$4,$5,$6,$7,$8,$9,"1",$10,1,1,1,"A"}' OFS="\t" $mainDir/original_data/public_scores/Meta-analysis_Wood_et_al+UKBiobank_2018.txt > HEIGHT-Yengo.txt &

# Wood
cp $mainDir/derived_data/2_Formatted/public/HEIGHT-Wood.txt.gz .

#-------------------------------------#

# MENARCHE
cp $mainDir/derived_data/2_Formatted/public/MENARCHE-Day.txt.gz .
cp $mainDir/derived_data/2_Formatted/public/MENARCHE-Perry.txt.gz .

#-------------------------------------#

# MORNING
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=449734;print $1,$2,$3,$4,$5,$6,$8,$9,$10,$7,N,1,1,$11,"A"}' OFS="\t" $mainDir/original_data/public_scores/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt > MORNING-Jones.txt &

#-------------------------------------#

# NEBwomen
cp $mainDir/derived_data/2_Formatted/public/NEBwomen-Barban.txt.gz . &
# Neale version above with other Neale Lab GWAS

#-------------------------------------#
# NEURO
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{print $1,$3,$4,$5,$6,$7,$10,$11,$12,$14,$13,1,1,1,"A"}' OFS="\t" $mainDir/original_data/public_scores/sumstats_neuro_sum_ctg_format.txt > NEURO-Nagel.txt &

#-------------------------------------#

# RISK
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=466571;print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,N,1,1,1,"A"}' OFS="\t" $mainDir/original_data/public_scores/RISK_GWAS_MA_UKB+replication.txt?dl=0 > RISK-Linner.txt &

#-------------------------------------#

# SWB
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=204978; print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,N,1,1,1,"A"}' OFS="\t" $mainDir/original_data/public_scores/SWB_Full.txt > SWB-Okbay.txt &

#-------------------------------------#


#------------------------------------------------------------------------------------------------------------#
# ----------------------------------- HANDLE SNPID FORMAT, EAF, SE ISSUES -----------------------------------#
#------------------------------------------------------------------------------------------------------------#

## EAF ISSUES
# Files with rsID and missing EAF
for file in EA-Rietveld.txt; do
	awk -F"\t" 'NR==FNR{ref[$1]=$3;raf[$1]=$5;next} \
	FNR==1{print;next} \
	FNR>1 && $1 in ref && ref[$1]==$4{$6=raf[$1];print;next}
	FNR>1 && $1 in ref && ref[$1]==$5{$6=1-raf[$1];print}' OFS="\t" $p2_HRCref tmp/$file > $file &
done
mv EA-Rietveld.txt tmp/EA-Rietveld.txt

#---------------------------------------#

## Missing SE
for file in EA-Rietveld.txt; do
	Rscript $mainDir/code/2_Formatting/2.3.1_Add_SE.R $mainDir/derived_data/2_Formatted/public_scores/tmp/$file $mainDir/derived_data/2_Formatted/public_scores/$file &
done
mv EA-Rietveld.txt tmp/EA-Rietveld.txt

#---------------------------------------#

## Missing CHR/BP
for file in EA-Rietveld.txt DEP-Howard.txt; do
	awk -F"\t" 'NR==FNR{ChrPosID[$1]=$2;next} \
		FNR==1{print;next} FNR>1 && $1 in ChrPosID {split(ChrPosID[$1],a,":");$2=a[1];$3=a[2];print}' OFS="\t" $p2_HRCref tmp/$file > $file &
done
wait

rm tmp/*
gzip *.txt



