dirIn="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/original_data/public_scores"
dirOut="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public_scores"
dirCode="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code/2_Formatting"
export R_LIBS=/homes/nber/aokbay/R/x86_64-redhat-linux-gnu-library/3.1/:$R_LIBS
cd $dirOut
mkdir -p tmp

rename bgz gz $dirIn/*bgz
mv $dirIn/SmokingInitiation.txt.gz?sequence=34 $dirIn/SmokingInitiation.txt.gz
mv $dirIn/DrinksPerWeek.txt.gz?sequence=32 $dirIn/DrinksPerWeek.txt.gz
mv $dirIn/CigarettesPerDay.txt.gz?sequence=31 $dirIn/CigarettesPerDay.txt.gz
mv $dirIn/GPC-1.BigFiveNEO.zip?dl=0 $dirIn/GPC-1.BigFiveNEO.zip
mv $dirIn/GWAS_CP_all.txt?dl=0 $dirIn/GWAS_CP_all.txt
mv $dirIn/DRINKS_PER_WEEK_GWAS.txt\?dl\=0 $dirIn/DRINKS_PER_WEEK_GWAS.txt
mv $dirIn/EVER_SMOKER_GWAS_MA_UKB+TAG.txt?dl=0 $dirIn/EVER_SMOKER_GWAS_MA_UKB+TAG.txt

gunzip $dirIn/*.gz

for file in *.zip; do
	unzip $file
done

rm *.zip

rm GPC-1.NEO-CONSCIENTIOUSNESS.full.txt  GPC-1.NEO-NEUROTICISM.full.txt GPC-1.NEO-AGREEABLENESS.full.txt  GPC-1.NEO-EXTRAVERSION.full.txt

##############################################################################################################################


#-------------------------------------#
# BMI and HEIGHT
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{print $3,$1,$2,$4,$5,$6,$7,$8,$9,"1",$10,1,1,1,"A"}' OFS="\t" $dirIn/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt > BMI-Yengo.txt &

cp /disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public/BMI-Locke.txt.gz BMI-Locke.txt.gz
cp /disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public/Height-Wood.txt.gz HEIGHT-Wood.txt.gz

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{print $3,$1,$2,$4,$5,$6,$7,$8,$9,"1",$10,1,1,1,"A"}' OFS="\t" $dirIn/Meta-analysis_Wood_et_al+UKBiobank_2018.txt > HEIGHT-Yengo.txt &


#-------------------------------------#

# Extraversion
cp /disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public/EXTRA-GPC2.txt.gz .

#-------------------------------------#

# Morning person
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=449734;print $1,$2,$3,$4,$5,$6,$8,$9,$10,$7,N,1,1,$11,"A"}' OFS="\t" $dirIn/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt > MORNING-Jones.txt &

#-------------------------------------#

# Risk
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{N=466571;print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,N,1,1,1,"A"}' OFS="\t" $dirIn/RISK_GWAS_MA_UKB+replication.txt?dl=0 > RISK-Linner.txt &

#-------------------------------------#

# Neuroticism
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
NR>1{print $1,$3,$4,$5,$6,$7,$10,$11,$12,$14,$13,1,1,1,"A"}' OFS="\t" $dirIn/sumstats_neuro_sum_ctg_format.txt > NEURO-Nagel.txt &

#-------------------------------------#

# Religiosity, Asthma, Hayfever, Nearsightedness, Cannabis
cp $dirIn/6160_3.gwas.imputed_v3.both_sexes.tsv tmp/RELIGATT-Neale.txt
cp $dirIn/20002_1111.gwas.imputed_v3.both_sexes.tsv tmp/ASTHMA-Neale.txt
cp $dirIn/6152_9.gwas.imputed_v3.both_sexes.tsv tmp/HAYFEVER-Neale.txt
cp $dirIn/20002_1387.gwas.imputed_v3.both_sexes.tsv tmp/Rhinitis-Neale.txt
cp $dirIn/6147_1.gwas.imputed_v3.both_sexes.v2.tsv tmp/NEARSIGHTED-Neale.txt
cp $dirIn/G43.gwas.imputed_v3.both_sexes.tsv tmp/MIGRAINE-Neale.txt
cp $dirIn/20453.gwas.imputed_v3.both_sexes.tsv tmp/CANNABIS-Neale.txt
cp $dirIn/2734.gwas.imputed_v3.female.tsv tmp/NEBwomen-Neale.txt
cp $dirIn/2178.gwas.imputed_v3.both_sexes.tsv tmp/SELFHEALTH-Neale.txt
cp $dirIn/4526.gwas.imputed_v3.both_sexes.tsv tmp/SWB-Neale.txt
cp $dirIn/4559.gwas.imputed_v3.both_sexes.tsv tmp/FAMSAT-Neale.txt
cp $dirIn/4570.gwas.imputed_v3.both_sexes.tsv tmp/FRIENDSAT-Neale.txt
cp $dirIn/2405.gwas.imputed_v3.male.tsv tmp/NEBmen-Neale.txt
cp $dirIn/4537.gwas.imputed_v3.both_sexes.tsv tmp/WORKSAT-Neale.txt
cp $dirIn/4581.gwas.imputed_v3.both_sexes.tsv tmp/FINSAT-Neale.txt
cp $dirIn/2020.gwas.imputed_v3.both_sexes.tsv tmp/LONELY-Neale.txt
cp $dirIn/22130.gwas.imputed_v3.both_sexes.tsv tmp/COPD-Neale.txt


for file in NEBmen LONELY COPD NEBwomen RELIGATT ASTHMA HAYFEVER NEARSIGHTED MIGRAINE; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR==FNR{chr[$1]=$2; bp[$1]=$3; ea[$1]=$5; oa[$1]=$4; rs[$1]=$6; info[$1]=$10; cr[$1]=$11; eaf[$1]=$13; hwe[$1]=$16;next} \
		(FNR>1){print rs[$1],chr[$1],bp[$1],ea[$1],oa[$1],eaf[$1],$9,$10,$12,info[$1],$6,1,cr[$1],hwe[$1],"A" }' OFS="\t" $dirIn/variants.tsv tmp/${file}-Neale.txt > ${file}-Neale.txt &
done

for file in WORKSAT FINSAT FRIENDSAT FAMSAT SWB SELFHEALTH CANNABIS; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR==FNR{chr[$1]=$2; bp[$1]=$3; ea[$1]=$5; oa[$1]=$4; rs[$1]=$6; info[$1]=$10; cr[$1]=$11; eaf[$1]=$13; hwe[$1]=$16;next} \
		(FNR>1){print rs[$1],chr[$1],bp[$1],ea[$1],oa[$1],eaf[$1],$8,$9,$11,info[$1],$5,1,cr[$1],hwe[$1],"A" }' OFS="\t" $dirIn/variants.tsv tmp/${file}-Neale.txt > ${file}-Neale.txt &
done

#-------------------------------------#

# Openness
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	{N=17375;print $1,"NA","NA",toupper($4),toupper($5),"NA",$6,$7,$8,$9,N,1,1,1,"A"}' OFS="\t" $dirIn/GPC-1.NEO-OPENNESS.full.txt > tmp/OPEN-deMoor.txt &


#-------------------------------------#

# EA
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=766345; print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,N,1,1,1,"A"}' OFS="\t" $dirIn/GWAS_EA_excl23andMe.txt?dl=0 > EA-Lee.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=328917; print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,N,1,1,1,"A"}' OFS="\t" $dirIn/EduYears_Main.txt > EA-Okbay.txt &

awk 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=126559; print $1,"NA","NA",toupper($2),toupper($3),"NA",$4,"NA",$5,1,N,1,1,1,"A"}' OFS="\t" $dirIn/MA_EA_1st_stage.txt > tmp/EA-Rietveld.txt &

#-------------------------------------#

# Intelligence
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{print $1,$3,$4,toupper($5),toupper($6),$7,$9,$10,$11,$13,$12,1,1,1,"A"}' OFS="\t" $dirIn/SavageJansen_2018_intelligence_metaanalysis.txt > CP-Savage.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=257841; print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,N,1,1,1,"A"}' OFS="\t" $dirIn/GWAS_CP_all.txt > CP-Lee.txt &



#-------------------------------------#

# SWB
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=204978; print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,N,1,1,1,"A"}' OFS="\t" $dirIn/SWB_Full.txt > SWB-Okbay.txt &


#-------------------------------------#

# Depression
# N=500,199 (170,756 cases and 329,443 controls)
awk 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=500199; print $1,"NA","NA",toupper($2),toupper($3),$4,$5,$6,$7,1,N,1,1,1,"A"}' OFS="\t" $dirIn/PGC_UKB_depression_genome-wide.txt?sequence=3 > tmp/DEP-Howard.txt &


#-------------------------------------#
 
# ADHD
cp  /disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public/ADHD-Demontis.txt.gz .

#-------------------------------------#

# Self-rated health
awk 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=111749; print $3,$1,$2,$4,$5,"NA",$7,"NA",$6,1,N,1,1,1,"A"}' OFS="\t" $dirIn/Harris2016_UKB_self_rated_health_summary_results_10112016.txt?sequence=2 > tmp/SELFHEALTH-Harris.txt
# Neale version above with other Neale Lab GWAS

#-------------------------------------#

# Ever-smoker, Cigarettes per day, Drinks per week
cp $dirIn/CigarettesPerDay.txt tmp/CPD-GSCAN.txt
cp $dirIn/SmokingInitiation.txt tmp/EVERSMOKE-GSCAN.txt
cp $dirIn/DrinksPerWeek.txt tmp/DPW-GSCAN.txt

for file in CPD-GSCAN.txt EVERSMOKE-GSCAN.txt DPW-GSCAN.txt; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR>1{print $3,$1,$2,$5,$4,$6,$9,$10,$8,1,$11,1,1,1,"A"}' OFS="\t" tmp/$file > $file &
done

# Drinks per week - Linner
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=414343; print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,N,1,1,1,"A"}' OFS="\t" $dirIn/DRINKS_PER_WEEK_GWAS.txt > DPW-Linner.txt &

# Ever-smoker - Linner
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=518633; print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,N,1,1,1,"A"}' OFS="\t" $dirIn/EVER_SMOKER_GWAS_MA_UKB+TAG.txt > EVERSMOKE-Linner.txt &
	
#-------------------------------------#

# Asthma/Eczema/Rhinitis
cp /disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public/ASTECZRHI-Ferreira.txt.gz . &

#-------------------------------------#

# NEB women, NEB pooled
cp /disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public/NEBwomen-Barban.txt.gz . &
cp /disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public/NEBpooled-Barban.txt.gz . &
cp /disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public/NEBmen-Barban.txt.gz . &
# Neale version above with other Neale Lab GWAS

#-------------------------------------#
# AFB 
cp /disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public/AFBpooled-Barban.txt.gz AFB-Barban.txt.gz &

#-------------------------------------#

# Menarche
gunzip -c /disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public/AgeFirstMenses-Day.txt.gz > tmp/MENARCHE-Day.txt
cp /disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/public/AgeFirstMenses-Perry.txt.gz MENARCHE-Perry.txt.gz



#--------------------------------------------------------------------------------------#
########################################################################################
#--------------------------------------------------------------------------------------#

## HANDLE EAF ISSUES
## Files with rsID and missing EAF
for file in EA-Rietveld.txt OPEN-deMoor.txt SELFHEALTH-Harris.txt; do
	awk -F"\t" 'NR==FNR{ref[$1]=$3;raf[$1]=$5;next} \
	FNR==1{print;next} \
	FNR>1 && $1 in ref && ref[$1]==$4{$6=raf[$1];print;next}
	FNR>1 && $1 in ref && ref[$1]==$5{$6=1-raf[$1];print}' OFS="\t" /disk/genetics/ukb/aokbay/reffiles/HRC/HRC_r1-1.GRCh37.wgs.mac5.maf001.rsID_cptid_alleles_raf tmp/$file > $file &
done
wait
mv SELFHEALTH-Harris.txt tmp/SELFHEALTH-Harris.txt
mv OPEN-GPC1.txt tmp/OPEN-GPC1.txt
mv EA-Rietveld.txt tmp/EA-Rietveld.txt

## Files with no SE 
for file in EA-Rietveld.txt SELFHEALTH-Harris.txt; do
	Rscript $dirCode/2.3.1_Add_SE.R $dirOut/tmp/$file $dirOut/$file &
done
wait
mv EA-Rietveld.txt tmp/EA-Rietveld.txt

## Files with no Chr BP
for file in EA-Rietveld.txt OPEN-deMoor.txt DEP-Howard.txt; do
	awk -F"\t" 'NR==FNR{ChrPosID[$1]=$2;next} \
		FNR==1{print;next} FNR>1 && $1 in ChrPosID {split(ChrPosID[$1],a,":");$2=a[1];$3=a[2];print}' OFS="\t" /disk/genetics/ukb/aokbay/reffiles/HRC/HRC_r1-1.GRCh37.wgs.mac5.maf001.rsID_cptid_alleles_raf tmp/$file > $file &
done

## Files with no rsID
awk -F"\t" 'NR==FNR{a[$2]=$1;next} \
NR>1 && $1 in a {$1=a[$1];print;next}{print}' OFS="\t" /disk/genetics/ukb/aokbay/reffiles/HRC/HRC_r1-1.GRCh37.wgs.mac5.maf001.rsID_cptid_alleles_raf tmp/MENARCHE-Day.txt > MENARCHE-Day.txt 

rm tmp/*

gzip *.txt



