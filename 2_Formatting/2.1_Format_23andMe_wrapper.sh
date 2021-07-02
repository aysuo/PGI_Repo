#!/bin/bash

cd $mainDir/derived_data/2_Formatted/23andMe

for version in 3.0 4.0 4.1 5.0 5.1 5.2 6.1 7.0; do
	gunzip $mainDir/original_data/23andMe/v${version}/*.dat.gz

	filesIn=$(ls -1 $mainDir/original_data/23andMe/v${version}/*dat | sed -z 's/\n/,/g' | sed 's/,$//g')
	
	filesOut=$(echo $filesIn | sed "s,$mainDir/original_data/23andMe/v${version}/,23andMe_,g")

	sh $mainDir/code/2_Formatting/2.1.1_Format_23andMe.sh \
		--annot_all "$mainDir/original_data/23andMe/v${version}/${version}_Annotation/all_snp_info-${version}.txt" \
		--annot_gt "$mainDir/original_data/23andMe/v${version}/${version}_Annotation/gt_snp_stat-${version}.txt" \
		--annot_im "$mainDir/original_data/23andMe/v${version}/${version}_Annotation/im_snp_stat-${version}.txt" \
		--annot_version $version \
		--gwas $filesIn \
		--out $filesOut &
done 

wait

#################################################################

## Files that came annotated, some with missing columns

cd $mainDir/original_data/23andMe/annotated/
gunzip *.gz

## Barban et al. NEB, AFB
## P-value column is -log_10(P)
for file in Barban*; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR>1{P=1/10^$9;print $1,$2,$3,$4,$5,$10,$7,$8,P,$16,$13,$14,$12,$11,"A"}' $file > $mainDir/derived_data/2_Formatted/23andMe/23andMe_$file &
done 

## Demontis et al. ADHD
for file in Demontis*; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR>1{N=$6+$7;gsub(/chr/,"",$2);print $1,$2,$3,$8,$9,$10,$11,$12,$13,$14,N,1,1,1,"A"}' $file > $mainDir/derived_data/2_Formatted/23andMe/23andMe_$file &
done

## Ferreira 2014. Hayfever
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=4230+10842;gsub(/chr/,"",$1);print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,N,1,1,1,"A"}' Ferreira_2014_asthma_with_hayfever-Hayfever.meta.dat > $mainDir/derived_data/2_Formatted/23andMe/23andMe_Ferreira_2014_asthma_with_hayfever-Hayfever.meta.dat &


## Hinds et al. Allergy
for i in cat dust_mites pollen; do
	awk -F"\t" -v N=46646 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{split($4,a,"/");gsub(/chr/,"",$2);print $1,$2,$3,a[1],a[2],$5,$8,$9,$10,$6,N,1,1,1,"A"}' Hinds_2013-${i}_allergy.dat > $mainDir/derived_data/2_Formatted/23andMe/23andMe_Hinds_2013-${i}_allergy.dat &
done

#################################################################

## POST-FORMATTING 

## Morning person re-format (im.num.0 and im.num.1 is NA for imputed SNPs, so EAF is set to dose.b.0 for imputed SNPs and N is set to NA)
## Take N for imputed/genotyped SNPs from html files instead
awk -F"\t" '($12==1 && $11=="NA") {$11=89283;print;next} \
	($12==0 && $11=="NA") {$11=91967;print;next}
	{print}' OFS="\t" $mainDir/derived_data/2_Formatted/23andMe/23andMe_Hu_2016-morning_person-3.0.dat > $mainDir/derived_data/2_Formatted/23andMe/TMP/tmp_morning
mv $mainDir/derived_data/2_Formatted/23andMe/TMP/tmp_morning $mainDir/derived_data/2_Formatted/23andMe/23andMe_Hu_2016-morning_person-3.0.dat


## v4.0 files re-format (no im.num.0 and im.num.1 columns, EAF is set to dose.b.0 for imputed SNPs)
## Take N from html file
N_agreeableness="59176"
N_conscientiousness="59176"‬
N_extraversion="59225‬"
N_neuroticism="59206‬"
N_openness="59176‬"
N_age_first_menses="76831"‬
N_age_voice_deepened="55871"

for path in $mainDir/derived_data/2_Formatted/23andMe/23andMe_Lo*;
do
	pheno=$(echo $path | rev | cut -f1 -d"_" | cut -d"-" -f2 | rev)
	eval N='$'{N_$pheno}
	awk -F"\t" -v N=$N '$11=="NA" {$11=N+0;print;next} {print}' OFS="\t" $path > $path.tmp &
done

for path in $mainDir/derived_data/2_Formatted/23andMe/23andMe_Day*;
do
	pheno=$(echo $path | cut -f2 -d"-")
	eval N='$'{N_$pheno}
	awk -F"\t" -v N=$N '$11=="NA" {$11=N+0;print;next} {print}' OFS="\t" $path > $path.tmp &
done
wait

for path in $mainDir/derived_data/2_Formatted/23andMe/*.tmp;
do
	out=$(echo $path | sed 's/.tmp//g')
	mv $path $out
done

##################################################################################

#Rename files
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Barban-23andMe.AFB.women.association-results.20130103.txt $mainDir/derived_data/2_Formatted/23andMe/AFB-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Barban-23andMe.NEB.men.association-results.20130126.txt $mainDir/derived_data/2_Formatted/23andMe/NEBmen-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Barban-23andMe.NEB.women.association-results.20130103.txt $mainDir/derived_data/2_Formatted/23andMe/NEBwomen-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Day_2015-age_first_menses-4.0.dat $mainDir/derived_data/2_Formatted/23andMe/MENARCHE-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Day_2015-age_voice_deepened-4.0.dat $mainDir/derived_data/2_Formatted/23andMe/VOICEDEEP-23andMe.txt
#mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Demontis_2017-23andMe_V2_adhd_20170720.txt $mainDir/derived_data/2_Formatted/23andMe/ADHDv2-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Demontis_2017-23andMe_V3_adhd_20170720.txt $mainDir/derived_data/2_Formatted/23andMe/ADHD-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Ferreira_2014_asthma_with_hayfever-Hayfever.meta.dat $mainDir/derived_data/2_Formatted/23andMe/HAYFEVER-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Ferreira_2017_Asthma_hayfever_eczema-any_asthma_eczema_rhinitis_4.1V2.dat $mainDir/derived_data/2_Formatted/23andMe/ASTECZRHI-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_general_health.dat $mainDir/derived_data/2_Formatted/23andMe/SELFHEALTH-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Hinds_2013-cat_allergy.dat $mainDir/derived_data/2_Formatted/23andMe/ALLERGYCAT-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Hinds_2013-dust_mites_allergy.dat $mainDir/derived_data/2_Formatted/23andMe/ALLERGYDUST-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Hinds_2013-pollen_allergy.dat $mainDir/derived_data/2_Formatted/23andMe/ALLERGYPOLLEN-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Hu_2016-morning_person-3.0.dat $mainDir/derived_data/2_Formatted/23andMe/MORNING-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Hyde_2016-depression-5.0.dat $mainDir/derived_data/2_Formatted/23andMe/DEP-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_iqb.adventurous.dat $mainDir/derived_data/2_Formatted/23andMe/ADVENTURE-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_iqb.age_started_reading.dat $mainDir/derived_data/2_Formatted/23andMe/READING-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_iqb.comfort_taking_risks.dat $mainDir/derived_data/2_Formatted/23andMe/RISK-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_iqb.left_out_social_activity.dat $mainDir/derived_data/2_Formatted/23andMe/LEFTOUT-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_iqb.life_satisfaction.dat $mainDir/derived_data/2_Formatted/23andMe/SWB-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_iqb.narcissism.dat $mainDir/derived_data/2_Formatted/23andMe/NARCIS-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_iqb.religious.dat $mainDir/derived_data/2_Formatted/23andMe/RELIGBLF-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Lee_2018-highest_math_over_25yo.dat $mainDir/derived_data/2_Formatted/23andMe/HIGHMATH-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Lee_2018-iqb.self_rated_math_ability.dat $mainDir/derived_data/2_Formatted/23andMe/SELFMATH-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Lo_2016-personality_agreeableness-4.0.dat $mainDir/derived_data/2_Formatted/23andMe/AGREE-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Lo_2016-personality_conscientiousness-4.0.dat $mainDir/derived_data/2_Formatted/23andMe/CONSC-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Lo_2016-personality_extraversion-4.0.dat $mainDir/derived_data/2_Formatted/23andMe/EXTRA-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Lo_2016-personality_neuroticism-4.0.dat $mainDir/derived_data/2_Formatted/23andMe/NEURO-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Lo_2016-personality_openness-4.0.dat $mainDir/derived_data/2_Formatted/23andMe/OPEN-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_new_mets_quant_norm.dat $mainDir/derived_data/2_Formatted/23andMe/ACTIVITY-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Pasman_2018-SUC_marijuana_DecisionMaking_Filtered.5.2.dat $mainDir/derived_data/2_Formatted/23andMe/CANNABIS-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Pickrell_2016-migraine_diagnosis-4.1.dat $mainDir/derived_data/2_Formatted/23andMe/MIGRAINE-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Pickrell_2016-nearsightedness-4.1.dat $mainDir/derived_data/2_Formatted/23andMe/NEARSIGHTED-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_recharge_by_socializing_qtl.dat $mainDir/derived_data/2_Formatted/23andMe/RECHARGE-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_SanchezRoige_2018-DelayDiscounting_MCQ_mean_Log10_K_DecisionMaking_Filtered.5.0.dat $mainDir/derived_data/2_Formatted/23andMe/DELAYDISC-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_SanchezRoige_2019-AUDIT_Log10_total_clean_DecisionMaking_Filtered.dat $mainDir/derived_data/2_Formatted/23andMe/AUDIT-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_Warrier_2018-empathy_qt-4.1.dat $mainDir/derived_data/2_Formatted/23andMe/COGEMP-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_ever_tobacco_user.dat $mainDir/derived_data/2_Formatted/23andMe/EVERSMOKE-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_packs_per_day.dat $mainDir/derived_data/2_Formatted/23andMe/CPD-23andMe.txt
mv $mainDir/derived_data/2_Formatted/23andMe/23andMe_GSCAN_drinking_per_week_tx.dat $mainDir/derived_data/2_Formatted/23andMe/DPW-23andMe.txt

gzip $mainDir/derived_data/2_Formatted/23andMe/* &
wait