#!/bin/bash

source paths2

cd ${p2_23andMeOut}

for version in 3.0 4.0 4.1 5.0 5.1 5.2 6.1 7.0; do
	gunzip ${p2_23andMeIn}/v${version}/*.dat.gz

	filesIn=$(ls -1 ${p2_23andMeIn}/v${version}/*dat | sed -z 's/\n/,/g' | sed 's/,$//g')
	
	filesOut=$(echo $filesIn | sed "s,${p2_23andMeIn}/v${version}/,23andMe_,g")

	sh $p2_code/2.1.1_Format_23andMe_slave.sh \
		--annot_all "${p2_23andMeIn}/v${version}/${version}_Annotation/all_snp_info-${version}.txt" \
		--annot_gt "${p2_23andMeIn}/v${version}/${version}_Annotation/gt_snp_stat-${version}.txt" \
		--annot_im "${p2_23andMeIn}/v${version}/${version}_Annotation/im_snp_stat-${version}.txt" \
		--annot_version $version \
		--gwas $filesIn \
		--out $filesOut &
done 

wait

#################################################################

## Files that came annotated, some with missing columns

cd ${p2_23andMeIn}/annotated/
gunzip *.gz

## Barban et al. NEB, AFB
## P-value column is -log_10(P)
for file in Barban*; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR>1{P=1/10^$9;print $1,$2,$3,$4,$5,$10,$7,$8,P,$16,$13,$14,$12,$11,"A"}' $file > ${p2_23andMeOut}/23andMe_$file &
done 

## Demontis et al. ADHD
for file in Demontis*; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR>1{N=$6+$7;gsub(/chr/,"",$2);print $1,$2,$3,$8,$9,$10,$11,$12,$13,$14,N,1,1,1,"A"}' $file > ${p2_23andMeOut}/23andMe_$file &
done

## Ferreira 2014. Hayfever
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=4230+10842;gsub(/chr/,"",$1);print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,N,1,1,1,"A"}' Ferreira_2014_asthma_with_hayfever-Hayfever.meta.dat > ${p2_23andMeOut}/23andMe_Ferreira_2014_asthma_with_hayfever-Hayfever.meta.dat &


## Hinds et al. Allergy
for i in cat dust_mites pollen; do
	awk -F"\t" -v N=46646 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{split($4,a,"/");gsub(/chr/,"",$2);print $1,$2,$3,a[1],a[2],$5,$8,$9,$10,$6,N,1,1,1,"A"}' Hinds_2013-${i}_allergy.dat > ${p2_23andMeOut}/23andMe_Hinds_2013-${i}_allergy.dat &
done

#################################################################

## POST-FORMATTING 

## Morning person re-format (im.num.0 and im.num.1 is NA for imputed SNPs, so EAF is set to dose.b.0 for imputed SNPs and N is set to NA)
## Take N for imputed/genotyped SNPs from html files instead
awk -F"\t" '($12==1 && $11=="NA") {$11=89283;print;next} \
	($12==0 && $11=="NA") {$11=91967;print;next}
	{print}' OFS="\t" $p2_23andMeOut/23andMe_Hu_2016-morning_person-3.0.dat > $p2_23andMeOut/TMP/tmp_morning
mv $p2_23andMeOut/TMP/tmp_morning $p2_23andMeOut/23andMe_Hu_2016-morning_person-3.0.dat


## v4.0 files re-format (no im.num.0 and im.num.1 columns, EAF is set to dose.b.0 for imputed SNPs)
## Take N from html file
N_agreeableness="59176"
N_conscientiousness="59176"‬
N_extraversion="59225‬"
N_neuroticism="59206‬"
N_openness="59176‬"
N_age_first_menses="76831"‬
N_age_voice_deepened="55871"

for path in $p2_23andMeOut/23andMe_Lo*;
do
	pheno=$(echo $path | rev | cut -f1 -d"_" | cut -d"-" -f2 | rev)
	eval N='$'{N_$pheno}
	awk -F"\t" -v N=$N '$11=="NA" {$11=N+0;print;next} {print}' OFS="\t" $path > $path.tmp &
done

for path in $p2_23andMeOut/23andMe_Day*;
do
	pheno=$(echo $path | cut -f2 -d"-")
	eval N='$'{N_$pheno}
	awk -F"\t" -v N=$N '$11=="NA" {$11=N+0;print;next} {print}' OFS="\t" $path > $path.tmp &
done
wait

for path in $p2_23andMeOut/*.tmp;
do
	out=$(echo $path | sed 's/.tmp//g')
	mv $path $out
done

##################################################################################

#Rename files
mv $p2_23andMeOut/23andMe_Barban-23andMe.AFB.women.association-results.20130103.txt $p2_23andMeOut/AFB-23andMe.txt
mv $p2_23andMeOut/23andMe_Barban-23andMe.NEB.men.association-results.20130126.txt $p2_23andMeOut/NEBmen-23andMe.txt
mv $p2_23andMeOut/23andMe_Barban-23andMe.NEB.women.association-results.20130103.txt $p2_23andMeOut/NEBwomen-23andMe.txt
mv $p2_23andMeOut/23andMe_Day_2015-age_first_menses-4.0.dat $p2_23andMeOut/MENARCHE-23andMe.txt
mv $p2_23andMeOut/23andMe_Day_2015-age_voice_deepened-4.0.dat $p2_23andMeOut/VOICEDEEP-23andMe.txt
#mv $p2_23andMeOut/23andMe_Demontis_2017-23andMe_V2_adhd_20170720.txt $p2_23andMeOut/ADHDv2-23andMe.txt
mv $p2_23andMeOut/23andMe_Demontis_2017-23andMe_V3_adhd_20170720.txt $p2_23andMeOut/ADHD-23andMe.txt
mv $p2_23andMeOut/23andMe_Ferreira_2014_asthma_with_hayfever-Hayfever.meta.dat $p2_23andMeOut/HAYFEVER-23andMe.txt
mv $p2_23andMeOut/23andMe_Ferreira_2017_Asthma_hayfever_eczema-any_asthma_eczema_rhinitis_4.1V2.dat $p2_23andMeOut/ASTECZRHI-23andMe.txt
mv $p2_23andMeOut/23andMe_general_health.dat $p2_23andMeOut/SELFHEALTH-23andMe.txt
mv $p2_23andMeOut/23andMe_Hinds_2013-cat_allergy.dat $p2_23andMeOut/ALLERGYCAT-23andMe.txt
mv $p2_23andMeOut/23andMe_Hinds_2013-dust_mites_allergy.dat $p2_23andMeOut/ALLERGYDUST-23andMe.txt
mv $p2_23andMeOut/23andMe_Hinds_2013-pollen_allergy.dat $p2_23andMeOut/ALLERGYPOLLEN-23andMe.txt
mv $p2_23andMeOut/23andMe_Hu_2016-morning_person-3.0.dat $p2_23andMeOut/MORNING-23andMe.txt
mv $p2_23andMeOut/23andMe_Hyde_2016-depression-5.0.dat $p2_23andMeOut/DEP-23andMe.txt
mv $p2_23andMeOut/23andMe_iqb.adventurous.dat $p2_23andMeOut/ADVENTURE-23andMe.txt
mv $p2_23andMeOut/23andMe_iqb.age_started_reading.dat $p2_23andMeOut/READING-23andMe.txt
mv $p2_23andMeOut/23andMe_iqb.comfort_taking_risks.dat $p2_23andMeOut/RISK-23andMe.txt
mv $p2_23andMeOut/23andMe_iqb.left_out_social_activity.dat $p2_23andMeOut/LEFTOUT-23andMe.txt
mv $p2_23andMeOut/23andMe_iqb.life_satisfaction.dat $p2_23andMeOut/SWB-23andMe.txt
mv $p2_23andMeOut/23andMe_iqb.narcissism.dat $p2_23andMeOut/NARCIS-23andMe.txt
mv $p2_23andMeOut/23andMe_iqb.religious.dat $p2_23andMeOut/RELIGBLF-23andMe.txt
mv $p2_23andMeOut/23andMe_Lee_2018-highest_math_over_25yo.dat $p2_23andMeOut/HIGHMATH-23andMe.txt
mv $p2_23andMeOut/23andMe_Lee_2018-iqb.self_rated_math_ability.dat $p2_23andMeOut/SELFMATH-23andMe.txt
mv $p2_23andMeOut/23andMe_Lo_2016-personality_agreeableness-4.0.dat $p2_23andMeOut/AGREE-23andMe.txt
mv $p2_23andMeOut/23andMe_Lo_2016-personality_conscientiousness-4.0.dat $p2_23andMeOut/CONSC-23andMe.txt
mv $p2_23andMeOut/23andMe_Lo_2016-personality_extraversion-4.0.dat $p2_23andMeOut/EXTRA-23andMe.txt
mv $p2_23andMeOut/23andMe_Lo_2016-personality_neuroticism-4.0.dat $p2_23andMeOut/NEURO-23andMe.txt
mv $p2_23andMeOut/23andMe_Lo_2016-personality_openness-4.0.dat $p2_23andMeOut/OPEN-23andMe.txt
mv $p2_23andMeOut/23andMe_new_mets_quant_norm.dat $p2_23andMeOut/ACTIVITY-23andMe.txt
mv $p2_23andMeOut/23andMe_Pasman_2018-SUC_marijuana_DecisionMaking_Filtered.5.2.dat $p2_23andMeOut/CANNABIS-23andMe.txt
mv $p2_23andMeOut/23andMe_Pickrell_2016-migraine_diagnosis-4.1.dat $p2_23andMeOut/MIGRAINE-23andMe.txt
mv $p2_23andMeOut/23andMe_Pickrell_2016-nearsightedness-4.1.dat $p2_23andMeOut/NEARSIGHTED-23andMe.txt
mv $p2_23andMeOut/23andMe_recharge_by_socializing_qtl.dat $p2_23andMeOut/RECHARGE-23andMe.txt
mv $p2_23andMeOut/23andMe_SanchezRoige_2018-DelayDiscounting_MCQ_mean_Log10_K_DecisionMaking_Filtered.5.0.dat $p2_23andMeOut/DELAYDISC-23andMe.txt
mv $p2_23andMeOut/23andMe_SanchezRoige_2019-AUDIT_Log10_total_clean_DecisionMaking_Filtered.dat $p2_23andMeOut/AUDIT-23andMe.txt
mv $p2_23andMeOut/23andMe_Warrier_2018-empathy_qt-4.1.dat $p2_23andMeOut/COGEMP-23andMe.txt
mv $p2_23andMeOut/23andMe_ever_tobacco_user.dat $p2_23andMeOut/EVERSMOKE-23andMe.txt
mv $p2_23andMeOut/23andMe_packs_per_day.dat $p2_23andMeOut/CPD-23andMe.txt
mv $p2_23andMeOut/23andMe_GSCAN_drinking_per_week_tx.dat $p2_23andMeOut/DPW-23andMe.txt

gzip $p2_23andMeOut/* &
wait