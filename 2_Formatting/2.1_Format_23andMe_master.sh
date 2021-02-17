#!/bin/bash

dirIn=/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/23andMe
dirOut=/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/2_Formatted/23andMe
dirCode=/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/code/2_Formatting
cd ${dirOut}

for version in 3.0 4.0 4.1 5.0 5.1 5.2 6.1 7.0; do
	gunzip ${dirIn}/v${version}/*.dat.gz

	filesIn=$(ls -1 ${dirIn}/v${version}/*dat | sed -z 's/\n/,/g' | sed 's/,$//g')
	
	filesOut=$(echo $filesIn | sed "s,${dirIn}/v${version}/,23andMe_,g")

	sh $dirCode/2.1.1_Format_23andMe_slave.sh \
		--annot_all "${dirIn}/v${version}/${version}_Annotation/all_snp_info-${version}.txt" \
		--annot_gt "${dirIn}/v${version}/${version}_Annotation/gt_snp_stat-${version}.txt" \
		--annot_im "${dirIn}/v${version}/${version}_Annotation/im_snp_stat-${version}.txt" \
		--annot_version $version \
		--gwas $filesIn \
		--out $filesOut &
done 

wait


#################################################################

## Files that came annotated
cd ${dirIn}/annotated/
gunzip *.gz

## Barban et al. NEB, AFB
## P-value column is -log_10(P)
for file in Barban*; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR>1{P=1/10^$9;print $1,$2,$3,$4,$5,$10,$7,$8,P,$16,$13,$14,$12,$11,"A"}' $file > ${dirOut}/23andMe_$file &
done 

## Demontis et al. ADHD
for file in Demontis*; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR>1{N=$6+$7;gsub(/chr/,"",$2);print $1,$2,$3,$8,$9,$10,$11,$12,$13,$14,N,1,1,1,"A"}' $file > ${dirOut}/23andMe_$file &
done

## Ferreira 2014. Hayfever
#4,230 cases, 10,842 controls
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=4230+10842;gsub(/chr/,"",$1);print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,N,1,1,1,"A"}' Ferreira_2014_asthma_with_hayfever-Hayfever.meta.dat > ${dirOut}/23andMe_Ferreira_2014_asthma_with_hayfever-Hayfever.meta.dat &


## Hinds et al. Allergy
## N=46646	

for i in cat dust_mites pollen; do
	awk -F"\t" -v N=46646 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{split($4,a,"/");gsub(/chr/,"",$2);print $1,$2,$3,a[1],a[2],$5,$8,$9,$10,$6,N,1,1,1,"A"}' Hinds_2013-${i}_allergy.dat > ${dirOut}/23andMe_Hinds_2013-${i}_allergy.dat &
done

#################################################################

## POST-FORMATTING 

## Morning person re-format (im.num.0 and im.num.1 is NA for imputed SNPs, so EAF is set to dose.b.0 for imputed SNPs and N is set to NA)
## Take N for imputed/genotyped SNPs from html files instead
awk -F"\t" '($12==1 && $11=="NA") {$11=89283;print;next} \
	($12==0 && $11=="NA") {$11=91967;print;next}
	{print}' OFS="\t" $dirOut/23andMe_Hu_2016-morning_person-3.0.dat > $dirOut/TMP/tmp_morning
mv $dirOut/TMP/tmp_morning $dirOut/23andMe_Hu_2016-morning_person-3.0.dat


## v4.0 files re-format (no im.num.0 and im.num.1 columns, EAF is set to dose.b.0 for imputed SNPs)
## Take N from html file
N_agreeableness="59176"
N_conscientiousness="59176"‬
N_extraversion="59225‬"
N_neuroticism="59206‬"
N_openness="59176‬"
N_age_first_menses="76831"‬
N_age_voice_deepened="55871"

for path in $dirOut/23andMe_Lo*;
do
	pheno=$(echo $path | rev | cut -f1 -d"_" | cut -d"-" -f2 | rev)
	eval N='$'{N_$pheno}
	awk -F"\t" -v N=$N '$11=="NA" {$11=N+0;print;next} {print}' OFS="\t" $path > $path.tmp &
done

for path in $dirOut/23andMe_Day*;
do
	pheno=$(echo $path | cut -f2 -d"-")
	eval N='$'{N_$pheno}
	awk -F"\t" -v N=$N '$11=="NA" {$11=N+0;print;next} {print}' OFS="\t" $path > $path.tmp &
done
wait

for path in $dirOut/*.tmp;
do
	out=$(echo $path | sed 's/.tmp//g')
	mv $path $out
done


##################################################################################

#Rename files
mv $dirOut/23andMe_Barban-23andMe.AFB.women.association-results.20130103.txt $dirOut/AFBwomen-23andMe.txt
mv $dirOut/23andMe_Barban-23andMe.NEB.men.association-results.20130126.txt $dirOut/NEBmen-23andMe.txt
mv $dirOut/23andMe_Barban-23andMe.NEB.women.association-results.20130103.txt $dirOut/NEBwomen-23andMe.txt
mv $dirOut/23andMe_Day_2015-age_first_menses-4.0.dat $dirOut/AgeFirstMenses-23andMe.txt
mv $dirOut/23andMe_Day_2015-age_voice_deepened-4.0.dat $dirOut/AgeVoiceDeepened-23andMe.txt
mv $dirOut/23andMe_Demontis_2017-23andMe_V2_adhd_20170720.txt $dirOut/ADHDv2-23andMe.txt
mv $dirOut/23andMe_Demontis_2017-23andMe_V3_adhd_20170720.txt $dirOut/ADHDv3-23andMe.txt
mv $dirOut/23andMe_Ferreira_2014_asthma_with_hayfever-Hayfever.meta.dat $dirOut/Hayfever-23andMe.txt
mv $dirOut/23andMe_Ferreira_2017_Asthma_hayfever_eczema-any_asthma_eczema_4.1V2.dat $dirOut/AsthEcz-23andMe.txt
mv $dirOut/23andMe_Ferreira_2017_Asthma_hayfever_eczema-any_asthma_eczema_rhinitis_4.1V2.dat $dirOut/AstEczRhi-23andMe.txt
mv $dirOut/23andMe_Ferreira_2017_Asthma_hayfever_eczema-any_asthma_rhinitis_4.1V2.dat $dirOut/AstRhi-23andMe.txt
mv $dirOut/23andMe_Ferreira_2017_Asthma_hayfever_eczema-any_rhinitis_eczema_4.1V2.dat $dirOut/RhiEcz-23andMe.txt
mv $dirOut/23andMe_general_health.dat $dirOut/SelfHealth-23andMe.txt
mv $dirOut/23andMe_Hinds_2013-cat_allergy.dat $dirOut/AllergyCat-23andMe.txt
mv $dirOut/23andMe_Hinds_2013-dust_mites_allergy.dat $dirOut/AllergyDust-23andMe.txt
mv $dirOut/23andMe_Hinds_2013-pollen_allergy.dat $dirOut/AllergyPollen-23andMe.txt
mv $dirOut/23andMe_Hu_2016-morning_person-3.0.dat $dirOut/Morning-23andMe.txt
mv $dirOut/23andMe_Hyde_2016-depression-5.0.dat $dirOut/DEP-23andMe.txt
mv $dirOut/23andMe_iqb.adventurous.dat $dirOut/Adventure-23andMe.txt
mv $dirOut/23andMe_iqb.age_started_reading.dat $dirOut/StartedReading-23andMe.txt
mv $dirOut/23andMe_iqb.comfort_taking_risks.dat $dirOut/Risk-23andMe.txt
mv $dirOut/23andMe_iqb.left_out_social_activity.dat $dirOut/LeftOut-23andMe.txt
mv $dirOut/23andMe_iqb.life_satisfaction.dat $dirOut/SWB-23andMe.txt
mv $dirOut/23andMe_iqb.narcissism.dat $dirOut/Narcissism-23andMe.txt
mv $dirOut/23andMe_iqb.religious.dat $dirOut/Religious-23andMe.txt
mv $dirOut/23andMe_Lee_2018-highest_math_over_25yo.dat $dirOut/HighMath-23andMe.txt
mv $dirOut/23andMe_Lee_2018-iqb.self_rated_math_ability.dat $dirOut/SelfMath-23andMe.txt
mv $dirOut/23andMe_Lo_2016-personality_agreeableness-4.0.dat $dirOut/Agreeable-23andMe.txt
mv $dirOut/23andMe_Lo_2016-personality_conscientiousness-4.0.dat $dirOut/Conscientious-23andMe.txt
mv $dirOut/23andMe_Lo_2016-personality_extraversion-4.0.dat $dirOut/Extraversion-23andMe.txt
mv $dirOut/23andMe_Lo_2016-personality_neuroticism-4.0.dat $dirOut/Neuro-23andMe.txt
mv $dirOut/23andMe_Lo_2016-personality_openness-4.0.dat $dirOut/Open-23andMe.txt
mv $dirOut/23andMe_new_mets_quant_norm.dat $dirOut/Activity-23andMe.txt
mv $dirOut/23andMe_Pasman_2018-SUC_marijuana_DecisionMaking_Filtered.5.2.dat $dirOut/Cannabis-23andMe.txt
mv $dirOut/23andMe_Pickrell_2016-migraine_diagnosis-4.1.dat $dirOut/Migraine-23andMe.txt
mv $dirOut/23andMe_Pickrell_2016-nearsightedness-4.1.dat $dirOut/Nearsighted-23andMe.txt
mv $dirOut/23andMe_recharge_by_socializing_qtl.dat $dirOut/Recharge-23andMe.txt
mv $dirOut/23andMe_SanchezRoige_2018-DelayDiscounting_MCQ_mean_Log10_K_DecisionMaking_Filtered.5.0.dat $dirOut/DelayDisc-23andMe.txt
mv $dirOut/23andMe_SanchezRoige_2019-AUDIT_Log10_total_clean_DecisionMaking_Filtered.dat $dirOut/AUDIT-23andMe.txt
mv $dirOut/23andMe_Warrier_2018-empathy_qt-4.1.dat $dirOut/Empathy-23andMe.txt
mv $dirOut/23andMe_Warrier_2018-empathy_qt_F-4.1.dat $dirOut/EmpathyWomen-23andMe.txt
mv $dirOut/23andMe_Warrier_2018-empathy_qt_M-4.1.dat $dirOut/EmpathyMen-23andMe.txt
mv $dirOut/23andMe_age_started_smoking.dat $dirOut/AgeStartedSmoking-23andMe.txt
mv $dirOut/23andMe_ever_tobacco_user.dat $dirOut/Eversmoke-23andMe.txt
mv $dirOut/23andMe_GSCAN_smoking_cessation.dat $dirOut/SmokingCessation-23andMe.txt
mv $dirOut/23andMe_packs_per_day.dat $dirOut/CPD-23andMe.txt
mv $dirOut/23andMe_GSCAN_drinking_per_week_tx.dat $dirOut/DPW-23andMe.txt

gzip $dirOut/* &
wait