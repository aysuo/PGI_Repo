#!/bin/bash

source $PGI_Repo/code/paths
export R_LIBS=$Rlib:$R_LIBS

cd $PGI_Repo/derived_data/2_Formatted/public
mkdir tmp

#------------------------------------------------------------------------------------------------------------#
# -------------------------------------------- PRE-PROCESSING -----------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

# Get openness from GPC1, remove rest 
unzip $PGI_Repo/original_data/public/GPC-1.BigFiveNEO.zip
mv GPC-1.NEO-OPENNESS.full.txt tmp/OPEN-deMoor.txt
mv ReadmeGPC-1.pdf $PGI_Repo/original_data/public/ReadMe/OPEN-deMoor_ReadMe.pdf
rm -r __MACOSX GPC-1.NEO-CONSCIENTIOUSNESS.full.txt  GPC-1.NEO-NEUROTICISM.full.txt GPC-1.NEO-AGREEABLENESS.full.txt  GPC-1.NEO-EXTRAVERSION.full.txt

# Unzip BRCA and convert to tab-delimited
unzip $PGI_Repo/original_data/public/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt.zip
sed 's/ /\t/g' icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt > tmp/BRCA-Zhang.txt
rm -r __MACOSX*

# Extract  EUR IBD sumstats, and ReadMe, remove rest (trans-ethnic file contains only 157k SNPs)
tar -xvf $PGI_Repo/original_data/public/iibdgc-trans-ancestry-filtered-summary-stats.tgz EUR.IBD.gwas_info03_filtered.assoc README
mv README $PGI_Repo/original_data/public/ReadMe/IBD.ReadMe.txt
mv EUR.IBD.gwas_info03_filtered.assoc tmp/IBD-Liu.txt

# Extract CAD-Nikpay
unzip $PGI_Repo/original_data/public/cad.additive.Oct2015.pub.zip 
mv cad.add.readme $PGI_Repo/original_data/public/ReadMe/
mv cad.add.160614.website.txt tmp/CAD-Nikpay.txt

# Extract MIGRAINE
tar -xvf $PGI_Repo/original_data/public/transfer_2674936_files_a3984b76.tar
mv README* $PGI_Repo/original_data/public/ReadMe/
mv migraine_ihgc2021_all_passed_variant_ids_with_cohorts3.txt.gz $PGI_Repo/original_data/public/
gunzip -c migraine_IHGC2021_no23andMe_gwama_v2.txt.gz > tmp/MIGRAINE-Hautakangas.txt
gunzip -c any_mig.no23andMe.gwama.out.isq75.nstud10.clean.gz > tmp/MIGRAINE-Gormley.txt &
rm migraine* any_mig*

# Extract AUDIT
tar -xzvf $PGI_Repo/original_data/public/AUDIT-Walters.tar.gz pgc_alcdep.eur_discovery.aug2018_release.txt.gz
gunzip -c pgc_alcdep.eur_discovery.aug2018_release.txt.gz > tmp/AUDIT-Walters.txt 
rm pgc_alcdep.eur_discovery.aug2018_release.txt.gz

#-----------------------------------#
## Decompress and rename
gunzip -c $PGI_Repo/original_data/public/Doherty-2018-NatureComms-overall-activity.csv.gz > tmp/ACTIVITY-Doherty.txt &
gunzip -c $PGI_Repo/original_data/public/PA_GWAS_ID1.tbl.gz > tmp/ACTIVITY-Meddens.txt &
gunzip -c $PGI_Repo/original_data/public/ADHD-adhd_eur_jun2017.gz > tmp/ADHD-Demontis.txt &
gunzip -c $PGI_Repo/original_data/public/AFBpooled-AgeFirstBirth_Pooled.txt.gz > tmp/AFB-Barban.txt &
gunzip -c $PGI_Repo/original_data/public/AFB-Mills.txt.gz > tmp/AFB-Mills.txt &
gunzip -c $PGI_Repo/original_data/public/GCST90000047_buildGRCh37.tsv.gz > tmp/AFS-Mills.txt &
gunzip -c $PGI_Repo/original_data/public/PGCALZ2sumstatsExcluding23andMe.txt.gz > tmp/ALZ-Wightman.txt &
gunzip -c $PGI_Repo/original_data/public/ANOREX-Duncan.txt.gz > tmp/ANOREX-Duncan.txt &
gunzip -c $PGI_Repo/original_data/public/ANOREX-Watson.txt.gz > tmp/ANOREX-Watson.txt & 
gunzip -c $PGI_Repo/original_data/public/ASD-Grove.txt.gz > tmp/ASD-Grove.txt &
gunzip -c $PGI_Repo/original_data/public/AgeofInitiation.txt.gz > tmp/ASI-Liu.txt
gunzip -c $PGI_Repo/original_data/public/AgeOfInitiation.WithoutUKB.txt.gz > tmp/ASI-LiuSansUKB.txt
gunzip -c $PGI_Repo/original_data/public/AstEczRhi-SHARE-without23andMe.LDSCORE-GC.SE-META.v0.gz > tmp/ASTECZRHI-Ferreira.txt &
gunzip -c $PGI_Repo/original_data/public/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz > tmp/ASTHMA-Tsuo.txt
gunzip -c $PGI_Repo/original_data/public/BIPOLAR-Mullins.txt.gz > tmp/BIPOLAR-Mullins.txt &
gunzip -c $PGI_Repo/original_data/public/BIPOLAR-Stahl.txt.gz > tmp/BIPOLAR-Stahl.txt &
gunzip -c $PGI_Repo/original_data/public/BMI-SNP_gwas_mc_merge_nogc.tbl.uniq.gz > tmp/BMI-Locke.txt &
gunzip -c $PGI_Repo/original_data/public/Evangelou_30224653_DBP.txt.gz > tmp/BPdia-Evangelou.txt &
gunzip -c $PGI_Repo/original_data/public/Evangelou_30224653_SBP.txt.gz > tmp/BPsys-Evangelou.txt &
gunzip -c $PGI_Repo/original_data/public/Evangelou_30224653_PP.txt.gz > tmp/BPpulse-Evangelou.txt &
gunzip -c $PGI_Repo/original_data/public/ICBP_DBP_02082017.txt.gz > tmp/BPdia-EvangelouSansUKB.txt &
gunzip -c $PGI_Repo/original_data/public/ICBP_SBP_02082017.txt.gz > tmp/BPsys-EvangelouSansUKB.txt &
gunzip -c $PGI_Repo/original_data/public/ICBP_PP_02082017.txt.gz > tmp/BPpulse-EvangelouSansUKB.txt &
gunzip -c $PGI_Repo/original_data/public/oncoarray_bcac_public_release_oct17.txt.gz > tmp/BRCA-Michailidou.txt &
gunzip -c $PGI_Repo/original_data/public/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz > tmp/CAD-Nelson.txt
gunzip -c $PGI_Repo/original_data/public/CHOHDL-SinnottArmstrong.txt.gz > tmp/CHOHDL-SinnottArmstrong.txt &
gunzip -c $PGI_Repo/original_data/public/CHOLDL-SinnottArmstrong.txt.gz > tmp/CHOLDL-SinnottArmstrong.txt &
gunzip -c $PGI_Repo/original_data/public/CHOTOT-SinnottArmstrong.txt.gz > tmp/CHOTOT-SinnottArmstrong.txt &
gunzip -c $PGI_Repo/original_data/public/COPD_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz > tmp/COPD-Tsuo.txt &
gunzip -c $PGI_Repo/original_data/public/Intelligence-cogent.hrc.meta.chr.bp.rsid.assoc.full.gz > tmp/CP-Trampush.txt &
gunzip -c $PGI_Repo/original_data/public/CPD-tag.cpd.tbl.gz > tmp/CPD-Furberg.txt &
gunzip -c $PGI_Repo/original_data/public/CigarettesPerDay.txt.gz?sequence=31 > tmp/CPD-Liu.txt.gz &
gunzip -c $PGI_Repo/original_data/public/CigarettesPerDay.WithoutUKB.txt.gz > tmp/CPD-LiuSansUKB.txt &
gunzip -c $PGI_Repo/original_data/public/KP_DEPR_BETA_EAF.txt.gz > tmp/DEP-GERA.txt &
gunzip -c $PGI_Repo/original_data/public/daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.gz > tmp/DEP-WraySansUKB.txt &
gunzip -c $PGI_Repo/original_data/public/DrinksPerWeek.txt.gz?sequence=32 > tmp/DPW-Liu.txt.gz &
gunzip -c $PGI_Repo/original_data/public/DrinksPerWeek.WithoutUKB.txt.gz > tmp/DPW-LiuSansUKB.txt &
gunzip -c $PGI_Repo/original_data/public/EVERSMOKE-tag.evrsmk.tbl.gz > tmp/EVERSMOKE-Furberg.txt &
gunzip -c $PGI_Repo/original_data/public/SmokingInitiation.txt.gz?sequence=34 > tmp/EVERSMOKE-Liu.txt.gz &
gunzip -c $PGI_Repo/original_data/public/SmokingInitiation.WithoutUKB.txt.gz > tmp/EVERSMOKE-LiuSansUKB.txt &
gunzip -c $PGI_Repo/original_data/public/HEIGHT-GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz > tmp/HEIGHT-Wood.txt &
gunzip -c $PGI_Repo/original_data/public/Insomnia_sumstats_Jansenetal.txt.gz > tmp/INSOMNIA-Jansen.txt
gunzip -c $PGI_Repo/original_data/public/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020.txt.gz > tmp/NEARSIGHTED-Hysi.txt &
gunzip -c $PGI_Repo/original_data/public/NEBmen-NumberChildrenEverBorn_Male.txt.gz > tmp/NEBmen-Barban.txt &
gunzip -c $PGI_Repo/original_data/public/NEBwomen-NumberChildrenEverBorn_Female.txt.gz > tmp/NEBwomen-Barban.txt &
unzip $PGI_Repo/original_data/public/meta_v3_onco_euro_overall_ChrAll_1_release.zip && mv meta_v3_onco_euro_overall_ChrAll_1_release.txt tmp/PRCA-Schumacher.txt &
gunzip -c $PGI_Repo/original_data/public/RISK_GWAS_MA_Nweighted_ID15_2017_08_06.tbl.gz > tmp/RISK-Linner.txt &
gunzip -c $PGI_Repo/original_data/public/SCZ-PGC3.txt.gz > tmp/SCZ-PGC3.txt
gunzip -c $PGI_Repo/original_data/public/SmokingCessation.txt.gz > tmp/SMCESS-Liu.txt
gunzip -c $PGI_Repo/original_data/public/SmokingCessation.WithoutUKB.txt.gz > tmp/SMCESS-LiuSansUKB.txt
gunzip -c $PGI_Repo/original_data/public/SWB_excl_PGSrepo_ldscGC.meta.gz > tmp/SWB-Okbay.txt &
unzip $PGI_Repo/original_data/public/Mahajan.NatGenet2018b.T2D.European.zip && mv Mahajan.NatGenet2018b.T2D.European.txt tmp/T2D-Mahajan.txt &
unzip $PGI_Repo/original_data/public/METAANALYSIS_DIAGRAM_SE1.zip && mv METAANALYSIS_DIAGRAM_SE1.txt tmp/T2D-Scott.txt &
gunzip -c $PGI_Repo/original_data/public/TRYGL-SinnottArmstrong.txt.gz > tmp/TRYGL-SinnottArmstrong.txt
wait
#-----------------------------------#

#-----------------------------------#
# Copy into tmp and rename
cp $PGI_Repo/original_data/public/AGEFIRSTMENSES-Menarche_1KG_NatGen2017_WebsiteUpload.txt tmp/MENARCHE-Day.txt &
cp $PGI_Repo/original_data/public/AGEFIRSTMENSES-Menarche_Nature2014_GWASMetaResults_17122014.txt tmp/MENARCHE-Perry.txt &
cp $PGI_Repo/original_data/public/Kunkle_etal_Stage1_results.txt?file=1 tmp/ALZ-Kunkle.txt &
cp $PGI_Repo/original_data/public/EA4_excl_PGIrepo_2021_12_27.meta tmp/EA-OkbayExclPGIrepo.txt &
cp $PGI_Repo/original_data/public/EA4_excl_UKB_2020_04_09.meta tmp/EA-OkbayExclUKB.txt &
cp $PGI_Repo/original_data/public/GPC-2.EXTRAVERSION.full.txt tmp/EXTRA-vandenBerg.txt &
cp $PGI_Repo/original_data/public/META_NEUROTICISM_ALLIwv_iwv_20150402_1.dat tmp/NEURO-deMoor.txt &
cp $PGI_Repo/original_data/public/CCI_discovery_MA_13_samples_07-18-2015_rs.txt.txt tmp/CANNABIS-Stringer.txt &
wait
#-----------------------------------#

#-----------------------------------#
# CANNABIS-Pasman
# Per-chromosome files, no header (get from readme)
# Header: SNP (rs-number, or CHR:BP position on build GRCh37 if no rs-number was available), Allele1 (effect allele), Allele2 (reference allele), MAF (minor allele frequency), Effect (beta regression coefficient from the meta-analysis), StdErr (standard error), P (p-value), Direction (direction of the effect per sample, order: ICC, [23andMe,] UKB), Chr (position on chromosome), Bp (position in basepairs), and N (sample size, depending on in how many samples this SNP was present).
awk 'BEGIN{OFS="\t";print "SNPID","EFFECT_ALLELE","OTHER_ALLELE","MAF","BETA","SE","P","CHR","BP","N"} \
{print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11}' $PGI_Repo/original_data/public/Pasman_cannabis/cannabis_icc_ukb_chr*.txt > tmp/CANNABIS-Pasman.txt
#-----------------------------------#

#-----------------------------------#
# CP-UKB
# Tar file, per chromosome
tar -xzf $PGI_Repo/original_data/public/Intelligence-IQ_BOLT_LMM_UKB_v2_BGEN.tar.gz
awk -F"\t" 'NR==1{print}NR>1 && $1!="SNP"{print}' OFS="\t" IQ_BOLT_LMM_UKB_v2_BGEN_Chr* > tmp/CP-UKB.txt
rm IQ_BOLT_LMM_UKB_v2_BGEN_Chr*
#-----------------------------------#

#-----------------------------------#
# AUDIT-Walters, IBD-Liu, BP*-Evangelou, ALZ-Kunkle, CANNABIS-Stringer, ASTECZRHI-Ferreira, CANNABIS-Pasman, CP-Trampush
# Convert to tab-delimited
sed -i 's/ /\t/g' tmp/AUDIT-Walters.txt tmp/IBD-Liu.txt tmp/BPdia-Evangelou.txt tmp/BPsys-Evangelou.txt tmp/BPpulse-Evangelou.txt tmp/ALZ-Kunkle.txt tmp/CANNABIS-Stringer.txt tmp/ASTECZRHI-Ferreira.txt tmp/CANNABIS-Pasman.txt tmp/CP-Trampush.txt
sed -i 's/,/\t/g' tmp/ACTIVITY-Doherty.txt
#-----------------------------------#




#------------------------------------------------------------------------------------------------------------#
# ---------------------- FORMAT EACH FILE, IGNORE SNPID FORMAT OR EAF ISSUES --------------------------------#
#------------------------------------------------------------------------------------------------------------#


# ACTIVITY
# Doherty
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=91105;print $1,$2,$3,toupper($5),toupper($6),$7,$11,$12,$16,$8,N,1,1,1,"A"}' tmp/ACTIVITY-Doherty.txt > ACTIVITY-Doherty.txt &

# Meddens
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{split($1,a,":"); SE=1/sqrt(2*$8*$4*(1-$4)); BETA=$9*SE} \
	NR>1{print $1,a[1],a[2],toupper($2),toupper($3),$4,BETA,SE,$10,1,$8,1,1,1,"A"}' tmp/ACTIVITY-Meddens.txt > tmp/tmp_ACTIVITY-Meddens.txt &

#---------------------------------------#

# ADHD: 20,183 cases, 35,191 controls
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=20183+35191;print $2,$1,$3,$4,$5,"NA",log($7),$8,$9,$6,N,1,1,1,"A"}' tmp/ADHD-Demontis.txt > tmp/tmp_ADHD-Demontis.txt &

#---------------------------------------#

# AFB
# Barban 
# all (N=251,151)-9,370=241781
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1 && $6>0 && $6<1{N=241781;print $1,$2,$3,$4,$5,$6,$7/sqrt(2*N*$6*(1-$6)),1/sqrt(2*N*$6*(1-$6)),$8,1,N,1,1,1,"A"}' tmp/AFB-Barban.txt > AFB-Barban.txt &

# Mills
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1 {N=397338;print $1,$7,$8,$2,$3,EAF,$4,$5,$6,1,N,1,1,1,"A"}' tmp/AFB-Mills.txt > tmp/tmp_AFB-Mills.txt &

#---------------------------------------#

# AFS
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1 {N=397338;print $1,$7,$8,$2,$3,EAF,$4,$5,$6,1,N,1,1,1,"A"}' tmp/AFS-Mills.txt > tmp/tmp_AFS-Mills.txt &

#---------------------------------------#

# ALZ
# IGAP-Kunkle - 21,982 cases + 41,944 controls 
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=21892+41944;print $3,$1,$2,$4,$5,"NA",$6,$7,$8,1,N,1,1,1,"A"}' tmp/ALZ-Kunkle.txt > tmp/tmp_ALZ-Kunkle.txt &

# Wightman
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{print $1":"$2,$1,$2,$3,$4,"NA",$5,"NA",$6,1,$7,1,1,1,"A"}' tmp/ALZ-Wightman.txt > tmp/tmp_ALZ-Wightman.txt &
	# Add EAF from HRC (below)
	# Add Beta and SE using HRC EAF
awk -F"\t" 'NR==1{print;next} \
	{SE=1/sqrt(2*$11*$6*(1-$6)); BETA=$7*SE; print $1,$2,$3,$4,$5,$6,BETA,SE,$9,$10,$11,$12,$13,$14,$15}' OFS="\t" tmp/tmp_ALZ-Wightman.txt > ALZ-Wightman.txt &
#---------------------------------------#

# ANOREX
# Watson
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	$1!="CHROM" && !/#/ {N=$12+$13;print $3,$1,$2,$4,$5,"NA",$6,$7,$8,$10,N,1,1,1,"A"}' tmp/ANOREX-Watson.txt > tmp/tmp_ANOREX-Watson.txt &

# Duncan
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=3495+10982;print $1":"$3,$1,$3,$4,$5,"NA",log($7),$8,$9,$6,N,1,1,1,"A"}' tmp/ANOREX-Duncan.txt > tmp/tmp_ANOREX-Duncan.txt &

#---------------------------------------#

# ASD
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=18381+27969;print $2,$1,$3,$4,$5,"NA",log($7),$8,$9,$6,N,1,1,1,"A"}' tmp/ASD-Grove.txt > tmp/tmp_ASD-Grove.txt &

#---------------------------------------#

# ASI
for file in ASI-Liu.txt ASI-LiuSansUKB.txt; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR>1{print $3,$1,$2,$5,$4,$6,$9,$10,$8,$12/$11,$11,1,1,1,"A"}' OFS="\t" tmp/$file > $file &
done

mv ASI-LiuSansUKB.txt tmp/tmp_ASI-LiuSansUKB.txt
# Add EAF below for ~1.5mil SNPs with missing EAF 

#---------------------------------------#

# ASTECZRHI								
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1 && toupper($4)==$11{N=$NF;print $10,$2,$3,toupper($4),toupper($5),$12,$6,$7,$8,1,N,1,1,1,"A"} \
	NR>1 && toupper($5)==$11{N=$NF;print $10,$2,$3,toupper($4),toupper($5),1-$12,$6,$7,$8,1,N,1,1,1,"A"}' tmp/ASTECZRHI-Ferreira.txt > ASTECZRHI-Ferreira.txt &

#---------------------------------------#

# ASTHMA
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1 && $16=="no" && $17=="no"{N=$12+$13;print $5,$1,$2,$4,$3,$6,$7,$8,$9,1,N,1,1,1,"A"}' tmp/ASTHMA-Tsuo.txt > ASTHMA-Tsuo.txt &

#---------------------------------------#

# AUDIT
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=11569+34999; print $1":"$3,$1,$3,$4,$5,"NA",$6,"NA",$7,1,N,1,1,1,"A"}' tmp/AUDIT-Walters.txt > tmp/tmp_AUDIT-Walters.txt &
	# Add EAF from HRC (below)
	# Add Beta and SE using HRC EAF
awk -F"\t" 'NR==1{print;next} \
	{SE=1/sqrt(2*$11*$6*(1-$6)); BETA=$7*SE; print $1,$2,$3,$4,$5,$6,BETA,SE,$9,$10,$11,$12,$13,$14,$15}' OFS="\t" tmp/tmp_AUDIT-Walters.txt > AUDIT-Walters.txt &

#---------------------------------------#

# BIPOLAR
# Mullins
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1 && !/#/ {N=$14+$15;print $3,$1,$2,$4,$5,$11,$6,$7,$8,$12,N,1,1,1,"A"}' tmp/BIPOLAR-Mullins.txt > BIPOLAR-Mullins.txt &

# Stahl
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=$17+$18;print $2,$1,$3,$4,$5,$7,log($9),$10,$11,$8,N,1,1,1,"A"}' tmp/BIPOLAR-Stahl.txt > BIPOLAR-Stahl.txt &

#---------------------------------------#

# BMI
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=$NF;print $1,"NA","NA",$2,$3,$4,$5,$6,$7,1,N,1,1,1,"A"}' tmp/BMI-Locke.txt > tmp/tmp_BMI-Locke.txt &

#---------------------------------------#

# BPsys, BPdia, BPpulse
# Evangelou
for file in BPsys-Evangelou.txt BPdia-Evangelou.txt BPpulse-Evangelou.txt; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR>1{split($1,a,":"); print a[1]":"a[2],a[1],a[2],toupper($2),toupper($3),$4,$5,$6,$7,$9/$8,$8,1,1,1,"A"}' tmp/$file > tmp/tmp_$file &
done

for file in BPsys-EvangelouSansUKB.txt BPdia-EvangelouSansUKB.txt BPpulse-EvangelouSansUKB.txt; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR>1{split($1,a,":"); print a[1]":"a[2],a[1],a[2],toupper($2),toupper($3),$4,$6,$7,$8,1,$15/$14,1,1,1,"A"}' tmp/$file > tmp/tmp_$file &
done

#---------------------------------------#

# BRCA
# Meta-analysis of ICOGS, GWAS and ONCOarray
# ICOGS : N_controls=37818, N=76167
# GWAS: N_controls=17588, N=32498
# ONCO: N_controls=58383, N=138508
# Total: N_controls=133384, Ncases=113789
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1 {N=133384+113789; EAF=((37818*$13)+(58383*$29)+(32498*$4))/(37818+58383+32498); if ($15<$31) {INFO=$15} else {INFO=$31}} \
	NR>1 {print $9":"$10,$9,$10,$40,$41,EAF,$42,sqrt($43),$46,INFO,N,1,1,1,"A"}' tmp/BRCA-Zhang.txt > tmp/tmp_BRCA-Zhang.txt &

#---------------------------------------#

# CAD
# Nelson
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$13,$11,1,1,1,"A"}' tmp/CAD-Nelson.txt > CAD-Nelson.txt &

# Nikpay 
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=60801+123504; print $1,$2,$3,$4,$5,$6,$9,$10,$11,$7,N,1,1,1,"A"}' tmp/CAD-Nikpay.txt > CAD-Nikpay.txt &

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

# CHOLESTEROL AND TRYGLYCERIDS
for pheno in CHOHDL CHOLDL CHOTOT CHOHDL TRYGL; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR>1{N=355891;print $1":"$2,$1,$2,$5,$4,"NA",$6,$7,$8,$15,N,1,1,1,"A"}' tmp/$pheno-SinnottArmstrong.txt > tmp/tmp_$pheno-SinnottArmstrong.txt &
done

#---------------------------------------#

# COPD 
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1 && $16=="no" && $17=="no"{N=$12+$13;print $5,$1,$2,$4,$3,$6,$7,$8,$9,1,N,1,1,1,"A"}' tmp/COPD-Tsuo.txt > COPD-Tsuo.txt &

#---------------------------------------#

# CP
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=$13;print $2,$3,$4,$6,$7,$8,$10,$11,$12,$9,N,1,1,1,"A"}' tmp/CP-Trampush.txt > CP-Trampush.txt

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=222543;print $1,$2,$3,$5,$6,$7,$9,$10,$11,$8,N,1,1,1,"A"}' tmp/CP-UKB.txt > CP-UKB.txt

#---------------------------------------#

# CPD
# Furberg 
# chr=hg18, FRQ_A=FRQ_U (only had "freq1"), INFO=1 for all (unavailable), OR is not an OR !!! It's the linear regression beta for the continuous variables and logistic regression beta for the discrete variables. 
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=38181;print $2,"NA","NA",$4,$5,$6,$9,$10,$11,$8,N,1,1,1,"A"}' tmp/CPD-Furberg.txt > tmp/tmp_CPD-Furberg.txt &

# Liu
for file in CPD-Liu.txt CPD-LiuSansUKB.txt; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR>1{print $3,$1,$2,$5,$4,$6,$9,$10,$8,$12/$11,$11,1,1,1,"A"}' OFS="\t" tmp/$file > $file &
done
 
mv CPD-LiuSansUKB.txt tmp/tmp_CPD-LiuSansUKB.txt
# Add EAF below for ~1.5mil SNPs with missing EAF 

#---------------------------------------#

# DEP
# Ntot= 56,368, Ncases = 7,231, Ncontrols = 49,137
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	{N=56368; split($2, a, ":", seps)}{if ($2~/^[1-9]/ ||  $2~/^X:/) $2=a[1]":"a[2] ; else $2=a[1]} NR>1{print $2,$1,$3,$4,$5,$14,$13,$11,$12,$8,N,1,1,1,"A"}' tmp/DEP-GERA.txt > DEP-GERA.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=$17+$18;EAF=(($6*45396)+($7*97250))/(45396+97250); print $2,$1,$3,$4,$5,EAF,log($9),$10,$11,$8,N,1,1,1,"A"}' tmp/DEP-WraySansUKB.txt > DEP-WraySansUKB.txt &

#---------------------------------------#

# DPW
for file in DPW-Liu.txt DPW-LiuSansUKB.txt; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR>1{print $3,$1,$2,$5,$4,$6,$9,$10,$8,$12/$11,$11,1,1,1,"A"}' OFS="\t" tmp/$file > $file &
done

mv DPW-LiuSansUKB.txt tmp/tmp_DPW-LiuSansUKB.txt
# Add EAF below for ~1.5mil SNPs with missing EAF 

#---------------------------------------#

# EA
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{print $3,$1,$2,$5,$6,$7,$17,$18,$16,1,$11,1,1,1,"A"}' tmp/EA-OkbayExclPGIrepo.txt > EA-OkbayExclPGIrepo.txt &

awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{print $3,$1,$2,$5,$6,$7,$17,$18,$16,1,$11,1,1,1,"A"}' tmp/EA-OkbayExclUKB.txt > EA-OkbayExclUKB.txt &

#---------------------------------------#

# EXTRA
# N=63,030
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=63030;print $1,$2,$3,toupper($4),toupper($5),"NA",$6,$7,$8,1,N,1,1,1,"A"}' tmp/EXTRA-vandenBerg.txt > tmp/tmp_EXTRA-vandenBerg.txt &

#---------------------------------------#

# EVERSMOKE
# Furberg 
# N=74,035, 41969 cases 32066 controls
# chr=hg18, FRQ_A=FRQ_U (only had "freq1"), INFO=1 for all (unavailable), OR is not an OR !!! It's the linear regression beta for the continuous variables and logistic regression beta for the discrete variables. 
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=41969+32066;print $2,"NA","NA",$4,$5,$6,$9,$10,$11,$8,N,1,1,1,"A"}' tmp/EVERSMOKE-Furberg.txt > tmp/tmp_EVERSMOKE-Furberg.txt &

# Liu
for file in EVERSMOKE-Liu.txt EVERSMOKE-LiuSansUKB.txt; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR>1{print $3,$1,$2,$5,$4,$6,$9,$10,$8,$12/$11,$11,1,1,1,"A"}' OFS="\t" tmp/$file > $file &
done

mv EVERSMOKE-LiuSansUKB.txt tmp/tmp_EVERSMOKE-LiuSansUKB.txt
# Add EAF below for ~1.5mil SNPs with missing EAF 
#---------------------------------------#

# HEIGHT 
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=$NF;print $1,"NA","NA",$2,$3,$4,$5,$6,$7,1,N,1,1,1,"A"}' tmp/HEIGHT-Wood.txt > tmp/tmp_HEIGHT-Wood.txt &

#---------------------------------------#

# IBD
# N_cases_EUR=12882, N_controls_EUR=21770
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=12882+21770; print $2,$1,$3,$4,$5,$7,log($9),$10,$11,$8,N,1,1,1,"A"}' OFS="\t" tmp/IBD-Liu.txt > IBD-Liu.txt &

#---------------------------------------#

# INSOMNIA
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{print $1,$3,$4,$5,$6,$7,log($8),$9,$10,$12,$11,1,1,1,"A"}' tmp/INSOMNIA-Jansen.txt > INSOMNIA-Jansen.txt &

#---------------------------------------#

# MENARCHE
# 329,345 = ReproGen consortium (N = 179,117) + 23andMe (N = 76,831) + UK Biobank (N = 73,397) studies 
# Npublic = 329,345 - 76,831
# Day
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=252514;split($1,a,":");print a[1]":"a[2],a[1],a[2],toupper($2),toupper($3),"NA",$4,"NA",$5,1,N,1,1,1,"A"}' tmp/MENARCHE-Day.txt | sed 's/chr//g' >  tmp/tmp_MENARCHE-Day.txt &

# Perry
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=132989;print $1,"NA","NA",toupper($2),toupper($3),$4,$5,"NA",$6,1,N,1,1,1,"A"}' tmp/MENARCHE-Perry.txt > tmp/tmp_MENARCHE-Perry.txt &

#---------------------------------------#

# MIGRAINE
# Gormley
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{print $3,$1,$2,$4,$5,$6,$7,$8,$12,1,$18,1,1,1,"A"}' tmp/MIGRAINE-Gormley.txt > MIGRAINE-Gormley.txt &

# Hautakangas
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{print $1,$2,$3,$4,$5,$6,$7,$8,$12,1,$18,1,1,1,"A"}' tmp/MIGRAINE-Hautakangas.txt > MIGRAINE-Hautakangas.txt &

#---------------------------------------#

# NEARSIGHTED
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{split($1,a,":"); SNPID=a[2]":"a[3]; SE=1/sqrt(2*$6*$4*(1-$4)); BETA=$7*SE}
	NR>1{print SNPID,a[2],a[3],toupper($2),toupper($3),$4,BETA,SE,$8,1,$6,1,1,1,"A"}' OFS="\t" tmp/NEARSIGHTED-Hysi.txt > tmp/tmp_NEARSIGHTED-Hysi.txt &

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

# PRCA
# N_cases = 46939+32255 = 79194, N_controls = 27910+33202 = 61112
# N_total = 140306
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1 {N=140306; if ($16!="NA") INFO=$16 ; else INFO=1; print $3,$4,$5,toupper($6),toupper($7),$8,$12,$13,$14,INFO,N,1,1,1,"A"}' OFS="\t" tmp/PRCA-Schumacher.txt > PRCA-Schumacher.txt &

#---------------------------------------#

# RISK 
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1 && $6>0 && $6<1{N=$10;split($3,a,":");print $2,a[1],a[2],toupper($4),toupper($5),$6,$11/sqrt(2*N*$6*(1-$6)),1/sqrt(2*N*$6*(1-$6)),$12,1,N,1,1,1,"A"}' tmp/RISK-Linner.txt > RISK-Linner.txt &

#---------------------------------------#

# SCZ
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1{N=$17+$18;print $2,$1,$3,$4,$5,$7,log($9),$10,$11,$8,N,1,1,1,"A"}' tmp/SCZ-PGC3.txt > SCZ-PGC3.txt &

#---------------------------------------#

# SMCESS
for file in SMCESS-Liu.txt SMCESS-LiuSansUKB.txt; do
	awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
		NR>1{print $3,$1,$2,$5,$4,$6,$9,$10,$8,$12/$11,$11,1,1,1,"A"}' OFS="\t" tmp/$file > $file &
done

mv SMCESS-LiuSansUKB.txt tmp/tmp_SMCESS-LiuSansUKB.txt
# Add EAF below for ~1.5mil SNPs with missing EAF 

#---------------------------------------#

# SWB
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1 && $6>0 && $6<1{N=$10;print $1,$2,$3,toupper($4),toupper($5),$6,$11/sqrt(2*N*$6*(1-$6)),1/sqrt(2*N*$6*(1-$6)),$12,1,N,1,1,1,"A"}' tmp/SWB-Okbay.txt > SWB-Okbay.txt &

#---------------------------------------#

# T2D
# Mahajan
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1 {N=74124+824006;print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,N,1,1,1,"A"}' tmp/T2D-Mahajan.txt > tmp/tmp_T2D-Mahajan.txt &

# Scott
awk -F"\t" 'BEGIN{OFS="\t"; print "SNPID","CHR","BP","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","P","INFO","N","IMPUTED","CALLRATE","HWE_PVAL","PLOIDY"} \
	NR>1 {split($1,a,":"); print $1,a[1],a[2],$2,$3,"NA",$4,$5,$6,1,$7,1,1,1,"A"}' tmp/T2D-Scott.txt > tmp/tmp_T2D-Scott.txt &


wait


#------------------------------------------------------------------------------------------------------------#
# ----------------------------------- HANDLE SNPID FORMAT, EAF, SE ISSUES -----------------------------------#
#------------------------------------------------------------------------------------------------------------#

## EAF ISSUES
# Files with rsID and missing EAF
for file in ASI-LiuSansUKB.txt EVERSMOKE-LiuSansUKB.txt DPW-LiuSansUKB.txt CPD-LiuSansUKB.txt SMCESS-LiuSansUKB.txt ANOREX-Watson.txt ASD-Grove.txt AFB-Mills.txt AFS-Mills.txt ALZ-Kunkle.txt ADHD-Demontis.txt EXTRA-vandenBerg.txt OPEN-deMoor.txt ; do
	awk -F"\t" 'NR==FNR{ref[$1]=$3;alt[$1]=$4;raf[$1]=$5;next}
	FNR==1{print;next}
	$6=="NA" && ref[$1]==$4 && alt[$1]==$5 {$6=raf[$1];print;next}
	$6=="NA" && ref[$1]==$5 && alt[$1]==$4 {$6=1-raf[$1];print;next}
	{print}' OFS="\t" $HRC_raf tmp/tmp_$file > $file &
done
wait
mv OPEN-deMoor.txt tmp/tmp_OPEN-deMoor.txt

## Files with ChrPosID and missing EAF (also adding rsID)
for file in CHOHDL-SinnottArmstrong.txt CHOLDL-SinnottArmstrong.txt CHOTOT-SinnottArmstrong.txt CHOHDL-SinnottArmstrong.txt TRYGL-SinnottArmstrong.txt AUDIT-Walters.txt ALZ-Wightman.txt ANOREX-Duncan.txt T2D-Scott.txt NEURO-deMoor.txt MENARCHE-Day.txt; do
	awk -F"\t" 'NR==FNR{rs[$2]=$1;ref[$2]=$3;alt[$2]=$4;raf[$2]=$5;next}
	FNR==1{print;next}
	$6=="NA" && ref[$1]==$4 && alt[$1]==$5 {$6=raf[$1];$1=rs[$1];print;next}
	$6=="NA" && ref[$1]==$5 && alt[$1]==$4 {$6=1-raf[$1];$1=rs[$1];print}
	{print}' OFS="\t" $HRC_raf tmp/tmp_$file > $file &
done 
wait
mv MENARCHE-Day.txt tmp/tmp_MENARCHE-Day.txt
mv ALZ-Wightman.txt tmp/tmp_ALZ-Wightman.txt
mv AUDIT-Walters.txt tmp/tmp_AUDIT-Walters.txt
#---------------------------------------#

## Files with ChrPosID (missing rsID)
for file in MIGRAINE-Hautakangas.txt ACTIVITY-Meddens.txt NEARSIGHTED-Hysi.txt BPsys-Evangelou.txt BPsys-EvangelouSansUKB.txt BPpulse-Evangelou.txt BPpulse-EvangelouSansUKB.txt BPdia-Evangelou.txt BPdia-EvangelouSansUKB.txt T2D-Mahajan.txt BRCA-Zhang.txt; do
	awk -F"\t" 'NR==FNR{rs[$2]=$1;next} \
		FNR==1{print;next} \
		FNR>1 && $1 in rs {$1=rs[$1];print;next} {print}' OFS="\t" $HRC_raf tmp/tmp_$file > $file &
done

#---------------------------------------#

## Missing CHR/BP
for file in BMI-Locke.txt MENARCHE-Perry.txt HEIGHT-Wood.txt EVERSMOKE-Furberg.txt CPD-Furberg.txt OPEN-deMoor.txt; do
	awk -F"\t" 'NR==FNR{ChrPosID[$1]=$2;next} \
		FNR==1{print;next} FNR>1 && $1 in ChrPosID {split(ChrPosID[$1],a,":");$2=a[1];$3=a[2];print}' OFS="\t" $HRC_raf tmp/tmp_$file > $file &
done
mv MENARCHE-Perry.txt tmp/tmp_MENARCHE-Perry.txt

#---------------------------------------#

## Missing SE (when Beta is given)
for file in MENARCHE-Day.txt MENARCHE-Perry.txt; do
	Rscript $PGI_Repo/code/2_Formatting/2.3.1_Add_SE.R $PGI_Repo/derived_data/2_Formatted/public/tmp/tmp_$file $PGI_Repo/derived_data/2_Formatted/public/$file &
done
wait


rm tmp/*
gzip *.txt
