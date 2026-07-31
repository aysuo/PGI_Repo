source $PGI_Repo/code/paths
mkdir -p $gf_dir_genotyped_STR/psych/preImputationQC/tmp $gf_dir_genotyped_STR/psych/preImputationQC/reffiles
cd $gf_dir_genotyped_STR/psych/preImputationQC

## Download support files for Infinium Global Screening Array v1.0 
cd reffiles
#Loci Name to rsID Conversion File (GRCh37) (Note: The array used is v1.0 but support files for versions <1.2 are not available)
wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/infinium-psycharray/v1-2/infinium-psycharray-24-v1-2-a1-b144-rsids.zip
unzip infinium-psycharray-24-v1-2-a1-b144-rsids.zip
awk -F"\t" '{split($2,a,","); print $1,a[1]}' OFS="\t" InfiniumPsychArray-24v1-2_A1_b144_rsids.txt > ../tmp/InfiniumPsychArray-24v1-2_A1_b144_rsids_single.txt

cd ..
# Per cohort SNP level QC, exclude SNPs removed by illumina, convert CATSS-PS to bed format
plink2 --bfile $gf_orig_genotyped_STRsalty --maf 0.0001 --geno 0.05 --hwe midp 1e-15 --make-bed --out tmp/STRsalty
plink2 --pfile $gf_orig_genotyped_STRcatssPS --maf 0.0001 --geno 0.05 --hwe midp 1e-10 --make-bed --out tmp/STRcatssPS

# Update SNP names, drop duplicated ID, duplicated BP, get allele freq 
for study in salty catssPS
do
	plink1.9 --bfile tmp/STR$study --update-name tmp/InfiniumPsychArray-24v1-2_A1_b144_rsids_single.txt --make-bed --out tmp/STR${study}_rsID
	plink2 --bfile tmp/STR${study}_rsID --set-missing-var-ids @:# --make-bed --out tmp/STR${study}_rsID_noMissID
	plink2 --bfile tmp/STR${study}_rsID_noMissID --rm-dup list force-first --make-bed --out tmp/STR${study}_rsID_noMissID_noDupID
	mv tmp/STR${study}_rsID_noMissID_noDupID.rmdup.list exclude_dupID_STR${study}.snps	
	plink1.9 --bfile tmp/STR${study}_rsID_noMissID_noDupID --list-duplicate-vars suppress-first ids-only --out tmp/STR${study}_rsID_noMissID_noDupID_BPdups
	sed 's/ /\n/g' tmp/STR${study}_rsID_noMissID_noDupID_BPdups.dupvar > exclude_dupBP_STR${study}.snps
	plink2 --bfile tmp/STR${study}_rsID_noMissID_noDupID --exclude exclude_dupBP_STR${study}.snps --make-bed --out tmp/STR${study}_rsID_noMissID_noDupID_noDupBP
	plink1.9 --bfile tmp/STR${study}_rsID_noMissID_noDupID_noDupBP --freq --out tmp/STR${study}_MAF
done

# Get Z scores for differences in MAF across cohorts
awk 'NR==FNR{maf[$2]=$5;n[$2]=$6;next}FNR==1{print "SNP","MAF_salty","N_salty","MAF_catss","N_catss";next} \
	($2 in maf && ($5!=0 || maf[$2]!=0)){print $2,$5,$6,maf[$2],n[$2]}' tmp/STRcatssPS_MAF.frq tmp/STRsalty_MAF.frq > tmp/tmp
awk 'NR==1{print $0, "Z"; next}
	{Z=($2-$4)/sqrt(($2*(1-$2)/$3)+($4*(1-$4)/$5)) ; $1=$1; print $0, Z}' OFS="\t" tmp/tmp > MAF_report.txt

# Get list of SNPs for which the P-value for H0: MAFdiff=0 is less than 0.001
awk -F"\t" '$6<-3.29 || $6>3.29 {print $1}' MAF_report.txt > exclude_MAFdiff_p001.snps

#--------------------------------------------------------------------#

## Merge SALTY and CATSS-PS

# Attempt merge and get list of problematic SNPs
plink1.9 --bfile tmp/STRsalty_rsID_noMissID_noDupID_noDupBP \
	--bmerge tmp/STRcatssPS_rsID_noMissID_noDupID_noDupBP \
	--make-bed \
	--out tmp/STRpsych
mv tmp/STRpsych-merge.missnp exclude_merge_missnp.snps

# Add merge-problematic SNPs to SNPs with significant MAF difference
cat exclude_MAFdiff_p001.snps exclude_merge_missnp.snps | sort | uniq > tmp/exclude_MAFdiff_p001_merge_missnp.snps

# Exclude problematic SNPs from each study
for study in salty catssPS
do
    plink1.9 --bfile tmp/STR${study}_rsID_noMissID_noDupID_noDupBP --exclude tmp/exclude_MAFdiff_p001_merge_missnp.snps --make-bed --out tmp/STR$study.nomissnp
done

# Merge
plink1.9 --bfile tmp/STRsalty.nomissnp \
	--bmerge tmp/STRcatssPS.nomissnp \
    --make-bed \
    --out tmp/STRpsych

#--------------------------------------------------------------------#

## Variant QC: call rate> 0.95 , maf>0.001, hwe 1e-10
plink2 --bfile tmp/STRpsych \
    --geno 0.05 \
    --maf 0.001 \
    --hwe midp 1e-10 \
    --make-bed \
    --out tmp/STRpsych_geno05_maf001_hwe1e-10


#--------------------------------------------------------------------#

## Sex check
plink1.9 --bfile tmp/STRpsych_geno05_maf001_hwe1e-10 \
    --chr X \
    --maf 0.01 \
    --indep-pairwise 1000kb 5 0.5 \
    --out tmp/STRpsych_geno05_maf001_hwe1e-10_pruned

plink1.9 --bfile tmp/STRpsych_geno05_maf001_hwe1e-10 \
    --extract tmp/STRpsych_geno05_maf001_hwe1e-10_pruned.prune.in \
    --check-sex \
    --out tmp/STRpsych_geno05_maf001_hwe1e-10_pruned_sexcheck

awk '$5!="OK"{print $1,$2}' tmp/STRpsych_geno05_maf001_hwe1e-10_pruned_sexcheck.sexcheck > exclude_sex_mismatch.samples


#--------------------------------------------------------------------#

## Get individuals with per chromosome subject-level missingness rate > 0.05
for chr in {1..22}
do
    plink2 -bfile tmp/STRpsych_geno05_maf001_hwe1e-10 \
        --chr $chr \
        --missing \
        --out tmp/STRpsych_geno05_maf001_hwe1e-10_chr$chr &
done
wait

for chr in {1..22}
do
    awk -F"\t" 'NR>1 && $5>0.05{print $1,$2}' tmp/STRpsych_geno05_maf001_hwe1e-10_chr$chr.smiss > tmp/chr${chr}_mind05.txt &
done
wait

cat tmp/chr* | sort | uniq > exclude_mind05_anychr.samples

#--------------------------------------------------------------------#

## Het/hom outliers
plink2 --bfile tmp/STRpsych_geno05_maf001_hwe1e-10 \
    --het \
    --out tmp/STRpsych_geno05_maf001_hwe1e-10_het

mv tmp/STRpsych_geno05_maf001_hwe1e-10_het.het hethom_report.txt 
awk -F"\t" '$6>0.06 || $6<-0.06{print $1,$2}' hethom_report.txt > exclude_het_hom_outlier.samples

# Merge individuals failing sex check and per chr missingness and drop them
cat exclude_sex_mismatch.samples exclude_mind05_anychr.samples exclude_het_hom_outlier.samples | sort | uniq > tmp/exclude.samples
plink2 --bfile tmp/STRpsych_geno05_maf001_hwe1e-10 \
    --remove tmp/exclude.samples \
    --make-bed \
    --out STRpsych_QC

################################################################################################

sh $PGI_Repo/code/7_Genotypes/7.0_HRCimputation_prep.sh $gf_dir_genotyped_STR/psych/preImputationQC/STRpsych_QC TOPmed | tee -a $PGI_Repo/code/7_Genotypes/7.0_HRCimputation_prep_STRpsych.sh.log

mkdir -p logs preImputationQC_info VCFfinal plink
mv *.vcf.gz VCFfinal
mv STRpsych_QC.bed STRpsych_QC.bim STRpsych_QC.fam STRpsych_QC_frq.afreq plink/
mv tmp/*.log *.log logs/
mv *.ref *.nonSnp *.mono *.geno *.dup *.af *.txt *.snps *.samples Run-* preImputationQC_info/
rm -r tmp *.vcf *updated* *.bed *.bim *.fam