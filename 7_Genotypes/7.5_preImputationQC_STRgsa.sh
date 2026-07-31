source $PGI_Repo/code/paths
mkdir -p $gf_dir_genotyped_STR/gsa/preImputationQC/tmp
cd $gf_dir_genotyped_STR/gsa/preImputationQC

## Download support files for Infinium Global Screening Array v1.0 
cd reffiles
#Loci Name to rsID Conversion File (GRCh37)
wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/global-screening-array-24/infinium-global-screening-array-24-v1-0-c1-b150-rsids.zip
unzip infinium-global-screening-array-24-v1-0-c2-b150-rsids.zip 
awk -F"\t" '{split($2,a,","); print $1,a[1]}' OFS="\t" GSA-24v1-0_C1_b150_rsids.txt > ../tmp/GSA-24v1-0_C1_b150_rsids_single.txt

# List of markers removed from manifest A1
wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/global-screening-array-24/infinium-global-screening-array-v1-0-a1-vs-c1-cut-list.zip
unzip infinium-global-screening-array-v1-0-a1-vs-c1-cut-list.zip
cd ..
#--------------------------------------------------------------------#

# Per cohort SNP level QC, exclude SNPs removed by illumina, convert CATSS-GSA to bed format
plink2 --pfile $gf_orig_genotyped_STRcatssGSA --maf 0.0001 --geno 0.05 --hwe midp 1e-10 --exclude reffiles/GSA-24v1-0_A1_vs_C1_CutList.txt --make-bed --out tmp/STRcatssGSA
plink2 --bfile $gf_orig_genotyped_STRyatss --maf 0.0001 --geno 0.05 --hwe midp 1e-10 --exclude reffiles/GSA-24v1-0_A1_vs_C1_CutList.txt --make-bed --out tmp/STRyatss
plink2 --bfile $gf_orig_genotyped_STRstage --maf 0.0001 --geno 0.05 --hwe midp 1e-15 --exclude reffiles/GSA-24v1-0_A1_vs_C1_CutList.txt --make-bed --out tmp/STRstage


# Update SNP names, drop duplicated ID, duplicated BP, get allele freq 
for study in stage catssGSA yatss 
do
	plink1.9 --bfile tmp/STR$study --update-name tmp/GSA-24v1-0_C1_b150_rsids_single.txt --make-bed --out tmp/STR${study}_rsID
	plink2 --bfile tmp/STR${study}_rsID --set-missing-var-ids @:# --make-bed --out tmp/STR${study}_rsID_noMissID
	plink2 --bfile tmp/STR${study}_rsID_noMissID --rm-dup list force-first --make-bed --out tmp/STR${study}_rsID_noMissID_noDupID
	mv tmp/STR${study}_rsID_noMissID_noDupID.rmdup.list exclude_dupID_STR${study}.snps	
	plink1.9 --bfile tmp/STR${study}_rsID_noMissID_noDupID --list-duplicate-vars suppress-first ids-only --out tmp/STR${study}_rsID_noMissID_noDupID_BPdups
	sed 's/ /\n/g' tmp/STR${study}_rsID_noMissID_noDupID_BPdups.dupvar > exclude_dupBP_STR${study}.snps
	plink2 --bfile tmp/STR${study}_rsID_noMissID_noDupID --exclude exclude_dupBP_STR${study}.snps --make-bed --out tmp/STR${study}_rsID_noMissID_noDupID_noDupBP
	plink1.9 --bfile tmp/STR${study}_rsID_noMissID_noDupID_noDupBP --freq --out tmp/STR${study}_MAF
done

# Get Z scores for differences in MAF across cohorts
awk 'NR==FNR{maf[$2]=$5; n[$2]=$6; next} FNR==1{print "SNP","MAF_yatss","N_yatss","MAF_stage","N_stage";next}
    ($2 in maf) {print $2, $5, $6, maf[$2], n[$2]}' tmp/STRstage_MAF.frq tmp/STRyatss_MAF.frq > tmp/MAF_report.txt
awk 'NR==FNR{maf[$2]=$5; n[$2]=$6; next} FNR==1{print "SNP","MAF_yatss","N_yatss","MAF_stage","N_stage","MAF_catss","N_catss";next}
    ($1 in maf) {print $1, $2, $3, $4, $5, maf[$1], n[$1]}' tmp/STRcatssGSA_MAF.frq tmp/MAF_report.txt > tmp/tmp
awk 'NR==1{print $0, "Zsy", "Zsc", "Zyc"; next}
   	{if ($2!=0 || $4!=0) Zys=($2-$4)/sqrt(($2*(1-$2)/$3)+($4*(1-$4)/$5)); else Zys=0}
	{if ($2!=0 || $6!=0) Zyc=($2-$6)/sqrt(($2*(1-$2)/$3)+($6*(1-$6)/$7)); else Zyc=0}
	{if ($4!=0 || $6!=0) Zsc=($4-$6)/sqrt(($4*(1-$4)/$5)+($6*(1-$6)/$7)); else Zsc=0}
	{$1=$1; print $0, Zys, Zyc, Zsc}' OFS="\t" tmp/tmp > MAF_report.txt

# Get list of SNPs for which the P-value for H0: MAFdiff=0 is less than 0.001
awk -F"\t" '$8<-3.29 || $8>3.29 || $9<-3.29 || $9>3.29 || $10<-3.29 || $10>3.29 {print $1}' MAF_report.txt > exclude_MAFdiff_p001.snps

#--------------------------------------------------------------------#

## Merge YATSS, STAGE and CATSS-GSA

# Create list of files to merge
echo "tmp/STRstage_rsID_noMissID_noDupID_noDupBP
tmp/STRcatssGSA_rsID_noMissID_noDupID_noDupBP
tmp/STRyatss_rsID_noMissID_noDupID_noDupBP" > tmp/mergelist

# Attempt merge and get list of problematic SNPs
plink1.9 --merge-list tmp/mergelist \
    --make-bed \
    --out tmp/STRgsa
mv tmp/STRgsa-merge.missnp exclude_merge_missnp.snps

# Add merge-problematic SNPs to SNPs with sign MAF difference
cat exclude_MAFdiff_p001.snps exclude_merge_missnp.snps | sort | uniq > tmp/exclude_MAFdiff_p001_merge_missnp.snps

# Exclude problematic SNPs from each study
for study in catssGSA yatss stage
do
    plink1.9 --bfile tmp/STR${study}_rsID_noMissID_noDupID_noDupBP --exclude tmp/exclude_MAFdiff_p001_merge_missnp.snps --make-bed --out tmp/STR$study.nomissnp
done

# Create list of files w/o problematic SNPs to merge
echo "tmp/STRcatssGSA.nomissnp
tmp/STRyatss.nomissnp
tmp/STRstage.nomissnp" > tmp/mergelist

# Merge
plink1.9 --merge-list tmp/mergelist \
    --make-bed \
    --out tmp/STRgsa

#--------------------------------------------------------------------#

## Variant QC: call rate> 0.95 , maf>0.001, hwe 1e-10
plink2 --bfile tmp/STRgsa \
    --geno 0.05 \
    --maf 0.001 \
    --hwe midp 1e-10 \
    --make-bed \
    --out tmp/STRgsa_geno05_maf001_hwe1e-10

#--------------------------------------------------------------------#

## Sex check
plink1.9 --bfile tmp/STRgsa_geno05_maf001_hwe1e-10 \
    --chr X \
    --maf 0.01 \
    --indep-pairwise 1000kb 5 0.5 \
    --out tmp/STRgsa_chrX_geno05_maf01_hwe1e-10_pruned

plink1.9 --bfile tmp/STRgsa_geno05_maf001_hwe1e-10 \
    --extract tmp/STRgsa_chrX_geno05_maf01_hwe1e-10_pruned.prune.in \
    --check-sex \
    --out tmp/STRgsa_chrX_geno05_maf01_hwe1e-10_pruned_sexcheck

awk '$5!="OK"{print $1,$2}' tmp/STRgsa_chrX_geno05_maf01_hwe1e-10_pruned_sexcheck.sexcheck > exclude_sex_mismatch.samples

#--------------------------------------------------------------------#

## Get individuals with per chromosome subject-level missingness rate > 0.05
for chr in {1..22}
do
    plink2 -bfile tmp/STRgsa_geno05_maf001_hwe1e-10 \
        --chr $chr \
        --missing \
        --out tmp/STRgsa_geno05_maf001_hwe1e-10_chr$chr &
done
wait

for chr in {1..22}
do
    awk -F"\t" 'NR>1 && $5>0.05{print $1,$2}' tmp/STRgsa_geno05_maf001_hwe1e-10_chr$chr.smiss > tmp/chr${chr}_mind05.txt &
done
wait

cat tmp/chr* | sort | uniq > exclude_mind05_anychr.samples

#--------------------------------------------------------------------#

## Het/hom outliers
plink2 --bfile tmp/STRgsa_geno05_maf001_hwe1e-10 \
    --het \
    --out tmp/STRgsa_geno05_maf001_hwe1e-10_het

mv tmp/STRgsa_geno05_maf001_hwe1e-10_het.het hethom_report.txt 
awk -F"\t" '$6>0.13 || $6<-0.13{print $1,$2}' hethom_report.txt > exclude_het_hom_outlier.samples

# Merge individuals failing sex check and per chr missingness and drop them
cat exclude_sex_mismatch.samples exclude_mind05_anychr.samples exclude_het_hom_outlier.samples | sort | uniq > tmp/exclude.samples
plink2 --bfile tmp/STRgsa_geno05_maf001_hwe1e-10 \
    --remove tmp/exclude.samples \
    --make-bed \
    --out STRgsa_QC



################################################################################################

sh $PGI_Repo/code/7_Genotypes/7.0_HRCimputation_prep.sh $gf_dir_genotyped_STR/gsa/preImputationQC/STRgsa_QC TOPmed | tee -a $PGI_Repo/code/7_Genotypes/7.0_HRCimputation_prep_STRgsa.sh.log

mkdir -p logs preImputationQC_info VCFfinal plink
mv *.vcf.gz VCFfinal
mv STRgsa_QC.bed STRgsa_QC.bim STRgsa_QC.fam STRgsa_QC_frq.afreq plink/
mv tmp/*.log *.log logs/
mv *.ref *.nonSnp *.mono *.geno *.dup *.af *.txt *.snps *.samples Run-* preImputationQC_info/
rm -r tmp *.vcf *updated* *.bed *.bim *.fam