source $PGI_Repo/code/paths
source $PGI_Repo/code/7_Genotypes/7.6.0_formatConversion.sh

gf_out=$gf_dir_MIDUS

# Create merged info files for PC construction later
for chr in {1..22}; do
    infoChrOmni11=$(echo "$info_orig_MIDUSomni11" | sed "s/\[1:22\]/$chr/g")
    infoChrOmni10=$(echo "$info_orig_MIDUSomni10" | sed "s/\[1:22\]/$chr/g")

    gunzip -c "$infoChrOmni10" > ${gf_out}/info/omni10_chr$chr.info
    gunzip -c "$infoChrOmni11" > ${gf_out}/info/omni11_chr$chr.info

    awk 'NR==FNR{info[$1]=$7;next}
        FNR==1{print "SNPID","info";next}
        {split($1,id,":"); chrposid=id[1]":"id[2]}
        ($1 in info && $7<info[$1]) || !($1 in info) {print chrposid,$7} 
        ($1 in info && $7>info[$1]){print chrposid,info[$1]}' OFS="\t" "${gf_out}/info/omni10_chr$chr.info" "${gf_out}/info/omni11_chr$chr.info" > ${gf_out}/info/MIDUS_chr${chr}_minINFO.txt
done


for subcohort in omni10 omni11; do

    eval info='$'info_orig_MIDUS${subcohort}
    eval gf='$'gf_orig_MIDUS${subcohort}
    eval gfDir='$'gf_dir_MIDUS${subcohort}
        
    for chr in {1..22}; do
        infoChr=$(echo "$info" | sed "s/\[1:22\]/$chr/g")
        gfChr=$(echo "$gf" | sed "s/\[1:22\]/$chr/g")
            
        zcat "$infoChr" | awk '$7>0.6{split($1,a,":"); print a[1],a[2]}' OFS="\t" > ${gf_out}/plink2/${subcohort}_chr${chr}_info60.snps &
    done
    wait

    for chr in {1..22}; do
        gfChr=$(echo "$gf" | sed "s/\[1:22\]/$chr/g")
        $vcftools --gzvcf "$gfChr" \
            --positions ${gf_out}/plink2/${subcohort}_chr${chr}_info60.snps \
            --keep $gfDir/sampleQC/MIDUS${subcohort}_EUR_FID_IID.txt \
            --min-alleles 2 \
            --max-alleles 2 \
            --maf 0.005 \
            --max-missing 0.1 \
            --hwe 1e-5 \
            --freq \
            --out ${gf_out}/plink2/${subcohort}_chr${chr}_filtered &
    done
    wait

    for chr in {1..22}; do
        awk -F"\t" 'NR>1{split($5,a,":"); split($6,b,":")} a[2]<b[2] {minAllele=a[1];MAF=a[2]} a[2]>b[2] {minAllele=b[1];MAF=b[2]}{print $1,$2,minAllele, MAF, $4}' OFS="\t" ${gf_out}/plink2/${subcohort}_chr${chr}_filtered.frq > ${gf_out}/plink2/tmp 
        mv ${gf_out}/plink2/tmp ${gf_out}/plink2/${subcohort}_chr${chr}_filtered.frq 
    done
done
            
for chr in {1..22}; do
    awk 'NR==FNR{id=$1 FS $2; minAllele[id]=$3; maf[id]=$4; n[id]=$5; next} FNR==1{print "CHR","BP","MinAllele","MAF_omni11","N_omni11","MAF_omni10","N_omni10";next}
        {id=$1 FS $2} id in minAllele && $3==minAllele[id] {print $1, $2, $3, $4, $5, maf[id], n[id]}' OFS="\t" ${gf_out}/plink2/omni10_chr${chr}_filtered.frq ${gf_out}/plink2/omni11_chr${chr}_filtered.frq >  ${gf_out}/plink2/tmp_chr${chr}

    awk 'NR==1{print $0, "Z"; next}
        {if ($4!=0 || $6!=0)  {p_pooled=(($4*$5)+($6*$7))/($5+$7);    Z=($4-$6)/sqrt( (p_pooled*(1-p_pooled))/$5 + (p_pooled*(1-p_pooled))/$7 ); } else Z=0}
        {$1=$1; print $0, Z}' OFS="\t" ${gf_out}/plink2/tmp_chr${chr} > ${gf_out}/plink2/MAF_report_chr${chr}.txt
    
    # Get list of SNPs for which the P-value for H0: MAFdiff=0 is less than 0.01
    awk -F"\t" '$8>-2.5758293 && $8<2.5758293 {print $1,$2}' OFS="\t" ${gf_out}/plink2/MAF_report_chr${chr}.txt > ${gf_out}/plink2/chr${chr}.snps
done

for subcohort in omni10 omni11; do
    eval gf='$'gf_orig_MIDUS${subcohort}
    for chr in {1..22}; do
        gfChr=$(echo "$gf" | sed "s/\[1:22\]/$chr/g")
            
        $vcftools --gzvcf "$gfChr" \
            --positions ${gf_out}/plink2/chr${chr}.snps \
            --recode \
            --stdout | bgzip -c > ${gf_out}/plink2/${subcohort}_chr${chr}_filtered.vcf.gz &
    done
    wait

    for chr in {1..22}; do
        $bcftools index ${gf_out}/plink2/${subcohort}_chr${chr}_filtered.vcf.gz &
    done
    wait
done

for chr in {1..22}; do
    $bcftools merge -m none \
        ${gf_out}/plink2/omni10_chr${chr}_filtered.vcf.gz \
        ${gf_out}/plink2/omni11_chr${chr}_filtered.vcf.gz \
        -O z -o ${gf_out}/plink2/MIDUS_chr$chr.vcf.gz 
done

mkdir -p ${gf_out}/plink2/logs 
mv ${gf_out}/plink2/*.log ${gf_out}/plink2/logs 
mv ${gf_out}/plink2/MAF* ${gf_out}/plink2/logs  
rm ${gf_out}/plink2/omni10_* ${gf_out}/plink2/omni11* ${gf_out}/plink2/tmp*  ${gf_out}/plink2/*.snps ${gf_out}/plink2/*.vcf.gz      


vcf2plink2 "${gf_out}/plink2/MIDUS_chr[1:22].vcf.gz" "${gf_out}/plink2/MIDUS" ChrPosIDhg19
subsetSNPs "${gf_out}/plink2/MIDUS_chr[1:22]" $HM3_ChrPosIDhg19 "${gf_out}/plink/HM3/MIDUS_HM3"
mergePlink "${gf_out}/plink/HM3/MIDUS_HM3_chr[1:22]" "${gf_out}/plink/HM3/MIDUS_HM3"
        