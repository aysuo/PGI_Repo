source $PGI_Repo/code/paths
source $PGI_Repo/code/7_Genotypes/7.6.0_formatConversion.sh

gf_out=$gf_dir_NCDS

awk -F"," '$2~"Infinium_HumanHap_550K_v3"{gsub(/"/,"",$1);print $1"_"$1}' $array_indicator_NCDS > ${gf_out}/plink2/Infinium_500kv3.id
awk -F"," '$2~"Infinium_HumanHap_550K_v1.1"{gsub(/"/,"",$1);print $1"_"$1}' $array_indicator_NCDS > ${gf_out}/plink2/Infinium_500kv11.id
awk -F"," '$2~"Illumina_Human_660Quad"{gsub(/"/,"",$1);print $1"_"$1}' $array_indicator_NCDS > ${gf_out}/plink2/Illumina_660quad.id
awk -F"," '$2~"Illumina_1.2M"{gsub(/"/,"",$1);print $1"_"$1}' $array_indicator_NCDS > ${gf_out}/plink2/Illumina_12.id
awk -F"," '$2~"Affymetrix_v6"{gsub(/"/,"",$1);print $1"_"$1}' $array_indicator_NCDS > ${gf_out}/plink2/affy_v6.id

for subcohort in Infinium_500kv3 Infinium_500kv11 Illumina_660quad Illumina_12 affy_v6; do

    eval info='$'info_orig_NCDS_${subcohort}
    eval gf='$'gf_orig_NCDS_${subcohort}
        
    for chr in {1..22}; do
        infoChr=$(echo "$info" | sed "s/\[1:22\]/$chr/g")
        gfChr=$(echo "$gf" | sed "s/\[1:22\]/$chr/g")
            
        zcat $infoChr | awk -F"\t" 'BEGIN{print "ChrPosID","rsID","REF","ALT","MAF","AF","R2"}
                                        !($1~"#"){ 
                                        gsub(/R2=/,"",$8);
                                        gsub(/MAF=/,"",$8);
                                        gsub(/AF=/,"",$8); 
                                        split($8,a,";"); 
                                        if (a[1]=="IMPUTED" || (a[1]=="TYPED" && a[2]!="IMPUTED")) {print $1":"$2,$3,$4,$5,a[3],a[2],a[5];} else {print $1":"$2,$3,$4,$5,a[4],a[3],a[6];}} ' OFS="\t" > ${gf_out}/info/${subcohort}_chr${chr}.info &
    done
    wait
done

#############################################
# Get min info for making PCs later
for chr in {1..22}; do
    awk -F"\t" 'NR==FNR{info[$1]=$7;next} 
        ($1 in info && $7<info[$1]) || !($1 in info) {print $1,$7} 
        ($1 in info && $7>info[$1]){print $1,info[$1]}' OFS="\t" ${gf_out}/info/Infinium_500kv3_chr${chr}.info ${gf_out}/info/Infinium_500kv11_chr${chr}.info > ${gf_out}/info/tmp
    awk -F"\t" 'NR==FNR{info[$1]=$2;next} 
        ($1 in info && $7<info[$1]) || !($1 in info) {print $1,$7} 
        ($1 in info && $7>info[$1]){print $1,info[$1]}' OFS="\t" ${gf_out}/info/tmp ${gf_out}/info/Illumina_660quad_chr${chr}.info > ${gf_out}/info/NCDS_minINFO_chr${chr}.txt
    awk -F"\t" 'NR==FNR{info[$1]=$2;next} 
        ($1 in info && $7<info[$1]) || !($1 in info) {print $1,$7} 
        ($1 in info && $7>info[$1]){print $1,info[$1]}' OFS="\t" ${gf_out}/info/NCDS_minINFO_chr${chr}.txt ${gf_out}/info/Illumina_12_chr${chr}.info > ${gf_out}/info/tmp
    awk -F"\t" 'NR==FNR{info[$1]=$2;next} 
        ($1 in info && $7<info[$1]) || !($1 in info) {print $1,$7} 
        ($1 in info && $7>info[$1]){print $1,info[$1]}' OFS="\t" ${gf_out}/info/tmp ${gf_out}/info/affy_v6_chr${chr}.info > ${gf_out}/info/NCDS_minINFO_chr${chr}.txt
done
#############################################


for subcohort in Infinium_500kv3 Infinium_500kv11 Illumina_660quad Illumina_12 affy_v6; do
    echo PROCESSING: "$subcohort"
    eval info='$'info_orig_NCDS_${subcohort}
    eval gf='$'gf_orig_NCDS_${subcohort}
    
    # echo "Obtaining list of SNPs with R2>0.6 and MAF>0.005"
    for chr in {1..22}; do
        awk '$5>0.005 && $7>0.6{split($1,a,":");print "chr"a[1],a[2]}' OFS="\t" ${gf_out}/info/${subcohort}_chr${chr}.info > ${gf_out}/plink2/${subcohort}_chr${chr}_maf005_info60.snps &
    done
    wait

    echo "Obtaining frequencies for biallelic SNPs with MAF>0.5%, HWE P>1e-6, per chr call rate > 0.9, R2>0.6.."
    for chr in {1..22}; do
        gfChr=$(echo "$gf" | sed "s/\[1:22\]/$chr/g")
        $vcftools --gzvcf "$gfChr" \
            --keep ${gf_out}/plink2/${subcohort}.id \
            --positions ${gf_out}/plink2/${subcohort}_chr${chr}_maf005_info60.snps \
            --min-alleles 2 \
            --max-alleles 2 \
            --max-missing 0.1 \
            --hwe 1e-6 \
            --freq \
            --out ${gf_out}/plink2/${subcohort}_chr${chr}_filtered &
            
    done
    wait

    echo "Re-formatting frequency files.."
    for chr in {1..22}; do
        awk -F"\t" 'NR>1{split($5,a,":"); split($6,b,":")} a[2]<b[2] {minAllele=a[1];MAF=a[2]} a[2]>b[2] {minAllele=b[1];MAF=b[2]}{print $1,$2,minAllele, MAF, $4}' OFS="\t" ${gf_out}/plink2/${subcohort}_chr${chr}_filtered.frq >  mv ${gf_out}/plink2/tmp
        mv ${gf_out}/plink2/tmp ${gf_out}/plink2/${subcohort}_chr${chr}_filtered.frq 
    done
done

echo "Comparing MAF between each subcohort and Infinium_500kv3 (largest subcohort), filtering out SNPs with MAF significantly different at 1%.."
for chr in {1..22}; do
    echo "Chromosome $chr"
    awk 'NR==FNR{id=$1 FS $2; minAllele[id]=$3; maf[id]=$4; n[id]=$5; next} FNR==1{print "CHR","BP","MinAllele","MAF_Illumina12","N_Illumina12","MAF_Infinium_500kv3","N_Infinium_500kv3";next}
        {id=$1 FS $2} id in minAllele && $3==minAllele[id] {print $1, $2, $3, $4, $5, maf[id], n[id]}' OFS="\t" ${gf_out}/plink2/Infinium_500kv3_chr${chr}_filtered.frq ${gf_out}/plink2/Illumina_12_chr${chr}_filtered.frq > ${gf_out}/plink2/MAF_report_chr${chr}.txt
    awk 'NR==FNR{id=$1 FS $2; minAllele[id]=$3; maf[id]=$4; n[id]=$5; next} FNR==1{print $0,"MAF_Illumina_660quad","N_Illumina_660quad";next}
        {id=$1 FS $2} id in minAllele && $3==minAllele[id] {print $1, $2, $3, $4, $5, $6, $7, maf[id], n[id]}' OFS="\t" ${gf_out}/plink2/Illumina_660quad_chr${chr}_filtered.frq ${gf_out}/plink2/MAF_report_chr${chr}.txt > ${gf_out}/plink2/tmp_chr${chr}
    awk 'NR==FNR{id=$1 FS $2; minAllele[id]=$3; maf[id]=$4; n[id]=$5; next} FNR==1{print $0,"MAF_Infinium_500kv11","N_Infinium_500kv11";next}
        {id=$1 FS $2} id in minAllele && $3==minAllele[id] {print $1, $2, $3, $4, $5, $6, $7, $8, $9, maf[id], n[id]}' OFS="\t" ${gf_out}/plink2/Infinium_500kv11_chr${chr}_filtered.frq ${gf_out}/plink2/tmp_chr${chr} > ${gf_out}/plink2/MAF_report_chr${chr}.txt
    awk 'NR==FNR{id=$1 FS $2; minAllele[id]=$3; maf[id]=$4; n[id]=$5; next} FNR==1{print $0,"MAF_affy_v6","N_affy_v6";next}
        {id=$1 FS $2} id in minAllele && $3==minAllele[id] {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, maf[id], n[id]}' OFS="\t" ${gf_out}/plink2/affy_v6_chr${chr}_filtered.frq ${gf_out}/plink2/MAF_report_chr${chr}.txt > ${gf_out}/plink2/tmp_chr${chr}

    awk 'NR==1{print $0, "Z1", "Z2", "Z3", "Z4"; next}
        {if ($4!=0 || $6!=0)  {p_pooled=(($4*$5)+($6*$7))/($5+$7);    Z1=($4-$6)/sqrt( (p_pooled*(1-p_pooled))/$5 + (p_pooled*(1-p_pooled))/$7 ); } else Z1=0}
        {if ($4!=0 || $8!=0)  {p_pooled=(($4*$5)+($8*$9))/($5+$9);    Z2=($4-$8)/sqrt( (p_pooled*(1-p_pooled))/$5 + (p_pooled*(1-p_pooled))/$9 ); } else Z2=0}
        {if ($4!=0 || $10!=0) {p_pooled=(($4*$5)+($10*$11))/($5+$11); Z3=($4-$10)/sqrt( (p_pooled*(1-p_pooled))/$5 + (p_pooled*(1-p_pooled))/$11 ); } else Z3=0}
        {if ($4!=0 || $12!=0) {p_pooled=(($4*$5)+($12*$13))/($5+$13); Z4=($4-$12)/sqrt( (p_pooled*(1-p_pooled))/$5 + (p_pooled*(1-p_pooled))/$13);} else Z4=0}
        {$1=$1; print $0, Z1, Z2, Z3, Z4}' OFS="\t" ${gf_out}/plink2/tmp_chr${chr} > ${gf_out}/plink2/MAF_report_chr${chr}.txt
    
    # Get list of SNPs for which the P-value for H0: MAFdiff=0 is less than 0.01
    awk -F"\t" '$14>-2.5758293 && $14<2.5758293 && $15>-2.5758293 && $15<2.5758293 && $16>-2.5758293 && $16<2.5758293 && $17>-2.5758293 && $17<2.5758293 {print $1,$2}' OFS="\t" ${gf_out}/plink2/MAF_report_chr${chr}.txt > ${gf_out}/plink2/chr${chr}.snps
done


for subcohort in Infinium_500kv3 Infinium_500kv11 Illumina_660quad Illumina_12 affy_v6; do
    eval gf='$'gf_orig_NCDS_${subcohort}
    echo "Extracting SNPs from vcf files after all filters for $subcohort and indexing.."
    for chr in {1..22}; do
        gfChr=$(echo "$gf" | sed "s/\[1:22\]/$chr/g")
            
        $vcftools --gzvcf "$gfChr" \
            --keep ${gf_out}/plink2/${subcohort}.id \
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

echo "Merging subcohorts.."
for chr in {1..22}; do
    $bcftools merge -m none \
        ${gf_out}/plink2/Infinium_500kv3_chr${chr}_filtered.vcf.gz \
        ${gf_out}/plink2/Infinium_500kv11_chr${chr}_filtered.vcf.gz \
        ${gf_out}/plink2/Illumina_660quad_chr${chr}_filtered.vcf.gz \
        ${gf_out}/plink2/Illumina_12_chr${chr}_filtered.vcf.gz \
        ${gf_out}/plink2/affy_v6_chr${chr}_filtered.vcf.gz \
        -O z -o ${gf_out}/plink2/NCDS_chr$chr.vcf.gz 
    echo "Chromosome $chr done."  
done

echo "Cleaning up.."
mkdir -p ${gf_out}/plink2/logs 
mv ${gf_out}/plink2/*.log ${gf_out}/plink2/logs
mv ${gf_out}/plink2/MAF* ${gf_out}/plink2/logs
rm ${gf_out}/plink2/*tmp* ${gf_out}/plink2/*.vcf.gz ${gf_out}/plink2/*.snps ${gf_out}/plink2/Illumina_* ${gf_out}/plink2/affy* ${gf_out}/plink2/Infinium*     

echo "Converting to plink2.."
vcf2plink2 "${gf_out}/plink2/NCDS_chr[1:22].vcf.gz" ${gf_out}/plink2/NCDS ChrPosIDhg38

echo "Extracting HM3 SNPs.."
subsetSNPs "${gf_out}/plink2/NCDS_chr[1:22]" $HM3_ChrPosIDhg38 "${gf_out}/plink/HM3/NCDS_HM3"

echo "Merging HM3 per chromosome files.."
mergePlink "${gf_out}/plink/HM3/NCDS_HM3_chr[1:22]" "${gf_out}/plink/HM3/NCDS_HM3"

echo "DONE."
        

