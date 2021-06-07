
#!/bin/bash

source paths1

cd $mainDir/derived_data/1_UKB_GWAS

###### GET LIST OF BRAIN MRI SAMPLE AND RELATED INDIVIDUALS FOR PARTITIONING ########
echo "Getting list of brain-imaged and up to 3rd degree related individuals.." 

### UKB brain imaging sample (file "k11425.f12188.R.tab" based on data field 12188)
awk -F"\t" 'NR>1 && !($2==0 && $3==0){print $1,"1"}' OFS="\t" $p1_UKBbrain > tmp/brain_sample

### UKB relatedness data - Contains relatives up to third degree
# Remove withdrawn individuals (negative ID), collapse, reverse sort on ID and kinship  
awk 'NR>1 && $1>0 && $2>0{print $1,$5;print $2,$5}' OFS="\t" $p1_UKBrelatedness | sort -k1,5 -r > tmp/tmp

# Remove duplicate ID's, keeping the one with larger kinship, then reverse sort on kinship
awk -F"\t" '!_[$1]++' OFS="\t" tmp/tmp | sort -k2 -g -r > tmp/UKBv2_rel_uniq_sorted

# Get list of individuals where first the brain sample is listed, then in reverse kinship order 
awk -F"\t" 'BEGIN{OFS="\t";print "n_eid","partition_order"}!_[$1]++{print $1,NR}' tmp/brain_sample tmp/UKBv2_rel_uniq_sorted > tmp/IDs_assignPartition_ordered.txt

echo "List created."
###########################################################################

echo "Running Stata to partition sample and obtain phenotype file for each partition.."

stata -b do $mainDir/code/1_UKB_GWAS/1.1_prep_UKB_phenos.do \
	$mainDir/derived_data/1_UKB_GWAS \
	$p1_UKBcrosswalk \
	$mainDir/derived_data/1_UKB_GWAS/tmp/IDs_assignPartition_ordered.txt \
	$p1_pheno_data_1 \
	$p1_pheno_data_2 \
	$p1_pheno_data_3 \
	$p1_covar_data \
	$mainDir/code/1_UKB_GWAS/1.1_prep_UKB_phenos_test.do.log

# Replace missing values (.) with NA
for file in input/*.pheno; do
	sed -i 's/ $/ NA/g' $file &
done 

wait

echo "Phenotype files created."

