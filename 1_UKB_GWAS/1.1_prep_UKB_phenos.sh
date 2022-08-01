
#!/bin/bash

source $PGI_Repo/code/paths

cd $PGI_Repo/derived_data/1_UKB_GWAS

###### GET LIST OF BRAIN MRI SAMPLE AND RELATED INDIVIDUALS FOR PARTITIONING ########
echo "Getting list of brain-imaged and up to 3rd degree related individuals.." 

### UKB brain imaging sample (file "k11425.f12188.R.tab" based on data field 12188)
awk -F"\t" 'NR>1 && !($2==0 && $3==0){print $1,"1"}' OFS="\t" $UKB_brain_sample > tmp/brain_sample

### UKB relatedness data - Contains relatives up to third degree
# Remove withdrawn individuals (negative ID), collapse, reverse sort on ID and kinship  
awk 'NR>1 && $1>0 && $2>0{print $1,$5;print $2,$5}' OFS="\t" $UKB_relatedness | sort -k1,5 -r > tmp/tmp

# Remove duplicate ID's, keeping the one with larger kinship, then reverse sort on kinship
awk -F"\t" '!_[$1]++' OFS="\t" tmp/tmp | sort -k2 -g -r > tmp/UKBv2_rel_uniq_sorted

# Get list of individuals where first the brain sample is listed, then in reverse kinship order 
awk -F"\t" 'BEGIN{OFS="\t";print "n_eid","partition_order"}!_[$1]++{print $1,NR}' tmp/brain_sample tmp/UKBv2_rel_uniq_sorted > tmp/IDs_assignPartition_ordered.txt

echo "List created."
###########################################################################

echo "Running Stata to partition sample and obtain phenotype file for each partition.."

stata -b do $PGI_Repo/code/1_UKB_GWAS/1.1_prep_UKB_phenos.do \
	$PGI_Repo/derived_data/1_UKB_GWAS \
	$UKBv3_crosswalk \
	$PGI_Repo/derived_data/1_UKB_GWAS/tmp/IDs_assignPartition_ordered.txt \
	$UKB_pheno_data_1 \
	$UKB_pheno_data_2 \
	$UKB_pheno_data_3 \
	$UKB_pheno_data_4 \
	$UKB_covar_data \
	$PGI_Repo/code/1_UKB_GWAS/1.1_prep_UKB_phenos.do.log \
	$UKB_withdrawn

# Replace missing values (.) with NA
for file in input/*.pheno; do
	sed -i 's/ $/ NA/g' $file &
done 

wait

echo "Phenotype files created."
