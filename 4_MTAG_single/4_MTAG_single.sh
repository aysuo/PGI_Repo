#!/bin/bash

source $PGI_Repo/code/4_MTAG_single/4.0_MTAG_single_functions.sh

cd $PGI_Repo/derived_data/4_MTAG_single

# Correct reverse-coded sumstats
# for study in CHILDLESS-Mathieson SWB-23andMe Risk-23andMe SMCESS-Liu SMCESS-LiuSansUKB SMCESS-SaundersEUR SMCESS-SaundersEURSansUKB; do
#   path=$(cut -f2 $PGI_Repo/code/4_MTAG_single/singleMTAG_input_filelist.txt | awk -F"," -v study=$study '{for(i=1;i<=NF;i++) if ($i~study) print $i}'  | sort | uniq)
#   eval path=$path
#   # Unzip and rename original file as *_revcoded
#   unzipped=$(echo $path | sed 's/\.gz//g')
#   if [[ $path == *.gz ]]
#   then
#     gunzip $path
#   fi 
#   mv $unzipped ${unzipped}_revcoded
#   # Reverse sign of effect 
#   awk -F"\t" 'NR==1{print}NR>1{$8=-$8;print}' OFS="\t" ${unzipped}_revcoded > ${unzipped}
#   gzip ${unzipped}
# done


# Check which phenotypes aren't done
checkStatusMTAG $PGI_Repo/code/4_MTAG_single/singleMTAG_input_filelist.txt

# Run MTAG for unfinished phenotypes, 20 at a time
if [[ $status == 1 ]]
then
  i=0
  while read row; do
    pheno=$(echo $row | cut -d" " -f1)
    sumstats=$(echo $row | cut -d" " -f2)
    eval sumstats=$sumstats
    mkdir -p $pheno

    phenos_force=("ANOREX3" "ECZEMA1" "ECZEMA2" "ECZEMA4" "SCZ1" "SCZ2" "SCZ3" "SCZ4" "SCZ5")
    phenos_noOverlap=("BL_CHOL3" "BL_CHOL4" "BL_CHOL5" "BL_HDL3" "BL_HDL4" "BL_HDL5" "BL_LDL3" "BL_LDL4" "BL_LDL5" "BL_nonHDL2" "BL_nonHDL3" "BL_nonHDL4" "BL_nonHDL5" "BMI1" "BMI2" "BMI3" "BMI4" "BMI5" "EA1" "EA2" "EA3" "EA4" "EA5" "EA6" "EA7" "EA8" "EA9" "EA10" "EA14" "EA15" "HEIGHT1" "HEIGHT2" "HEIGHT3" "HEIGHT4" "HEIGHT5" "HEIGHT6")
    phenos_metaformat=("ALZ2" "ALZ3" "ALZ4" "ALZ5") #"ALZ1"
    phenos_cptid=("PRCA1" "PRCA2" "PRCA3" "PRCA4" "ALZ1" "ALZ3" "ALZ5")

    if [[ ${phenos_force[@]} =~ $pheno  ]]
    then
      force=1
    else
      force=0
    fi

    if [[ ${phenos_noOverlap[@]} =~ $pheno  ]]
    then
      nooverlap=1
    else
      nooverlap=0
    fi
    
    if [[ ${phenos_metaformat[@]} =~ $pheno  ]]
    then
      metaformat=1
    else
      metaformat=0
    fi

    if [[ ${phenos_cptid[@]} =~ $pheno  ]]
    then
      cptid=1
    else
      cptid=0
    fi

    MTAG_single $pheno $sumstats $pheno/$pheno $force $nooverlap $metaformat $cptid & 
      
    let i+=1
    
    if [[ $i == 20 ]]
    then
      wait
      i=0
    fi
  done < $PGI_Repo/code/4_MTAG_single/singleMTAG_input_filelist.txt.rerun
  wait
fi

# Check if everything ran successfully
checkStatusMTAG $PGI_Repo/code/4_MTAG_single/singleMTAG_input_filelist.txt

# # Check if SNPs have been dropped due to P=0 in MTAG, should be 0.
# rm $PGI_Repo/derived_data/4_MTAG_single/SNPs_dropped_Pval0_input.txt
# for pheno in *
# do 
#   echo -n -e "$pheno\t" >> $PGI_Repo/derived_data/4_MTAG_single/SNPs_dropped_Pval0_input.txt
#   SNPsdropped=$(grep "SNPs with out-of-bounds p-values." $PGI_Repo/derived_data/4_MTAG_single/$pheno/$pheno.log | cut -d" " -f2)
#   echo $SNPsdropped >> $PGI_Repo/derived_data/4_MTAG_single/SNPs_dropped_Pval0_input.txt
# done

# # Check how many SNPs with P=0 there are in the MTAG output
# for pheno in *
# do 
#   for file in $pheno/${pheno}_trait*_SBayesR.txt
#   do
#     N0=$(awk 'NR==1{x=0;next}$7==0{x++}END{print x}' $file)
#     echo -e "$pheno\t$file\t$N0" >> NSNPs_Pval0_in_output.txt
#   done
# done
