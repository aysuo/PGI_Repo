#!/bin/bash
source $PGI_Repo/code/paths

usage() {
echo "Usage: 
      format_23andMe.sh 
       --annot_all <path to all SNP annotation file>
       --annot_im <path to imputed SNP annotation file>
       --annot_gt <path to genotyped SNP annotation file>
       --annot_version <version of the annotation files, e.g. 7.0>
       --gwas <path to gwas results>
       --out <output file name>"
echo ""
echo "Note: order of options is not important"
}

echo ""
echo "-------------------------------------------------------------"
echo ""

echo -n "Script started on "
date
start=$(date +%s)

#######################################################################
############################ PARSE ARGUMENTS ##########################
#######################################################################

ARGUMENT_LIST=(
    "annot_all"
    "annot_im"
    "annot_gt"
    "annot_version"
    "gwas"
    "out"
)


# Read arguments
opts=$(getopt \
    --longoptions "$(printf "%s:," "${ARGUMENT_LIST[@]}")" \
    --name "$(basename "$0")" \
    --options "" \
    -- "$@"
)

eval set --$opts

# Assign arguments to variables
while [[ $# -gt 0 ]]; do
    case "$1" in
        --annot_all)
            annot_all=$2
            shift 2
            ;;
        --annot_im)
            annot_im=$2
            shift 2
            ;;
        --annot_gt)
            annot_gt=$2
            shift 2
            ;;
        --annot_version)
            annot_version=$2
            shift 2
            ;;
        --gwas)
            gwas=$2
            shift 2
            ;;
        --out)
            out=$2
            shift 2
            ;;
        *)
            usage
            break
            ;;
    esac
done


#######################################################################

echo ""
echo "-------------------------------------------------------------"
echo ""
echo "Checking required arguments.."
echo ""

# Check if required arguments have been supplied"

if [[ -f $annot_all ]]
then
  echo "All SNP annotation file: $annot_all"
elif [[ -d $annot_all ]]
  then
    echo "Error: The path you have specified for the all SNP annotation file is a directory."
      exit 1
    else
      echo "Error: Please specify a valid path for the all SNP annotation file."
    exit 1
fi


if [[ -f $annot_gt ]]
then
  echo "Genotyped SNP annotation file: $annot_gt"
elif [[ -d $annot_gt ]]
  then
    echo "Error: The path you have specified for the genotyped SNP annotation file is a directory."
      exit 1
    else
      echo "Error: Please specify a valid path for the genotyped SNP annotation file."
      exit 1
fi


if [[ -f $annot_im ]]
then
  echo "Imputed SNP annotation file: $annot_im"
elif [[ -d $annot_im ]]
  then
    echo "Error: The path you have specified for the imputed SNP annotation file is a directory."
    exit 1
    else
      echo "Error: Please specify a valid path for the imputed SNP annotation file."
      exit 1
fi

if [[ -z $annot_version ]]
then
  echo "Error: You have not specified the version of annotation files."
  exit 1
else
  echo "Version of annotation files is ${annot_version}."
fi

  
if [[ -z $gwas ]]
then
  echo "Error: No GWAS file(s) have been suppied." 
  exit 1
else
  GWAS=$(echo $gwas | sed 's/,/ /g')
  declare -a GWAS=$(echo "("$GWAS")")
  
  i=1
  for j in ${GWAS[*]}
  do
    if [[ -f $j ]]
    then
      echo "GWAS file $i: $j"
      let i+=1
    elif [[ -d $j ]]
      then
        echo "Error: The path you have specified for GWAS file $i is a directory."
        exit 1
      else
        echo "Error: Please specify a valid path for the GWAS file."
        exit 1
    fi
  done  
fi

if [[ -z $out ]]
then
  echo "Error: No output name(s) have been supplied. Please specify an output name for each GWAS."
    exit 1
else
  OUT=$(echo $out | sed 's/,/ /g')
  declare -a OUT=$(echo "("$OUT")")
  
  i=1
  for j in ${OUT[*]}
  do
    if [[ -s $j ]]; then
      echo -n "Output file $j already exists. Do you want to rewrite? : "
      read $answer
      
      if [[ $answer="Y" ]]; then
        echo "Output name $i: $j"
        let i+=1
      else
        echo -n "Please specify a new output name for GWAS $i: "
        read ${new_path}
        OUT[$((i-1))]=${new_path}  
        echo "Output name $i: $j"
        let i+=1 
      fi
    elif [[ -d $j ]]; then
      echo "Error: The output name you have specified for GWAS file $i is a directory. Please specify a file name."
      exit 1
    else
      echo "Output name $i: $j"
      let i+=1
    fi
  done
fi



echo ""
echo "-------------------------------------------------------------"
echo ""

#######################################################################
##################### CREATE WORKING DIRECTORIES ######################
#######################################################################

echo "Creating working directories."
echo ""

mkdir TMP || echo "Using the existing TMP directory."
#mkdir OUTPUT || echo "Using the existing OUTPUT directory."
mkdir LOG || echo "Using the existing LOG directory."

echo ""
echo "-------------------------------------------------------------"
echo ""

#######################################################################
########################## DEFINE FUNCTIONS ###########################
#######################################################################

## Function to obtain the following columns from GWAS results: "all.data.id","IMPUTED", "PVALUE", "EFFECT", "SE", "EAF", "N"
## N is obtained as sum of genotype counts if the results are for genotyped SNP (src=G), otherwise (src=I) it is the sum of imputed cases and controls.
## EAF is obtained using genotype counts if the results are for genotyped SNP (src=G), otherwise (src=I) it is obtained using average dose in cases and controls.
## EAF is set to NA if N=0
## SNPs that fail 23andMe QC are dropped

gwasCols(){
  gwas_path=$1
  gwas_file=$(echo ${gwas_path} | sed 's,.\+/,,g')
  awk -F"\t" '
    NR==1 {
      for (i=1; i<=NF; i++){ ix[$i] = i }
      { print "all.data.id","IMPUTED", "P", "BETA", "SE", "EAF", "N" }
    }
    NR>1 { 
      if ("im.num.0" in ix && "im.num.1" in ix && $ix["im.num.0"]!="NA") {
        N_im=$ix["im.num.0"]+$ix["im.num.1"]
        if (N_im>0) {
          EAF_im=($ix["im.num.0"]*$ix["dose.b.0"]+$ix["im.num.1"]*$ix["dose.b.1"])/($ix["im.num.0"]+$ix["im.num.1"])
        }
        else {
        EAF_im="NA"
        }
      }
      else {
        N_im="NA"
        EAF_im=$ix["dose.b.0"]
      }
        
      N_gt=$ix["AA.0"]+$ix["AB.0"]+$ix["BB.0"]+$ix["AA.1"]+$ix["AB.1"]+$ix["BB.1"] 
      
      if (N_gt>0) {
        EAF_gt=(2*($ix["BB.0"]+$ix["BB.1"])+$ix["AB.0"]+$ix["AB.1"])/(N_gt*2)
      }
      else {
        EAF_gt="NA"
      }
    } 
    NR>1 && $ix["pass"]=="Y" && $ix["src"]=="G"{ 
      print $ix["all.data.id"] , 0 , $ix["pvalue"] , $ix["effect"] , $ix["stderr"] , EAF_gt, N_gt 
    } 
    NR>1 && $ix["pass"]=="Y" && $ix["src"]=="I"{ 
      print $ix["all.data.id"] , 1 , $ix["pvalue"] , $ix["effect"] , $ix["stderr"] , EAF_im, N_im 
    }' OFS="\t" ${gwas_path} > TMP/tmp_${gwas_file}

}


#############################################################


main(){
  if ! [[ -f TMP/SNPstats_v${annot_version}.txt ]]
  then
    
    if [ "${annot_version}" = "7.0" || "${annot_version}" = "6.1" ]
    then
      im_cols="im.data.id,rsq"
    else
      im_cols="im.data.id,avg.rsqr"
    fi

    # Extract required columns from annotation files
    echo "Extracting necessary columns from annotation files"

    sh $PGI_Repo/code/2_Formatting/2.1.1.1_ExtractCols.sh \
    --file ${annot_all} \
    --cols "all.data.id,gt.data.id,im.data.id,assay.name,scaffold,position,alleles,ploidy" \
    --out TMP/tmp_all_v${annot_version} > LOG/tmp_all_v${annot_version}.log &

    sh $PGI_Repo/code/2_Formatting/2.1.1.1_ExtractCols.sh \
    --file ${annot_im} \
    --cols ${im_cols} \
    --out TMP/tmp_im_v${annot_version} > LOG/tmp_im_v${annot_version}.log &

    sh $PGI_Repo/code/2_Formatting/2.1.1.1_ExtractCols.sh \
    --file ${annot_gt} \
    --cols "gt.data.id,gt.rate,hw.p.value" \
    --out TMP/tmp_gt_v${annot_version} > LOG/tmp_gt_v${annot_version}.log &
    
    wait

    ## Merge annotation files
    echo -e "Merging annotation files"

    sh $PGI_Repo/code/2_Formatting/2.1.1.2_Merge.sh \
    --mergeType L \
    --mergeCol1 3 \
    --mergeCol2 1 \
    --file1 TMP/tmp_all_v${annot_version} \
    --file2 TMP/tmp_im_v${annot_version} \
    --out TMP/tmp_all_im_v${annot_version} > LOG/tmp_all_im_v${annot_version}.log 
    
    sh $PGI_Repo/code/2_Formatting/2.1.1.2_Merge.sh \
    --mergeType L \
    --mergeCol1 2 \
    --mergeCol2 1 \
    --file1 TMP/tmp_all_im_v${annot_version} \
    --file2 TMP/tmp_gt_v${annot_version} \
    --out TMP/tmp_all_im_gt_v${annot_version} > LOG/tmp_all_im_gt_v${annot_version}.log
    
    echo "Removing 'chr' string from chromosome column and obtaining alleles."
    
    awk -F"\t" '
      NR==1 {
        print "all.data.id","SNPID","CHR","BP","OTHER_ALLELE","EFFECT_ALLELE","PLOIDY","INFO","CALLRATE","HWE_PVAL"
      }
      NR>1 {
        gsub("chr","",$5);gsub("/","\t",$7);print $1,$4,$5,$6,$7,$8,$10,$12,$13
      }' OFS="\t" TMP/tmp_all_im_gt_v${annot_version} > TMP/SNPstats_v${annot_version}.txt &
  
    rm TMP/tmp_all_v${annot_version} TMP/tmp_im_v${annot_version} TMP/tmp_gt_v${annot_version} TMP/tmp_all_im_v${annot_version} TMP/tmp_all_im_gt_v${annot_version} 

  fi


  ## Get all necessary columns from GWAS files. 
  echo "Obtaining necessary columns from GWAS file(s)."
  for gwas_path in ${GWAS[*]}
  do
    gwasCols ${gwas_path} &
  done
  
  wait
  
  ## Merge GWASs with annotation data.
  echo "Merging GWAS results with annotation data." 
  for gwas_path in ${GWAS[*]}
  do
    gwas_file=$(echo ${gwas_path} | sed 's,.\+/,,g')
    sh $PGI_Repo/code/2_Formatting/2.1.1.2_Merge.sh \
    --mergeType R \
    --mergeCol1 1 \
    --mergeCol2 1 \
    --file1 TMP/SNPstats_v${annot_version}.txt \
    --file2 TMP/tmp_${gwas_file} \
    --out TMP/tmp_${gwas_file}_annot > LOG/tmp_${gwas_file}_annot.log &
  done 
  wait

  ## Remove unnecessary columns
  i=0
  for gwas_path in ${GWAS[*]}
  do
    gwas_file=$(echo ${gwas_path} | sed 's,.\+/,,g')
    sh $PGI_Repo/code/2_Formatting/2.1.1.1_ExtractCols.sh \
    --file TMP/tmp_${gwas_file}_annot \
    --cols "SNPID,CHR,BP,EFFECT_ALLELE,OTHER_ALLELE,EAF,BETA,SE,P,INFO,N,IMPUTED,CALLRATE,HWE_PVAL,PLOIDY" \
    --out ${OUT[$i]} > LOG/${OUT[$i]}.log &

    let i+=1

  done
  wait

  ## Clean TMP directory
  for gwas_path in ${GWAS[*]}
  do
    gwas_file=$(echo ${gwas_path} | sed 's,.\+/,,g')
    rm TMP/tmp_${gwas_file} TMP/tmp_${gwas_file}_annot
  done

  ## Remove flipped 1kG imputation SNPs from versions 5.0, 5.1, 5.2
  i=0
  for gwas_path in ${GWAS[*]}
  do
    if [[ ${annot_version} == "5"* ]]
    then
      gwas_file=$(echo ${gwas_path} | sed 's,.\+/,,g')
      preN=$(wc -l ${OUT[$i]} | cut -d" " -f1) 
      awk 'NR==FNR{a[$1]=$1;next}{ChrPosID=$2":"$3}!(ChrPosID in a){print $0}' OFS="\t" ${flippedSNPs_1kG} ${OUT[$i]} > TMP/tmp_${gwas_file}
      mv TMP/tmp_${gwas_file} ${OUT[$i]}  
      postN=$(wc -l ${OUT[$i]} | cut -d" " -f1) 
      echo "$(($preN - $postN)) flipped SNPs were dropped from ${OUT[$i]}."
    fi
    let i+=1
  done
 
}

main

echo -n "Script finished running on "
date

end=$(date +%s)
echo "Execution time was $(($end - $start)) seconds."
