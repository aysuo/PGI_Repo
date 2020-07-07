#!/bin/bash

## To do:
## Add option to provide list of SNPs to be included in the score rather than a bim file.
## Optional: Add function to make P&T scores


usage() {
echo "Usage: 
      bash make_PGS.sh 
       --valgf <validation data file path; enter per chr files as e.g. HRS_chr[1:22].gen.gz>
       --valgfFormat <valgf format: plink, plinkchr, gen, bgen, vcf >
       --sample <sample file path for bgen/gen/gprobs files (optional)>
       --genoFormat <dosage or hardcall>
       --weights <path to weights>
       --weightCols <column numbers for SNPID, effect allele and effect (comma-separated)>
       --sampleKeep <list of individuals to make scores for>
       --out <output file prefix>"

echo "note: order of options is not important"
}

echo ""
echo "-------------------------------------------------------------"
echo ""


echo -n "Script started on "
date
start=$(date +%s)


snpid_map="/disk/genetics3/PGS/PGS_Repo/data/REFFILES/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.rsid_map"


#######################################################################
############################ PARSE ARGUMENTS ##########################
#######################################################################

ARGUMENT_LIST=(
    "valgf"
    "sample"
    "valgfFormat"
    "genoFormat"
    "weights"
    "weightCols"
    "sampleKeep"
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
        --valgf)
            valgf=$2
            shift 2
            ;;
        --valgfFormat)
            valgfFormat=$2
            shift 2
            ;;
        --sample)
            sample=$2
            shift 2
            ;;
        --genoFormat)
            genoFormat=$2
            shift 2
            ;;
        --weights)
            weights=$2
            shift 2
            ;;
        --weightCols)
            weightCols=$2
            shift 2
            ;;
        --sampleKeep)
            sampleKeep=$2
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

if [[ -z $valgf ]]; then
  echo "Error: No validation genotype file path has been specified."
  exit 1
fi

case "$valgfFormat" in
    plink)
        echo "Validation genotype data is in plink format: $valgf"
        ;;
    plinkChr)
        echo "Validation genotype data is in per chromosome plink format: $valgf"
        ;;
    gen)
        echo -e "Validation genotype data is in per chromosome gen format: $valgf"
        ;;           
    bgen)
        echo "Validation genotype data is in bgen format: $valgf"
        ;;
    vcf)
        echo "Validation genotype data is in vcf format: $valgf"
        ;;
    *)
        echo "Error: Please specify a valid genotype data format."    
        usage
        exit 1
        ;;
esac

if [[ $valgfFormat == "gen" || $valgfFormat == "bgen" ]]; then
  if [[ -z $sample ]]; then
    echo "Error: No sample file path for gen data has been specified."
    exit 1
  else
    echo "Sample file path for gen data: $sample."
  fi
fi


if [[ $genoFormat == "dosage" ]]; then
  if [[ $valgfFormat == "gen" ]]; then
    echo "PGS will be made using dosages".
  else
    echo "Error: Plink cannot use dosages from $valgfFormat data."
  fi
else
  echo "PGS will be made using hard calls."
fi

if [[ -z $weights ]]
  then
    echo "Error: No file path for the weights has been specified."
    exit 1
  else
    echo "Weights: $weights."
fi

if [[ -z $weightCols ]]
  then
    echo "No column numbers for SNPID, effect allele and effect have been specified, assuming 1,2,3."
    weightCols="1 2 3"
  else
    echo "SNPID, effect allele and effect columns: $weightCols."
    weightCols=$(echo $weightCols | sed 's/,/ /g')
fi

if [[ -z $sampleKeep ]]
  then
    echo "No list of individuals to restrict PGS sample has been specified. PGS will be made for the whole sample."
    sampleKeep="NA"
  else
    echo "PGS will be made for the subsample of individuals specified in $sampleKeep."
fi


if [[ -z $out ]]
  then
    echo "Error: No output prefix has been specified."
    exit 1
  else
    echo "Output file prefix: $out."
fi

echo ""
echo "-------------------------------------------------------------"
echo ""


#######################################################################################################


mergeChr(){
    out=$1
    genoFormat=$2
    valgfFormat=$3

    paste scores/PGS_${out}_chr*.profile > tmp/PGS_${out}

    if [[ $valgfFormat == "bgen" ]]; then
      awk -v out=${out} 'BEGIN{OFS="\t"; print "FID","IID","CNT","CNT2","PGS_"out} \
        NR>1{PGS=0; for(i=1;i<=NF;i++){if (i%4 == 0) PGS+=$i}; print $1,$2,"NA","NA",PGS}' OFS="\t" tmp/PGS_${out} > scores/PGS_${out}.txt
      elif ! [[ $genoFormat == "dosage" ]]; then
        awk -v out=${out} 'BEGIN{OFS="\t"; print "FID","IID","CNT","CNT2","PGS_"out} \
          NR>1{CNT=0;CNT2=0;PGS=0; for(i=1;i<=NF;i++){if (i%6 == 4) CNT+=$i ; if (i%6 == 5) CNT2+=$i ; if (i%6 == 0) PGS+=$i}; print $1,$2,CNT,CNT2,PGS}' OFS="\t" tmp/PGS_${out} > scores/PGS_${out}.txt
      else
        awk -v out=${out} 'BEGIN{OFS="\t"; print "FID","IID","CNT","CNT2","PGS_"out} \
          NR>1{CNT=0;PGS=0; for(i=1;i<=NF;i++){if (i%5 == 4) CNT+=$i ; if (i%5 == 0) PGS+=$i}; print $1,$2,CNT,"NA",PGS}' OFS="\t" tmp/PGS_${out} > scores/PGS_${out}.txt
    fi
    
    cat scores/PGS_${out}_chr*.log > logs/PGS_${out}_make_PGS.log 
    rm tmp/PGS_${out} scores/PGS_${out}_chr*
}


make_PGS_plink(){
  valgf=$1
  sampleKeep=$2
  out=$3
  weights=$4
  weightCols=$5
  
  if [[ $sampleKeep == "NA" ]]
  then
    sampleKeep=${valgf}.fam
  fi
    
  plink1.9 --bfile $valgf \
     --keep $sampleKeep \
     --score $weights $weightCols \
     --out scores/PGS_${out}

  mv scores/PGS_${out}.profile scores/PGS_${out}.txt
  sed -i 's/ \+/\t/g' scores/PGS_${out}.txt
  sed -i 's/^\t//g' scores/PGS_${out}.txt
}


make_PGS_plink_chr(){
  valgf=$1
  sampleKeep=$2
  out=$3
  genoFormat=$4
  weights=$5
  weightCols=$6

  valplink=$(echo $valgf | sed 's/\[1:22\]/\${chr}/g')

  if [[ $sampleKeep == "NA" ]]
  then
    fam1=$(echo $valgf | sed 's/\[1:22\]/1/g')
    sampleKeep=${fam1}.fam
  fi

  for chr in {1..22}
  do
      eval valplinkchr=$valplink
      plink1.9 --bfile ${valplinkchr} \
        --keep $sampleKeep \
        --score $weights $weightCols \
        --out scores/PGS_${out}_chr$chr &
  done  
  wait

  mergeChr $out $genoFormat
}

make_PGS_gen(){
  valgf=$1
  sampleKeep=$2
  out=$3
  genoFormat=$4
  weights=$5
  weightCols=$6
  sample=$7

  valgen=$(echo "$valgf" | sed 's/\[1:22\]/\${chr}/g')
  
  if [[ $sample == *.gz ]]; then
    gunzip -c $sample > tmp/${out}.sample
    sample="tmp/${out}.sample"
  fi

  awk 'NR>2{print $1,$2,0,0,0,-9}' ${sample} > tmp/${out}.fam

  if [[ $sampleKeep == "NA" ]]
  then
    sampleKeep="tmp/${out}.fam"
  fi
  
  if [[ $genoFormat == "dosage" ]]
  then
    for chr in {1..22}; do
      eval valgenchr=$valgen
      plink1.9 --dosage ${valgenchr} skip0=1 skip1=1 format=3 noheader \
        --fam tmp/${out}.fam \
        --keep $sampleKeep \
        --score $weights $weightCols include-cnt \
        --out scores/PGS_${out}_chr$chr &
    done
  else
    for chr in {1..22}; do
      eval valgenchr=$valgen
      plink1.9 --gen ${valgenchr} \
        --sample ${sample} \
        --keep $sampleKeep \
        --score $weights $weightCols include-cnt \
        --out scores/PGS_${out}_chr$chr &
    done
  fi

  wait

  mergeChr $out $genoFormat
}


make_PGS_bgen(){
  valgf=$1
  sampleKeep=$2
  out=$3
  genoFormat=$4
  weights=$5
  weightCols=$6
  sample=$7

  valbgen=$(echo "$valgf" | sed 's/\[1:22\]/\${chr}/g')

  if [[ $sampleKeep == "NA" ]]
  then
    for chr in {1..22}; do
      eval valbgenchr=$valbgen
      plink2 --bgen ${valbgenchr} \
        --score $weights $weightCols \
        --sample $sample \
        --out scores/PGS_${out}_chr$chr &
    done
  else
    for chr in {1..22}; do
      eval valbgenchr=$valbgen
      plink2 --bgen ${valbgenchr} \
        --keep $sampleKeep \
        --score $weights $weightCols \
        --sample $sample \
        --out scores/PGS_${out}_chr$chr &
    done
  fi
  wait

  #mergeChr $out $genoFormat
}


make_PGS_vcf(){
  valgf=$1
  sampleKeep=$2
  out=$3
  genoFormat=$4
  weights=$5
  weightCols=$6

  valvcf=$(echo $valgf | sed 's/\[1:22\]/\${chr}/g')

  if [[ $sampleKeep == "NA" ]]
  then
    for chr in {1..22}; do
      eval valvcfchr=$valvcf
      plink1.9 --vcf ${valvcfchr} \
        --score $weights $weightCols \
        --include-cnt \
        --out scores/PGS_${out}_chr$chr &
    done
  else
    for chr in {1..22}; do
      eval valvcfchr=$valvcf
      plink1.9 --bgen ${valvcfchr} \
        --keep $sampleKeep \
        --score $weights $weightCols \
        --include-cnt \
        --out scores/PGS_${out}_chr$chr &
    done
  fi
  wait

  mergeChr $out $genoFormat
}


main(){
  case "$valgfFormat" in
    plink)
        make_PGS_plink "$valgf" $sampleKeep $out $weights "$weightCols"
        ;;
    plinkChr)
        make_PGS_plink_chr "$valgf" $sampleKeep $out $genoFormat $weights "$weightCols"
        ;;
    gen)
        make_PGS_gen "$valgf" $sampleKeep $out $genoFormat $weights "$weightCols" $sample
        ;;           
    bgen)
        make_PGS_bgen "$valgf" $sampleKeep $out $genoFormat $weights "$weightCols" $sample
        ;;
    vcf)
        make_PGS_vcf "$valgf" $sampleKeep $out $genoFormat $weights "$weightCols"
        ;;
  esac
}


main