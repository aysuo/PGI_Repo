#!/bin/bash

## To do:
## Add option to provide list of SNPs to be included in the score rather than a bim file.
## Add function to make P&T scores


usage() {
echo "Usage: 
      bash make_PGS.sh 
       --valgf <validation data file path in plink2 per chr format>
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


#######################################################################
############################ PARSE ARGUMENTS ##########################
#######################################################################

ARGUMENT_LIST=(
    "valgf"
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

    paste scores/PGS_${out}_chr*.sscore > tmp/PGS_${out}
    awk -v out=${out} 'BEGIN{OFS="\t"; print "FID","IID","ALLELE_CT","NAMED_ALLELE_DOSAGE_SUM","PGS_"out} \
      NR>1{nallele=0;dosagesum=0;PGS=0; for(i=1;i<=NF;i++){if (i%5 == 3) nallele+=$i; if (i%5 == 4) dosagesum+=$i ; if (i%5 == 0) PGS+=$i}; print $1,$2,nallele,dosagesum,PGS/nallele}' OFS="\t" tmp/PGS_${out} > scores/PGS_${out}.txt
     
    cat scores/PGS_${out}_chr*.log > logs/PGS_${out}_make_PGS.log 
    rm tmp/PGS_${out} scores/PGS_${out}_chr*
}

makePGS(){
  valgf=$1
  sampleKeep=$2
  out=$3
  weights=$4
  weightCols=$5

  if [[ $sampleKeep == "NA" ]]; then
    chr1=$(echo $valgf | sed 's/\[1:22\]/1/g')
    sampleKeep=${chr1}.psam
  fi

  i=0
  for chr in {1..22}; do
    valgfChr=$(echo $valgf | sed "s/\[1:22\]/$chr/g")
    plink2 --pfile ${valgfChr} \
      --keep $sampleKeep \
      --rm-dup force-first \
      --score $weights $weightCols cols=fid,nallele,dosagesum,scoresums \
      --out scores/PGS_${out}_chr$chr &

    let i+=1

    if [[ $i == 5 ]]; then
		  wait
		  i=0
	  fi
  done  
  
  wait

  mergeChr $out
}


main(){
  makePGS "$valgf" $sampleKeep $out $weights "$weightCols"
}

main
