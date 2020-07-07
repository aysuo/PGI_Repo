#!/bin/bash

## To do:
## Add option to provide list of SNPs to be included in the score rather than a bim file.

usage() {
echo "Usage: 
      bash LDpred.sh 
       --efftype <effect type:LINREG,OR,LOGOR,BLUP>
       --sumstats <path to sumstats>
       --LDgf <LD reference plink file path>
       --Valbim <Bim file prefix with list of SNPs to include in scores>
       --P <List of fraction of causal SNPs for LDpred>
       --hm3 <1 if restrict to hm3>
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
    "efftype"
    "sumstats"
    "LDgf"
    "Valbim"
    "P"
    "hm3"
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
        --efftype)
            efftype=$2
            shift 2
            ;;
        --sumstats)
            sumstats=$2
            shift 2
            ;;
        --LDgf)
            LDgf=$2
            shift 2
            ;;
        --Valbim)
            Valbim=$2
            shift 2
            ;;
        --P)
            P=$2
            shift 2
            ;;
        --hm3)
            hm3=$2
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

if [[ -z $efftype ]]
  then
    efftype="LINREG"
    echo "No effect type has been specified, assuming the default (LINREG)."
  else 
    echo "Effect type: $efftype"
fi 

if [[ -z $sumstats ]]
  then
    echo "Error: No GWAS summary statistics file path has been supplied."
    exit 1
  else
    echo "GWAS sumstats file path: $sumstats."
fi

if [[ -z $LDgf ]]
  then
    echo "Error: No LD reference genotype file path has been supplied."
    exit 1
  else
    echo "LD reference genotype file path: $LDgf."
fi

if [[ -z $Valbim ]]
  then
    echo "Error: No validation bim file path has been supplied."
    exit 1
  else
    echo "Validation genotype bim file path: $Valbim.bim."
fi


if [[ -z $P ]]
  then
    echo "Error: No fraction of causal SNPs has been supplied."
    exit 1
  else
    echo "P: $P."
fi

if [[ -z $hm3 || $hm3==0 ]]
  then
    hm3=0
  else
    echo "Restricting the analysis to HapMap3 SNPs."
fi

if [[ -z $out ]]
  then
    $out=PGS
    echo "No output file prefix has been supplied. Output file prefix is set to PGS."
  else
    echo "Output file prefix: $out."
fi

echo ""
echo "-------------------------------------------------------------"
echo ""

#######################################################################
##################### CREATE WORKING DIRECTORIES ######################
#######################################################################

echo "Creating working directories."
echo ""

if [ -a ./coord ]
  then
    echo "'coord' directory already exists. Using the existing directory."
  else
    mkdir coord
fi

if [ -a ./pickled ]
  then
    echo "'pickled' directory already exists. Using the existing directory."
  else
    mkdir pickled
fi

if [ -a ./weights ]
  then
    echo "'weights' directory already exists. Using the existing directory."
  else
    mkdir weights
fi

if [ -a ./scores ]
  then
    echo "'scores' directory already exists. Using the existing directory."
  else
    mkdir scores
fi

if [ -a ./logs ]
  then
    echo "'logs' directory already exists. Using the existing directory."
  else
    mkdir logs
fi

if [ -a ./tmp ]
  then
    echo "'tmp' directory already exists. Using the existing directory."
  else
    mkdir tmp
fi

echo ""
echo "-------------------------------------------------------------"
echo ""

#######################################################################
########################## DEFINE FUNCTIONS ###########################
#######################################################################

# Function to format GWAS summary statistics as LDpred input


#######################################################################

# Function to coordinate genotypes and GWAS summary statistics
ldpred_coord(){
  if [[ $hm3 == 1 ]]; then
    echo "Coordinating genotypes and summary statistics for HapMap3 SNPs.."
    ldpred coord \
    --gf=$1 \
    --ssf=$2 \
    --vgf=$3 \
    --ssf-format LDPRED \
    --match-genomic-pos \
    --max-freq-discrep 0.1 \
    --eff_type $4 \
    --z-from-se \
    --only-hm3 \
    --out=coord/$5.coord > tmp/$5_coord.log \
    && echo "Coordination completed" \
    || (echo "Error: Genotypes and summary statistics could not be coordinated." && exit 1)
  else
    echo "Coordinating genotypes and summary statistics.."
    ldpred coord \
    --gf=$1 \
    --ssf=$2 \
    --vgf=$3 \
    --ssf-format LDPRED \
    --match-genomic-pos \
    --max-freq-discrep 0.1 \
    --eff_type $4 \
    --z-from-se \
    --out=coord/$5.coord > tmp/$5_coord.log \
    && echo "Coordination completed" \
    || (echo "Error: Genotypes and summary statistics could not be coordinated." && exit 1)
  fi
}

#######################################################################

# Function to obtain LD radius
get_LDradius(){
  echo "Obtaining LD radius (#SNPs/3000)"
  common_SNPs=$(sed -n '/SNPs retained after filtering:/p' tmp/$1_coord.log | cut -d":" -f2 | sed 's/ //g') \
  && ld_radius=$(( $common_SNPs/3000 )) \
  && echo "LD radius will be set to $ld_radius" \
  || (echo "Error: LD radius could not be obtained." && exit 1)
}

#######################################################################

# Function to calculate LD-adjusted weights
ldpred_weights(){
  echo "Calculating LD-adjusted SNP weights with LD radius=$ld_radius , N=$median_N and P=$P"
  ldpred gibbs \
  --cf=coord/$1.coord \
  --ldr=$2 \
  --ldf=pickled/$1 \
  --f=$3 \
  --out=weights/$1_weights \
  && echo "LD adjusted weights calculated." \
  || (echo "Error: LDpred weights couldn't be obtained." && exit 1)
}

#######################################################################


#######################################################################
####################### DEFINE AND EXECUTE MAIN #######################
#######################################################################


main(){
  ldpred_coord $LDgf $sumstats $Valbim $efftype $out $hm3
  
  get_LDradius $out 
  
  ldpred_weights $out $ld_radius $P 
}


(main && echo "Script executed successfully.") || echo "Error: Script could not be executed successfully." 

# Script over.

echo -n "Script finished running on "
date

end=$(date +%s)
echo "Execution time was $((($end - $start)/3600)) hours $(((($end - $start)%3600)/60)) minutes."
