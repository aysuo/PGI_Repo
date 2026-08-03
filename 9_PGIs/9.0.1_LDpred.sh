#!/bin/bash

## To do:
## Add option to provide list of SNPs to be included in the score rather than a bim file.

usage() {
echo "Usage: 
      bash LDpred.sh 
       --efftype <effect type:LINREG,OR,LOGOR,BLUP>
       --chr <chr col>
       --A1 <Effect allele col>
       --A2 <Other allele col>
       --reffreq <freq A1 col>
       --info <info col>
       --snpid <snpid col>
       --pval <pval col>
       --eff <effect col>
       --ncol <N col>
       --pos <bp col>
       --se <se col>
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
    "chr"
    "A1"
    "A2"
    "reffreq"
    "info"
    "snpid"
    "pval"
    "eff"
    "ncol"
    "pos"
    "se"
    "efftype"
    "sumstats"
    "LDgf"
    "Valbim"
    "P"
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
        --chr)
          chr=$2
          shift 2
          ;;
        --A1)
          A1=$2
          shift 2
          ;;
        --A2)
          A2=$2
          shift 2
          ;;
        --reffreq)
          reffreq=$2
          shift 2
          ;;
        --info)
          info=$2
          shift 2
          ;;
        --snpid)
          snpid=$2
          shift 2
          ;;
        --pval)
          pval=$2
          shift 2
          ;;
        --eff)
          eff=$2
          shift 2
          ;;
        --ncol)
          ncol=$2
          shift 2
          ;;
        --pos)
          pos=$2
          shift 2
          ;;
        --se)
          se=$2
          shift 2
          ;;
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
if [[ -z $chr ]]
  then
    chr="CHR"
    echo "No chromosome column has been specified, assuming the default (CHR)."
  else 
    echo "Chr column: $chr"
fi 

if [[ -z $A1 ]]
  then
    A1="EFFECT_ALLELE"
    echo "No effect allele column has been specified, assuming the default (EFFECT_ALLELE)."
  else 
    echo "Effect allele column: $A1"
fi 

if [[ -z $A2 ]]
  then
    A2="OTHER_ALLELE"
    echo "No other allele column has been specified, assuming the default (OTHER_ALLELE)."
  else 
    echo "Other allele column: $A2"
fi 

if [[ -z $reffreq ]]
  then
    reffreq="EAF"
    echo "No effect allele frequency column has been specified, assuming the default (EAF)."
  else 
    echo "Effect allele frequency column: $reffreq"
fi 

#if [[ -z $info ]]
#  then
#    info="INFO"
#    echo "No imputation accuracy column has been specified, assuming the default (INFO)."
#  else 
#    echo "Imputation accuracy column: $info"
#fi 

if [[ -z $snpid ]]
  then
    snpid="SNPID"
    echo "No snp identifier column has been specified, assuming the default (SNPID)."
  else 
    echo "SNP identifier column: $snpid"
fi 

if [[ -z $pval ]]
  then
    pval="PVALUE"
    echo "No P-value column has been specified, assuming the default (PVALUE)."
  else 
    echo "P-value column: $pval"
fi 

if [[ -z $eff ]]
  then
    eff="EFFECT"
    echo "No effect size column has been specified, assuming the default (EFFECT)."
  else 
    echo "Effect size column: $eff"
fi 

if [[ -z $ncol ]]
  then
    ncol="N"
    echo "No sample size column has been specified, assuming the default (N)."
  else 
    echo "Sample size column: $ncol"
fi 

if [[ -z $pos ]]
  then
    pos="BP"
    echo "No base pair position column has been specified, assuming the default (BP)."
  else 
    echo "Base pair position column: $pos"
fi 

if [[ -z $se ]]
  then
    se="SE"
    echo "No standard error column has been specified, assuming the default (SE)."
  else 
    echo "Standard error column: $se"
fi 

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

#if [ -a ./scores ]
#  then
#    echo "'scores' directory already exists. Using the existing directory."
#  else
#    mkdir scores
#fi

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

# Function to coordinate genotypes and GWAS summary statistics
ldpred_coord(){
    echo "Coordinating genotypes and summary statistics.."
    ldpred coord \
    --gf $1 \
    --ssf $2 \
    --vgf $3 \
    --chr $4 \
    --A1 $5 \
    --A2 $6 \
    --reffreq $7 \
    --se $8 \
    --rs $9 \
    --pval ${10} \
    --eff ${11} \
    --ncol ${12} \
    --pos ${13} \
    --eff_type ${14} \
    --ssf-format CUSTOM \
    --match-genomic-pos \
    --max-freq-discrep 0.1 \
    --z-from-se \
    --out=coord/${15}.coord > tmp/${15}_coord.log \
    && echo "Coordination completed" \
    || (echo "Error: Genotypes and summary statistics could not be coordinated." && exit 1)
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
  ldpred_coord $LDgf $sumstats $Valbim $chr $A1 $A2 $reffreq $se $snpid $pval $eff $ncol $pos $efftype $out
  
  get_LDradius $out 
  
  ldpred_weights $out $ld_radius $P 
}

(main && echo "Script executed successfully.") || echo "Error: Script could not be executed successfully." 

# Script over.

echo -n "Script finished running on "
date

end=$(date +%s)
echo "Execution time was $((($end - $start)/3600)) hours $(((($end - $start)%3600)/60)) minutes."

