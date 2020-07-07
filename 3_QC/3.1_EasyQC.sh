#!/bin/bash
export R_LIBS=/homes/nber/aokbay/R/x86_64-redhat-linux-gnu-library/3.1/:$R_LIBS

#-------------------------------------------------------------------------#
# Creates and runs EasyQC script
# Date Updated: 03/01/2020
# Author: Aysu Okbay
#-------------------------------------------------------------------------#


## To do: Add option for other effect types (e.g. OR, logOR, etc)

usage() {
echo "Usage: 
      bash EasyQC.sh 
       --fileIn <input file path>
       --sep <input file delimiter, default = WHITESPACE>
       --miss <input file missing value coding, default = NA>
       --pathOut <output file path>
       --tag <output file name tag>
       --SNPID <SNP id column name>
       --SNPIDtype <SNP id type (rs or chrpos)>
       --EA <effect allele column name>
       --OA <alternative allele column name>
       --EAF <effect allele frequency column name>
       --EFFECT <effect column name>
       --SE <SE column name>
       --P <p-value column name>
       --N <sample size column name>
       --EFFECTtype <effect type, default = BETA>
       --CHR <chromosome column name, optional>
       --BP <base pair column name, optional>
       --IMPUTED <imputed SNP identifier column equal to 1 if imputed, 0 if genotyped>
       --INFO <imputation accuracy column name, optional>
       --HWE <HWE p-value column name, optional>
       --CALLRATE <Callrate column name, optional>
       --cutoff_MAF <MAF cutoff, default = 0.01>
       --cutoff_MAC <MAC cutoff, optional>
       --cutoff_R2 <R2 cutoff, optional>
       --cutoff_INFO <INFO cutoff, default = 0.7>
       --cutoff_HWE <HWE pval cutoff, optional>
       --cutoff_CALLRATE <call rate cutoff, optional>
       --cutoff_SE <SE ratio cutoff>
       --cutoff_N <N cutoff>
       --XY <1 if sex chromosomes are to be kept, 0 otherwise, default = 0>
       --INDEL <1 if INDELs are to be kept, 0 otherwise, default = 0>
       --cptref <cptid reference file path, default = /disk/genetics/ukb/aokbay/reffiles/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.rsid_map>
       --cptref_marker <marker name in cptid reference file, default = rsid>
       --cptref_chr <chr name in cptid reference file, default = chr>
       --cptref_bp <bp name in cptid reference file, default = pos>
       --afref <allele freq reference file path, default = /disk/genetics/ukb/aokbay/reffiles/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.cptid.maf001.gz>
       --afref_marker <marker name in allele freq reference file, default = cptid>
       --afref_ref  <reference allele name in allele freq reference file, default = ref>
       --afref_alt  <alternative allele name in allele freq reference file, default = alt>
       --afref_raf <reference allele frequency column in allele freq reference file, default = raf>
       --snpStd <1 if genotypes are standardized, 0 otherwise, default=0>
       --SDy <SD of phenotype, required if phenoStd=0>
       --dom <temporary option to handle dominance r2 filter>"

echo "note: order of options is not important"
}

echo ""
echo "-------------------------------------------------------------"
echo ""

#######################################################################
############################ PARSE ARGUMENTS ##########################
#######################################################################

ARGUMENT_LIST=(
   "fileIn"
   "sep"
   "miss"
   "pathOut"
   "tag"
   "SNPID"
   "SNPIDtype"
   "EA"
   "OA"
   "EAF"
   "EFFECT"
   "SE"
   "P"
   "N"
   "EFFECTtype"
   "CHR"
   "BP"
   "IMPUTED"
   "INFO"
   "HWE"
   "CALLRATE"
   "cutoff_MAF"
   "cutoff_MAC"
   "cutoff_R2"
   "cutoff_INFO"
   "cutoff_HWE"
   "cutoff_CALLRATE"
   "cutoff_SE"
   "cutoff_N"
   "XY"
   "INDEL"
   "cptref"
   "cptref_marker"
   "cptref_chr"
   "cptref_bp"
   "afref"
   "afref_marker"
   "afref_ref"
   "afref_alt"
   "afref_raf"
   "snpStd"
   "SDy"
   "dom"
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
        --fileIn)
            fileIn=$2
            shift 2
            ;;
        --sep)
            sep=$2
            shift 2
            ;;
        --miss)
            miss=$2
            shift 2
            ;;
        --pathOut)
            pathOut=$2
            shift 2
            ;;
        --tag)
            tag=$2
            shift 2
            ;;
        --SNPID)
            SNPID=$2
            shift 2
            ;;
        --SNPIDtype)
            SNPIDtype=$2
            shift 2
            ;;
        --EA)
            EA=$2
            shift 2
            ;;
        --OA)
            OA=$2
            shift 2
            ;;        
        --EAF)
            EAF=$2
            shift 2
            ;;
        --EFFECT)
            EFFECT=$2
            shift 2
            ;;
        --SE)
            SE=$2
            shift 2
            ;;
        --P)
            P=$2
            shift 2
            ;;
        --N)
            N=$2
            shift 2
            ;;
        --EFFECTtype)
            EFFECTtype=$2
            shift 2
            ;;
        --CHR)
            CHR=$2
            shift 2
            ;;
        --BP)
            BP=$2
            shift 2
            ;;
        --IMPUTED)
            IMPUTED=$2
            shift 2
            ;;
        --INFO)
            INFO=$2
            shift 2
            ;;
        --HWE)
            HWE=$2
            shift 2
            ;;
        --CALLRATE)
            CALLRATE=$2
            shift 2
            ;;
        --cutoff_MAF)
            cutoff_MAF=$2
            shift 2
            ;;
        --cutoff_MAC)
            cutoff_MAC=$2
            shift 2
            ;;
        --cutoff_R2)
            cutoff_R2=$2
            shift 2
            ;;            
        --cutoff_INFO)
            cutoff_INFO=$2
            shift 2
            ;;
        --cutoff_HWE)
            cutoff_HWE=$2
            shift 2
            ;;
        --cutoff_CALLRATE)
            cutoff_CALLRATE=$2
            shift 2
            ;;
        --cutoff_SE)
            cutoff_SE=$2
            shift 2
            ;;
        --cutoff_N)
            cutoff_N=$2
            shift 2
            ;;
        --XY)
            XY=$2
            shift 2
            ;;
        --INDEL)
            INDEL=$2
            shift 2
            ;;
        --cptref)
            cptref=$2
            shift 2
            ;;
        --cptref_marker)
            cptref_marker=$2
            shift 2
            ;;
        --cptref_chr)
            cptref_chr=$2
            shift 2
            ;;
        --cptref_bp)
            cptref_bp=$2
            shift 2
            ;;
        --afref)
            afref=$2
            shift 2
            ;;
        --afref_marker)
            afref_marker=$2
            shift 2
            ;;
        --afref_ref)
            afref_ref=$2
            shift 2
            ;;
        --afref_alt)
            afref_alt=$2
            shift 2
            ;;
        --afref_raf)
            afref_raf=$2
            shift 2
            ;;
        --snpStd)
            snpStd=$2
            shift 2
            ;;
        --SDy)
            SDy=$2
            shift 2
            ;;
        --dom)
            dom=$2
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

if [[ -z $fileIn ]]; then
    echo "Error: Input file path not specified."
    exit 1
else
  echo "Input file: $fileIn ."
fi

if [[ -z $sep ]]; then
  sep="WHITESPACE"
fi
echo "Input file separator: $sep ."

if [[ -z $miss ]]; then
  miss="NA"
fi
echo "Missing value coding in input file: $miss ."

if [[ -z $pathOut ]]; then
    echo "Error: Output path not specified."
    exit 1
else
  echo "Output file: $pathOut ."
fi

if [[ -z $tag ]]; then
    echo "Error: Output file tag not specified."
    exit 1
else
  echo "Output file tag: $tag ."
fi

if [[ -z $SNPID ]]; then
    echo "Error: SNP identifier column name not specified."
    exit 1
else
  echo "SNP identifier column: $SNPID ."
fi
 
if [[ -z $SNPIDtype ]]; then
    echo "Error: SNP identifier type not specified."
    exit 1
else
  echo "SNP identifier type: $SNPIDtype ."
fi

if [[ -z $EA ]]; then
    echo "Error: Effect allele column not specified."
    exit 1
else
  echo "Effect allele column: $EA ."
fi

if [[ -z $OA ]]; then
    echo "Error: Non-effect allele column not specified."
    exit 1
else
  echo "Non-effect allele column: $OA ."
fi

if [[ -z $EAF ]]; then
    echo "Error: Effect allele frequency column not specified."
    exit 1
else
  echo "Effect allele frequency column: $EAF ."
fi

if [[ -z $EFFECT ]]; then
    echo "Error: Effect column not specified."
    exit 1
else
  echo "Effect column: $EFFECT ."
fi

if [[ -z $EFFECTtype ]]; then
    EFFECTtype="BETA"
fi
echo "Effect type: $EFFECTtype ."

if [[ -z $SE ]]; then
    echo "Error: SE column not specified."
    exit 1
else
  echo "SE column: $SE ."
fi
  
if [[ -z $P ]]; then
    echo "Error: P-value column not specified."
    exit 1
else
  echo "P-value column: $P ."
fi

if [[ -z $N ]]; then
    echo "Error: Sample size column not specified."
    exit 1
else
  echo "Sample size column: $N ."
fi

if [[ ! -z $CHR ]]; then
    echo "Chr column: $CHR ."
fi

if [[ ! -z $BP ]]; then
    echo "BP column: $BP ."
fi

if [[ ! -z $IMPUTED ]]; then
    echo "Imputed SNP identifier column: $IMPUTED ."
fi

if [[ ! -z $INFO ]]; then
    echo "Imputation accuracy column: $INFO ."
fi

if [[ ! -z $HWE ]]; then
    echo "HWE P-value column: $HWE ."
fi

if [[ ! -z $CALLRATE ]]; then
    echo "Call rate column: $CALLRATE ."
fi

if [[ -z $cutoff_MAF ]]; then
    cutoff_MAF=0.01
fi
echo "MAF cutoff: $cutoff_MAF ."

if [[ ! -z $cutoff_MAC ]]; then
    echo "MAC cutoff = $cutoff_MAC ."
fi

if [[ ! -z $cutoff_R2 ]]; then
    echo "R2 cutoff: $cutoff_R2 ."  
fi


if [[ ! -z $cutoff_SE ]]; then
    echo "SE ratio cutoff: $cutoff_SE ."
fi


if [[ -z $cutoff_INFO ]]; then
    cutoff_INFO=0.7
fi
echo "INFO cutoff: $cutoff_INFO ."

if [[ ! -z $cutoff_HWE ]]; then
  echo "HWE P-value cutoff = $cutoff_HWE ."
fi

if [[ ! -z $cutoff_CALLRATE ]]; then
  echo "Callrate cutoff = $cutoff_CALLRATE ."
fi

if [[ -z $XY || $XY == 0 ]]; then
  XY=0
  echo "Sex chromosomes will be dropped."
else
  echo "Sex chromosomes will be kept."
fi

if [[ -z $INDEL || $INDEL == 0 ]]; then
  INDEL=0
  echo "INDELs will be dropped."
else
  echo "INDELs will be kept."
fi

if [[ -z $cptref ]]; then
  cptref="/disk/genetics/ukb/aokbay/reffiles/1kG_phase1/rsmid_map.1000G_ALL_p1v3.merged_mach_impute.v1.txt"
fi
echo "Reference file to obtain cptid's: $cptref ."

if [[ -z $cptref_marker ]]; then
  cptref_marker="rsid"
fi
echo "Marker name in cptid reference file: $cptref_marker ."

if [[ -z $cptref_chr ]]; then
  cptref_chr="chr"
fi
echo "Chromosome name in cptid reference file: $cptref_chr ."

if [[ -z $cptref_bp ]]; then
  cptref_bp="pos"
fi
echo "Base pair name in cptid reference file: $cptref_bp ."

if [[ -z $afref ]]; then
  afref="/disk/genetics/ukb/aokbay/reffiles/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.cptid.maf001.gz"
fi
echo "Allele frequency reference file: $afref ."

if [[ -z $afref_marker ]]; then
  afref_marker="cptid"
fi
echo "Marker name in allele frequency reference file: $afref_marker ."

if [[ -z $afref_ref ]]; then
  afref_ref="ref"
fi
echo "Reference allele name in allele frequency reference file: $afref_ref ."

if [[ -z $afref_alt ]]; then
  afref_alt="alt"
fi
echo "ALternative allele name in allele frequency reference file: $afref_alt ."

if [[ -z $afref_raf ]]; then
  afref_raf="raf"
fi
echo "Allele frequency name in allele frequency reference file: $afref_raf ."

if [[ -z $snpStd || $snpStd -eq 0 ]]; then
  snpStd=0
  echo "Genotypes are not standardized."
else
  echo "Genotypes are standardized."
fi

if [[ -z $SDy ]]; then
  echo "Error: Please provide SD of the phenotype."
  exit 1
else
  echo "SD of phenotype: $SDy ."
fi


echo ""
echo "-------------------------------------------------------------"
echo ""

if [ -a ./$pathOut/QC_${tag}_$(date +"%Y_%m_%d") ]
  then
    echo "Output directory already exists. Using the existing directory."
  else
    mkdir $pathOut/QC_${tag}_$(date +"%Y_%m_%d")
fi


## DEFINE OUTPUT COLUMN NAMES AND CLASSES
SNPID_out=SNPID ; SNPID_class=character
EA_out=EFFECT_ALLELE ; EA_class=character
OA_out=OTHER_ALLELE ; OA_class=character
EAF_out=EAF ; EAF_class=numeric
EFFECT_out=EFFECT ; EFFECT_class=numeric
SE_out=SE ; SE_class=numeric
P_out=PVALUE ; P_class=numeric
N_out=N ; N_class=numeric
CHR_out=CHR ; CHR_class=character
BP_out=BP ; BP_class=integer
IMPUTED_out=IMPUTED ; IMPUTED_class=integer
INFO_out=INFO ; INFO_class=numeric
HWE_out=HWE_PVAL ; HWE_class=numeric 
CALLRATE_out=CALLRATE ; CALLRATE_class=numeric

# Required columns 
cols_in=$(echo -n "${SNPID};${EA};${OA};${EAF};${EFFECT};${SE};${P};${N}")
cols_out=$(echo -n "${SNPID_out};${EA_out};${OA_out};${EAF_out};${EFFECT_out};${SE_out};${P_out};${N_out}")
cols_class=$(echo -n "${SNPID_class};${EA_class};${OA_class};${EAF_class};${EFFECT_class};${SE_class};${P_class};${N_class}")

# Optional columns
for i in CHR BP IMPUTED INFO HWE CALLRATE; do
  eval col='$'{$i}
  if [[ ! -z $col ]]; then
    cols_in=$(echo "$cols_in;$col")
    eval col_out='$'{${i}_out}
    cols_out=$(echo "$cols_out;$col_out")
    eval col_class='$'{${i}_class}
    cols_class=$(echo "$cols_class;${col_class}")
  fi
done

# Start writing ecf file
echo "DEFINE  --pathOut $pathOut/QC_${tag}_$(date +"%Y_%m_%d")
  --strMissing $miss
  --strSeparator $sep
  --acolIn $cols_in
  --acolInClasses $cols_class
  --acolNewName $cols_out
" > $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf

filelist=$(echo $fileIn | sed 's/,/ /g')
taglist=$(echo $tag | sed 's/,/ /g')
declare -a files=$(echo "($filelist)")
declare -a tags=$(echo "($taglist)")
N_files=${#files[@]}


for (( i=0; i<${N_files}; i++ )); do
  echo "EASYIN  --fileIn ${files[$i]} --fileInShortName ${tags[$i]}" >>  $pathOut/QC_${tags[$i]}_$(date +"%Y_%m_%d")/QC_${tags[$i]}.ecf
done

echo "
START EASYQC

####################
## 1. Sanity checks: 

CLEAN --rcdClean is.na(${EA_out}) | is.na(${OA_out}) --strCleanName numDrop_Missing_Alleles --blnWriteCleaned 0
CLEAN --rcdClean is.na(${P_out}) --strCleanName numDrop_Missing_P --blnWriteCleaned 0
CLEAN --rcdClean is.na(${EFFECT_out}) --strCleanName numDrop_Missing_EFFECT --blnWriteCleaned 0
CLEAN --rcdClean is.na(${SE_out}) --strCleanName numDrop_Missing_SE --blnWriteCleaned 0
CLEAN --rcdClean is.na(${EAF_out}) --strCleanName numDrop_Missing_EAF --blnWriteCleaned 0
CLEAN --rcdClean is.na(${N_out}) --strCleanName numDrop_Missing_N --blnWriteCleaned 0" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf

if [[ ! -z $IMPUTED ]]; then
  echo "CLEAN --rcdClean is.na(${IMPUTED_out}) --strCleanName numDrop_Missing_IMPUTED --blnWriteCleaned 0" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
  echo "CLEAN --rcdClean ${IMPUTED_out}==1 & is.na(${INFO_out}) --strCleanName numDrop_Missing_INFO --blnWriteCleaned 0" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
  echo "CLEAN --rcdClean ${IMPUTED_out}!=0 & ${IMPUTED_out}!=1 --strCleanName numDrop_invalid_IMPUTED --blnWriteCleaned 0" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
fi

echo "
CLEAN --rcdClean ${P_out}<0 | ${P_out}>1 --strCleanName numDrop_invalid_PVAL --blnWriteCleaned 0
CLEAN --rcdClean ${SE_out}<=0 | ${SE_out}==Inf --strCleanName numDrop_invalid_SE --blnWriteCleaned 0
CLEAN --rcdClean ${EAF_out}<0 | ${EAF_out}>1 --strCleanName numDrop_invalid_EAF --blnWriteCleaned 0


####################
## 2. Prepare files for filtering and apply minimum thresholds: 

CLEAN --rcdClean ${EAF_out}==0 | ${EAF_out}==1 --strCleanName numDrop_Monomorph --blnWriteCleaned 0" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf

if [[ ! -z $IMPUTED && ! -z $cutoff_HWE ]]; then
  echo "CLEAN --rcdClean ${IMPUTED_out}==0 & ${HWE_out}<${cutoff_HWE} --strCleanName numDrop_HWE_${cutoff_HWE} --blnWriteCleaned 0" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
fi

if [[ ! -z $IMPUTED && ! -z $cutoff_CALLRATE ]]; then
  echo "CLEAN --rcdClean ${IMPUTED_out}==0 & ${CALLRATE_out}<${cutoff_CALLRATE} --strCleanName numDrop_CALLRATE_${cutoff_CALLRATE} --blnWriteCleaned 0" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
fi

echo "CLEAN --rcdClean ${EAF_out}<${cutoff_MAF} | ${EAF_out}>1-${cutoff_MAF} --strCleanName numDrop_MAF_${cutoff_MAF} --blnWriteCleaned 0" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf

if [[ ! -z $cutoff_MAC ]]; then
  echo "ADDCOL --rcdAddCol signif(2*pmin(${EAF_out},1-${EAF_out})*${N_out},4) --colOut MAC" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
  echo "CLEAN --rcdClean MAC<${cutoff_MAC} --strCleanName numDrop_MAC_${cutoff_MAC} --blnWriteCleaned 0" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
fi

## Exclude imputed SNPs with missing or low Imputation Quality 
if [[ ! -z $IMPUTED ]]; then
  echo "GETNUM  --rcdGetNum ${IMPUTED_out}==1 --strGetNumName num_Imputed" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
elif [[ ! -z $INFO ]]; then
  echo "GETNUM  --rcdGetNum ${INFO_out}<1 --strGetNumName num_Imputed" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
fi

if [[ ! -z ${IMPUTED} && ! -z ${INFO} ]]; then
  echo "CLEAN --rcdClean ${IMPUTED_out}==1 & ${INFO_out}<${cutoff_INFO} --strCleanName numDrop_Imputed_lowImpQual --blnWriteCleaned 0" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
elif [[ -z ${IMPUTED} && ! -z ${INFO} ]]; then
  echo "CLEAN --rcdClean ${INFO_out}<${cutoff_INFO} --strCleanName numDrop_Imputed_lowImpQual --blnWriteCleaned 0" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
fi


if [[ ! -z ${cutoff_N} ]]; then
   echo "CLEAN --rcdClean ${N_out} < ${cutoff_N}  --strCleanName numDrop_N --blnWriteCleaned 0" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
fi

## To do: Add option for other effect types (e.g. OR, logOR, etc)
if [[ $dom != 1 ]]; then
  if [[ ${EFFECTtype} == "BETA" ]]; then
    if [[ ${snpStd} == 0 ]]; then
      echo "ADDCOL --rcdAddCol ${EFFECT_out}*${EFFECT_out}*2*${EAF_out}*(1-${EAF_out})/(${SDy}*${SDy}) --colOut R2" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
      echo "ADDCOL --rcdAddCol ${SDy}/sqrt(2*${N_out}*${EAF_out}*(1-${EAF_out})) --colOut SEpred" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
    elif [[ ${snpStd} == 1 ]]; then
      echo "ADDCOL --rcdAddCol ${EFFECT_out}*${EFFECT_out}/(${SDy}*${SDy}) --colOut R2" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
      echo "ADDCOL --rcdAddCol ${SDy}/sqrt(${N_out}) --colOut SEpred" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
    fi
  fi
else
  echo "ADDCOL --rcdAddCol ${EFFECT_out}*${EFFECT_out}*2*${EAF_out}*(1-${EAF_out})/(${SDy}*${SDy}*(1-2*${EAF_out}*(1-${EAF_out}))) --colOut R2" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
  echo "ADDCOL --rcdAddCol ${SDy}*sqrt(1-(2*${EAF_out}*(1-${EAF_out})))/sqrt(2*${N_out}*${EAF_out}*(1-${EAF_out})) --colOut SEpred" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
fi

if [[ ! -z ${cutoff_R2} ]]; then
   echo "CLEAN --rcdClean R2 > ${cutoff_R2}  --strCleanName numDrop_R2 --blnWriteCleaned 0" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
fi

if [[ ! -z ${cutoff_SE} ]]; then
   echo "CLEAN --rcdClean ${SE_out}/SEpred > ${cutoff_SE} | SEpred/${SE_out} > ${cutoff_SE} --strCleanName numDrop_SEratio --blnWriteCleaned 1" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
fi



echo "
####################
#### 3. Harmonization of allele coding (I/D)
## The aim of this step is to compile uniform allele codes A/C/G/T or I/D from different versions of given alleles

HARMONIZEALLELES --colInA1 ${EA_out} --colInA2 ${OA_out}
" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf

## Remove INDELs if specified
if [[ ${INDEL} == 0 ]]; then
  echo "CLEAN --rcdClean (${EA_out}%in%c('I','D')) | (${OA_out}%in%c('I','D'))   --strCleanName numDrop_INDEL --blnWriteCleaned 1" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
fi

echo "
####################
## 4. Harmonization of marker names (compile 'cptid')
" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf

if [[ -z $CHR && -z $BP ]]; then
  echo "CREATECPTID --fileMap ${cptref}
  --colMapMarker ${cptref_marker}
  --colMapChr ${cptref_chr}
  --colMapPos ${cptref_bp}
  --colInMarker ${SNPID_out}
  --colInA1 ${EA_out}
  --colInA2 ${OA_out}" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
else
  echo "CREATECPTID --fileMap ${cptref}
  --colMapMarker ${cptref_marker}
  --colMapChr ${cptref_chr}
  --colMapPos ${cptref_bp}
  --colInMarker ${SNPID_out}
  --colInA1 ${EA_out}
  --colInA2 ${OA_out}
  --colInChr ${CHR_out}
  --colInPos ${BP_out}" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
fi

if [[ -z $CHR && -z $BP ]]; then
  echo -e "\nSTRSPLITCOL --colSplit cptid --strSplit : --numSplitIdx 1 --colOut ${CHR_out}" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
  echo "STRSPLITCOL --colSplit cptid --strSplit : --numSplitIdx 2 --colOut ${BP_out}" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
fi

if [[ ${XY} == 0 ]]; then
  echo -e "\nCLEAN  --rcdClean !${CHR_out}%in%c(1:22,NA) --strCleanName numDropSNP_ChrXY --blnWriteCleaned 1" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
fi

#if [[ ${XY} == 0 ]]; then
#  echo -e "\nCLEAN  --rcdClean !${CHR_out}%in%c(1:22) --strCleanName numDropSNP_ChrXY --blnWriteCleaned 1" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
#fi

echo " 
####################
## 5.Filter duplicate SNPs

CLEANDUPLICATES --colInMarker cptid --strMode removeall
        

####################
## 6. AF Checks

MERGE --colInMarker cptid
  --fileRef ${afref}
  --acolIn ${afref_marker};${afref_ref};${afref_alt};${afref_raf} 
  --acolInClasses character;character;character;numeric
  --strRefSuffix .ref
  --colRefMarker cptid
  --blnWriteNotInRef 1
  --blnInAll 0
" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf


echo "ADJUSTALLELES --colInA1 ${EA_out} 
    --colInA2 ${OA_out} 
    --colInFreq ${EAF_out} 
    --colInBeta ${EFFECT_out}
    --colRefA1 ${afref_ref}.ref
    --colRefA2 ${afref_alt}.ref
    --blnMetalUseStrand 1
    --blnRemoveMismatch 1
    --blnRemoveInvalid 1
    --blnWriteMismatch 1
    --blnRemoveMismatch 1
    --blnWriteInvalid 1" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf
## All mismatches will be removed (e.g. A/T in input, A/C in reference) 

echo "
AFCHECK --colInFreq ${EAF_out}
  --colRefFreq ${afref_raf}.ref
  --numLimOutlier 0.2
  --blnPlotAll 0
" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf


echo "
####################
## 7. Rearrange columns and Write CLEANED output" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf

cols_out=$(echo -n "cptid;${SNPID_out};${CHR_out};${BP_out};${EA_out};${OA_out};${EAF_out};${EFFECT_out};${SE_out};${P_out};${N_out}")

#for col in "IMPUTED" "INFO" "HWE" "CALLRATE"; do
#  if [[ ! -z $col ]]; then
#    eval col_out='$'{${col}_out}
#    cols_out=$(echo "${cols_out};$col_out")
#  fi
#done

for i in IMPUTED INFO HWE CALLRATE; do
  eval col='$'{$i}
  if [[ ! -z $col ]]; then
    eval col_out='$'{${i}_out}
    cols_out=$(echo "$cols_out;$col_out")
  fi
done

echo -e "\nGETCOLS --acolOut ${cols_out}" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf

echo -e "\nWRITE  --strPrefix CLEANED. --strMissing NA --strMode gz" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf

echo "
####################
## 8.  Plot Z versus P

PZPLOT --colBeta ${EFFECT_out} --colSe ${SE_out} --colPval ${P_out}" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf


echo "
####################
## 9.  QQ plot

QQPLOT  --acolQQPlot ${P_out} --numPvalOffset 0.05 --strMode subplot" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf


echo "
####################
## 10. Summary Stats post-QC

CALCULATE --rcdCalc max(${N_out},na.rm=T) --strCalcName N_max

GC --colPval ${P_out} --blnSuppressCorrection 1

####################
    
STOP EASYQC" >> $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf


#################################################################################################################
#################################################################################################################

## Create R script to call EasyQC

echo "
#!/usr/bin/env Rscript
library(EasyQC)
EasyQC('$pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.ecf')" > $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.R

## Run EasyQC
Rscript --vanilla $pathOut/QC_${tag}_$(date +"%Y_%m_%d")/QC_${tag}.R