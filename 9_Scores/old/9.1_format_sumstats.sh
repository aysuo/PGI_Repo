#!/bin/bash

## To do:
## Add option to provide list of SNPs to be included in the score rather than a bim file.
## Optional: Add function to make P&T scores


usage() {
echo "Usage: 
      bash format_sumstats_LDpred.sh 
       --snpidtype <rs or chrpos>
       --rsid <rsid column name>
       --chrposid <chr:bp id column name>
       --chr <chr column name>
       --bp <base pair column name>
       --effall <effect allele column name>
       --altall <alternative allele column name>
       --eaf <effect allele frequency column name>
       --effect <effect column name>
       --se <se column name>
       --efftype <effect type:LINREG,OR,LOGOR,BLUP>
       --pval <p-value column name>
       --zscore <zscore column name>
       --info <imputation accuracy column name>
       --N <sample size column name>
       --sumstats <path to sumstats>
       --out <output file prefix>"

echo "note: order of options is not important"
}

echo ""
echo "-------------------------------------------------------------"
echo ""

start=$(date +%s)

#######################################################################
############################ PARSE ARGUMENTS ##########################
#######################################################################

ARGUMENT_LIST=(
    "snpidtype"
    "rsid"
    "chrposid"
    "chr"
    "bp"
    "effall"
    "altall"
    "eaf"
    "effect"
    "se"
    "efftype"
    "pval"
    "zscore"
    "info"
    "N"
    "sumstats"
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
        --snpidtype)
            snpidtype=$2
            shift 2
            ;;
        --rsid)
            rsid=$2
            shift 2
            ;;
        --chrposid)
            chrposid=$2
            shift 2
            ;;
        --chr)
            chr=$2
            shift 2
            ;;
        --bp)
            bp=$2
            shift 2
            ;;
        --effall)
            effall=$2
            shift 2
            ;;
        --altall)
            altall=$2
            shift 2
            ;;        
        --eaf)
            eaf=$2
            shift 2
            ;;
        --effect)
            effect=$2
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
        --pval)
            pval=$2
            shift 2
            ;;
        --zscore)
            zscore=$2
            shift 2
            ;;
        --info)
            info=$2
            shift 2
            ;;
        --N)
            N=$2
            shift 2
            ;;
        --sumstats)
            sumstats=$2
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

if [[ -z $snpidtype || $snpidtype == "rs" ]]
  then
    echo "SNP identifiers are rs-numbers."
    snpidtype="rs"
  elif [[ $snpidtype == "chrpos" ]]
    then
      echo "SNP identifiers are of the type chr:pos."
    else
      echo "Error: --snpidtype cannot take values other than 'rs' and 'chrpos'."
      exit 1
fi

if [[ (-z $rsid || $rsid = "NA") && (-z $chrposid || $chrposid = "NA") ]]
  then
    echo "Error: Neither rs-id nor ChrPosID column names have been supplied."
    exit 1
  elif [[ $snpidtype == "chrpos" && ! -z $chrposid && $chrposid != "NA" ]]
    then
      echo "SNP identifier column: $chrposid."
      snpid=$chrposid
    elif [[ $snpidtype == "chrpos" && (-z $chrposid || $chrposid = "NA") ]]
      then  
        echo "$rsid column will be mapped to ChrPosID's using the reference file $snpid_map."
        snpid=$rsid
      elif [[ $snpidtype != "chrpos" && ! -z $rsid && $rsid != "NA" ]]
        then
          echo "SNP identifier column: $rsid."
          snpid=$rsid
        else 
          echo "$chrposid column will be mapped to rs-id's using the reference file $snpid_map."
          snpid=$chrposid
fi


if [[ (-z $chr || $chr = "NA") && ! -z $rsid ]]
  then
    echo "Chromosome numbers not specified. They will be obtained using rs-id's and reference file."
    chr=$snpid
  elif [[ (-z $chr || $chr = "NA") && -z $rsid ]]
    then
      echo "Chromosome numbers not specified. They will be obtained from ChrPosID's and reference file."
      chr=$snpid
    else 
      echo "Chromosome: $chr".
fi


if [[ (-z $bp || $bp = "NA") && ! -z $rsid ]]
  then
    echo "Base pair positions not specified. They will be obtained using rs-id's and reference file."
    bp=$snpid
  elif [[ (-z $bp || $bp = "NA") && -z $rsid ]]
    then
      echo "Base pair positions not specified. They will be obtained using ChrPosID's and reference file."
      bp=$snpid
    else
      echo "Base pair position: $bp".
fi


if [[ -z $effall ]]
  then
    echo "Error: No effect allele column has been supplied."
    exit 1
  else 
    echo "Effect allele: $effall."
fi

if [[ -z $altall ]]
  then
    echo "Error: No alternative allele column has been supplied."
    exit 1
  else
    echo "Alternative allele: $altall."
fi

if [[ -z $eaf ]]
  then
    echo "Error: No effect allele frequency column has been supplied."
    exit 1
  else 
    echo "Effect allele frequency: $eaf."
fi

if [[ -z $efftype ]]
  then
    efftype="LINREG"
    echo "No effect type has been specified, assuming the default (LINREG)."
  else 
    echo "Effect type: $efftype"
fi 


if [[ $efftype = "LINREG" ]]
then 
  if [[ (-z $effect || $effect = "NA") && (-z $zscore || $zscore = "NA") ]]
  then
    echo "Error: No linear regression effect or Z-score column has been supplied."
    exit 1
  elif [[ (-z $effect || $effect = "NA") && ! -z $zscore && $zscore != "NA" ]]
    then
      effect=NA
      echo "Effect will be obtained from Z-score and EAF."
  else
    zscore=NA       
    echo "Effect: $effect"
  fi
elif [[ $efftype="OR" || $efftype="LOGOR" || $efftype="BLUP" ]]
  then
    if [[ -z $effect ]]
    then
      echo "Error: No OR/LOGOR/BLUP effect column has been supplied."
      exit 1
    else
      echo "Effect: $effect"
    fi
fi

if [[ -z $se || $se = "NA" ]]
then
    se="NA"
    echo "SE will be obtained from EAF and N"
  else      
    echo "SE: $se"
fi

if [[ -z $pval ]]
  then
    echo "Error: No p-value column has been supplied."
    exit 1
  else
    echo "P-value: $pval."
fi

if [[ -z $info ]]
  then
    info="NA"
  else
    echo "Imputation accuracy: $info."
fi

if [[ -z $N ]]
  then
    echo "Error: No sample size column has been supplied."
    exit 1
  else
    echo "Sample size: $N."
fi

if [[ -z $sumstats ]]
  then
    echo "Error: No GWAS summary statistics file path has been supplied."
    exit 1
  else
    echo "GWAS sumstats file path: $sumstats."
fi

if [[ -z $out ]]
  then
    echo "Error: No output prefix has been supplied."
    exit 1
  else
    echo "Output file prefix: $out."
fi

if [[ -z $snpid_map ]]
  then
    snpid_map="/disk/genetics3/PGS/PGS_Repo/data/REFFILES/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.rsid_map"
  else
    snpid_map=$snpid_map
fi
echo "SNPID mapping file: $snpid_map"

echo ""
echo "-------------------------------------------------------------"
echo ""



#######################################################################
##################### CREATE WORKING DIRECTORIES ######################
#######################################################################

echo "Creating working directories."
echo ""

if [ -a ./sumstats ]
  then
    echo "'sumstats' directory already exists. Using the existing directory."
  else
    mkdir sumstats
fi

echo ""
echo "-------------------------------------------------------------"
echo ""


if [[ $sumstats == *.gz ]]; then
  gunzip -c $sumstats > sumstats/tmp_sumstats_${out}.txt
  sumstats="sumstats/tmp_sumstats_${out}.txt"
fi

echo "Formatting GWAS summary statistics.."
awk -F"\t" \
  -v chr=$chr \
  -v bp=$bp \
  -v altall=$altall \
  -v effall=$effall \
  -v eaf=$eaf \
  -v snpid=$snpid \
  -v pval=$pval \
  -v zscore=$zscore \
  -v effect=$effect \
  -v N=$N \
  -v se=$se \
  'NR==1{ for (i=1; i<=NF; i++) {ix[$i] = i}; print "CHR","POS","SNP_ID","REF","ALT","REF_FRQ","PVAL","BETA","SE","N"} \
  NR>1 && effect!="NA"{ print "chr"$ix[chr], $ix[bp], $ix[snpid], toupper($ix[effall]), toupper($ix[altall]), $ix[eaf], $ix[pval], $ix[effect], $ix[se], $ix[N]} \
  NR>1 && effect=="NA" && $ix[eaf]>0 && $ix[eaf]<1 { print "chr"$ix[chr], $ix[bp], $ix[snpid], toupper($ix[effall]), toupper($ix[altall]), $ix[eaf], $ix[pval], $ix[zscore]/sqrt(2*$ix[N]*$ix[eaf]*(1-$ix[eaf])), 1/sqrt(2*$ix[N]*$ix[eaf]*(1-$ix[eaf])), $ix[N]}' OFS="\t" ${sumstats} > sumstats/sumstats_${out}.txt  \
  && echo "GWAS summary statistics formatted as LDpred input." || (echo "Error: Summary statistics could not be formatted." && exit 1)


if [[ $snpidtype == "chrpos" && $snpid == $rsid ]]
then
  echo "Obtaining ChrPosID from rs-id and reference file."
  awk -F"\t" 'NR==FNR {a[$1]=$4;next}FNR==1{print;next}FNR>1{$3=a[$3];print}' OFS="\t" $snpid_map sumstats/sumstats_${out}.txt > sumstats/tmp_sumstats_${out}.txt \
    && mv sumstats/tmp_sumstats_${out}.txt sumstats/sumstats_${out}.txt \
    && echo "ChrPosID's obtained." \
    || (echo "Error: ChrPosID's could not be obtained." && exit 1)
fi

if [[ $snpidtype == "rs" && $snpid == $chrposid ]]
then
  echo "Obtaining rsID's from ChrPosID and reference file."
  awk -F"\t" 'NR==FNR {a[$4]=$1;next}FNR==1{print;next}FNR>1{$3=a[$3];print}' OFS="\t" $snpid_map sumstats/sumstats_${out}.txt > sumstats/tmp_sumstats_${out}.txt \
    && mv sumstats/tmp_sumstats_${out}.txt sumstats/sumstats_${out}.txt \
    && echo "rsID's obtained." \
    || (echo "Error: rsID's could not be obtained." && exit 1)
fi

if [[ $chr == $snpid || $bp == $snpid ]]
then
  if [[ $snpidtype == "chrpos" ]]
  then 
    echo "Obtaining chromosome and base pair positions using ChrPosID's and reference file."
    awk -F"\t" -v CHR=$chr -v BP=$bp -v SNP=$snpid \
      'NR==FNR {chr[$4]=$2; bp[$4]=$3; next} 
      FNR==1 {print;next}\
      FNR>1 && CHR==SNP && BP==SNP && $3 in chr {$1="chr"chr[$3];$2=bp[$3];print;next} \
      FNR>1 && CHR!=SNP && BP==SNP && $3 in chr {$2=bp[$3];print;next} \
      FNR>1 && CHR==SNP && BP!=SNP && $3 in chr {$1="chr"chr[$3];print;next}' OFS="\t" $snpid_map sumstats/sumstats_${out}.txt > sumstats/tmp_sumstats_${out}.txt \
      && mv sumstats/tmp_sumstats_${out}.txt sumstats/sumstats_${out}.txt \
      && echo "Chromosome and/or base pair positions obtained." \
      || (echo "Error: Chromosome and/or base pair positions could not be obtained" && exit 1)
  else 
    echo "Obtaining chromosome and base pair positions using rs-id's and reference file."
    awk -F"\t" -v CHR=$chr -v BP=$bp -v SNP=$snpid \
      'NR==FNR {chr[$1]=$2; bp[$1]=$3; next} 
      FNR==1 {print;next}\
      FNR>1 && CHR==SNP && BP==SNP && $3 in chr {$1="chr"chr[$3];$2=bp[$3];print;next} \
      FNR>1 && CHR!=SNP && BP==SNP && $3 in chr {$2=bp[$3];print;next} \
      FNR>1 && CHR==SNP && BP!=SNP && $3 in chr {$1="chr"chr[$3];print;next}' OFS="\t" $snpid_map sumstats/sumstats_${out}.txt > sumstats/tmp_sumstats_${out}.txt \
      && mv sumstats/tmp_sumstats_${out}.txt sumstats/sumstats_${out}.txt \
      && echo "Chromosome and/or base pair positions obtained." \
      || (echo "Error: Chromosome and/or base pair positions could not be obtained" && exit 1)
  fi
fi

rm -f sumstats/tmp_sumstats_${out}.txt
echo "The sumstats look like this:"
head sumstats/sumstats_${out}.txt


