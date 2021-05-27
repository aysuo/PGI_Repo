#!/bin/bash

## Script to merge 2 tab-delimited files   

usage() {
echo "Usage: 
      merge.sh 
       --mergeType <'L' for left merge, 'R' for right merge, 'I' for inner merge, 'O' for outer merge >
       --mergeCol1 <Column number to merge on in file 1>
       --mergeCol2 <Column number to merge on in file 2>
       --file1 <path to file 1>
       --file2 <path to file 2>
       --out <output file path>"
echo ""
echo "Note: order of options is not important"
}

echo ""
echo "-------------------------------------------------------------"
echo ""


echo -n "Merge started on "
date
start=$(date +%s)

#######################################################################
############################ PARSE ARGUMENTS ##########################
#######################################################################

ARGUMENT_LIST=(
  "mergeType"
  "mergeCol1"
  "mergeCol2"
  "file1"
  "file2"
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
        --mergeType)
            mergeType=$2
            shift 2
            ;;
        --mergeCol1)
            mergeCol1=$2
            shift 2
            ;;
        --mergeCol2)
            mergeCol2=$2
            shift 2
            ;;
        --file1)
            file1=$2
            shift 2
            ;;
        --file2)
            file2=$2
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

if [[ $mergeType != "I" && $mergeType != "O" && $mergeType != "L" && $mergeType != "R" ]]; then
  echo -e "Error: Invalid merge type parameter. Please choose on of the following: \n \
  'I' for inner merge, 'L' for left merge, 'R' for right merge, 'O' for outer merge." 1>&2
  exit 1
fi


if ! [[ $mergeCol1 =~ ^[0-9]+$ && $mergeCol2 =~ ^[0-9]+$ ]]; then
  echo "Error: Please specify column numbers to merge on, not column names." 1>&2
  exit 1
fi


if ! [[ -f $file1 ]]; then
    echo "Error: Could not locate $file1" 1>&2
    exit 1
fi


if ! [[ -f $file2 ]]; then
    echo "Error: Could not locate $file2" 1>&2
    exit 1
fi


if [[ -z $out ]]; then
  echo "Error: No output path have been supplied."
  exit 1
elif [[ -s $out ]]; then
  echo -n "Output file $out already exists. Do you want to rewrite? (Y/N) : "
  read $answer
  if [[ $answer="Y" ]]; then
    echo "Output path: $out"
  else
    echo -n "Please specify a new output path:"
    read ${new_path}
    out=${new_path}   
  fi
elif [[ -d $j ]]; then
    echo "Error: The output path you have specified is a directory. Please specify a full file path."
    exit 1
else
  echo "Output path: $out"
fi


echo ""
echo "-------------------------------------------------------------"
echo ""


case "$mergeType" in
  I)
    awk -F"\t" -v m1=$mergeCol1 -v m2=$mergeCol2 'NR==FNR{a[$m1]=$0;next}($m2 in a){print a[$m2],$0}' OFS="\t" $file1 $file2 | \
    awk -F"\t" 'NR==1{Q=NF;print} NR>1{for(i=1;i<=Q;i++){if(!$i && $i!=0){$i="NA"}};print}' OFS="\t" > $out
    ;;
  L)
    awk -F"\t" -v m1=$mergeCol1 -v m2=$mergeCol2 'NR==FNR{a[$m2]=$0;next}{print $0, a[$m1]}' OFS="\t" $file2 $file1 | \
    awk -F"\t" 'NR==1{Q=NF;print} NR>1{for(i=1;i<=Q;i++){if(!$i && $i!=0){$i="NA"}};print}' OFS="\t" > $out
    ;;
  R)
    awk -F"\t" -v m1=$mergeCol1 -v m2=$mergeCol2 'NR==FNR{a[$m1]=$0;next}{print a[$m2],$0}' OFS="\t" $file1 $file2 | \
    awk -F"\t" 'NR==1{Q=NF;print} NR>1{for(i=1;i<=Q;i++){if(!$i && $i!=0){$i="NA"}};print}' OFS="\t" > $out
    ;;
  O)
    awk -F"\t" -v m1=$mergeCol1 -v m2=$mergeCol2 'NR==FNR{a[$m2]=$0;next}{print $0, a[$m1]}' OFS="\t" $file2 $file1 > TMP/mrgO1
    awk -F"\t" -v m1=$mergeCol1 -v m2=$mergeCol2 'NR==FNR{a[$m1]=$0;next}FNR>1 && !($m2 in a){print a[$m2], $0}' OFS="\t" $file1 $file2 > TMP/mrgO2
    cat TMP/mrgO1 TMP/mrgO2 | awk -F"\t" 'NR==1{Q=NF;print} NR>1{for(i=1;i<=Q;i++){if(!$i && $i!=0){$i="NA"}};print}' OFS="\t" > $out
    rm TMP/mrgO1 TMP/mrgO
    ;;
esac

echo ""
echo "-------------------------------------------------------------"
echo ""

echo -n "Merge finished running on "
date

end=$(date +%s)
echo "Execution time was $(($end - $start)) seconds."