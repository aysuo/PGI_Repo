#!/bin/bash

## Script to extract columns by name        
## Assumes tab-delimited file               
## Note: If there are multiple columns with the same name, it takes the first one.

usage(){
echo "Usage: 
      extractCols.sh 
       --file <path to file, do not specify if input comes from stdin>
       --cols <columns to keep, separated by comma> 
       --out <output file path>"
echo ""
echo "Note: order of options is not important"
}

echo ""
echo "-------------------------------------------------------------"
echo ""


echo -n "extractCols started on "
date
start=$(date +%s)

#######################################################################
############################ PARSE ARGUMENTS ##########################
#######################################################################

ARGUMENT_LIST=(
  "file"
  "cols"
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
        --file)
            file=$2
            shift 2
            ;;
        --cols)
            cols=$2
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

# Check if required arguments have been supplied"


echo ""
echo "-------------------------------------------------------------"
echo ""
echo "Checking required arguments.."
echo ""


if [[ -z $file ]]
then
  file=/dev/stdin
  if [[ -s $file ]]
  then
    echo "No input data has been supplied." 1>&2 && usage 
    exit 1
  else
    echo "Input data redirected from stdin."
  fi
elif ! [[ -f $file ]]
  then
    echo "Error: Could not locate $file" 1>&2 && usage
    exit 1
  else
    echo "Input file: $file"
fi


if [[ -z $cols ]]
then
   echo "Error: No column names have been specified." 1>&2 && usage
   exit 1
else
  colKeep=$(echo $cols | sed 's/,/ /g')
  echo "Columns to keep: $colKeep"
fi


if [[ -z $out ]]
then
  echo "Error: No output path have been supplied." 1>&2 && usage
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
    echo "Error: The output path you have specified is a directory. Please specify a full file path." 1>&2 && usage
    exit 1
else
  echo "Output path: $out"
fi


#######################################################################

echo ""
echo "-------------------------------------------------------------"
echo ""



declare -a Acols=$(echo "("$colKeep")")
ncols=${#Acols[@]}

awk -F"\t" -v cols="${Acols[*]}" 'BEGIN{ncols=split(cols,c," ")} \
    NR==1 { for (i=1; i<=NF; i++) { ix[$i] = i } } \
    NR>=1 { for (i=1; i<=ncols; i++) {ORS=i%ncols?"\t":"\n" ;  print $ix[c[i]] } }' $file > $out || echo "extractCols has not been executed successfully." exit 1

echo -n "extractCols finished running on "
date

end=$(date +%s)
echo "Execution time was $(($end - $start)) seconds."