#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Check if IMPUTE subset worked and create dummy file if not
# ----------------------------------------------------------------------------------------

CHR=$1
PART=$2

echo "Checking part $PART of chr$CHR..."

PREFIX=results/impute/AGRHUM_EASTERN_100x267251.1M.$CHR.1000g.gen.impute2_pt$PART

# Check to see if IMPUTE run worked

if [[ -s $PREFIX ]]; then
    echo -e "\tIMPUTE output exists and is non-zero"
    exit
else
    if [[ -e $PREFIX ]]; then
        echo -e "\tIMPUTE output exists BUT is non-zero. Created by this script?"
        exit
    fi
fi

if [[ -s ${PREFIX}_summary ]]; then

    NO_SNP_COUNT=`grep -c -e "not be any SNPs" \
        -e "no SNPs in the imputation" \
        -e "no type 2 SNPs" ${PREFIX}_summary`

    if [[ $NO_SNP_COUNT -gt 0 ]]; then
        echo -e "\tIMPUTE ran, but no SNPs present"
        echo -e "\tCreating empty output file"
        touch $PREFIX
    else
        echo -e "\tIMPUTE ran, but failed for a reason other than no SNPs"
    fi

else
    echo -e "\tIMPUTE did not run"
    exit
fi
