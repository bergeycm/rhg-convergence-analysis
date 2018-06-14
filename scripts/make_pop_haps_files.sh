#!/bin/sh

# ========================================================================================
# --- Split hap file by population
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Input files
# ----------------------------------------------------------------------------------------

chr=$1

IN_HAP=results/impute/AGRHUM_EASTERN_100x267251.1M.${chr}.1000g.gen.impute2_haps
SAMP_LIST=results/impute/AGRHUM_EASTERN_100x267251.1M.${chr}.sample

# Lists of individuals to include comes from script:
# scripts/make_pop_lists_for_iHS.R

# Figure out line numbers for samples we want to include
KIGA_INDS=(`grep -n -f data/eAGR_IDs.iHS.txt $SAMP_LIST | cut -d':' -f1`)
TWA_INDS=(` grep -n -f data/eRHG_IDs.iHS.txt $SAMP_LIST | cut -d':' -f1`)

# Subtract 3 to get to zero-indexed position in list of samples
# (so line 3 represents the 0th sample in the list)
# Then figure out nth pair of columns that goes with this nth individual, now 1-indexed
# But then add 5 to get to the 1-based column number in the hap file
# (E.g. 0th column is actually the sixth column counting starting with 1 in the haps file)

KIGA_COLS=()
for (( i = 0 ; i < ${#KIGA_INDS[@]} ; i++ )); do
	# Get zero-indexed position in sample list
    INDEX_IN_LIST=$((${KIGA_INDS[$i]} - 3))
    KIGA_COL1=$(($(($(($INDEX_IN_LIST + 1)) * 2)) - 1))
    KIGA_COL2=$(($KIGA_COL1 + 1))
    KIGA_COL1=$(($KIGA_COL1 + 5))
    KIGA_COL2=$(($KIGA_COL2 + 5))
    KIGA_COLS=( "${KIGA_COLS[@]}" $KIGA_COL1 $KIGA_COL2 )
done

# Smoosh values into comma-separated string
KIGA_COL_STR=$(printf ",%s" "${KIGA_COLS[@]}")
KIGA_COL_STR=${KIGA_COL_STR:1}

TWA_COLS=()
for (( i = 0 ; i < ${#TWA_INDS[@]} ; i++ )); do
	# Get zero-indexed position in sample list
    INDEX_IN_LIST=$((${TWA_INDS[$i]} - 3))
    TWA_COL1=$(($(($(($INDEX_IN_LIST + 1)) * 2)) - 1))
    TWA_COL2=$(($TWA_COL1 + 1))
    TWA_COL1=$(($TWA_COL1 + 5))
    TWA_COL2=$(($TWA_COL2 + 5))
    TWA_COLS=( "${TWA_COLS[@]}" $TWA_COL1 $TWA_COL2 )
done

# Smoosh values into comma-separated string
TWA_COL_STR=$(printf ",%s" "${TWA_COLS[@]}")
TWA_COL_STR=${TWA_COL_STR:1}

# Reduce haps file to only columns we want
cut -d' ' -f 1-5,$KIGA_COL_STR $IN_HAP > ${IN_HAP/_haps/_kiga_haps}
cut -d' ' -f 1-5,$TWA_COL_STR  $IN_HAP > ${IN_HAP/_haps/_twa_haps}
