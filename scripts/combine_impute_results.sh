#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Combine IMPUTE2 results (processed in hunks)
# ----------------------------------------------------------------------------------------

chr=$1

PREFIX=results/impute/AGRHUM_EASTERN_100x267251.1M.$chr.1000g.gen.impute2

# --- Combine main IMPUTE output
IMPUTE_OUT_FILES=($(ls -v ${PREFIX}_pt* | grep "_pt[^_]*$"))
cat `echo ${IMPUTE_OUT_FILES[*]}` > $PREFIX

# --- Combine haps files
HAPS_OUT_FILES=($(ls -v ${PREFIX}_pt*haps | grep "_pt[^_]*_haps$"))
cat `echo ${HAPS_OUT_FILES[*]}` > ${PREFIX}_haps
