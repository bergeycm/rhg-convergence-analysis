#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download 1000 Genomes data (SNPs only in PED format to get individual info)
# ----------------------------------------------------------------------------------------

PREFIX_1KG=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502

mkdir -p data/1000genomes

# --- Get population information

# Get all individuals with relatedness info
POP_URL=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502
POP_URL=${POP_URL}/integrated_call_samples.20130502.ALL.ped

wget $POP_URL -O data/1000genomes/integrated_call_samples.20130502.ALL.ped

# Get just individuals in panel
POP_URL=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502
POP_URL=${POP_URL}/integrated_call_samples_v3.20130502.ALL.panel

wget $POP_URL -O data/1000genomes/integrated_call_samples.20130502.ALL.panel

# Merge them to exclude individuals not sequenced
TMP1=`mktemp`
TMP2=`mktemp`

sort -t $'\t' -k 1 data/1000genomes/integrated_call_samples.20130502.ALL.ped   > $TMP1
sort -t $'\t' -k 1 data/1000genomes/integrated_call_samples.20130502.ALL.panel > $TMP2

join --nocheck-order -t $'\t' -j1 \
    $TMP1 $TMP2 > data/1000genomes/integrated_call_samples.20130502.seq.ped

exit
