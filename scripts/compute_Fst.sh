#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Compute Fst between two populations
# ----------------------------------------------------------------------------------------

# E.g. data/AGRHUM_EASTERN_100x267251.vcf.gz
IN_VCF_GZ=$1

# Figure out output suffix
SUFFIX=$(echo $IN_VCF_GZ | grep -o "\.[^\.]*syn")

# --- Compute Fst on a per-SNP basis

FST_CMD="vcftools --gzvcf $IN_VCF_GZ \
    --weir-fst-pop data/eAGR_IDs.txt \
    --weir-fst-pop data/eRHG_IDs.txt \
    --out results/eAGR_eRHG${SUFFIX}"

`$FST_CMD`

# --- Compute windowed Fst

WIN_SIZE=100000

`$FST_CMD --fst-window-size $WIN_SIZE --fst-window-step $WIN_SIZE`

exit
