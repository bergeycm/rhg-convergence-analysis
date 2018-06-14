#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Compute various stats by population
# ----------------------------------------------------------------------------------------

WIN_SIZE=100000

# E.g. data/AGRHUM_EASTERN_100x267251.vcf.gz
IN_VCF_GZ=$1

# Figure out output suffix
SUFFIX=$(echo $IN_VCF_GZ | grep -o "\.[^\.]*syn")

# --- Compute windowed pi

for pop in eAGR eRHG; do
    vcftools --gzvcf $IN_VCF_GZ \
        --window-pi ${WIN_SIZE} \
        --window-pi-step $WIN_SIZE \
        --keep data/${pop}_IDs.txt \
        --out results/${pop}${SUFFIX}
done

# --- Compute Tajima's D

for pop in eAGR eRHG; do
    vcftools --gzvcf $IN_VCF_GZ \
        --TajimaD $WIN_SIZE \
        --keep data/${pop}_IDs.txt \
        --out results/${pop}${SUFFIX}
done

exit
