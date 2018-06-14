#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Get random subset of unrelated individuals from 1000G comparison populations
# ----------------------------------------------------------------------------------------

mkdir -p data/1000genomes/BAMs

# British in England and Scotland
awk '{ if ($8 == "unrel" && $9 == 0 && $7 == "GBR") print $2 }' \
    data/1000genomes/integrated_call_samples.20130502.seq.ped | shuf | \
    head -n 30 > data/1000genomes/subset_GBR.txt

# Luhya in Webuye, Kenya
awk '{ if ($8 == "unrel" && $9 == 0 && $7 == "LWK") print $2 }' \
    data/1000genomes/integrated_call_samples.20130502.seq.ped | shuf | \
    head -n 30 > data/1000genomes/subset_LWK.txt

# --- Make file for each individual so Snakemake easily knows what needs processing

while read ind; do
    touch data/1000genomes/BAMs/${ind}.txt
done < data/1000genomes/subset_GBR.txt

while read ind; do
    touch data/1000genomes/BAMs/${ind}.txt
done < data/1000genomes/subset_LWK.txt

exit
