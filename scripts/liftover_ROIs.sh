#!/bin/bash

# ========================================================================================
# --- LiftOver regions of interest from hg18 to hg19 coordinates
# ========================================================================================

HG18_TO_HG19_CHAIN=~/liftover_chain_files/hg18ToHg19.over.chain

# ----------------------------------------------------------------------------------------
# --- Perry et al 2014
# ----------------------------------------------------------------------------------------

PERRY_PREFIX=data/Perry_etal_2014/

IN_BEDS=(
    GWAS_assoc_regions/Perry_etal_2014-S2-batwa_assoc.bed
    BayeScan_outlier_regions/Perry_etal_2014-S3-baka-vs-nzebi-bzime_BayeScan_outliers.bed
    BayeScan_outlier_regions/Perry_etal_2014-S3-batwa-vs-bakiga_BayeScan_outliers.bed
    iHS_outlier_regions/Perry_etal_2014-S3-baka_iHS_outliers.bed
    iHS_outlier_regions/Perry_etal_2014-S3-bakiga_iHS_outliers.bed
    iHS_outlier_regions/Perry_etal_2014-S3-batwa_iHS_outliers.bed
    iHS_outlier_regions/Perry_etal_2014-S3-nzebi-nzime_iHS_outliers.bed
)

for IN_BED in ${IN_BEDS[*]}; do
    ~/bin/liftOver \
        ${PERRY_PREFIX}${IN_BED} \
        $HG18_TO_HG19_CHAIN \
        ${PERRY_PREFIX}${IN_BED/.bed/.hg19.bed} \
        ${PERRY_PREFIX}${IN_BED/.bed/.hg19.unmapped.bed}
done
