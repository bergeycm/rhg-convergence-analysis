#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Combine all iHS output files and normalize iHS
# ----------------------------------------------------------------------------------------

iHS_PREFIX=results/selscan/AGRHUM_EASTERN_100x267251.1M
NORM=/storage/home/cxb585/bin/selscan/bin/linux/norm

for pop in twa kiga; do
    cat `ls -v ${iHS_PREFIX}.*.1000g.${pop}.selscan.ihs.out` > \
        results/ihs.${pop}.txt
    $NORM --ihs --files results/ihs.${pop}.txt
    # Replace original file with noramlized one
    mv results/ihs.${pop}.txt.100bins.norm results/ihs.${pop}.txt
done
