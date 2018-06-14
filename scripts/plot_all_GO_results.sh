#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Plot all GO results
# ----------------------------------------------------------------------------------------

module load R/3.3.0

# --- PBS and Bayenv

STATS=( bayenv_BF_pass
        pbs.Anda.LWK
        pbs.Bakiga.GBR
        pbs.Batwa.GBR
        pbs.BR1.LWK
      )

for STAT in ${STATS[*]}; do
    for ONT in BP MF CC; do
        Rscript scripts/plot_GO_results.R \
            results/$STAT.pvals.$ONT.overrep.results.txt \
            reports/$STAT.pvals.$ONT.GO.overrep.pdf
    done
done
