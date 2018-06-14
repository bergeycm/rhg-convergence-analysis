#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Make all gene lists - N vs S lists
# ----------------------------------------------------------------------------------------

stats_files=(results/eAGR_eRHG.weir.fst
             results/ihs.twa.txt
             results/ihs.kiga.txt
             results/pbs.Batwa.GBR
             results/pbs.Bakiga.GBR
             results/pbs.Anda.LWK
             results/pbs.BIR.LWK
             results/pbs.BR1.LWK
             results/bayenv_BF_pass
            )

for stat_file in ${stats_files[@]}; do
	for group in syn nonsyn; do
        BED=${stat_file}.anno.${group}.bed
        echo "--- Processing $stat_file ---"
        perl scripts/make_gene_lists.pl $BED
    done
done

exit
