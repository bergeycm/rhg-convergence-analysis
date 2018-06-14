#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Combine BP and MF results into single files
# ----------------------------------------------------------------------------------------

for cor_type in "" "_sizeCor" "_mafCor"; do

    # Overrep results, pop-specific

    for pop in Batwa.GBR Bakiga.GBR Anda.LWK BR1.LWK; do
        cp results${cor_type}/pbs.$pop.pvals.BP.overrep.results.txt \
            results${cor_type}/pbs.$pop.pvals.BPMF.overrep.results.txt
        cat results${cor_type}/pbs.$pop.pvals.MF.overrep.results.txt | \
            grep -v "^GO\.ID" >> \
            results${cor_type}/pbs.$pop.pvals.BPMF.overrep.results.txt
    done

    # Overrep results, convergent

    for pair in Batwa.GBR-Anda.LWK Bakiga.GBR-BR1.LWK; do

        cp results${cor_type}/convergent_permuted_pval_$pair.BP.overrep.txt \
            results${cor_type}/convergent_permuted_pval_$pair.BPMF.overrep.txt
        cat results${cor_type}/convergent_permuted_pval_$pair.MF.overrep.txt | \
            grep -v "^\"GO\.ID" >> \
            results${cor_type}/convergent_permuted_pval_$pair.BPMF.overrep.txt
    done

    # ------------------------------------------------------------------------------------

    # Enrichment results, pop-specific

    for pop in Batwa.GBR Bakiga.GBR Anda.LWK BR1.LWK; do
        cp results${cor_type}/pbs.$pop.pvals.BP.enrich.results.txt \
            results${cor_type}/pbs.$pop.pvals.BPMF.enrich.results.txt
        cat results${cor_type}/pbs.$pop.pvals.MF.enrich.results.txt | \
            grep -v "^GO\.ID" >> \
            results${cor_type}/pbs.$pop.pvals.BPMF.enrich.results.txt
    done

    # Enrichment results, convergent

    for pair in Batwa.GBR-Anda.LWK Bakiga.GBR-BR1.LWK; do

        cp results${cor_type}/convergent_permuted_pval_$pair.BP.enrich.txt \
            results${cor_type}/convergent_permuted_pval_$pair.BPMF.enrich.txt
        cat results${cor_type}/convergent_permuted_pval_$pair.MF.enrich.txt | \
            grep -v "^\"GO\.ID" >> \
            results${cor_type}/convergent_permuted_pval_$pair.BPMF.enrich.txt
    done

done
