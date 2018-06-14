#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Plot genotype frequencies for all Bayenv outliers
# ----------------------------------------------------------------------------------------

module load python/3.3.2
module load R

LOCI=($(sort -n -r -k 4 results/bayenv_BF_pass.txt | head -n 10 | cut -f1-2 | \
    sed -e "s/chr\(.*\)\t/\1_/"))

for LOC in ${LOCI[*]}; do

    echo "Plotting genotype frequency for chr${LOC/_/:}..."

    Rscript scripts/plot_SNP_genotype_freqs_by_pop.R ${LOC/_/ } \
        data/all.pops.merged.clean.recode.vcf data/all.pop.list.txt
done

# Touch sentinel file to show Makefile that this script was run
touch results/bayenv_outlier_plotting_sentinel.txt
