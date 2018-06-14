#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Figure out missingness by population for merged file that Bayenv uses
# ----------------------------------------------------------------------------------------

IN_VCF=data/all.pops.merged.clean.recode.vcf

POP_LISTS=(data/Mondal_indiv_list_Anda.txt
           data/Mondal_indiv_list_BR1.txt
           data/eRHG_IDs.txt
           data/eAGR_IDs.txt)

OUTPUT_MISSING=results/all.pops.merged.clean.overly.missing.txt
rm -f $OUTPUT_MISSING

for pop in ${POP_LISTS[*]}; do
    echo "Processing $pop..."
    OUT_PREFIX=${pop/data\//results/missing_report.}
    vcftools --vcf $IN_VCF --keep $pop --missing-site --out $OUT_PREFIX

    awk 'BEGIN {OFS="\t"} { if ($6 > 0.5) print $1,$2 }' $OUT_PREFIX.lmiss | \
        grep -v "^CHR" >> $OUTPUT_MISSING
done

sort $OUTPUT_MISSING | uniq > ${OUTPUT_MISSING/.txt/.uniq.txt}
