#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Add rs IDs to SNPs that show convergently high PBS in multiple populations
# ----------------------------------------------------------------------------------------

sed -e "s/\"//g" results/joint_SNP_pvalues.txt | \
    awk -F'[ ]' '{ if ($10 < 0.0001) print $0 }' | \
    while read line; do
        SNP=`echo $line | cut -d' ' -f 2`
        CHR=`echo $SNP  | cut -d":" -f 1`
        POS=`echo $SNP  | cut -d":" -f 2`

        RSID=`grep -P "^$CHR\t$POS\t" data/AGRHUM_EASTERN_100x267251.vcf | cut -f 3`

        echo -e "$line\t$RSID"
    done
