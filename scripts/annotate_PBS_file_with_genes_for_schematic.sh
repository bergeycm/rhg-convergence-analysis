#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Annotate Fst, T, and PBS file with genes for PBS schematic generation
# ----------------------------------------------------------------------------------------

PBS_GBR=results/GBR.FST_for_PBS.txt
PBS_LWK=results/LWK.FST_for_PBS.txt

# --- Annotate BED with genes

module load bedtools/2.26.0

for BED in $PBS_GBR $PBS_LWK; do

    sed -e "s/^chr//g" $BED | sort -k1,1n -k2,2n | sed -e "s/^/chr/g" > $BED.sort

    mapBed -a $BED.sort \
        -b refGene/refGene.sort.simple.justGenes.gtf \
        -c 9 -o collapse | \
    awk '{ if ($18 != ".") print $0 }' > $BED.anno.bed

    rm $BED.sort
done
