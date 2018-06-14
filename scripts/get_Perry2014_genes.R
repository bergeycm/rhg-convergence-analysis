#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Find genes in Perry et al 2014 phenotype-associated regions
# ----------------------------------------------------------------------------------------

module load bedtools

REFGENE=refGene/refGene.sort.simple.justGenes.gtf

PERRY2014=data/Perry_etal_2014/GWAS_assoc_regions/Perry_etal_2014-S2-batwa_assoc.hg19.bed

bedtools intersect \
    -a $REFGENE -b $PERRY2014 \
    > ${PERRY2014/.bed/.genes.bed}

cut -f 9 ${PERRY2014/.bed/.genes.bed} | sort | uniq > ${PERRY2014/.bed/.justGenes.txt}

exit
