#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Make combined set of growth-associated genes
# ----------------------------------------------------------------------------------------

grep -v "^#" data/Wood_etal_2014/OMIM_height_genes.txt > data/all_growth_genes.tmp.wood.txt
wc -l data/all_growth_genes.tmp.wood.txt

grep "MP:0005378" data/MGI/HMD_HumanPhenotype.rpt | cut -f 1 > data/all_growth_genes.tmp.mgi.txt
wc -l data/all_growth_genes.tmp.mgi.txt

grep -v "^#" data/Lui_etal_2012/growth_plate_expressed.txt > data/all_growth_genes.tmp.lui.txt
wc -l data/all_growth_genes.tmp.lui.txt

grep -v "^#" data/GO_growth_genes/GO_0040007_genes.txt | cut -f 3 > data/all_growth_genes.tmp.gogrowth.txt
wc -l data/all_growth_genes.tmp.gogrowth.txt

cat data/all_growth_genes.tmp.wood.txt \
    data/all_growth_genes.tmp.mgi.txt \
    data/all_growth_genes.tmp.lui.txt \
    data/all_growth_genes.tmp.gogrowth.txt | \
    sort | uniq > data/all_growth_genes.txt
wc -l data/all_growth_genes.txt

rm data/all_growth_genes.tmp.*.txt

# --- Also count genes in Perry et al 2014 regions

wc -l data/Perry_etal_2014/GWAS_assoc_regions/Perry_etal_2014-S2-batwa_assoc.hg19.justGenes.txt
