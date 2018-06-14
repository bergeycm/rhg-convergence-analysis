#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download and convert gene set annotations
# ----------------------------------------------------------------------------------------

mkdir -p refGene
cd refGene

# Download refGene file
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz
gzip -d refGene.txt.gz

# Get rid of bin column using genePredToGtf
cut -f 2- refGene.txt | genePredToGtf file stdin refGene.gtf

# Sort (and remove non-autosome stuff)
grep -v "^chr[0-9]*_" refGene.gtf | \
    grep -v "^chr[XYMU]" | sed -e "s/^chr//" | \
    sort -k1,1n -k4,4n | sed -e "s/^/chr/" > refGene.sort.gtf

# Create simple version of refGene that just has gene and exon info
awk '{ if ($3 == "exon") print $0}' refGene.sort.gtf | \
    sed -e 's/gene_id "\([^"]*\)"; transcript_id "[^"]*"; exon_number "\([^"]*\)".*$/\1_exon\2/' > \
    refGene.sort.simple.gtf

# Make simple refGene BED with gene rather than exon
sed -e "s/_exon[0-9]*//" refGene.sort.simple.gtf > \
    refGene.sort.simple.justGenes.gtf

cd ..

exit
