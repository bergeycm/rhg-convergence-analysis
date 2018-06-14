#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download reference genome used by 1000 Genomes project, make sequence dicationary,
# --- and make a ROD file of Batwa/Bakiga SNP locations
# ----------------------------------------------------------------------------------------

mkdir -p genomes/hs37d5

URL=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference
URL=$URL/phase2_reference_assembly_sequence/hs37d5.fa.gz

wget $URL -O genomes/hs37d5/hs37d5.fa.gz

gunzip genomes/hs37d5/hs37d5.fa.gz

# ----------------------------------------------------------------------------------------

# Make sequence dictionary...

module load picard/1.136

PICARD_PATH=/usr/global/picard-tools/1.99

java -jar $PICARD_PATH/CreateSequenceDictionary.jar \
    R= genomes/hs37d5/hs37d5.fa O= genomes/hs37d5/hs37d5.dict

# ...and make index

module load samtools

samtools faidx genomes/hs37d5/hs37d5.fa

# ----------------------------------------------------------------------------------------

# While we're at it, make VCF file with Batwa/Bakiga SNP locations
# to serve as ROD file for GATK
gunzip -c data/AGRHUM_EASTERN_100x267251.vcf.gz | \
    sed -e "s/^chr//" | \
    cut -f1-10 > data/AGRHUM_EASTERN_100x267251.SNP.locs.vcf

exit
