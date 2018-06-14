#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Make genome sequence dictionary
# ----------------------------------------------------------------------------------------

module load picard/1.136

PICARD_PATH=/usr/global/picard-tools/1.99

java -jar $PICARD_PATH/CreateSequenceDictionary.jar \
    R= genomes/hg19/hg19.fa O= genomes/hg19/hg19.dict

# ----------------------------------------------------------------------------------------

module load samtools

samtools faidx genomes/hg19/hg19.fa

exit
