#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download 1000 Genomes data (BAM files, aligned reads)
# ----------------------------------------------------------------------------------------

# --- Download BAM file for an individual

ind=$1

echo "--- Downloading $ind ---"
URL_PATH=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data
URL_PATH=$URL_PATH/$ind/exome_alignment/

# Figure out individual's population
pop=`grep "^$ind" data/1000genomes/integrated_call_samples.20130502.seq.ped | cut -f7`

bam=data/1000genomes/BAMs/${ind}.mapped.ILLUMINA.bwa.${pop}.exome.bam

RETRY_PARAMS="--retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 20"

# Remove files if they already exist
rm -f $bam
rm -f $bam.bai

wget $RETRY_PARAMS -r --no-parent -A '*.mapped.ILLUMINA.bwa.*.exome.*.bam' \
    $URL_PATH -O $bam
wget $RETRY_PARAMS -r --no-parent -A '*.mapped.ILLUMINA.bwa.*.exome.*.bam.bai' \
    $URL_PATH -O $bam.bai

exit
