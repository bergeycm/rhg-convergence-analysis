#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Concatenate 1000 Genomes SNPs from separate chromosomes
# ----------------------------------------------------------------------------------------

module load perl
export PERL5LIB=$HOME/bin/vcftools-vcftools-490848f/src/perl/

module load tabix

pop=$1

vcf-concat `ls -v data/1000genomes/hs37d5_snps/$pop.chr[0-9]*.pass.snp.vcf` | \
    sed -e "s/^\([0-9]\)/chr\1/" | \
    bgzip -c > data/1000genomes/hs37d5_snps/hs37d5.$pop.pass.snp.vcf.gz

exit
