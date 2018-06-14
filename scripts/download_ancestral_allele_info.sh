#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download SNP ancestral allele info
# ----------------------------------------------------------------------------------------

AA_PATH=ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/database/organism_data
AA_FILE=SNPAncestralAllele.bcp.gz
AA_URL=$AA_PATH/$AA_FILE

wget $AA_URL -O data/$AA_FILE

wget ftp://ftp.ncbi.nlm.nih.gov/snp/database/shared_data/Allele.bcp.gz \
    -O data/Allele.bcp.gz

gunzip data/$AA_FILE
gunzip data/Allele.bcp.gz

# rsID is stored as an integer in the first column of the SNPAncestralAllele.bcp table.
# Corresponding ancestral allele ID is located in the second column
# You can then look in the Allele.bcp table for the allele letter that corresponds to
# that Ancestral allele ID

# NOT USED: 1000G ancestral allele info:
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/\
#    ancestral_alignments/human_ancestor_GRCh37_e59.tar.bz2 \
#    -O data/human_ancestor_GRCh37_e59.tar.bz2

exit
