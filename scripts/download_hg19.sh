#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download hg19
# ----------------------------------------------------------------------------------------

mkdir -p genomes/hg19
cd genomes/hg19
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar -zxvf chromFa.tar.gz
rm *_*.fa chrM.fa
cat `ls -v chr*.fa` > hg19.fa
cd ../..

sh scripts/make_seq_dict.sh

exit
