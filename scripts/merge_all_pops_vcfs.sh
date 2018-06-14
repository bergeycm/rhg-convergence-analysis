#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Merge Batwa/Bakiga and Anda/BR1 SNP sets and convert to Bayenv format
# ----------------------------------------------------------------------------------------

BATWA_ETAL_VCF=data/AGRHUM_EASTERN_100x267251.vcf.gz
ANDA_ETAL_VCF=data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.vcf.gz

module load vcftools

# --- Merge VCFs

vcf-merge $BATWA_ETAL_VCF $ANDA_ETAL_VCF | bgzip -c > data/all.pops.merged.vcf.gz

# --- Filter for missingness, biallelic only

vcftools --gzvcf data/all.pops.merged.vcf.gz \
    --min-alleles 2 --max-alleles 2 \
    --max-missing 0.5 \
    --recode --out data/all.pops.merged.clean

# head -n 100000 data/all.pops.merged.clean.recode.vcf > mini.vcf
# vcftools --vcf mini.vcf --max-missing 0.5 --recode --out mini.clean

# --- Prepare to convert to Bayenv format

wget https://raw.githubusercontent.com/chrishah/genomisc/master/popogeno/vcf_2_div.py \
    -O scripts/vcf_2_div.py

module load python/2.7.8
export PYTHONPATH=$PYTHONPATH:/storage/home/cxb585/local_python/lib/python2.7/site-packages/

# --- Make population list

POP=data/all.pop.list.txt

awk 'BEGIN { OFS="\t" } { print $1,"ANDA" }'   data/Mondal_indiv_list_Anda.txt > $POP
awk 'BEGIN { OFS="\t" } { print $1,"BR1" }'    data/Mondal_indiv_list_BR1.txt >> $POP
awk 'BEGIN { OFS="\t" } { print $1,"BATWA" }'  data/eRHG_IDs.txt >> $POP
awk 'BEGIN { OFS="\t" } { print $1,"BAKIGA" }' data/eAGR_IDs.txt >> $POP

# --- Do conversion

python scripts/vcf_2_div.py \
    -m $POP \
    --exclude_singletons \
    -v -b \
    -o data/all.pops.merged.cleaned \
    data/all.pops.merged.clean.recode.vcf

# Writes
# data/all.pops.merged.cleaned.bayenv.SNPfile
# data/all.pops.merged.cleaned.bayenv.SNPfile

#  -m <FILE>, --popmap <FILE>
#                        Tab delimited text file to assign individuals to
#                        populations. col1 = population ID, col2 = individual
#                        (as in vcf).
#  --exclude_singletons  exclude any singleton loci, i.e. minor allele observed
#                        only once
#  -b, --bayenv          output bayenv format (default)
#  -o <STR>, --out <STR>
#                        prefix for output files (default: out)

# --- Create thinned file for covariation matrix

vcftools --vcf data/all.pops.merged.clean.recode.vcf \
    --thin 500000 --recode --out data/all.pops.merged.clean.thin

python scripts/vcf_2_div.py \
    -m $POP \
    --exclude_singletons \
    -v -b \
    -o data/all.pops.merged.cleaned.thin \
    data/all.pops.merged.clean.thin.recode.vcf

exit
