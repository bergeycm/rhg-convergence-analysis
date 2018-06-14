#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Run ADMIXTURE on African and Asian datasets
# ----------------------------------------------------------------------------------------

module load vcftools
module load tabix
module load plink/1.90b3.40
module load admixture/1.3.0

# ----------------------------------------------------------------------------------------

VCF_AFR=data/AGRHUM_EASTERN_100x267251.GBR.vcf
VCF_ASN=data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.LWK.vcf

mkdir -p results/admixture

for VCF in $VCF_AFR $VCF_ASN; do

    echo "Processing [$VCF]..."

    bgzip -c $VCF > $VCF.gz
    tabix -p vcf $VCF.gz

    plink --vcf $VCF.gz \
        --geno 0.999 \
        --make-bed --out $VCF

    for k in `seq 2 6`; do
        echo "Running ADMIXTURE with k=$k..."
        admixture -j8 $VCF.bed $k
        for let in P Q; do
            mv `basename $VCF`.$k.$let results/admixture/
        done
    done
done
