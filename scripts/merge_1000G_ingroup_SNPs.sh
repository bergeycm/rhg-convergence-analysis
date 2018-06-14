#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Merge GBR SNPs with Batwa/Bakiga SNPs (for GBR) or Indian/Andamanese SNPs (for LWK)
# ----------------------------------------------------------------------------------------

module load perl
export PERL5LIB=$HOME/bin/vcftools-vcftools-490848f/src/perl/
module load tabix

pop=$1

VCF_1KG=data/1000genomes/hs37d5_snps/hs37d5.$pop.pass.snp.vcf.gz

if [ "$pop" == "GBR" ]; then
    VCF_INGROUP=data/AGRHUM_EASTERN_100x267251.vcf.gz
elif [ "$pop" == "LWK" ]; then
    VCF_INGROUP=data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.vcf.gz
else
    echo "ERROR: Population must be either GBR or LWK. Aborting."
    exit 1
fi

tabix -p vcf $VCF_1KG

vcf-merge \
    $VCF_INGROUP \
    $VCF_1KG > \
    ${VCF_INGROUP/.vcf.gz}.$pop.TMP.vcf

# --- Remove SNPs that don't have at least 10 1kG individuals (last 30 individuals)

awk '/^#/ {
        print $0
    }
    /^[^#]/ {
        missing=0;
        for (x = NF-29; x <= NF; x++) {
            if ($x == ".") {
                missing++;
            }
        };
        if (missing <= 20) {
            print $0
        }
    }' ${VCF_INGROUP/.vcf.gz}.$pop.TMP.vcf > ${VCF_INGROUP/.vcf.gz}.$pop.vcf

rm ${VCF_INGROUP/.vcf.gz}.$pop.TMP.vcf

exit
