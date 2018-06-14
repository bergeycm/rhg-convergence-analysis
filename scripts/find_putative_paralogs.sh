#!/bin/bash

# ========================================================================================
# --- Find possible paralogs (gene where all individuals heterozygous)
# ========================================================================================

module load vcftools
module load R
module load bedtools

# ----------------------------------------------------------------------------------------
# --- Annotate VCF
# ----------------------------------------------------------------------------------------

ANNO_VCF=results/AGRHUM_EASTERN_100x267251.hg19_multianno.vcf

# ----------------------------------------------------------------------------------------
# --- Find SNPs where all individuals are heterozygous (there are none)
# ----------------------------------------------------------------------------------------

SNPSIFT_JAR=$HOME/bin/snpEff/SnpSift.jar

echo "# --- Sites where all Batwa individuals are HET:"
vcftools --vcf $ANNO_VCF --keep data/eRHG_IDs.iHS.VCF.txt --recode --stdout | \
    java -jar $SNPSIFT_JAR filter "( countHom() == 0 )" | \
    grep -v "^#" > results/all_het_sites.eRHG.txt

echo "# --- Sites where all Bakiga individuals are HET:"
vcftools --vcf $ANNO_VCF --keep data/eAGR_IDs.iHS.VCF.txt --recode --stdout | \
    java -jar $SNPSIFT_JAR filter "( countHom() == 0 )" | \
    grep -v "^#" > results/all_het_sites.eAGR.txt

# ----------------------------------------------------------------------------------------
# --- Compute heterozygosity and inbreeding coefficient and filter based on them
# ----------------------------------------------------------------------------------------

for POP in AGR RHG; do
    vcftools --vcf $ANNO_VCF --keep data/e${POP}_IDs.iHS.VCF.txt \
        --hardy --out results/HWE_info.e${POP}
    Rscript scripts/filter_SNPs_for_het.R results/HWE_info.e${POP}.hwe
done

# ----------------------------------------------------------------------------------------
# --- Find genes containing SNPs flagged as putative paralogs
# ----------------------------------------------------------------------------------------

REFGENE=refGene/refGene.sort.simple.justGenes.gtf

for POP in AGR RHG; do

    PUT_PARA_SNPS=results/HWE_info.e${POP}.hwe.filtered.bed

    bedtools intersect \
        -a $REFGENE -b $PUT_PARA_SNPS \
        > ${PUT_PARA_SNPS/.bed/.genes.bed}

    cut -f 9 ${PUT_PARA_SNPS/.bed/.genes.bed} | \
        sort | uniq > ${PUT_PARA_SNPS/.bed/.justGenes.txt}
done

exit
