#!/bin/bash

# ========================================================================================
# --- Split per-SNP stats into nonsynonymous and synonymous datasets
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Make BEDs of synonymous and nonsynonymous SNP sites
# ----------------------------------------------------------------------------------------

# Also increment positions by 1 in BED file created by BEDOPS to make sure they match
# the coordinates in the SNP-based stats files (and the original VCF). Essentially,
# undo the automatic conversion to 0-based BED that BEDOPS does:

#  "[vcf2bed] converts 1-based, closed [a, b] VCF v4 data from standard input
#  into 0-based, half-open [a-1, b) extended BED, sent to standard output."

module load bedops/2.4.28

for class in syn nonsyn; do
    for PREFIX in AGRHUM_EASTERN_100x267251 GreatApe_Indiafinal_Austosome.sm.recode; do
        VCF=results/$PREFIX.$class.vcf.gz
        gunzip -c $VCF | vcf2bed --snvs | \
            awk -v OFS='\t' '{ $2++; $3++; print $0 }' > ${VCF/\.vcf\.gz/\.bed}
    done
done

# ----------------------------------------------------------------------------------------
# --- Split SNP-based stats files into nonsynonymous and synonymous sites
# ----------------------------------------------------------------------------------------

module load bedtools/2.26.0

# Batwa and Bakiga datasets

stats_files_afr=(
                 results/eAGR_eRHG.weir.fst
                 results/ihs.twa.txt
                 results/ihs.kiga.txt
                 results/pbs.Batwa.GBR
                 results/pbs.Bakiga.GBR
                )

for stat_file in ${stats_files_afr[@]}; do
    BED=$stat_file.anno.bed
    echo "--- Processing $BED ---"

    for class in syn nonsyn; do
        NS_BED=results/AGRHUM_EASTERN_100x267251.${class}.bed
        OUT_BED=${BED/.bed}.$class.bed
        bedtools intersect -a $BED -b $NS_BED > $OUT_BED
    done
done

# ----------------------------------------------------------------------------------------

# Mondal et al stats

stats_files_asn=(
                 results/pbs.Anda.LWK
                 results/pbs.BIR.LWK
                 results/pbs.BR1.LWK
                )

for stat_file in ${stats_files_asn[@]}; do
    BED=$stat_file.anno.bed
    echo "--- Processing $BED ---"

    for class in syn nonsyn; do
        NS_BED=results/GreatApe_Indiafinal_Austosome.sm.recode.${class}.bed
        OUT_BED=${BED/.bed}.$class.bed
        bedtools intersect -a $BED -b $NS_BED > $OUT_BED
    done
done

# ----------------------------------------------------------------------------------------

# Both together

stats_files_both=(
                 results/bayenv_BF_pass
                )

for stat_file in ${stats_files_both[@]}; do
    BED=$stat_file.anno.bed
    echo "--- Processing $BED ---"

    for class in syn nonsyn; do
        cat results/AGRHUM_EASTERN_100x267251.${class}.bed \
            results/GreatApe_Indiafinal_Austosome.sm.recode.${class}.bed | \
            cut -f 1-3 | sort | uniq > results/temp.both.$class.bed

        OUT_BED=${BED/.bed}.$class.bed
        bedtools intersect -a $BED -b results/temp.both.$class.bed > $OUT_BED
        rm results/temp.both.$class.bed
    done
done
