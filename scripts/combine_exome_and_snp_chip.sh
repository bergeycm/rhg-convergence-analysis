#!/bin/bash

# ========================================================================================
# --- Combine SNPs from exomes and 1M SNP chip
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Convert exome VCF to PLINK format
# ----------------------------------------------------------------------------------------

PREFIX_EX=results/AGRHUM_EASTERN_100x267251

plink --vcf data/AGRHUM_EASTERN_100x267251.vcf.gz \
    --recode --out ${PREFIX_EX}

# ----------------------------------------------------------------------------------------
# --- Convert 1M SNPs from binary PED to PED and then LiftOver to hg19 coordinates
# ----------------------------------------------------------------------------------------

PREFIX_1M=data/Perry_etal_2014-PNAS-1M/Batwa_Kiga.913651pos.230samples.PNAS2014

# Convert from BED to PED
plink --bfile $PREFIX_1M --recode --out $PREFIX_1M

# LiftOver to hg19 coordinates
 python scripts/LiftMap.py \
     -m $PREFIX_1M.map -p $PREFIX_1M.ped -o $PREFIX_1M.hg19 \
     -c /storage/home/cxb585/liftover_chain_files/hg18ToHg19.over.chain

# Clean-up
rm $PREFIX_1M.hg19.bed*

# ----------------------------------------------------------------------------------------
# --- Remove rs IDs to ease merging process
# ----------------------------------------------------------------------------------------

# Backup map file and then replace SNP ID with fake ID that redundantly includes position
mv $PREFIX_1M.hg19.map{,_BACKUP}
awk 'BEGIN {OFS="\t"} { print $1,$1":"$4,$3,$4}' $PREFIX_1M.hg19.map_BACKUP > \
    $PREFIX_1M.hg19.map

# Do the same for the exome data
mv ${PREFIX_EX}.map{,_BACKUP}
awk 'BEGIN {OFS="\t"} { print $1,$1":"$4,$3,$4}' ${PREFIX_EX}.map_BACKUP > \
    ${PREFIX_EX}.map

# ----------------------------------------------------------------------------------------
# --- Convert to binary format
# ----------------------------------------------------------------------------------------

plink --file $PREFIX_1M.hg19 --make-bed --out $PREFIX_1M.hg19
plink --file $PREFIX_EX      --make-bed --out $PREFIX_EX

# ----------------------------------------------------------------------------------------
# --- Try merge to find "multiallelic" sites that may be strand errors
# ----------------------------------------------------------------------------------------

plink --bfile $PREFIX_EX \
    --bmerge $PREFIX_1M.hg19 \
    --merge-mode 1 \
    --out $PREFIX_EX.1M

# ----------------------------------------------------------------------------------------
# --- Flip to see if multiallelic variants are just strand errors
# ----------------------------------------------------------------------------------------

plink --file $PREFIX_EX --flip $PREFIX_EX.1M.missnp --make-bed --out $PREFIX_EX.flip

# ----------------------------------------------------------------------------------------
# --- Try merge again now that multiallelic SNPs have been flipped
# ----------------------------------------------------------------------------------------

plink --bfile $PREFIX_EX.flip \
    --bmerge $PREFIX_1M.hg19 \
    --merge-mode 1 \
    --out $PREFIX_EX.1M

# ----------------------------------------------------------------------------------------
# --- Remove "true" multiallelic sites, edit sample IDs, and finally do merge
# ----------------------------------------------------------------------------------------

plink --bfile $PREFIX_EX.flip --exclude $PREFIX_EX.1M.missnp \
    --make-bed --out $PREFIX_EX.flip_tmp
plink --bfile $PREFIX_1M.hg19 --exclude $PREFIX_EX.1M.missnp \
    --make-bed --out $PREFIX_1M.hg19_tmp

mv $PREFIX_EX.flip_tmp.fam{,_BACKUP}
awk '{ $1=toupper($1); gsub(/.*TWA/, "twa", $1); \
    gsub(/.*KIGA/, "kiga", $1); print $1,$1,$3,$4,$5,$6 }' \
    $PREFIX_EX.flip_tmp.fam_BACKUP > $PREFIX_EX.flip_tmp.fam

# Finally do merge
plink --bfile $PREFIX_EX.flip_tmp --bmerge $PREFIX_1M.hg19_tmp \
    --make-bed --out $PREFIX_EX.1M

rm $PREFIX_EX.flip_tmp.*
rm $PREFIX_1M.hg19_tmp.*

# ----------------------------------------------------------------------------------------
# --- Convert to input format for IMPUTE2 and split by chromosome
# ----------------------------------------------------------------------------------------

mkdir -p results/impute
for chr in `seq 1 22`; do
    OUT_FILE=results/impute/AGRHUM_EASTERN_100x267251.1M.${chr}
    plink --bfile $PREFIX_EX.1M --chr $chr --recode oxford --out $OUT_FILE
    # Clean up
    rm $OUT_FILE.log $OUT_FILE.nosex
done

#    # ----------------------------------------------------------------------------------------
#    # --- Add back in SNP name (rs ID)
#    # ----------------------------------------------------------------------------------------
#
#    # See: https://www.biostars.org/p/171557/
#
#    # Get info (chrom, chromEnd, name) by selecting fields from hg19.snp147 table of UCSC
#    # Settings are:
#    #  clade:    Mammal
#    #  genome:   Human
#    #  assembly: hg19
#    #  group:    Variation
#    #  track:    All SNPs (147)
#    #  table:    All SNPs
#    #  region:   Genome
#
#    UCSC_DAT=data/UCSC-hg19-snp147.txt.gz
#
#    # Reformat for PLINK
#    gunzip -c $UCSC_DAT | grep -v "^#" | \
#        awk 'BEGIN {OFS="\t"} { sub(/chr/,"",$1); print $1":"$2,$3 }' > \
#        ${UCSC_DAT/.txt.gz/.simple.txt}
#
#    # Find SNPs with multiple rs IDs
#    cut -f1 ${UCSC_DAT/.txt.gz/.simple.txt} | sort | uniq -c | \
#        grep -v '^ *1 '> results/duplicate_rs_SNPs.txt
#
#    # Change SNP name
#    plink --bfile $PREFIX_EX.1M --update-name ${UCSC_DAT/.txt.gz/.simple.txt} \
#        --make-bed --out $PREFIX_EX.1M.rs
