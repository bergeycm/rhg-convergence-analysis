#!/bin/bash

# ========================================================================================
# --- Download Andamanese data and outside populations
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Download Andamanese data from Mondal et al 2016 (only available as VCF)
# ----------------------------------------------------------------------------------------

# Two Andamanese populations: Jarawa (4 individuals) and Onge (6 individuals)
# ~15x coverage

mkdir -p data/Mondal_etal_2016

# Skipping:
#wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ126/ERZ126015/Andamanese.vcf \
#    -O data/Mondal_etal_2016/Andamanese.vcf

# ----------------------------------------------------------------------------------------
# --- Download data from Mondal rebuttal paper (VCF via Dropbox)
# ----------------------------------------------------------------------------------------

ANDA_VCF=data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.vcf.gz

wget https://www.dropbox.com/s/rk6jjjocefhwns1/GreatApe_Indiafinal_Austosome.vcf.gz \
    -O $ANDA_VCF

# From the Simons Genome Diversity Project:
# China:          DAI
# France:         FRN
# Chinese:        HAN
# Mandenka:       MAD
# Mbuti:          MBT
# Papuan:         PAP
# San:            SAN
# Sardinia:       SAR
# Yoruba:         YRI

# Indian populations:
# Birhor:         BIR-08 BIR-11 BIR-12 BIR-13 BIR-14 BIR-15 BIR-16 BIR-18 BIR-22
# Brahmin  (?):   BR1-1 BR1-10 BR1-12 BR1-13 BR1-16 BR1-18 BR1-3 BR1-4 BR1-5 BR1-6
# Irula    (?):   IL-01 IL-04 IL-07 IL-08 IL-09 IL-10 IL-11 IL-12 IL-15 IL-16
# Punjabi  (???): PRBB-1 PRBB-2
# Rajput   (???): RAI-11 RAI-12 RAI-13 RAI-15 RAI-2 RAI-3 RAI-4 RAI-5 RAI-8 RAI-9
# Riang:          RIA-12 RIA-14 RIA-19 RIA-20 RIA-23 RIA-24 RIA-41 RIA-42 RIA-45 RIA-63
# Vellalar (???): VV1-01 VV1-03 VV1-04 VV1-05 VV1-06 VV1-07 VV1-14 VV1-16 VV1-20

# Andamanese populations:
# Jarawa:         JAR-27 JAR-32 JAR-54 JAR-61
# Onge:           ONG-1 ONG-12 ONG-14 ONG-4 ONG-8 ONG-9

# ----------------------------------------------------------------------------------------
# --- Reduce to just pops of interest and SNPs
# ----------------------------------------------------------------------------------------

IND_LIST=data/Mondal_indiv_list.txt
echo -e "BIR-08\nBIR-11\nBIR-12\nBIR-13\nBIR-14"      > $IND_LIST
echo -e "BIR-15\nBIR-16\nBIR-18\nBIR-22"             >> $IND_LIST
echo -e "BR1-1\nBR1-10\nBR1-12\nBR1-13\nBR1-16"      >> $IND_LIST
echo -e "BR1-18\nBR1-3\nBR1-4\nBR1-5\nBR1-6"         >> $IND_LIST
echo -e "JAR-27\nJAR-32\nJAR-54\nJAR-61"             >> $IND_LIST
echo -e "ONG-1\nONG-12\nONG-14\nONG-4\nONG-8\nONG-9" >> $IND_LIST

grep "^BIR" $IND_LIST > data/Mondal_indiv_list_BIR.txt
grep "^BR1" $IND_LIST > data/Mondal_indiv_list_BR1.txt
grep "^JAR" $IND_LIST > data/Mondal_indiv_list_Anda.txt
grep "^ONG" $IND_LIST >> data/Mondal_indiv_list_Anda.txt

vcftools --gzvcf $ANDA_VCF \
    --recode \
    --out ${ANDA_VCF/\.vcf\.gz/.sm} \
    --keep $IND_LIST --mac 1

module load tabix

bgzip -c data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.vcf \
    > data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.vcf.gz

tabix -p vcf data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.vcf.gz

# ----------------------------------------------------------------------------------------
# --- Make ROD file with SNP locations
# ----------------------------------------------------------------------------------------

# While we're at it, make VCF file with Jarawa/Onge/Birhor SNP locations
# to serve as ROD file for GATK
sed -e "s/^chr//" ${ANDA_VCF/\.vcf\.gz/.sm.recode.vcf} | \
    cut -f1-10 > ${ANDA_VCF/\.vcf\.gz/.sm.recode.SNP.locs.vcf}
