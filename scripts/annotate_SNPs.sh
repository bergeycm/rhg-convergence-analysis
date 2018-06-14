#!/bin/bash

# ========================================================================================
# --- Annotate VCF with ANNOVAR
# ========================================================================================

module load tabix

# --- ANNOVAR databases built with these commands:

#  cd /gpfs/cyberstar/ghp3/Bergey/annovar
#
#  perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
#  perl annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
#  perl annotate_variation.pl -buildver hg19 -downdb genomicSuperDups humandb/
#  perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all humandb/
#  perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2014oct humandb/
#  perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar snp138 humandb/
#  perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ljb26_all humandb/

# ----------------------------------------------------------------------------------------
# --- Annotate VCF
# ----------------------------------------------------------------------------------------

ANNO_DIR=$GRP/annovar

for VCF_PRE in data/AGRHUM_EASTERN_100x267251 data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode; do

    VCF_RES=results/`basename $VCF_PRE`

    # Temporarily decompress VCF
    if [[ ! -e ${VCF_PRE}.vcf ]]; then
        gunzip -c ${VCF_PRE}.vcf.gz > ${VCF_PRE}.vcf
    fi

    PROTOCOL_STR="refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,"
    PROTOCOL_STR="${PROTOCOL_STR}1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,"
    PROTOCOL_STR="${PROTOCOL_STR}snp138,ljb26_all"

    perl ${ANNO_DIR}/table_annovar.pl \
        --thread 8 \
        ${VCF_PRE}.vcf \
        ${ANNO_DIR}/humandb/ -buildver hg19 \
        -out ${VCF_RES} -remove \
        -protocol $PROTOCOL_STR \
        -operation g,r,r,f,f,f,f,f,f,f \
        -nastring . -vcfinput

    ANNO_VCF=$VCF_RES.hg19_multianno.vcf

    # ----------------------------------------------------------------------------------------
    # --- Tally up SNPs in each class and by Polyphen prediction
    # ----------------------------------------------------------------------------------------

    grep -o "Func.refGene=[^;]\+" $ANNO_VCF | sed -e "s/Func\.refGene=//" | \
        grep -v "\." | sort | uniq -c > reports/`basename $VCF_PRE`_SNP_class_counts.txt

    grep -o "Polyphen2_HDIV_pred=[^;]\+" $ANNO_VCF | sed -e "s/Polyphen2_HDIV_pred=//" | \
        grep -v "\." | sort | uniq -c > reports/`basename $VCF_PRE`_Polyphen_prediction_counts.txt

    # ----------------------------------------------------------------------------------------
    # --- Create separate files of synonymous and nonsynonymous SNPs
    # ----------------------------------------------------------------------------------------

    SNPSIFT_JAR=$HOME/bin/snpEff/SnpSift.jar

    # Tally up exonic SNPs by their function
    cat $ANNO_VCF | java -jar $SNPSIFT_JAR filter "( Func.refGene = 'exonic')" | \
        grep -v "^#" | awk -F';' '{ print $6 }' | sort | uniq -c

    # Make separate VCFs of synonymous and nonsynonymous SNPs
    for func in syn nonsyn; do
        cat $ANNO_VCF | java -jar $SNPSIFT_JAR filter \
            "( Func.refGene = 'exonic' & ExonicFunc.refGene = '${func}onymous_SNV')" |
            bgzip > $VCF_RES.$func.vcf.gz
    done

    # Make list of genes with stopgain, just for kicks
    cat $ANNO_VCF | java -jar $SNPSIFT_JAR filter "( ExonicFunc.refGene = 'stopgain')" | \
        grep -v "^#" | awk -F';' '{ print $4}' | cut -d'=' -f 2 > results/`basename $VCF_PRE`_stopgain.genes.txt

    # ----------------------------------------------------------------------------------------
    # --- Compress and clean up
    # ----------------------------------------------------------------------------------------

    rm ${VCF_PRE}.vcf
    rm $VCF_RES.avinput
    rm $VCF_RES.hg19_multianno.txt

done
