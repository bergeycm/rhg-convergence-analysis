#!/bin/bash

# ========================================================================================
# --- Call and filter 1000G alleles at Batwa/Bakiga SNP locations
# ========================================================================================

CHROM=$1
POP_1KG=$2  # GBR or LWK

GATK=$HOME/bin/
GENOME_FA=genomes/hs37d5/hs37d5.fa
GENOME_NAME=hs37d5

OUT_DIR=data/1000genomes/${GENOME_NAME}_snps
PREFIX=$OUT_DIR/${POP_1KG}.chr${CHROM}

if [ "$POP_1KG" == "GBR" ]; then
    ALLELES=data/AGRHUM_EASTERN_100x267251.SNP.locs.vcf
elif [ "$POP_1KG" == "LWK" ]; then
    ALLELES=data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.SNP.locs.vcf
else
    echo "ERROR: Population must be either GBR or LWK. Aborting."
    exit 1
fi

# ----------------------------------------------------------------------------------------
# --- Call SNPs
# ----------------------------------------------------------------------------------------

BAMS=(`ls data/1000genomes/BAMs/*.mapped.ILLUMINA.bwa.${POP_1KG}.exome.bam`)

count=0
for b in ${BAMS[*]}; do
    BAMS[$count]="-I "$b" "
    count=`expr $count + 1`
done

# Make output directory
mkdir -p $OUT_DIR

java -jar ${GATK}/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R ${GENOME_FA} \
    ${BAMS[*]} \
    -stand_call_conf 50.0 \
    -stand_emit_conf 10.0 \
    -o $PREFIX.raw.snps.indels.ROD.vcf \
    -nct 1 \
    -nt 8 \
    -L ${CHROM} \
    --genotyping_mode GENOTYPE_GIVEN_ALLELES \
    --alleles $ALLELES \
    --output_mode EMIT_ALL_SITES

# ----------------------------------------------------------------------------------------
# --- Filter SNPs
# ----------------------------------------------------------------------------------------

java -Xmx2g -jar ${GATK}/GenomeAnalysisTK.jar \
	-R ${GENOME_FA} \
	-T VariantFiltration \
    -o $PREFIX.flt.vcf \
	--variant $PREFIX.raw.snps.indels.ROD.vcf \
	--filterExpression "QD < 2.0" \
	--filterName "QDfilter" \
	--filterExpression "MQ < 40.0" \
	--filterName "MQfilter" \
	--filterExpression "FS > 60.0" \
	--filterName "FSfilter" \
	--filterExpression "HaplotypeScore > 13.0" \
	--filterName "HAPSCfilter" \
	--filterExpression "MQRankSum < -12.5" \
	--filterName "MQRSfilter" \
	--filterExpression "ReadPosRankSum < -8.0" \
	--filterName "RPRSfilter" \
	--missingValuesInExpressionsShouldEvaluateAsFailing

# Select variants with "FILTER=PASS"
java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R ${GENOME_FA} \
	--variant $PREFIX.flt.vcf \
	-select "vc.isNotFiltered()" \
	-o $PREFIX.pass.snp.vcf

# ----------------------------------------------------------------------------------------
# --- Clean up
# ----------------------------------------------------------------------------------------

rm $PREFIX.raw.snps.indels.ROD.vcf
rm $PREFIX.raw.snps.indels.ROD.vcf.idx
rm $PREFIX.flt.vcf
rm $PREFIX.flt.vcf.idx

exit
