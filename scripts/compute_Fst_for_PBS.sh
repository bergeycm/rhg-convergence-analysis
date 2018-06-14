#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Compute Fst between populations in 1000 Genomes data for PBS
# ----------------------------------------------------------------------------------------

outpop=$1

# --- Compute Fst on a per-SNP basis

function compute_fst {

    IN_VCF=$1
    POP_OF_INTEREST=$2
    POP_SISTER_GROUP=$3
    INDS_OF_INTEREST=$4
    INDS_SISTER_GROUP=$5

    # Group of interest vs. sister group
    POPS_STR=`echo -e "$POP_OF_INTEREST\n$POP_SISTER_GROUP" | sort | paste -s -d "_"`
    FST_CMD="vcftools --vcf $IN_VCF \
        --weir-fst-pop $INDS_OF_INTEREST \
        --weir-fst-pop $INDS_SISTER_GROUP \
        --out results/PBS_$POPS_STR"

    `$FST_CMD`

    # Sister group vs. outgroup
    POPS_STR=`echo -e "$POP_SISTER_GROUP\n$outpop" | sort | paste -s -d "_"`
    FST_CMD="vcftools --vcf $IN_VCF \
        --weir-fst-pop $INDS_SISTER_GROUP \
        --weir-fst-pop $INDS_OUT_1000G \
        --out results/PBS_$POPS_STR"

    `$FST_CMD`

    # Group of interest vs. outgroup
    POPS_STR=`echo -e "$POP_OF_INTEREST\n$outpop" | sort | paste -s -d "_"`
    FST_CMD="vcftools --vcf $IN_VCF \
        --weir-fst-pop $INDS_OF_INTEREST \
        --weir-fst-pop $INDS_OUT_1000G \
        --out results/PBS_$POPS_STR"

    `$FST_CMD`
}

# ----------------------------------------------------------------------------------------

INDS_OUT_1000G=data/1000genomes/subset_${outpop}.txt

if [ "$outpop" == "GBR" ]; then

    IN_VCF=data/AGRHUM_EASTERN_100x267251.${outpop}.vcf

    POP_OF_INTEREST=Batwa
    INDS_OF_INTEREST=data/eRHG_IDs.txt

    POP_SISTER_GROUP=Bakiga
    INDS_SISTER_GROUP=data/eAGR_IDs.txt

    compute_fst $IN_VCF $POP_OF_INTEREST $POP_SISTER_GROUP \
        $INDS_OF_INTEREST $INDS_SISTER_GROUP

elif [ "$outpop" == "LWK" ]; then

    IN_VCF=data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.${outpop}.vcf

    POP_OF_INTEREST=Anda
    INDS_OF_INTEREST=data/Mondal_indiv_list_Anda.txt

    for POP_SISTER_GROUP in BIR BR1; do
        INDS_SISTER_GROUP=data/Mondal_indiv_list_${POP_SISTER_GROUP}.txt
        compute_fst $IN_VCF $POP_OF_INTEREST $POP_SISTER_GROUP \
            $INDS_OF_INTEREST $INDS_SISTER_GROUP
    done

else
    echo "ERROR: Population must be either GBR or LWK. Aborting."
    exit 1
fi

exit
