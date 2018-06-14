#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Compute pi within populations in 1000 Genomes data
# ----------------------------------------------------------------------------------------

declare -A pop_ind_lists
pop_ind_lists=( ["GBR"]="data/1000genomes/subset_GBR.txt"
                ["LWK"]="data/1000genomes/subset_LWK.txt"
                ["Batwa"]="data/eRHG_IDs.txt"
                ["Bakiga"]="data/eRHG_IDs.txt"
                ["Anda"]="data/Mondal_indiv_list_Anda.txt"
                ["BIR"]="data/Mondal_indiv_list_BIR.txt"
                ["BR1"]="data/Mondal_indiv_list_BR1.txt" )

LWK_VCF=data/Mondal_etal_2016/GreatApe_Indiafinal_Austosome.sm.recode.LWK.vcf
GBR_VCF=data/AGRHUM_EASTERN_100x267251.GBR.vcf

declare -A pop_vcf_files
pop_vcf_files=( ["GBR"]=$GBR_VCF
                ["LWK"]=$LWK_VCF
                ["Batwa"]=$GBR_VCF
                ["Bakiga"]=$GBR_VCF
                ["Anda"]=$LWK_VCF
                ["BIR"]=$LWK_VCF
                ["BR1"]=$LWK_VCF )

for pop in ${!pop_ind_lists[@]}; do

    # Compute pi within pop

    IN_VCF=${pop_vcf_files[$pop]}
    IND_LIST=${pop_ind_lists[$pop]}

    PI_CMD="vcftools --vcf $IN_VCF \
        --keep $IND_LIST \
        --site-pi \
        --out results/pi_$pop"

    `$PI_CMD`

done

# --- Compute pi for full triad (pop of interest, sister pop, and outgroup)

# Bakiga, Batwa, GBR:

TMP_IND_LIST=`mktemp tmp.ind_list.XXXXX`

cat ${pop_ind_lists["Batwa"]} ${pop_ind_lists["Bakiga"]} ${pop_ind_lists["GBR"]} > \
    $TMP_IND_LIST

PI_CMD="vcftools --vcf ${pop_vcf_files["GBR"]} \
        --keep $TMP_IND_LIST \
        --site-pi \
        --out results/pi_GBR_all_pops"
`$PI_CMD`

rm $TMP_IND_LIST

# Anda, BIR, LWK:

TMP_IND_LIST=`mktemp tmp.ind_list.XXXXX`

cat ${pop_ind_lists["Anda"]} ${pop_ind_lists["BIR"]} ${pop_ind_lists["LWK"]} > \
    $TMP_IND_LIST

PI_CMD="vcftools --vcf ${pop_vcf_files["LWK"]} \
        --keep $TMP_IND_LIST \
        --site-pi \
        --out results/pi_LWK_all_pops_BIR"
`$PI_CMD`

rm $TMP_IND_LIST

# Anda, BR1, LWK:

TMP_IND_LIST=`mktemp tmp.ind_list.XXXXX`

cat ${pop_ind_lists["Anda"]} ${pop_ind_lists["BR1"]} ${pop_ind_lists["LWK"]} > \
    $TMP_IND_LIST

PI_CMD="vcftools --vcf ${pop_vcf_files["LWK"]} \
        --keep $TMP_IND_LIST \
        --site-pi \
        --out results/pi_LWK_all_pops_BR1"
`$PI_CMD`

rm $TMP_IND_LIST


exit
