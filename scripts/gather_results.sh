#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Gather together results for supplements
# ----------------------------------------------------------------------------------------

rm -r reports/supplementary_tables
mkdir -p reports/supplementary_tables

POP_STRS=(Batwa.GBR Bakiga.GBR Anda.LWK BR1.LWK)

# --- Table S1 - Per-Population PBS results by SNP

OUT_DIR=reports/supplementary_tables/Table_S1_per_pop_PBS_SNP
mkdir -p $OUT_DIR
for pop in ${POP_STRS[*]}; do
    cp results/pbs.$pop.bed $OUT_DIR/Table_S1_per_pop_PBS_SNP.$pop.bed
done

# --- Table S2 - Per-Population PBS results by gene

OUT_DIR=reports/supplementary_tables/Table_S2_per_pop_PBS_gene
mkdir -p $OUT_DIR
for pop in ${POP_STRS[*]}; do
    cp results/pbs.$pop.pvals.txt $OUT_DIR/Table_S2_per_pop_PBS_gene-$pop.txt
done

# --- Table S3 - Per-Population PBS results by GO term

OUT_DIR=reports/supplementary_tables/Table_S3_per_pop_PBS_GO
mkdir -p $OUT_DIR
for ont in BP MF CC; do
    for pop in ${POP_STRS[*]}; do
        cp results/pbs.$pop.pvals.$ont.results.txt \
            $OUT_DIR/Table_S3_per_pop_PBS_GO-$pop.$ont.txt
    done
done

# --- Table S4 - Per-Population PBS "N vs. S" results

# Cut

# ----------------------------------------------------------------------------------------

# --- Table S5 - Convergence by SNP

OUT_DIR=reports/supplementary_tables/Table_S5_convergence_PBS_SNP
mkdir -p $OUT_DIR
OUT_FILE=$OUT_DIR/Table_S5_convergence_PBS_SNP.joint_pvalues.txt
head -n1 results/joint_SNP_pvalues.txt | \
    sed -e "s/\.x/.pop1/g" -e "s/\.y/.pop2/g" -e "s/stat/stat.pop/g" > $OUT_FILE
sed -e "s/^\"[0-9]\+\" //" results/joint_SNP_pvalues.txt | tail -n +2 >> $OUT_FILE

# --- Table S6 - Convergence by gene

OUT_DIR=reports/supplementary_tables/Table_S6_convergence_PBS_gene
mkdir -p $OUT_DIR
OUT_FILE=$OUT_DIR/Table_S6_convergence_PBS_gene.joint_pvalues.txt
head -n1 results/joint_gene_pvalues.txt | \
    sed -e "s/\.x/.pop1/g" -e "s/\.y/.pop2/g" -e "s/stat/stat.pop/g" > $OUT_FILE
sed -e "s/^\"[0-9]\+\" //" results/joint_gene_pvalues.txt | tail -n +2 >> $OUT_FILE

# --- Table S7 - Convergence by GO term

OUT_DIR=reports/supplementary_tables/Table_S7_convergence_PBS_GO
mkdir -p $OUT_DIR
OUT_FILE=$OUT_DIR/Table_S7_convergence_PBS_GO.joint_pvalues.txt
head -n1 results/joint_GO_pvalues.txt | \
    sed -e "s/\.x/.pop1/g" -e "s/\.y/.pop2/g" -e "s/stat/stat.pop/g" > $OUT_FILE
sed -e "s/^\"[0-9]\+\" //" results/joint_GO_pvalues.txt | tail -n +2 | \
    grep "pbs" | grep -v "NvsS" >> $OUT_FILE

# --- Table S8 - Convergence in "N vs. S" statistic

# Cut

# ----------------------------------------------------------------------------------------

# --- Table S9 - Bayenv results by SNP

OUT_DIR=reports/supplementary_tables/Table_S9_Bayenv_BF_SNP
mkdir -p $OUT_DIR
cp results/bayenv_BF_pass.bed $OUT_DIR/Table_S9_Bayenv_BF_SNP.bed

# --- Table S10 - Bayenv results by gene

OUT_DIR=reports/supplementary_tables/Table_S10_Bayenv_BF_gene
mkdir -p $OUT_DIR
cp results/bayenv_BF_pass.pvals.txt $OUT_DIR/Table_S10_Bayenv_BF_gene.txt

# --- Table S11 - Bayenv results by GO

OUT_DIR=reports/supplementary_tables/Table_S11_Bayenv_BF_GO
mkdir -p $OUT_DIR
for ont in BP MF CC; do
    cp results/bayenv_BF_pass.pvals.$ont.results.txt \
        $OUT_DIR/Table_S11_Bayenv_BF_GO-$ont.txt
done
