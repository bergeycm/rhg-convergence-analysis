#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Write tables with SNP-based results
# ----------------------------------------------------------------------------------------

mkdir -p reports/latex
cd reports/latex

# Copy in main LaTeX file
cp ../../scripts/latex/extreme_SNP_results_tables.tex .

IN_TEX_FILES=("eAGR_eRHG.weir.fst.anno.extreme.vals"
              "eAGR_eRHG.weir.fst.anno.gwas.region.extreme.vals"
              "eAGR_eRHG.weir.fst.anno.deleterious.extreme.vals"
              "eAGR_eRHG.weir.fst.anno.chisq.overrep"
              "pbs.Bakiga.GBR.anno.extreme.vals"
              "pbs.Bakiga.GBR.anno.gwas.region.extreme.vals"
              "pbs.Bakiga.GBR.anno.deleterious.extreme.vals"
              "pbs.Bakiga.GBR.anno.chisq.overrep"
              "pbs.Batwa.GBR.anno.extreme.vals"
              "pbs.Batwa.GBR.anno.gwas.region.extreme.vals"
              "pbs.Batwa.GBR.anno.deleterious.extreme.vals"
              "pbs.Batwa.GBR.anno.chisq.overrep"
              "pbs.Anda.LWK.anno.extreme.vals"
              "pbs.Anda.LWK.anno.gwas.region.extreme.vals"
              "pbs.Anda.LWK.anno.deleterious.extreme.vals"
              "pbs.Anda.LWK.anno.chisq.overrep"
              "pbs.BIR.LWK.anno.extreme.vals"
              "pbs.BIR.LWK.anno.gwas.region.extreme.vals"
              "pbs.BIR.LWK.anno.deleterious.extreme.vals"
              "pbs.BIR.LWK.anno.chisq.overrep"
              "pbs.BR1.LWK.anno.extreme.vals"
              "pbs.BR1.LWK.anno.gwas.region.extreme.vals"
              "pbs.BR1.LWK.anno.deleterious.extreme.vals"
              "pbs.BR1.LWK.anno.chisq.overrep"
              "ihs.twa.txt.anno.extreme.vals"
              "ihs.twa.txt.anno.gwas.region.extreme.vals"
              "ihs.twa.txt.anno.deleterious.extreme.vals"
              "ihs.twa.txt.anno.chisq.overrep"
              "ihs.kiga.txt.anno.extreme.vals"
              "ihs.kiga.txt.anno.gwas.region.extreme.vals"
              "ihs.kiga.txt.anno.deleterious.extreme.vals"
              "ihs.kiga.txt.anno.chisq.overrep")

for TEX in ${IN_TEX_FILES[*]}; do
    cp ../$TEX.tex ${TEX//./_}.tex
done

rm -f extreme_SNP_results_tables.aux

for x in 1 2; do
    pdflatex extreme_SNP_results_tables.tex
done

for TEX in ${IN_TEX_FILES[*]}; do
    rm ${TEX//./_}.tex
done

# Move file to main reports folder
mv extreme_SNP_results_tables.pdf ..

cd ../..

rm -r reports/latex
