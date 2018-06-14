#!/bin/bash

# ========================================================================================
# --- Pull out interesting results from a priori gene list analysis
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Outlier-based test results
# ----------------------------------------------------------------------------------------

for POP in Batwa.GBR Anda.LWK Bakiga.GBR BR1.LWK; do

    echo "--------------------------------------------------"
    echo "--- $POP over-representation results:"
    echo "--------------------------------------------------"

    echo "--- All gene-based overrep results (0.01 cutoff):"
    awk 'BEGIN {FS="\t"; OFS="\t"} { if ($2 == 0.01) print $1,$2,$7}' \
        reports/pbs.$POP.pvals.fisher.overrep.txt

    echo "--- Significant gene-based overrep results (0.01 cutoff):"
    awk 'BEGIN {FS="\t"; OFS="\t"} { if ($2 == 0.01 && $7 < 0.1) print $1,$2,$7}' \
        reports/pbs.$POP.pvals.fisher.overrep.txt

    echo "--- Text on significant gene-based overrep results (0.01 cutoff):"
    awk 'BEGIN {FS="\t"} {
        if ($2 == 0.01 && $7 < 0.1) {
            print "\t",$1,"[",$2,"] : ",
                $3,"out of",$5+$3,"(",100 * $3 / ($5+$3),"%) high PBS genes in list, and",
                $4-$3,"out of",$4+$6-$3-$5,"(", 100 * ($4-$3) / (($4+$6)-($3+$5)),
                "%) other genes in list. Fisher p =", $7
        }
    }' reports/pbs.$POP.pvals.fisher.overrep.txt

    echo "--- Significant gene-based overrep results (all cutoffs):"
    awk 'BEGIN {FS="\t"; OFS="\t"} { if ($7 < 0.1) print $1,$2,$7}' \
        reports/pbs.$POP.pvals.fisher.overrep.txt

done

# ----------------------------------------------------------------------------------------
# --- Enrichment-based test results
# ----------------------------------------------------------------------------------------

echo "==========================================================================="

for POP in Batwa.GBR Anda.LWK Bakiga.GBR BR1.LWK; do

    echo "--------------------------------------------------"
    echo "--- $POP enrichment results:"
    echo "--------------------------------------------------"

    echo "--- All gene-based enrich results:"
    cat reports/pbs.$POP.pvals.wilcoxshift.txt

    echo "--- Significant gene-based enrich results:"
    awk 'BEGIN {FS="\t"} { if ($2 < 0.1) print $0 }' \
        reports/pbs.$POP.pvals.wilcoxshift.txt
done

# ----------------------------------------------------------------------------------------
# --- Bayenv test results - overrep and enrichment
# ----------------------------------------------------------------------------------------

echo "==========================================================================="

echo "--------------------------------------------------"
echo "--- Bayenv over-representation results:"
echo "--------------------------------------------------"

echo "--- All gene-based overrep results (0.01 cutoff):"
awk 'BEGIN {FS="\t"; OFS="\t"} { if ($2 == 0.01) print $1,$2,$7}' \
    reports/bayenv_BF_pass.pvals.fisher.overrep.txt

echo "--- Significant gene-based overrep results (0.01 cutoff):"
awk 'BEGIN {FS="\t"; OFS="\t"} { if ($2 == 0.01 && $7 < 0.1) print $1,$2,$7}' \
    reports/bayenv_BF_pass.pvals.fisher.overrep.txt

echo "--- Text on significant gene-based overrep results (0.01 cutoff):"
awk 'BEGIN {FS="\t"} {
    if ($2 == 0.01 && $7 < 0.1) {
        print "\t",$1,"[",$2,"] : ",
            $3,"out of",$5+$3,"(",100 * $3 / ($5+$3),"%) high PBS genes in list, and",
            $4-$3,"out of",$4+$6-$3-$5,"(", 100 * ($4-$3) / (($4+$6)-($3+$5)),
            "%) other genes in list. Fisher p =", $7
    }
}' reports/bayenv_BF_pass.pvals.fisher.overrep.txt

echo "--- Significant gene-based overrep results (all cutoffs):"
awk 'BEGIN {FS="\t"; OFS="\t"} { if ($7 < 0.1) print $1,$2,$7}' \
    reports/bayenv_BF_pass.pvals.fisher.overrep.txt

echo "--------------------------------------------------"
echo "--- Bayenv enrichment results:"
echo "--------------------------------------------------"

echo "--- All gene-based enrich results:"
cat reports/bayenv_BF_pass.pvals.wilcoxshift.txt
echo "--- Significant gene-based enrich results:"
awk 'BEGIN {FS="\t"} { if ($2 < 0.1) print $0 }' \
    reports/bayenv_BF_pass.pvals.wilcoxshift.txt

exit
