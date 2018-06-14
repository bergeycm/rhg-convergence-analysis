#!/bin/bash

module load r/3.4

# ========================================================================================
# --- Run overrep and enrichment tests on p-values (single and joint)
# ========================================================================================

pval_file=$1

if [[ `echo $pval_file | grep -c "jointp"` -eq 0 ]]; then

    # ------------------------------------------------------------------------------------
    # --- Single population p-values
    # ------------------------------------------------------------------------------------

    if [[ $pval_file =~ .*weir\.fst.* ]]; then
        label="Fst p-value"
    elif [[ $pval_file =~ .*ihs.* ]]; then
        label="iHS p-value"
    elif [[ $pval_file =~ .*pbs.* ]]; then
        label="PBS p-value"
    fi

    log_file=${pval_file/txt/overrep.log}

    Rscript scripts/do_enrich_overrep_tests.R \
        $pval_file FALSE "$label" 0.99 TRUE \
        data/GO_of_interest.txt &> reports/`basename $log_file`

else

    # ------------------------------------------------------------------------------------
    # --- Joint p-values
    # ------------------------------------------------------------------------------------

    label="Joint PBS p-value"

    log_file=${pval_file/txt/overrep.log}

    Rscript scripts/do_enrich_overrep_tests.R \
        $pval_file FALSE "$label" 0.99 TRUE \
        data/GO_of_interest.txt &> reports/`basename $log_file`

fi

exit
