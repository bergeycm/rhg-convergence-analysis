#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Compile overrep and enrichment test LaTeX
# ----------------------------------------------------------------------------------------

# Run pdflatex twice, in case table widths change

# Stuff in reports/gene_lists
for X in 1 2; do
    cd reports/gene_lists/
    for tex in $( ls *.tex ); do
        pdflatex $tex
    done
    cd ../..
done

# Stuff in reports
for X in 1 2; do
    cd reports/
    for tex in $( ls *.tex ); do
        # See if this is a full LaTeX document, rather than a snippet
        if grep "document" $tex; then
            pdflatex $tex
        fi
    done
    cd ..
done

# Stuff in results
for X in 1 2; do
    cd results/
    for tex in $( ls *.tex ); do
        # See if this is a full LaTeX document, rather than a snippet
        if grep "document" $tex; then
            pdflatex $tex
        fi
    done
    cd ..
done
