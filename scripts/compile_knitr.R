#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Combine knitr documents
# ----------------------------------------------------------------------------------------

require(knitr)
require(markdown)

#knit('reports/SOM_tables.Rmd', 'reports/SOM_tables.md')
#markdownToHTML('reports/SOM_tables.md', 'reports/SOM_tables.html')

knit2html('reports/SOM_tables.Rmd', 'reports/SOM_tables.html',
    stylesheet='reports/tables.css')
