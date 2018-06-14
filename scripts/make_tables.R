#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ========================================================================================
# --- Make tables for paper
# ========================================================================================

source("http://bioconductor.org/biocLite.R")
library("GO.db")
library(xtable)

# ----------------------------------------------------------------------------------------
# --- Make table with pop-specific over-represented GO term results
# ----------------------------------------------------------------------------------------

pops = c("Batwa.GBR", "Anda.LWK", "Bakiga.GBR", "BR1.LWK")
onts = c("BP", "MF", "CC")

overrep.results = do.call(rbind, lapply(pops, function(pop) {
    do.call(rbind, lapply(onts, function(ont) {
        df = read.table(paste0("results/pbs.", pop, ".pvals.", ont, ".results.txt"),
            header=TRUE, sep="\t", quote="~")
        df$ont = ont
        df$pop = pop
        df = df[order(df$classicFisher),]
        df = df[1:3,]
        df
    }))
}))

# # Reduce to just significant
# overrep.results = overrep.results[overrep.results$classicFisher < 0.005,]

# Get full GO term description
overrep.results$desc.full = sapply(overrep.results$GO.ID, function (x) {
    go.term = GOTERM[[x]]
    if (is.null(go.term)) {
        return(NA)
    } else {
        return(Term(go.term))
    }
 })

# Make table

cleaned = overrep.results[,c(8,7,9,3:6)]

cleaned$pop = gsub("\\..*", "", cleaned$pop)
cleaned$pop = gsub("Anda", "Andamanese", cleaned$pop)
cleaned$pop = gsub("BR1",  "Brahmin",    cleaned$pop)

names(cleaned) = c("Population", "", "GO", "Annotated", "Obs.", "Exp.", "$p$")

xt = xtable(cleaned, digits=c(0,0,0,0,0,0,2,-3), align='llll|rrrr')

out.file = "reports/per_pop_GO_results.tex"

sink(out.file)

cat("\\documentclass[a4paper,landscape]{article}",
    "\\usepackage{graphicx}",
    "\\usepackage{longtable}",
    "\\DeclareGraphicsExtensions{.pdf}",
    "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}",
    "\\usepackage{caption}",
    "\\captionsetup[table]{labelformat=empty}",
    "\\begin{document}", sep="\n")

print.xtable(xt,
        include.rownames = FALSE,
        size="scriptsize",
        tabular.environment = 'longtable', floating = FALSE,
        caption.placement = "top")

cat("\\end{document}", sep="\n")

sink()

# ----------------------------------------------------------------------------------------
# --- Make table with Batwa-Anday convergent GO term results
# ----------------------------------------------------------------------------------------

conv = read.table(file="results/joint_GO_pvalues.txt")

# Just Batwa and Anda
conv = conv[grepl("batwa", conv$stat1) & grepl("anda", conv$stat2),]

# Get full GO term description
conv$Term.x = sapply(conv$GO.ID, function (x) {
    go.term = GOTERM[[x]]
    if (is.null(go.term)) {
        return(NA)
    } else {
        return(Term(go.term))
    }
 })

conv$ont = gsub("go\\.([^\\.]+)\\..*", "\\1", conv$stat1)

conv = conv[,c(13,2:12)]

# Remove CC
conv = conv[conv$ont != "CC",]

conv$stat1 = gsub("go\\.([^\\.]+)\\.", "", conv$stat1)
conv$stat2 = gsub("go\\.([^\\.]+)\\.", "", conv$stat2)

# --- Write to LaTeX table

conv$stat1 = gsub("pbs.", "$PBS$ ", conv$stat1)
conv$stat2 = gsub("pbs.", "$PBS$ ", conv$stat2)

conv$stat1 = gsub("batwa", "Batwa", conv$stat1)
conv$stat2 = gsub("batwa", "Batwa", conv$stat2)

conv$stat1 = gsub("bakiga", "Bakiga", conv$stat1)
conv$stat2 = gsub("bakiga", "Bakiga", conv$stat2)

conv$stat1 = gsub("anda", "Andamanese", conv$stat1)
conv$stat2 = gsub("anda", "Andamanese", conv$stat2)

conv$stat1 = gsub("bir", "Birhor", conv$stat1)
conv$stat2 = gsub("bir", "Birhor", conv$stat2)

conv$stat1 = gsub("br1", "Brahmin", conv$stat1)
conv$stat2 = gsub("br1", "Brahmin", conv$stat2)

conv$stat1 = gsub("go.bp.", "GO ", conv$stat1)
conv$stat2 = gsub("go.bp.", "GO ", conv$stat2)

names(conv) = c("", "GO",
    "Stat 1", "Sig.", "Exp.", "$p$",
    "Stat 2", "Sig.", "Exp.", "$p$",
    "Joint $p$ Edington", "Joint $p$ Fisher")

# Reorder to PBS and then N vs. S
conv = rbind(conv[grep("PBS", conv[["Stat 1"]]),],
             conv[grep("PBS", conv[["Stat 1"]], invert=TRUE),])

horiz.lines = unique(c(
    which(c(conv[["Stat 1"]], NA) != c(NA, conv[["Stat 1"]])),
    which(c(conv[["Stat 2"]], NA) != c(NA, conv[["Stat 2"]]))
)) - 1

# Remove stats
conv = conv[,-c(3,7)]

# Remove all numbers but p-value
conv = conv[,-c(3,4,6,7)]
names(conv)[3:4] = c("$p$ Batwa", "$p$ Andamanese")

xt = xtable(conv, digits=c(0,0,0,-3,-3,-3,-3), align='lll|ll|rr')

out.file = "reports/joint_GO_pvalues_subset.tex"
sink(out.file)

cat("\\documentclass[a4paper,landscape]{article}",
    "\\usepackage{graphicx}",
    "\\usepackage{longtable}",
    "\\DeclareGraphicsExtensions{.pdf}",
    "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}",
    "\\usepackage{caption}",
    "\\captionsetup[table]{labelformat=empty}",
    "\\begin{document}", sep="\n")

print(xt, type="latex",
    include.rownames=FALSE,
    tabular.environment="longtable", floating = FALSE,
    size="scriptsize",
    sanitize.colnames.function = identity,
    sanitize.text.function = identity,
    hline.after = horiz.lines)

cat("\\end{document}", sep="\n")

sink()
