#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Find convergent genes (e.g. high PBS in both Batwa and Andamanese)
# ----------------------------------------------------------------------------------------

library(xtable)
library(metap)

# --- Read in gene p-values

stats = list()

stats[['pbs.batwa']] =  read.table("reports/pbs.Batwa.GBR.pvals.extreme.vals.all.txt",
                            header=TRUE, fill=TRUE)
stats[['pbs.bakiga']] = read.table("reports/pbs.Bakiga.GBR.pvals.extreme.vals.all.txt",
                            header=TRUE, fill=TRUE)
stats[['pbs.anda']] =   read.table("reports/pbs.Anda.LWK.pvals.extreme.vals.all.txt",
                            header=TRUE, fill=TRUE)
stats[['pbs.br1']] =    read.table("reports/pbs.BR1.LWK.pvals.extreme.vals.all.txt",
                            header=TRUE, fill=TRUE)

# --- Function to compute joint p-value using sum of p-values method (Edgington's method)
# --- as well as Fisher's method

compute.joint.pval = function (stat1, stat2) {
    stats.both = merge(stats[[stat1]], stats[[stat2]],
        by="gene")
    stats.both$stat1 = stat1
    stats.both$stat2 = stat2

    # Remove genes with missing p-values
    stats.both = stats.both[!(is.na(stats.both$p.val.x) | is.na(stats.both$p.val.y)),]

    edington.p = apply(cbind(stats.both$p.val.x, stats.both$p.val.y), 1, sump)

    # Set zero values to something slightly more than zero since log(0) = -Inf
    if (sum(stats.both$p.val.x == 0) > 0) {
        stats.both[stats.both$p.val.x == 0,]$p.val.x = 1e-30
    }
    if (sum(stats.both$p.val.y == 0) > 0) {
        stats.both[stats.both$p.val.y == 0,]$p.val.y = 1e-30
    }

    fisher.p = apply(cbind(stats.both$p.val.x, stats.both$p.val.y), 1, sumlog)

    stats.both$joint.p.edington = unlist(lapply(edington.p, function(x) { x$p }))
    stats.both$joint.p.fisher   = unlist(lapply(fisher.p,   function(x) { x$p }))

    combined = stats.both[order(stats.both$joint.p.edington, decreasing=FALSE),]
}

joint = list()

# --- Batwa PBS and Andamanese PBS

joint[['pbs.batwa.pbs.anda']]  = compute.joint.pval ('pbs.batwa', 'pbs.anda')

# --- Bakiga PBS and Andamanese PBS

joint[['pbs.bakiga.pbs.anda']] = compute.joint.pval ('pbs.bakiga', 'pbs.anda')

# --- Batwa PBS and BR1 PBS

joint[['pbs.batwa.pbs.br1']]   = compute.joint.pval ('pbs.batwa', 'pbs.br1')

# --- Bakiga PBS and BR1 PBS

joint[['pbs.bakiga.pbs.br1']]  = compute.joint.pval ('pbs.bakiga', 'pbs.br1')

# --- Write results to file

cols.to.keep = c(8:9,1:2,4:5,7,10:11)

combine.joint.info = function (joint.stats) {

    full = joint[[joint.stats]][,cols.to.keep]

    full = full[order(full$joint.p.edington),]

}

full = do.call(rbind, lapply(names(joint), combine.joint.info))

write.table(full, file="results/joint_gene_pvalues.txt")

cleaned = full[which(full$joint.p.edington < 0.0001),]

# --- Write to LaTeX table

write.to.latex = function (df, out.file) {

    horiz.lines = unique(c(
        which(c(df$stat1, NA) != c(NA, df$stat1)),
        which(c(df$stat2, NA) != c(NA, df$stat2))
    )) - 1

    df$stat1 = gsub("pbs.", "$PBS$ ", df$stat1)
    df$stat2 = gsub("pbs.", "$PBS$ ", df$stat2)

    df$stat1 = gsub("batwa", "Batwa", df$stat1)
    df$stat2 = gsub("batwa", "Batwa", df$stat2)

    df$stat1 = gsub("bakiga", "Bakiga", df$stat1)
    df$stat2 = gsub("bakiga", "Bakiga", df$stat2)

    df$stat1 = gsub("anda", "Andamanese", df$stat1)
    df$stat2 = gsub("anda", "Andamanese", df$stat2)

    df$stat1 = gsub("br1", "Brahmin", df$stat1)
    df$stat2 = gsub("br1", "Brahmin", df$stat2)

    names(df) = c(
        "Stat 1", "Stat 2", "Gene",
        "Stat 1 $p$", " ",
        "Stat 2 $p$", " ",
        "$p$ Edgington", "$p$ Fisher")

    xt = xtable(df, digits=c(0,0,0,0,-3,0,-3,0,-3,-3), align='llll|rr|rr|rr')

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
        size="\\fontsize{9pt}{10pt}\\selectfont",
        sanitize.colnames.function = identity,
        sanitize.text.function = identity,
        hline.after = horiz.lines)

    cat("\\end{document}", sep="\n")

    sink()

}

write.to.latex(df = cleaned, out.file = "reports/joint_gene_pvalues.tex")
