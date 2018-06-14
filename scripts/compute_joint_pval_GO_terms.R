#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Find convergent GO terms
# ----------------------------------------------------------------------------------------

library(metap)
library(xtable)

onts = c("BP", "MF", "CC")

test.types = c("overrep", "enrich")

in.files = list()

# --- Files with GO enrichment data for PBS p-values

for (tt in test.types) {

    in.files[[tt]] = list()

    for (ont in onts) {
        in.files[[tt]][[paste0('go.', ont, '.', tt, '.pbs.batwa')]]  =
            paste0("results/pbs.Batwa.GBR.pvals.",  ont, '.', tt, ".results.txt")
        in.files[[tt]][[paste0('go.', ont, '.', tt, '.pbs.bakiga')]] =
            paste0("results/pbs.Bakiga.GBR.pvals.", ont, '.', tt, ".results.txt")
        in.files[[tt]][[paste0('go.', ont, '.', tt, '.pbs.anda')]]   =
            paste0("results/pbs.Anda.LWK.pvals.",   ont, '.', tt, ".results.txt")
        in.files[[tt]][[paste0('go.', ont, '.', tt, '.pbs.br1')]]    =
            paste0("results/pbs.BR1.LWK.pvals.",    ont, '.', tt, ".results.txt")
    }
}

# --- Files with GO enrichment data for PBS p-values - Size or MAF corrected results

for (cor.type in c("sizeCor", "mafCor")) {

    for (tt in test.types) {

        for (ont in onts) {
            in.files[[tt]][[paste0('go.', ont, '.', tt, '.pbs.batwa.', cor.type)]]  =
                paste0("results_", cor.type, "/pbs.Batwa.GBR.pvals.",  ont, '.', tt, ".results.txt")
            in.files[[tt]][[paste0('go.', ont, '.', tt, '.pbs.bakiga.', cor.type)]] =
                paste0("results_", cor.type, "/pbs.Bakiga.GBR.pvals.", ont, '.', tt, ".results.txt")
            in.files[[tt]][[paste0('go.', ont, '.', tt, '.pbs.anda.', cor.type)]]   =
                paste0("results_", cor.type, "/pbs.Anda.LWK.pvals.",   ont, '.', tt, ".results.txt")
            in.files[[tt]][[paste0('go.', ont, '.', tt, '.pbs.br1.', cor.type)]]    =
                paste0("results_", cor.type, "/pbs.BR1.LWK.pvals.",    ont, '.', tt, ".results.txt")
        }
    }
}

# ----------------------------------------------------------------------------------------

for (tt in test.types) {

    # --- Read in all files

    stats = sapply(names(in.files[[tt]]), function (file.desc) {
        this.stat = list(read.table(in.files[[tt]][[file.desc]],
                                    header=TRUE, fill=TRUE, sep="\t", quote="~"))
        return(this.stat)
    })

    if (tt == "overrep") {
        stats.p  = lapply(stats, function (x) { x$p = as.numeric(x$classicFisher); x })
        stats.ln = lapply(stats, function (x) { x$ln = log(as.numeric(x$classicFisher)); x })
    }

    if (tt == "enrich") {
        stats.p  = lapply(stats, function (x) { x$p = as.numeric(x$classicKS); x })
        stats.ln = lapply(stats, function (x) { x$ln = log(as.numeric(x$classicKS)); x })
    }

    # --- Function to compute joint p-value using sum of p-values method (Edgington's method)
    # --- as well as Fisher's method

    compute.joint.pval = function (stat1, stat2) {
        stats.p.both = merge(stats.p[[stat1]], stats.p[[stat2]],
            by=names(stats.ln[[stat1]][1]))
        stats.p.both$stat1 = stat1
        stats.p.both$stat2 = stat2

        edington.p = apply(cbind(stats.p.both$p.x, stats.p.both$p.y), 1, sump)
        fisher.p   = apply(cbind(stats.p.both$p.x, stats.p.both$p.y), 1, sumlog)

        stats.p.both$joint.p.edington = unlist(lapply(edington.p, function(x) { x$p }))
        stats.p.both$joint.p.fisher   = unlist(lapply(fisher.p,   function(x) { x$p }))
        combined = stats.p.both[order(stats.p.both$joint.p.edington, decreasing=FALSE),]
    }

    pops = c("batwa", "bakiga", "anda", "br1")
    compares = t(combn(pops, 2))

    if (tt == "enrich") {
        compares = cbind(compares, 'pbs', 'pbs')
    }

    # Triple pop pairs and add ontology column
    compares = rbind(cbind(compares, 'BP'), cbind(compares, 'MF'), cbind(compares, 'CC'))

    # Triple it all again and
    # add correction type column (empty for first third)
    compares = rbind(cbind(compares, ''),
                     cbind(compares, 'sizeCor'),
                     cbind(compares, 'mafCor'))

    # Compute joint p-values for all comparisons

    joint = lapply(1:nrow(compares), function (row.num) {
        x = compares[row.num,]
        pop1      = x[1]
        pop2      = x[2]
        analysis1 = x[3]
        analysis2 = x[4]
        ontology  = x[5]
        cor.type  = x[6]

        analysis.name = paste(tt, ontology, analysis1, pop1, analysis2, pop2, sep=".")

        if (cor.type != "") {
            analysis.name = paste(analysis.name, cor.type, sep=".")
        }

        stat1 = paste("go", ontology, tt, analysis1, pop1, sep=".")
        stat2 = paste("go", ontology, tt, analysis2, pop2, sep=".")

        if (cor.type != "") {
            stat1 = paste(stat1, cor.type, sep=".")
            stat2 = paste(stat2, cor.type, sep=".")
        }

        result = list(compute.joint.pval(stat1, stat2))
        names(result) = analysis.name
        result
    })

    # Go from (unnamed) list of (named) lists to (named) list.
    joint = unlist(joint, recursive=FALSE)

    # --- Write results to file

    if (tt == "overrep") {
        cols.to.keep = c("stat1", "Significant.x", "Expected.x", "classicFisher.x",
                         "stat2", "Significant.y", "Expected.y", "classicFisher.y",
                         "joint.p.edington", "joint.p.fisher")
    }

    if (tt == "enrich") {
        cols.to.keep = c("stat1", "classicKS.x",
                         "stat2", "classicKS.y",
                         "joint.p.edington", "joint.p.fisher")
    }

    combine.joint.info = function (joint.stats) {

        GO.term = joint[[joint.stats]][,1:2]

        cleaned = cbind(GO.term, joint[[joint.stats]][,cols.to.keep])

        cleaned = cleaned[which(cleaned$joint.p.edington < 0.01),]

        cleaned = cleaned[order(cleaned$joint.p.edington),]

    }

    cleaned = do.call(rbind, lapply(names(joint), combine.joint.info))

    write.table(cleaned, file=paste0("results/joint_GO_pvalues.", tt, ".txt"))

    # --- Write to LaTeX table

    horiz.lines = unique(c(
        which(c(cleaned$stat1, NA) != c(NA, cleaned$stat1)),
        which(c(cleaned$stat2, NA) != c(NA, cleaned$stat2))
    )) - 1

    cleaned$stat1 = gsub("overrep.", "", cleaned$stat1)
    cleaned$stat2 = gsub("overrep.", "", cleaned$stat2)

    cleaned$stat1 = gsub("enrich.", "", cleaned$stat1)
    cleaned$stat2 = gsub("enrich.", "", cleaned$stat2)

    cleaned$stat1 = gsub("go.", "", cleaned$stat1)
    cleaned$stat2 = gsub("go.", "", cleaned$stat2)

    cleaned$stat1 = gsub("pbs.", "$PBS$ ", cleaned$stat1)
    cleaned$stat2 = gsub("pbs.", "$PBS$ ", cleaned$stat2)

    cleaned$stat1 = gsub("batwa", "Batwa", cleaned$stat1)
    cleaned$stat2 = gsub("batwa", "Batwa", cleaned$stat2)

    cleaned$stat1 = gsub("bakiga", "Bakiga", cleaned$stat1)
    cleaned$stat2 = gsub("bakiga", "Bakiga", cleaned$stat2)

    cleaned$stat1 = gsub("anda", "Anda.", cleaned$stat1)
    cleaned$stat2 = gsub("anda", "Anda.", cleaned$stat2)

    cleaned$stat1 = gsub("bir", "Birhor", cleaned$stat1)
    cleaned$stat2 = gsub("bir", "Birhor", cleaned$stat2)

    cleaned$stat1 = gsub("br1", "Brahmin", cleaned$stat1)
    cleaned$stat2 = gsub("br1", "Brahmin", cleaned$stat2)

    cleaned$stat1 = gsub("bp.", "BP ", cleaned$stat1, ignore.case=TRUE)
    cleaned$stat2 = gsub("bp.", "BP ", cleaned$stat2, ignore.case=TRUE)
    cleaned$stat1 = gsub("mf.", "MF ", cleaned$stat1, ignore.case=TRUE)
    cleaned$stat2 = gsub("mf.", "MF ", cleaned$stat2, ignore.case=TRUE)
    cleaned$stat1 = gsub("cc.", "CC ", cleaned$stat1, ignore.case=TRUE)
    cleaned$stat2 = gsub("cc.", "CC ", cleaned$stat2, ignore.case=TRUE)

    if (tt == "overrep") {
        names(cleaned) = c("GO", "Description",
            "Stat 1", "Sig.", "Exp.", "$p$",
            "Stat 2", "Sig.", "Exp.", "$p$",
            "$p$ Edington", "$p$ Fisher")
    }
    if (tt == "enrich") {
        names(cleaned) = c("GO", "Description",
            "Stat 1", "$p$",
            "Stat 2", "$p$",
            "$p$ Edington", "$p$ Fisher")
    }

    if (tt == "overrep") {
        xt = xtable(cleaned, digits=c(0,0,0,0,0,2,-3,0,0,2,-3,-3,-3),
            align='lll|lrrl|lrrl|rr')
    }

    if (tt == "enrich") {
        xt = xtable(cleaned, digits=c(0,0,0,0,-3,0,-3,-3,-3),
            align='lll|ll|ll|rr')
    }

    out.file = paste0("results/joint_GO_pvalues.", tt, ".tex")
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
