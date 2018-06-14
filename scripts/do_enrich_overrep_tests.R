#!/usr/bin/env Rscript

# module load R/3.3.0

options(stringsAsFactors=FALSE)

# ========================================================================================
# --- Do overrep and enrichment tests with topGO
# ========================================================================================

source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("topGO")
# biocLite("org.Hs.eg.db")

library(topGO)
library(xtable)
library(GOSemSim)

# Input is gene or SNP list with values

args = commandArgs(trailingOnly=TRUE)
in.file       = args[1]
has.header    = as.logical(args[2])  # Does data have header line?
description   = args[3]
cutoff        = as.numeric(args[4])  # Cutoff quantile for overrep test
is.pval       = as.logical(args[5])  # Should scores be flipped and treated as p-value?
ecdf.plots    = args[6]              # File of GO terms to plot,
                                     #  e.g. data/GO_of_interest.txt

if (grepl("sizeCor", in.file)) {
    dir.create("reports_sizeCor/gene_lists", showWarnings=FALSE)
} else if (grepl("mafCor", in.file)) {
    dir.create("reports_mafCor/gene_lists", showWarnings=FALSE)
} else {
    dir.create("reports/gene_lists", showWarnings=FALSE)
}

# If TRUE (default), gene mode. Otherwise, SNP mode.
if (length(args) == 6) {
    gene.mode = TRUE
} else {
    gene.mode = args[7]
}

cat("--- Input parameters: ---\n")
cat(paste("in.file:      ", in.file      ), "\n")
cat(paste("has.header:   ", has.header   ), "\n")
cat(paste("description:  ", description  ), "\n")
cat(paste("cutoff:       ", cutoff       ), "\n")
cat(paste("is.pval:      ", is.pval      ), "\n")
cat(paste("ecdf.plots:   ", ecdf.plots   ), "\n")
cat(paste("gene.mode:    ", gene.mode    ), "\n")

input = read.table(in.file, header=has.header)

# Collapse SNP info into one column
if (gene.mode == FALSE) {

    input = input[,c(5,4)]

    input = do.call(rbind, lapply(1:nrow(input), function (x) {
        this.row = input[x,]
        genes = strsplit(this.row$V5, split=",")
        genes.u = unique(genes[[1]])
        res = cbind(genes.u, this.row[rep(1, length(genes.u)),2])
        return(res)
    }))
    input = data.frame(input)
    input$V2 = as.numeric(input$V2)
}

# Lose column names
names(input) = paste0("V", 1:ncol(input))

# Lose NAs
input = input[!is.na(input$V2),]

# Lose Inf values by setting them to non-infinity maximum
if (sum(input$V2 == Inf, na.rm=TRUE) > 0) {
    input[input$V2 == Inf,]$V2 = max(input$V2[input$V2 != Inf])
}

if (sum(input$V2 == -Inf, na.rm=TRUE) > 0) {
    input[input$V2 == -Inf,]$V2 = min(input$V2[input$V2 != -Inf])
}

# ----------------------------------------------------------------------------------------

#refgene = read.table("refGene/refGene.sort.simple.genes.gtf")

#all.genes = unique(sort(refgene$V9))

if (ecdf.plots != "") {
    GO.of.interest = read.table(ecdf.plots, sep="\t")

    if (grepl("sizeCor", in.file)) {
        dir.create("reports_sizeCor/enrichment_ecdf_plots_sizeCor", showWarnings=FALSE)
    } else if (grepl("mafCor", in.file)) {
        dir.create("reports_mafCor/enrichment_ecdf_plots_mafCor", showWarnings=FALSE)
    } else {
        dir.create("reports/enrichment_ecdf_plots", showWarnings=FALSE)
    }
}

# ----------------------------------------------------------------------------------------
# --- Function to adjust p-values
# ----------------------------------------------------------------------------------------

cluster.go.and.adjust.p = function (allRes.full, ont) {

    hsGO = godata('org.Hs.eg.db', ont=ont)

    go.dists = mgoSim(allRes.full$GO.ID, allRes.full$GO.ID,
        semData=hsGO, measure="Jiang", combine=NULL)

    # Ward Hierarchical Clustering
    go.dists.d = as.dist(1 - go.dists)
    go.dists.d[is.na(go.dists.d)] = 1

    fit = hclust(go.dists.d, method="complete")

    #   pdf("go_dist_tree.pdf", height=5, width=25)
    #   plot(fit) # display dendogram
    #   dev.off()

    groups = cutree(fit, h=0.5)

    # For each group, find the member with the lowest p-value
    member.indices = lapply(1:max(groups), function (grp) {
        this.grp = which(groups == grp)
        best.grp.member = this.grp[1]
        random.member = sample(this.grp, 1)
        grp.members.to.cut = this.grp[-c(1)]
        return(list(go.idx.keep = best.grp.member,
                    go.idx.lose = grp.members.to.cut,
                    random.member = random.member))
    })

    go.idx.keep = do.call(c, lapply(member.indices, function (x) {
        return(x[[1]])
    }))

    allRes.sm = allRes.full[go.idx.keep,]

    p.col = "classicFisher"
    if ("classicKS" %in% names(allRes.sm)) {
        p.col = "classicKS"
    }

    allRes.sm$p.adj = p.adjust(as.numeric(allRes.sm[[p.col]]),
        method="fdr")

    return(allRes.sm)
}

# ----------------------------------------------------------------------------------------
# --- Do overrep and enrich tests and output ECDF plots
# ----------------------------------------------------------------------------------------

for (ont in c('BP', 'MF', 'CC')) {

    subset = input$V1

    # Keep values
    orig.vals = input
    all.genes = as.vector(t(orig.vals$V2))
    names(all.genes) = orig.vals$V1

    # Lose NA values
    all.genes = all.genes[!is.na(all.genes)]

    top.interest.genes = function (val) {
        if (is.pval) {
            # If this is a p-value, we don't compute a quantile to determine significance
            return(val <= 1 - cutoff)
        } else {
            return(val >= quantile(orig.vals$V2, cutoff, na.rm=TRUE))
        }
    }

    GOdata = new("topGOdata",
        description = description,
        ontology = ont,
        allGenes = all.genes,
        geneSel = top.interest.genes,
        annot = annFUN.org,
        ID = "alias",
        mapping = "org.Hs.eg",
        nodeSize = 50)

    # ------------------------------------------------------------------------------------
    # --- Do overrep tests
    # ------------------------------------------------------------------------------------

    result.fisher.classic  = runTest(GOdata, "classic",  "fisher")

    allRes.overrep = GenTable(GOdata,
        classicFisher = result.fisher.classic,
        orderBy = "classicFisher", topNodes = 50)

    # Also get all p-values for writing to simple table
    allGO = usedGO(object = GOdata)

    allRes.overrep.full = GenTable(GOdata,
        classicFisher  = result.fisher.classic,
        orderBy = "classicFisher",
        topNodes = length(allGO))

    # Adjust p-values with FDR using the full dataset
    allRes.overrep.full$p.adj = p.adjust(allRes.overrep.full$classicFisher, method="fdr")

    # Adjust p-values after pruning GO terms based on their semantic similarity
    allRes.overrep.sm = cluster.go.and.adjust.p(allRes.overrep.full, ont)

    # Add OR and its CI to table
    add.or.ci = function (resTable.row, fisher.res, GOdata.obj) {

        go.sig = resTable.row$Significant
        notgo.sig = sum(geneScore(GOdata.obj, whichGenes = genes(GOdata.obj)) <= 0.01) -
                resTable.row$Significant
        go.notsig = resTable.row$Annotated - resTable.row$Significant
        notgo.notsig = numGenes(GOdata.obj) -
                sum(geneScore(GOdata.obj, whichGenes = genes(GOdata.obj)) <= 0.01) -
                (resTable.row$Annotated - resTable.row$Significant)

        fish.m = matrix(c(
            go.sig,
            notgo.sig,
            go.notsig,
            notgo.notsig), nrow=2)

        fishy = fisher.test(fish.m)

        return(c(fishy$estimate, fishy$conf.int[[1]], fishy$conf.int[[2]]))
    }

    # Add OR and CI to small table...
    or.ci = data.frame(t(sapply(1:nrow(allRes.overrep.sm),
        function (x) {
            add.or.ci(allRes.overrep.sm[x,], result.fisher.classic, GOdata)
        })))
    names(or.ci) = c("odds.ratio", "CI.95.low", "CI.95.up")

    allRes.overrep.sm = cbind(allRes.overrep.sm, or.ci)

    # ...and to full table
    or.ci = data.frame(t(sapply(1:nrow(allRes.overrep.full),
        function (x) {
            add.or.ci(allRes.overrep.full[x,], result.fisher.classic, GOdata)
        })))
    names(or.ci) = c("odds.ratio", "CI.95.low", "CI.95.up")

    allRes.overrep.full = cbind(allRes.overrep.full, or.ci)

    # Permute gene-value associations and rerun GO overrep to generate
    # empirical distribution of GO overrep p-values.
    source("scripts/GO_iteration_functions.R")

    bs.p = bs.to.get.GO.pval.dist.overrep(num.iter = 1000)

    get.bs.p = function(go.term) {
        go.bs.pvals = bs.p[bs.p$X1 == go.term, 2:ncol(bs.p)]
        go.real.p = allRes.overrep.full[allRes.overrep.full$GO.ID == go.term,]$classicFisher
        less.extreme = sum(go.bs.pvals >= go.real.p)
        tot.trials = length(go.bs.pvals)
        return(1 - (less.extreme / tot.trials))
    }

    allRes.overrep.full$bs.p = sapply(allRes.overrep.full$GO.ID, function (go) {
        get.bs.p(go)
    })

    # Write table of empirical GO p-values
    bs.out.file = gsub("txt$", paste0(ont, ".bootstrapped.GO.p-values.overrep.txt"),
        gsub("results", "reports", in.file))
    write.table(bs.p, file=bs.out.file,
        sep="\t", quote=TRUE, row.names=FALSE, col.names=FALSE)

    # Do similar thing, but this time permute gene-GO relationship

    rm(list=c("bs.p", "get.bs.p"))

    bs.p = bs.to.get.GO.pval.dist.overrep.GOgene(num.iter = 1000)

    get.bs.p = function(go.term) {
        if (go.term %in% bs.p$X1) {
            go.bs.pvals = bs.p[bs.p$X1 == go.term, 2:ncol(bs.p)]
            go.real.p = allRes.overrep.full[allRes.overrep.full$GO.ID == go.term,]$classicFisher
            less.extreme = sum(go.bs.pvals >= go.real.p)
            tot.trials = length(go.bs.pvals)
            res = 1 - (less.extreme / tot.trials)
        } else {
            res = NA
        }

        return (res)
    }

    allRes.overrep.full$bs.p.geneGO = sapply(allRes.overrep.full$GO.ID, function (go) {
        get.bs.p(go)
    })

    # Write table of empirical GO p-values
    bs.out.file = gsub("txt$", paste0(ont, ".bootstrapped-geneGO.GO.p-values.overrep.txt"),
        gsub("results", "reports", in.file))
    write.table(bs.p, file=bs.out.file,
        sep="\t", quote=TRUE, row.names=FALSE, col.names=FALSE)

    # ------------------------------------------------------------------------------------

    # Write LaTeX table

    xt = xtable(head(allRes.overrep.sm, n=50))

    if (gene.mode) {
        out.file = gsub("txt$", paste0(ont, ".overrep.tex"),
            gsub("results", "reports", in.file))
    } else {
        out.file = gsub("bed$", paste0(ont, ".overrep.bySNP.tex"),
            gsub("results", "reports", in.file))
    }

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

    # ------------------------------------------------------------------------------------
    # --- Do enrich tests
    # ------------------------------------------------------------------------------------

    # From topGO source:
    #   scoreOrder =  TRUE (decreasing)  the max(score) is considered the best score
    #   scoreOrder = FALSE (increasing)  the min(score) is considered the best score

    if (is.pval) {
        scoreOrder = "increasing"
    } else {
        scoreOrder = "decreasing"
    }

    result.ks.classic  = runTest(GOdata, "classic",  "ks", scoreOrder = scoreOrder)
    #result.ks.elim     = runTest(GOdata, "elim",     "ks", scoreOrder = scoreOrder)
    #result.ks.weight01 = runTest(GOdata, "weight01", "ks", scoreOrder = scoreOrder)

    #allRes.enrich = GenTable(GOdata,
    #    classicKS  = result.ks.classic,
    #    elimKS     = result.ks.elim,
    #    weight01KS = result.ks.weight01,
    #    orderBy = "classicKS", topNodes = 50)

    allRes.enrich = GenTable(GOdata,
        classicKS  = result.ks.classic,
        orderBy = "classicKS", topNodes = 50)

    # Also get all p-values for writing to simple table
    allRes.enrich.full = GenTable(GOdata,
        classicKS = result.ks.classic,
        orderBy = "classicKS",
        topNodes = length(allGO))

    # Convert character "< 1e-30" to numeric 1e-30
    allRes.enrich.full$classicKS[allRes.enrich.full$classicKS == "< 1e-30"] = 1e-30
    allRes.enrich.full$classicKS = as.numeric(allRes.enrich.full$classicKS)

    # Adjust p-values with FDR using the full dataset
    allRes.enrich.full$p.adj = p.adjust(allRes.enrich.full$classicKS, method="fdr")

    # Adjust p-values after pruning GO terms based on their semantic similarity
    allRes.enrich.sm = cluster.go.and.adjust.p(allRes.enrich.full, ont)

    # Add SE to table
    add.se = function (resTable.row, ks.res, GOdata.obj) {

        vals.in.go = geneScore(GOdata.obj,
            whichGenes=genesInTerm(GOdata, resTable.row$GO.ID)[[1]])
        vals.all = geneScore(GOdata.obj, use.names=TRUE)
        vals.notin.go = vals.all[! names(vals.all) %in% names(vals.in.go)]

        if (length(vals.notin.go) == 0) {
            return(NA)
        }

        this.ks.res = ks.test(vals.in.go, vals.notin.go, alternative="greater")

        jack.p.vals = sapply(1:100, function (x) {
            vals.in.go.samp = sample(vals.in.go, size=length(vals.in.go) - 1,
                                     replace=TRUE)
            ks.test(vals.in.go.samp, vals.notin.go, alternative="greater")$p
        })

        this.SE = sd(jack.p.vals) / sqrt(length(jack.p.vals))

        return(this.SE)
    }

    row.count = 1

    # First for full table...
    all.se = data.frame(t(sapply(1:nrow(allRes.enrich.full),
        function (x) {
            cat(paste0("Processing row ", x, "...\n"))
            row.count = row.count + 1
            add.se(allRes.enrich.full[x,], result.ks.classic, GOdata)
        })))

    allRes.enrich.full = cbind(allRes.enrich.full, t(all.se))

    names(allRes.enrich.full)[ncol(allRes.enrich.full)] = "std.err"

    # ...and for small table
    all.se = data.frame(t(sapply(1:nrow(allRes.enrich.sm),
        function (x) {
            cat(paste0("Processing row ", x, "...\n"))
            row.count = row.count + 1
            add.se(allRes.enrich.sm[x,], result.ks.classic, GOdata)
        })))

    allRes.enrich.sm = cbind(allRes.enrich.sm, t(all.se))

    names(allRes.enrich.sm)[ncol(allRes.enrich.sm)] = "std.err"

    # ------------------------------------------------------------------------------------

    # Permute gene-value associations and rerun GO enrich to generate
    # empirical distirbution of GO enrich p-values.
    bs.p = bs.to.get.GO.pval.dist.enrich(num.iter = 1000)

    get.bs.p = function(go.term) {
        go.bs.pvals = bs.p[bs.p$X1 == go.term, 2:ncol(bs.p)]
        go.real.p = allRes.enrich.full[allRes.enrich.full$GO.ID == go.term,]$classicKS
        less.extreme = sum(go.bs.pvals >= go.real.p)
        tot.trials = length(go.bs.pvals)
        return(1 - (less.extreme / tot.trials))
    }

    allRes.enrich.full$bs.p = sapply(allRes.enrich.full$GO.ID, function (go) {
        get.bs.p(go)
    })

    # Write table of empirical GO p-values
    bs.out.file = gsub("txt$", paste0(ont, ".bootstrapped.GO.p-values.enrich.txt"),
        gsub("results", "reports", in.file))
    write.table(bs.p, file=bs.out.file,
        sep="\t", quote=TRUE, row.names=FALSE, col.names=FALSE)

    # Do similar thing, but this time permute gene-GO relationship

    rm(list=c("bs.p", "get.bs.p"))

    bs.p = bs.to.get.GO.pval.dist.enrich(num.iter = 1000)

    get.bs.p = function(go.term) {
        if (go.term %in% bs.p$X1) {
            go.bs.pvals = bs.p[bs.p$X1 == go.term, 2:ncol(bs.p)]
            go.real.p = allRes.enrich.full[allRes.enrich.full$GO.ID == go.term,]$classicKS
            less.extreme = sum(go.bs.pvals >= go.real.p)
            tot.trials = length(go.bs.pvals)
            res = 1 - (less.extreme / tot.trials)
        } else {
            res = NA
        }

        return(res)
    }

    allRes.enrich.full$bs.p.geneGO = sapply(allRes.enrich.full$GO.ID, function (go) {
        get.bs.p(go)
    })

    # Write table of empirical GO p-values
    bs.out.file = gsub("txt$", paste0(ont, ".bootstrapped-geneGO.GO.p-values.enrich.txt"),
        gsub("results", "reports", in.file))
    write.table(bs.p, file=bs.out.file,
        sep="\t", quote=TRUE, row.names=FALSE, col.names=FALSE)

    # ------------------------------------------------------------------------------------

    xt = xtable(head(allRes.enrich.sm, n=50))

    if (gene.mode) {
        out.file = gsub("txt$", paste0(ont, ".enrich.tex"),
            gsub("results", "reports", in.file))
    } else {
        out.file = gsub("bed$", paste0(ont, ".enrich.bySNP.tex"),
            gsub("results", "reports", in.file))
    }

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

    # ------------------------------------------------------------------------------------
    # --- Output simple text files
    # ------------------------------------------------------------------------------------

    if (gene.mode) {
        out.file.txt.overrep = gsub("txt$", paste0(ont, ".overrep.results.txt"), in.file)
        out.file.txt.enrich  = gsub("txt$", paste0(ont, ".enrich.results.txt"),  in.file)
    } else {
        out.file.txt.overrep = gsub("bed$", paste0(ont, ".overrep.results.txt"), in.file)
        out.file.txt.enrich  = gsub("bed$", paste0(ont, ".enrich.results.txt"),  in.file)
    }

    # --- Overrep

    write.table(allRes.overrep.full, out.file.txt.overrep,
        quote=FALSE, row.names=FALSE, sep="\t")

    write.table(allRes.overrep.sm, gsub(".txt$", ".sm.txt", out.file.txt.overrep),
        quote=FALSE, row.names=FALSE, sep="\t")

    # --- Enrich

    write.table(allRes.enrich.full, out.file.txt.enrich,
        quote=FALSE, row.names=FALSE, sep="\t")

    write.table(allRes.enrich.sm, gsub(".txt$", ".sm.txt", out.file.txt.enrich),
        quote=FALSE, row.names=FALSE, sep="\t")

    # --------------------------------------------------------------------------------
    # --- Plot cummulative density (ECDF) for GO terms of interest and
    # --- make output file with genes in plots of interest
    # --------------------------------------------------------------------------------

    if (ecdf.plots != "") {

        this.out.sig.genes.txt = gsub(".txt$",
            paste0(".sig.GO.sig.genes.", ont, ".txt"), in.file)
        this.out.sig.genes.tex = gsub("results", "reports",
            gsub("txt$", "tex", this.out.sig.genes.txt))

        # Touch both output files so they are created even if we have nothing to write
        system(paste0("touch ", this.out.sig.genes.txt))
        system(paste0("touch ", this.out.sig.genes.tex))

        go.terms = GO.of.interest[,1]
        go.genes = genesInTerm(GOdata, go.terms)

        if (length(go.genes)) {

            # Keyed by GO ID, but later smooshed into data.frame
            nifty.pathway.genes = list()

            for (i in 1:length(go.genes)) {
                this.term = names(go.genes[i])
                this.term.title = GO.of.interest[GO.of.interest$V1 == this.term,]$V2[1]
                this.term.genes = go.genes[this.term][[1]]
                this.term.values = input[input$V1 %in% this.term.genes, 2]
                matches.these.goes = allRes.overrep.full$GO.ID == this.term
                this.term.overrep.pval =
                    allRes.overrep.full[matches.these.goes,]$classicFisher

                # --- Gather info on genes in these regions of interest
                nifty.pathway.genes[[this.term]] = data.frame(cbind(this.term,
                    this.term.title, this.term.genes, this.term.values,
                    this.term.overrep.pval))

                # --- Plot ECDF

                this.out.ecdf.pdf = gsub("txt$",
                    paste0(ont, ".", gsub(":", "", this.term), ".pdf"),
                    gsub("results", "reports/enrichment_ecdf_plots", in.file))

                if (grepl("sizeCor", in.file)) {
                    this.out.ecdf.pdf = gsub("reports", "reports_sizeCor",
                        this.out.ecdf.pdf)
                } else if (grepl("mafCor", in.file)) {
                    this.out.ecdf.pdf = gsub("reports", "reports_mafCor",
                        this.out.ecdf.pdf)
                }

                pdf(this.out.ecdf.pdf)

                    all.ecdf = ecdf(input$V2)

                    plot(all.ecdf,
                            verticals = FALSE, do.points = TRUE,
                            main=paste0(description, ' - ', ont, "\n",
                                this.term, ' - ', this.term.title),
                            xlab="Estimate", ylab="Fraction", cex=0.25)

                    n.matches = length(this.term.values)
                    anno.ecdf = ecdf(this.term.values)

                    plot(anno.ecdf, add=TRUE, cex=0.5, col='red')

                dev.off()
            }

            # --- Write out significant genes in significant pathways
            nifty.pathway.genes.df = do.call(rbind, nifty.pathway.genes)
            names(nifty.pathway.genes.df) = c("GO_ID", "GO_title",
                "gene", "value", "GO.overrep.p")

            which.doubly.sig = top.interest.genes(nifty.pathway.genes.df$value) &
                               nifty.pathway.genes.df$GO.overrep.p < 0.05

            if (sum(which.doubly.sig) > 0) {

                sig.genes.sig.paths = nifty.pathway.genes.df[which.doubly.sig,]

                write.table(sig.genes.sig.paths, file=this.out.sig.genes.txt,
                    quote=FALSE, row.names=FALSE, sep="\t")

                # Write LaTex table

                sig.genes.sig.paths.clean = sig.genes.sig.paths[,c(1:2,5,3:4)]
                names(sig.genes.sig.paths.clean) = c("GO Term", "Description",
                    "GO overrep $p$", "Gene", "Gene Value")

                horiz.lines = c(1, which(c(sig.genes.sig.paths$GO_ID, NA) !=
                    c(NA, sig.genes.sig.paths$GO_ID)) - 1)

                xt = xtable(sig.genes.sig.paths.clean)

                sink(this.out.sig.genes.tex)

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
                        caption.placement = "top",
                        digits=c(0,0,0,-3,0,-3), align='lll|r||lr',
                        sanitize.colnames.function = identity,
                        sanitize.text.function = identity,
                        hline.after = horiz.lines)

                cat("\\end{document}", sep="\n")

                sink()
            }
        }
    }

    # ------------------------------------------------------------------------------------
    # --- Do SNP-based overrep and enrichment test
    # ----------------------------------------------------------------------------------------

    if (gene.mode == FALSE) {
        stop("SNP mode no longer implemented")
    }
}
