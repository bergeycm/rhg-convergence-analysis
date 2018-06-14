#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ========================================================================================
# --- Find extreme genes
# ========================================================================================

# Examples:
#  Rscript scripts/find_outlier_genes.R results/pbs.Batwa.GBR.pvals.txt FALSE

library(ggplot2)
library(xtable)

args = commandArgs(trailingOnly=TRUE)
stat.in = args[1]                          # e.g. results/pbs.Batwa.GBR.pvals.txt
fltr.for.paralogs = as.logical(args[2])    # Should filtering for paralogs be turned on?

flt.suffix = ""
if (fltr.for.paralogs) {
    flt.suffix = ".flt"
}

out.file.outliers = gsub(".txt$", paste0(".extreme.vals", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

out.file.gwas = gsub(".txt$", paste0(".gwas.region.extreme.vals", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

out.file.wood = gsub(".txt$", paste0(".wood.region.extreme.vals", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

out.file.mouse = gsub(".txt$", paste0(".MGI.region.extreme.vals", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

out.file.lui = gsub(".txt$", paste0(".lui.extreme.vals", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

out.file.growth = gsub(".txt$", paste0(".growthgo.extreme.vals", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

out.file.ophid_igf1 = gsub(".txt$", paste0(".ophid_igf1.extreme.vals", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

out.file.ophid_gh1 = gsub(".txt$", paste0(".ophid_gh1.extreme.vals", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

out.file.fisher = gsub(".txt$", paste0(".fisher.overrep", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

stat = read.table(stat.in)
names(stat) = c("gene", "p.val")

# ----------------------------------------------------------------------------------------
# --- Bring in filtered data from paralog-finding script, if desired
# ----------------------------------------------------------------------------------------

if (fltr.for.paralogs) {

    para.agr = read.table("results/HWE_info.eAGR.hwe.filtered.justGenes.txt")
    para.rhg = read.table("results/HWE_info.eRHG.hwe.filtered.justGenes.txt")

    para.agr$is.flt.agr = TRUE
    para.rhg$is.flt.rhg = TRUE

    para.all = merge(para.agr, para.rhg, all=TRUE)

    para.all$is.flt.agr[is.na(para.all$is.flt.agr)] = FALSE
    para.all$is.flt.rhg[is.na(para.all$is.flt.rhg)] = FALSE

    names(para.all)[1] = names(stat)[1]

    stat = merge(stat, para.all, all.x=TRUE, sort=FALSE)

    stat$is.flt = FALSE
    stat[which(stat$is.flt.agr | stat$is.flt.rhg),]$is.flt = TRUE
    # table(stat$is.flt)

    p = ggplot(stat, aes(is.flt, p.val)) +
        geom_boxplot() +
        xlab("Removed during paralog filtering?") + ylab("Statistic")

    ggsave(plot=p,
        filename=gsub(".txt$", ".paralog.flt.pdf", gsub("results", "reports", stat.in)),
        width=7, height=7)

    # Reduce to just SNPs that pass filter
    stat = stat[!stat$is.flt, 1:2]
}

# ----------------------------------------------------------------------------------------
# --- Remove NA values
# ----------------------------------------------------------------------------------------

stat = stat[!is.na(stat$p.val),]

# ----------------------------------------------------------------------------------------
# --- Adjust empirical p-value
# ----------------------------------------------------------------------------------------

# Adjust p-value
stat$p.adj = p.adjust(stat$p.val, method = "fdr")

stat$anno = ""

# Stats are written after addition of annotations for overlap with other gene sets, etc.

# ----------------------------------------------------------------------------------------
# --- Find extreme stats that overlap height-associated regions, other height gene lists
# ----------------------------------------------------------------------------------------

# --- Perry et al 2014

gwas = read.table(paste0("data/Perry_etal_2014/GWAS_assoc_regions/",
    "Perry_etal_2014-S2-batwa_assoc.hg19.justGenes.txt"))

matches.gwas = sapply(stat$gene, function (x) {
    genes = strsplit(x, split=",")[[1]]
    return (sum(genes %in% gwas$V1) > 0)
})

# Add annotation to original stat data.frame
stat[matches.gwas,]$anno = paste0(stat[matches.gwas,]$anno, "P14;")

gwas = stat[matches.gwas,]

# Find extreme values of stat
gwas.sig = gwas[gwas$p.val < 0.01,]
gwas.sig.sort = gwas.sig[order(gwas.sig$p.val, decreasing=FALSE),]

# Write Perry GWAS overlap results to file
write.table(gwas.sig.sort, file=out.file.gwas, quote=FALSE, sep="\t", row.names=FALSE)

# --- Wood et al 2014

wood = read.table("data/Wood_etal_2014/OMIM_height_genes.txt")

matches.wood = sapply(stat$gene, function (x) {
    genes = strsplit(x, split=",")[[1]]
    return (sum(genes %in% wood$V1) > 0)
})

# Add annotation to original stat data.frame
stat[matches.wood,]$anno = paste0(stat[matches.wood,]$anno, "W14;")

wood = stat[matches.wood,]

# Annotate stat with whether it overlaps Wood
wood.sig = wood[wood$p.val < 0.01,]
wood.sig.sort = wood.sig[order(wood.sig$p.val, decreasing=FALSE),]

write.table(wood.sig.sort, file=out.file.wood, quote=FALSE, sep="\t", row.names=FALSE)

# --- MGI - Mouse/Human Orthology with Phenotype Annotations (tab-delimited)

mgi = read.table("data/MGI/HMD_HumanPhenotype.rpt", fill=TRUE)

mgi$annos = sapply(1:nrow(mgi), function(x) {
    this=mgi[x,];
    gsub(";+", ";", paste(this[7:length(this)], collapse=";"))
})

# growth/size phenotype
growth.mp = "MP:0005378"
mgi$growth = sapply(1:nrow(mgi), function(x) {
    this=mgi[x,];
    growth.mp %in% this[7:length(this)]
})

mgi.growth.genes = mgi[mgi$growth,]$V1

matches.mgi.growth = sapply(stat$gene, function (x) {
    genes = strsplit(x, split=",")[[1]]
    return (sum(genes %in% mgi.growth.genes) > 0)
})

# Add annotation to original stat data.frame
stat[matches.mgi.growth,]$anno = paste0(stat[matches.mgi.growth,]$anno, "MGI;")

mouse = stat[matches.mgi.growth,]

mouse.sig = mouse[mouse$p.val < 0.01,]
mouse.sig.sort = mouse.sig[order(mouse.sig$p.val, decreasing=FALSE),]

write.table(mouse.sig.sort, file=out.file.mouse, quote=FALSE, sep="\t", row.names=FALSE)

# --- Lui et al 2012 - growth plate expressed from table S2

lui = read.table("data/Lui_etal_2012/growth_plate_expressed.txt")
lui$V1 = toupper(lui$V1)

matches.lui = sapply(stat$gene, function (x) {
    genes = strsplit(x, split=",")[[1]]
    return (sum(genes %in% lui$V1) > 0)
})

# Add annotation to original stat data.frame
stat[matches.lui,]$anno = paste0(stat[matches.lui,]$anno, "L12;")

lui = stat[matches.lui,]

# Annotate stat with whether it overlaps Lui
lui.sig = lui[lui$p.val < 0.01,]
lui.sig.sort = lui.sig[order(lui.sig$p.val, decreasing=FALSE),]

write.table(lui.sig.sort, file=out.file.lui, quote=FALSE, sep="\t", row.names=FALSE)

# --- Growth GO term (GO:0040007)

growth.go = read.table("data/GO_growth_genes/GO_0040007_genes.txt", sep="\t")

matches.growth = sapply(stat$gene, function (x) {
    genes = strsplit(x, split=",")[[1]]
    return (sum(genes %in% growth.go$V3) > 0)
})

# Add annotation to original stat data.frame
stat[matches.growth,]$anno = paste0(stat[matches.growth,]$anno, "GR;")

growth = stat[matches.growth,]

# Annotate stat with whether it overlaps growth GO genes
growth.sig = growth[growth$p.val < 0.01,]
growth.sig.sort = growth.sig[order(growth.sig$p.val, decreasing=FALSE),]

write.table(growth.sig.sort, file=out.file.growth, quote=FALSE, sep="\t", row.names=FALSE)

# --- OPHID genes associated with IGF1

ophid_igf1 = read.table("data/OPHID_IGF1.txt", sep="\t", header=TRUE)

matches.ophid_igf1 = sapply(stat$gene, function (x) {
    genes = strsplit(x, split=",")[[1]]
    return (sum(genes %in% ophid_igf1$Partner.Symbol) > 0)
})

# Add annotation to original stat data.frame
stat[matches.ophid_igf1,]$anno = paste0(stat[matches.ophid_igf1,]$anno, "IGF1;")

ophid_igf1 = stat[matches.ophid_igf1,]

# Annotate stat with whether it overlaps OPHID IGF1 GO genes
ophid_igf1.sig = ophid_igf1[ophid_igf1$p.val < 0.01,]
ophid_igf1.sig.sort = ophid_igf1.sig[order(ophid_igf1.sig$p.val, decreasing=FALSE),]

write.table(ophid_igf1.sig.sort, file=out.file.ophid_igf1,
    quote=FALSE, sep="\t", row.names=FALSE)

# --- OPHID genes associated with GH1

ophid_gh1 = read.table("data/OPHID_GH1.txt", sep="\t", header=TRUE)

matches.ophid_gh1 = sapply(stat$gene, function (x) {
    genes = strsplit(x, split=",")[[1]]
    return (sum(genes %in% ophid_gh1$Partner.Symbol) > 0)
})

# Add annotation to original stat data.frame
stat[matches.ophid_gh1,]$anno = paste0(stat[matches.ophid_gh1,]$anno, "GH1;")

ophid_gh1 = stat[matches.ophid_gh1,]

# Annotate stat with whether it overlaps OPHID GH1 genes
ophid_gh1.sig = ophid_gh1[ophid_gh1$p.val < 0.01,]
ophid_gh1.sig.sort = ophid_gh1.sig[order(ophid_gh1.sig$p.val, decreasing=FALSE),]

write.table(ophid_gh1.sig.sort, file=out.file.ophid_gh1,
    quote=FALSE, sep="\t", row.names=FALSE)

# ----------------------------------------------------------------------------------------
# --- Write significant stats
# ----------------------------------------------------------------------------------------

stat.sort = stat[order(stat$p.val, decreasing=FALSE),]

# Write table of all SNPs, not just significant ones
write.table(stat.sort, file=gsub(".txt$", ".all.txt", out.file.outliers),
    quote=FALSE, sep="\t", row.names=FALSE)

# Find and write out significant extreme values
stat.sig = stat[stat$p.val < 0.0001,]

stat.sig.sort = stat.sig[order(stat.sig$p.val, decreasing=FALSE),]

write.table(stat.sig.sort, file=out.file.outliers,
    quote=FALSE, sep="\t", row.names=FALSE)

# ----------------------------------------------------------------------------------------
# --- Do Fisher's test to see if significant high values have more annotations
# --- than expected
# ----------------------------------------------------------------------------------------

p.val.cutoffs = c(0.05, 0.01, 0.001)

tot.count = nrow(stat)

# Perry et al 2014
fisher.gwas = lapply(p.val.cutoffs, function (x) {
    p.val.cutoff = x

    sig.count = nrow(stat[stat$p.val <= p.val.cutoff,])

    sig.count.gwas = nrow(gwas[gwas$p.val <= p.val.cutoff,])
    tot.count.gwas = nrow(gwas)
    test.mat = matrix(c(sig.count.gwas,             tot.count.gwas,
                        sig.count - sig.count.gwas, tot.count - tot.count.gwas), ncol=2)
    list(c(test.mat), fisher.test(test.mat))
})
names(fisher.gwas) = paste0("p-cutoff_", p.val.cutoffs)

# Wood et al 2014
fisher.wood = lapply(p.val.cutoffs, function (x) {
    p.val.cutoff = x

    sig.count = nrow(stat[stat$p.val <= p.val.cutoff,])

    sig.count.wood = nrow(wood[wood$p.val <= p.val.cutoff,])
    tot.count.wood = nrow(wood)
    test.mat = matrix(c(sig.count.wood,             tot.count.wood,
                        sig.count - sig.count.wood, tot.count - tot.count.wood), ncol=2)
    list(c(test.mat), fisher.test(test.mat))

})
names(fisher.wood) = paste0("p-cutoff_", p.val.cutoffs)

# MGI
fisher.mgi = lapply(p.val.cutoffs, function (x) {
    p.val.cutoff = x

    sig.count = nrow(stat[stat$p.val <= p.val.cutoff,])

    sig.count.mgi = nrow(mouse[mouse$p.val <= p.val.cutoff,])
    tot.count.mgi = nrow(mouse)
    test.mat = matrix(c(sig.count.mgi,             tot.count.mgi,
                        sig.count - sig.count.mgi, tot.count - tot.count.mgi), ncol=2)
    list(c(test.mat), fisher.test(test.mat))
})
names(fisher.mgi) = paste0("p-cutoff_", p.val.cutoffs)

# Lui et al 2014
fisher.lui = lapply(p.val.cutoffs, function (x) {
    p.val.cutoff = x

    sig.count = nrow(stat[stat$p.val <= p.val.cutoff,])

    sig.count.lui = nrow(lui[lui$p.val <= p.val.cutoff,])
    tot.count.lui = nrow(lui)
    test.mat = matrix(c(sig.count.lui,             tot.count.lui,
                        sig.count - sig.count.lui, tot.count - tot.count.lui), ncol=2)
    list(c(test.mat), fisher.test(test.mat))
})
names(fisher.lui) = paste0("p-cutoff_", p.val.cutoffs)

# Growth GO
fisher.growth = lapply(p.val.cutoffs, function (x) {
    p.val.cutoff = x

    sig.count = nrow(stat[stat$p.val <= p.val.cutoff,])

    sig.count.growth = nrow(growth[growth$p.val <= p.val.cutoff,])
    tot.count.growth = nrow(growth)
    test.mat = matrix(c(sig.count.growth,             tot.count.growth,
                        sig.count - sig.count.growth, tot.count - tot.count.growth), ncol=2)
    list(c(test.mat), fisher.test(test.mat))
})
names(fisher.growth) = paste0("p-cutoff_", p.val.cutoffs)

# OPHID IGF1
fisher.ophid_igf1 = lapply(p.val.cutoffs, function (x) {
    p.val.cutoff = x

    sig.count = nrow(stat[stat$p.val <= p.val.cutoff,])

    sig.count.ophid_igf1 = nrow(ophid_igf1[ophid_igf1$p.val <= p.val.cutoff,])
    tot.count.ophid_igf1 = nrow(ophid_igf1)
    test.mat = matrix(c(sig.count.ophid_igf1,             tot.count.ophid_igf1,
                        sig.count - sig.count.ophid_igf1, tot.count - tot.count.ophid_igf1), ncol=2)
    list(c(test.mat), fisher.test(test.mat))
})
names(fisher.ophid_igf1) = paste0("p-cutoff_", p.val.cutoffs)


# OPHID GH1
fisher.ophid_gh1 = lapply(p.val.cutoffs, function (x) {
    p.val.cutoff = x

    sig.count = nrow(stat[stat$p.val <= p.val.cutoff,])

    sig.count.ophid_gh1 = nrow(ophid_gh1[ophid_gh1$p.val <= p.val.cutoff,])
    tot.count.ophid_gh1 = nrow(ophid_gh1)
    test.mat = matrix(c(sig.count.ophid_gh1,             tot.count.ophid_gh1,
                        sig.count - sig.count.ophid_gh1, tot.count - tot.count.ophid_gh1), ncol=2)
    list(c(test.mat), fisher.test(test.mat))
})
names(fisher.ophid_gh1) = paste0("p-cutoff_", p.val.cutoffs)

# --- Write Fisher results to file

fisher.gwas.df = cbind("Perry et. al 2014", p.val.cutoffs,
    data.frame(t(sapply(fisher.gwas, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fisher.wood.df = cbind("Wood et. al 2014", p.val.cutoffs,
    data.frame(t(sapply(fisher.wood, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fisher.mgi.df = cbind("MGI Phenotypes", p.val.cutoffs,
    data.frame(t(sapply(fisher.mgi, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fisher.lui.df = cbind("Lui et. al 2012", p.val.cutoffs,
    data.frame(t(sapply(fisher.lui, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fisher.growth.df = cbind("Growth GO", p.val.cutoffs,
    data.frame(t(sapply(fisher.growth, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fisher.ophid_igf1.df = cbind("OPHID IGF1", p.val.cutoffs,
    data.frame(t(sapply(fisher.ophid_igf1, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fisher.ophid_gh1.df = cbind("OPHID GH1", p.val.cutoffs,
    data.frame(t(sapply(fisher.ophid_gh1, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fish.names = c("dataset", "p.val.cutoffs", "sig.count.this", "tot.count.this",
              "sig.count.notthis", "tot.count.notthis",
              "pvalue", "conf.int", "estimate", "null.value")
names(fisher.gwas.df) = names(fisher.wood.df) =
    names(fisher.mgi.df) = names(fisher.lui.df) = names(fisher.growth.df) =
    names(fisher.ophid_igf1.df) = names(fisher.ophid_gh1.df) = fish.names

fisher.df = rbind(fisher.gwas.df, fisher.wood.df, fisher.mgi.df,
                 fisher.lui.df, fisher.growth.df,
                 fisher.ophid_igf1.df, fisher.ophid_gh1.df)

fisher.df$sig.count.this = unlist(fisher.df$sig.count.this)
fisher.df$tot.count.this = unlist(fisher.df$tot.count.this)
fisher.df$sig.count.notthis = unlist(fisher.df$sig.count.notthis)
fisher.df$tot.count.notthis = unlist(fisher.df$tot.count.notthis)
fisher.df$statistic = unlist(fisher.df$statistic)
fisher.df$parameter = unlist(fisher.df$parameter)
fisher.df$pvalue = unlist(fisher.df$pvalue)
fisher.df$method = unlist(fisher.df$method)

write.table(fisher.df[,-c(8:10)], file=out.file.fisher,
    quote=FALSE, sep="\t", row.names=FALSE)

# --- Write to LaTeX too

xt.fisher = xtable(fisher.df[,-c(8:10)], digits=c(0,0,3,0,0,0,0,-3))

print(xt.fisher, type="latex", file=gsub(".txt$", ".tex",  out.file.fisher),
    include.rownames=FALSE, tabular.environment="longtable",
    floating=FALSE,
    size="\\fontsize{9pt}{10pt}\\selectfont")

# ----------------------------------------------------------------------------------------
# --- Write LaTeX tables
# ----------------------------------------------------------------------------------------

first.col.names = c("Gene", "p", "adj. p")

stat.sig.sort.clean = stat.sig.sort
names(stat.sig.sort.clean) = c(first.col.names, "")

xt.stat  = xtable(stat.sig.sort.clean, digits=c(0,0,-3,-3,0))

print(xt.stat, type="html",  file=gsub(".txt$", ".html", out.file.outliers),
    include.rownames=FALSE)
print(xt.stat, type="latex", file=gsub(".txt$", ".tex",  out.file.outliers),
    include.rownames=FALSE, tabular.environment="longtable",
    floating=FALSE,
    size="\\fontsize{9pt}{10pt}\\selectfont")

# ----------------------------------------------------------------------------------------
# --- Test for shift in distribution
# ----------------------------------------------------------------------------------------

# Uses K-S Test to see if growth genes have shifted PBS in the focus group
gwas.KS       = ks.test(stat[matches.gwas,]$p.val,
                            stat[!matches.gwas,]$p.val,
                            alternative="greater")
wood.KS       = ks.test(stat[matches.wood,]$p.val,
                            stat[!matches.wood,]$p.val,
                            alternative="greater")
mgi.KS        = ks.test(stat[matches.mgi.growth,]$p.val,
                            stat[!matches.mgi.growth,]$p.val,
                            alternative="greater")
lui.KS        = ks.test(stat[matches.lui,]$p.val,
                            stat[!matches.lui,]$p.val,
                            alternative="greater")
growth.KS     = ks.test(stat[matches.growth,]$p.val,
                            stat[!matches.growth,]$p.val,
                            alternative="greater")
ophid_igf1.KS = ks.test(stat[matches.ophid_igf1,]$p.val,
                            stat[!matches.ophid_igf1,]$p.val,
                            alternative="greater")
ophid_gh1.KS  = ks.test(stat[matches.ophid_gh1,]$p.val,
                            stat[!matches.ophid_gh1,]$p.val,
                            alternative="greater")

out.file.KS = gsub(".txt$", paste0(".KSshift", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

all.KS = cbind(unique(fisher.df$dataset),
    c(gwas.KS$p.value, wood.KS$p.value, mgi.KS$p.value,
      lui.KS$p.value, growth.KS$p.value,
      ophid_igf1.KS$p.value, ophid_gh1.KS$p.value))

write.table(all.KS, out.file.KS, quote=FALSE,
    sep="\t", col.names=FALSE, row.names=FALSE)
