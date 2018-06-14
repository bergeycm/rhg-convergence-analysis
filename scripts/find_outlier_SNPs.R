#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ========================================================================================
# --- Find extreme values in SNP-based statistics and highlight those in pygmy-associated
# --- GWAS regions or with high predicted impact
# ========================================================================================

# Examples:
#  Rscript scripts/find_outlier_SNPs.R results/pbs.Batwa.GBR.anno.bed TRUE
#  Rscript scripts/find_outlier_SNPs.R results/pbs.Anda.LWK.anno.bed TRUE
#  Rscript scripts/find_outlier_SNPs.R results/pbs.Bakiga.GBR.anno.bed TRUE
#  Rscript scripts/find_outlier_SNPs.R results/pbs.BR1.LWK.anno.bed TRUE

library(ggplot2)
library(xtable)

args = commandArgs(trailingOnly=TRUE)
stat.in   = args[1]                        # e.g. results/pbs.Batwa.GBR.anno.bed
fltr.for.paralogs = as.logical(args[2])    # Should filtering for paralogs be turned on?

# Figure out sister file
if (grepl("Batwa",   stat.in)) sister.in = gsub("Batwa",  "Bakiga", stat.in)
if (grepl("Bakiga",  stat.in)) sister.in = gsub("Bakiga", "Batwa",  stat.in)
if (grepl("Anda",    stat.in)) sister.in = gsub("Anda",   "BR1",    stat.in)
if (grepl("BR1",     stat.in)) sister.in = gsub("BR1",    "Anda",   stat.in)
if (grepl("\\.twa",  stat.in)) sister.in = gsub("twa",    "kiga",   stat.in)
if (grepl("\\.kiga", stat.in)) sister.in = gsub("kiga",   "twa",    stat.in)

# Hack to let this run on Bayenv and Fst results (which do not have sister files):
# Run self as sister
if (grepl("bayenv",    stat.in)) sister.in = stat.in
if (grepl("eAGR_eRHG", stat.in)) sister.in = stat.in

flt.suffix = ""
if (fltr.for.paralogs) {
    flt.suffix = ".flt"
}

out.file.outliers = gsub(".bed$", paste0(".extreme.vals", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

out.file.gwas = gsub(".bed$", paste0(".gwas.region.extreme.vals", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

out.file.wood = gsub(".bed$", paste0(".wood.region.extreme.vals", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

out.file.mouse = gsub(".bed$", paste0(".MGI.region.extreme.vals", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

out.file.lui = gsub(".bed$", paste0(".lui.extreme.vals", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

out.file.growth = gsub(".bed$", paste0(".growthgo.extreme.vals", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

out.file.delet = gsub(".bed$", paste0(".deleterious.extreme.vals", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

out.file.fisher = gsub(".bed$", paste0(".fisher.overrep", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

stat = read.table(stat.in)
sister = read.table(sister.in)
names(stat) = names(sister) = c("chr", "start", "end", "stat.value", "genes")

# Reduce to just SNPs in genes
stat = stat[stat$genes != ".",]
sister = sister[sister$genes != ".",]

# Get rid of redundant genes in string
reduce.to.one.gene.only = function (x) {
    paste(unique(strsplit(x, ',')[[1]]), sep=",", collapse=",")
}

stat$genes   = do.call(rbind, lapply(stat$genes,   reduce.to.one.gene.only))
sister$genes = do.call(rbind, lapply(sister$genes, reduce.to.one.gene.only))

if (grepl("ihs", stat.in)) {
    stat$stat.value   = abs(stat$stat.value)
    sister$stat.value = abs(sister$stat.value)
}

# Remove infinite values
stat   = stat[abs(stat$stat.value)     != Inf,]
sister = sister[abs(sister$stat.value) != Inf,]

# ----------------------------------------------------------------------------------------
# --- Bring in filtered data from paralog-finding script, if desired
# ----------------------------------------------------------------------------------------

if (fltr.for.paralogs) {

    para.agr = read.table("results/HWE_info.eAGR.hwe.filtered.bed")
    para.rhg = read.table("results/HWE_info.eRHG.hwe.filtered.bed")

    para.agr$is.flt.agr = TRUE
    para.rhg$is.flt.rhg = TRUE

    para.all = merge(para.agr, para.rhg, all=TRUE)

    para.all$is.flt.agr[is.na(para.all$is.flt.agr)] = FALSE
    para.all$is.flt.rhg[is.na(para.all$is.flt.rhg)] = FALSE

    names(para.all)[1:3] = names(stat)[1:3]

    stat   = merge(stat,   para.all, all.x=TRUE, sort=FALSE)
    sister = merge(sister, para.all, all.x=TRUE, sort=FALSE)

    stat$is.flt = FALSE
    stat[which(stat$is.flt.agr | stat$is.flt.rhg),]$is.flt = TRUE
    sister$is.flt = FALSE
    sister[which(sister$is.flt.agr | sister$is.flt.rhg),]$is.flt = TRUE
    # table(stat$is.flt)

    p = ggplot(stat, aes(is.flt, stat.value)) +
        geom_boxplot() +
        xlab("Removed during paralog filtering?") + ylab("Statistic")

    ggsave(plot=p,
        filename=gsub(".bed$", ".paralog.flt.pdf", gsub("results", "reports", stat.in)),
        width=7, height=7)

    # Reduce to just SNPs that pass filter
    stat   = stat[!stat$is.flt, 1:5]
    sister = sister[!sister$is.flt, 1:5]
}

# ----------------------------------------------------------------------------------------
# --- Remove NA values
# ----------------------------------------------------------------------------------------

stat   = stat[!is.na(stat$stat.value),]
sister = sister[!is.na(sister$stat.value),]

# ----------------------------------------------------------------------------------------
# --- Compute empirical p-value
# ----------------------------------------------------------------------------------------

stat.ecdf   = ecdf(stat$stat.value)
sister.ecdf = ecdf(sister$stat.value)

stat$p.val   = 1 - sapply(stat$stat.value,   stat.ecdf)
sister$p.val = 1 - sapply(sister$stat.value, sister.ecdf)

# Adjust p-value
stat$p.adj   = p.adjust(stat$p.val,   method = "fdr")
sister$p.adj = p.adjust(sister$p.val, method = "fdr")

stat$anno   = ""
sister$anno = ""

# Stats are written after addition of annotations for overlap with other gene sets, etc.

# ----------------------------------------------------------------------------------------
# --- Find extreme stats that overlap height-associated regions, other height gene lists
# ----------------------------------------------------------------------------------------

# --- Perry et al 2014

gwas = read.table(paste0("data/Perry_etal_2014/GWAS_assoc_regions/",
    "Perry_etal_2014-S2-batwa_assoc.hg19.bed"))

overlap = do.call(rbind, lapply(1:nrow(gwas), function (x) {
    this = gwas[x,]
    stat[which(stat$chr == this$V1 & stat$start >= this$V2 & stat$end <= this$V3),]
}))

overlap.sister = do.call(rbind, lapply(1:nrow(gwas), function (x) {
    this = gwas[x,]
    sister[which(sister$chr == this$V1 & sister$start >= this$V2 & sister$end <= this$V3),]
}))

# Add annotation to original stat data.frame
stat[rownames(overlap),]$anno          = "P14;"
sister[rownames(overlap.sister),]$anno = "P14;"

# Find extreme values of stat
gwas.sig = overlap[overlap$p.val < 0.01,]
gwas.sig.sort = gwas.sig[order(gwas.sig$p.val, decreasing=FALSE),]

# Write Perry GWAS overlap results to file
write.table(gwas.sig.sort, file=out.file.gwas, quote=FALSE, sep="\t", row.names=FALSE)

# --- Wood et al 2014

wood.df = read.table("data/Wood_etal_2014/OMIM_height_genes.txt")

matches.wood = sapply(stat$genes, function (x) {
    genes = strsplit(x, split=",")[[1]]
    return (sum(genes %in% wood.df$V1) > 0)
})

matches.wood.sister = sapply(sister$genes, function (x) {
    genes = strsplit(x, split=",")[[1]]
    return (sum(genes %in% wood.df$V1) > 0)
})

# Add annotation to original stat data.frame
stat[matches.wood,]$anno          = paste0(stat[matches.wood,]$anno,          "W14;")
sister[matches.wood.sister,]$anno = paste0(sister[matches.wood.sister,]$anno, "W14;")

wood        = stat[matches.wood,]
wood.sister = sister[matches.wood.sister,]

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

matches.mgi.growth = sapply(stat$genes, function (x) {
    genes = strsplit(x, split=",")[[1]]
    return (sum(genes %in% mgi.growth.genes) > 0)
})

matches.mgi.growth.sister = sapply(sister$genes, function (x) {
    genes = strsplit(x, split=",")[[1]]
    return (sum(genes %in% mgi.growth.genes) > 0)
})

# Add annotation to original stat data.frame
stat[matches.mgi.growth,]$anno          = paste0(stat[matches.mgi.growth,]$anno,
                                                 "MGI;")
sister[matches.mgi.growth.sister,]$anno = paste0(sister[matches.mgi.growth.sister,]$anno,
                                                 "MGI;")

mouse        = stat[matches.mgi.growth,]
mouse.sister = sister[matches.mgi.growth.sister,]

mouse.sig = mouse[mouse$p.val < 0.01,]
mouse.sig.sort = mouse.sig[order(mouse.sig$p.val, decreasing=FALSE),]

write.table(mouse.sig.sort, file=out.file.mouse, quote=FALSE, sep="\t", row.names=FALSE)

# --- Lui et al 2012 - growth plate expressed from table S2

lui = read.table("data/Lui_etal_2012/growth_plate_expressed.txt")
lui$V1 = toupper(lui$V1)

matches.lui = sapply(stat$genes, function (x) {
    genes = strsplit(x, split=",")[[1]]
    return (sum(genes %in% lui$V1) > 0)
})

matches.lui.sister = sapply(sister$genes, function (x) {
    genes = strsplit(x, split=",")[[1]]
    return (sum(genes %in% lui$V1) > 0)
})

# Add annotation to original stat data.frame
stat[matches.lui,]$anno          = paste0(stat[matches.lui,]$anno,          "L12;")
sister[matches.lui.sister,]$anno = paste0(sister[matches.lui.sister,]$anno, "L12;")

lui        = stat[matches.lui,]
lui.sister = sister[matches.lui.sister,]

# Annotate stat with whether it overlaps Lui
lui.sig = lui[lui$p.val < 0.01,]
lui.sig.sort = lui.sig[order(lui.sig$p.val, decreasing=FALSE),]

write.table(lui.sig.sort, file=out.file.lui, quote=FALSE, sep="\t", row.names=FALSE)

# --- Growth GO term (GO:0040007)

growth.go = read.table("data/GO_growth_genes/GO_0040007_genes.txt", sep="\t")

matches.growth = sapply(stat$genes, function (x) {
    genes = strsplit(x, split=",")[[1]]
    return (sum(genes %in% growth.go$V3) > 0)
})

matches.growth.sister = sapply(sister$genes, function (x) {
    genes = strsplit(x, split=",")[[1]]
    return (sum(genes %in% growth.go$V3) > 0)
})

# Add annotation to original stat data.frame
stat[matches.growth,]$anno          = paste0(stat[matches.growth,]$anno,          "GR;")
sister[matches.growth.sister,]$anno = paste0(sister[matches.growth.sister,]$anno, "GR;")

growth        = stat[matches.growth,]
growth.sister = sister[matches.growth.sister,]

# Annotate stat with whether it overlaps growth GO
growth.sig = growth[growth$p.val < 0.01,]
growth.sig.sort = growth.sig[order(growth.sig$p.val, decreasing=FALSE),]

write.table(growth.sig.sort, file=out.file.growth, quote=FALSE, sep="\t", row.names=FALSE)

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
# --- Do Fisher's test to see if significant high values have
# --- more annotations than expected
# ----------------------------------------------------------------------------------------

p.val.cutoffs = c(0.1, 0.05, 0.01, 0.001)

tot.count        = nrow(stat)
tot.count.sister = nrow(sister)

# ----------------------------------------------------------------------------------------

# --- Perry et al 2014

fisher.gwas = lapply(p.val.cutoffs, function (x) {
    p.val.cutoff = x

    sig.count = nrow(stat[stat$p.val <= p.val.cutoff,])

    sig.count.gwas = nrow(overlap[overlap$p.val <= p.val.cutoff,])
    tot.count.gwas = nrow(overlap)
    test.mat = matrix(c(sig.count.gwas,             tot.count.gwas,
                        sig.count - sig.count.gwas, tot.count - tot.count.gwas), ncol=2)
    list(c(test.mat), fisher.test(test.mat))
})
names(fisher.gwas) = paste0("p-cutoff_", p.val.cutoffs)

fisher.gwas.sister = lapply(p.val.cutoffs, function (x) {
    p.val.cutoff = x

    sig.count = nrow(sister[sister$p.val <= p.val.cutoff,])

    sig.count.gwas = nrow(overlap.sister[overlap.sister$p.val <= p.val.cutoff,])
    tot.count.gwas = nrow(overlap.sister)
    test.mat = matrix(c(sig.count.gwas,             tot.count.gwas,
                        sig.count - sig.count.gwas, tot.count - tot.count.gwas), ncol=2)
    list(c(test.mat), fisher.test(test.mat))
})
names(fisher.gwas.sister) = paste0("p-cutoff_", p.val.cutoffs)

# ----------------------------------------------------------------------------------------

# --- Wood et al 2014

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

fisher.wood.sister = lapply(p.val.cutoffs, function (x) {
    p.val.cutoff = x

    sig.count = nrow(sister[sister$p.val <= p.val.cutoff,])

    sig.count.wood = nrow(wood.sister[wood.sister$p.val <= p.val.cutoff,])
    tot.count.wood = nrow(wood.sister)
    test.mat = matrix(c(sig.count.wood,             tot.count.wood,
                        sig.count - sig.count.wood, tot.count - tot.count.wood), ncol=2)
    list(c(test.mat), fisher.test(test.mat))
})
names(fisher.wood.sister) = paste0("p-cutoff_", p.val.cutoffs)

# ----------------------------------------------------------------------------------------

# --- MGI

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

fisher.mgi.sister = lapply(p.val.cutoffs, function (x) {
    p.val.cutoff = x

    sig.count = nrow(sister[sister$p.val <= p.val.cutoff,])

    sig.count.mgi = nrow(mouse.sister[mouse.sister$p.val <= p.val.cutoff,])
    tot.count.mgi = nrow(mouse.sister)
    test.mat = matrix(c(sig.count.mgi,             tot.count.mgi,
                        sig.count - sig.count.mgi, tot.count - tot.count.mgi), ncol=2)
    list(c(test.mat), fisher.test(test.mat))
})
names(fisher.mgi.sister) = paste0("p-cutoff_", p.val.cutoffs)

# ----------------------------------------------------------------------------------------

# --- Lui et al 2014

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

fisher.lui.sister = lapply(p.val.cutoffs, function (x) {
    p.val.cutoff = x

    sig.count = nrow(sister[sister$p.val <= p.val.cutoff,])

    sig.count.lui = nrow(lui.sister[lui.sister$p.val <= p.val.cutoff,])
    tot.count.lui = nrow(lui.sister)
    test.mat = matrix(c(sig.count.lui,             tot.count.lui,
                        sig.count - sig.count.lui, tot.count - tot.count.lui), ncol=2)
    list(c(test.mat), fisher.test(test.mat))
})
names(fisher.lui.sister) = paste0("p-cutoff_", p.val.cutoffs)

# ----------------------------------------------------------------------------------------

# --- GO Growth

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

fisher.growth.sister = lapply(p.val.cutoffs, function (x) {
    p.val.cutoff = x

    sig.count = nrow(sister[sister$p.val <= p.val.cutoff,])

    sig.count.growth = nrow(growth.sister[growth.sister$p.val <= p.val.cutoff,])
    tot.count.growth = nrow(growth.sister)
    test.mat = matrix(c(sig.count.growth,             tot.count.growth,
                        sig.count - sig.count.growth, tot.count - tot.count.growth), ncol=2)
    list(c(test.mat), fisher.test(test.mat))
})
names(fisher.growth.sister) = paste0("p-cutoff_", p.val.cutoffs)

# ----------------------------------------------------------------------------------------

# --- Write Fisher results to file

fisher.gwas.df = cbind("Perry et. al 2014", p.val.cutoffs,
    data.frame(t(sapply(fisher.gwas, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fisher.wood.df = cbind("Wood et. al 2014", p.val.cutoffs,
    data.frame(t(sapply(fisher.wood, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fisher.mgi.df = cbind("MGI Phenotypes", p.val.cutoffs,
    data.frame(t(sapply(fisher.mgi, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fisher.lui.df = cbind("Lui et. al 2012", p.val.cutoffs,
    data.frame(t(sapply(fisher.lui, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fisher.growth.df = cbind("GO growth", p.val.cutoffs,
    data.frame(t(sapply(fisher.growth, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fish.names = c("dataset", "p.val.cutoffs", "sig.count.this", "tot.count.this",
              "sig.count.notthis", "tot.count.notthis",
              "pvalue", "conf.int", "estimate", "null.value")

names(fisher.gwas.df) = names(fisher.wood.df) = names(fisher.mgi.df) =
    names(fisher.lui.df) = names(fisher.growth.df) = fish.names
fisher.df = rbind(fisher.gwas.df, fisher.wood.df, fisher.mgi.df,
                 fisher.lui.df, fisher.growth.df)

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

# --- Do same for sister Fisher results, but do not write to file

fisher.gwas.df.sis = cbind("Perry et. al 2014", p.val.cutoffs,
    data.frame(t(sapply(fisher.gwas.sister, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fisher.wood.df.sis = cbind("Wood et. al 2014", p.val.cutoffs,
    data.frame(t(sapply(fisher.wood.sister, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fisher.mgi.df.sis = cbind("MGI Phenotypes", p.val.cutoffs,
    data.frame(t(sapply(fisher.mgi.sister, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fisher.lui.df.sis = cbind("Lui et. al 2012", p.val.cutoffs,
    data.frame(t(sapply(fisher.lui.sister, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fisher.growth.df.sis = cbind("GO growth", p.val.cutoffs,
    data.frame(t(sapply(fisher.growth.sister, function(x) { c(x[1:3][[1]], x[1:3][[2]]) }))))[,1:10]

fish.names = c("dataset", "p.val.cutoffs", "sig.count.this", "tot.count.this",
              "sig.count.notthis", "tot.count.notthis",
              "pvalue", "conf.int", "estimate", "null.value")
names(fisher.gwas.df.sis) = names(fisher.wood.df.sis) = names(fisher.mgi.df.sis) =
    names(fisher.lui.df.sis) = names(fisher.growth.df.sis) = fish.names
fisher.df.sis = rbind(fisher.gwas.df.sis, fisher.wood.df.sis, fisher.mgi.df.sis,
                     fisher.lui.df.sis, fisher.growth.df.sis)

fisher.df.sis$sig.count.this = unlist(fisher.df.sis$sig.count.this)
fisher.df.sis$tot.count.this = unlist(fisher.df.sis$tot.count.this)
fisher.df.sis$sig.count.notthis = unlist(fisher.df.sis$sig.count.notthis)
fisher.df.sis$tot.count.notthis = unlist(fisher.df.sis$tot.count.notthis)
fisher.df.sis$statistic = unlist(fisher.df.sis$statistic)
fisher.df.sis$parameter = unlist(fisher.df.sis$parameter)
fisher.df.sis$pvalue = unlist(fisher.df.sis$pvalue)
fisher.df.sis$method = unlist(fisher.df.sis$method)

# --- Write to LaTeX too

xt.fisher = xtable(fisher.df[,-c(8,10)], digits=c(0,0,3,0,0,0,0,-3,-3))

print(xt.fisher, type="latex", file=gsub(".txt$", ".tex",  out.file.fisher),
    include.rownames=FALSE, tabular.environment="longtable",
    floating=FALSE,
    size="\\fontsize{9pt}{10pt}\\selectfont")

# --- Compare proportions of focus and sister groups

df = cbind(fisher.df[,1:4], fisher.df.sis[,3:4])

prop.pvals = sapply(1:nrow(df), function (row.idx) {
    prop.test(x = as.numeric(df[row.idx,c(3,5)]),
              n = as.numeric(df[row.idx,c(4,6)]),
              alternative="greater")$p.value
})

df.res = cbind(df, prop.pvals)

out.file.proptest = gsub(".bed$", paste0(".proptest", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

write.table(df.res, out.file.proptest)

# ----------------------------------------------------------------------------------------
# --- Test for shift in disribution
# ----------------------------------------------------------------------------------------

# Uses independent 2-group Mann-Whitney U Test to see if growth genes have shifted
# PBS in the focus group
gwas.wilcox   = wilcox.test(overlap$p.val,    overlap.sister$p.val, alternative="less")
wood.wilcox   = wilcox.test(wood$p.val,       wood.sister$p.val,    alternative="less")
mgi.wilcox    = wilcox.test(mouse$p.val,      mouse.sister$p.val,   alternative="less")
lui.wilcox    = wilcox.test(lui$p.val,        lui.sister$p.val,     alternative="less")
growth.wilcox = wilcox.test(growth$p.val,     growth.sister$p.val,  alternative="less")

out.file.wilcox = gsub(".bed$", paste0(".wilcoxshift", flt.suffix, ".txt"),
    gsub("results", "reports", stat.in))

all.wilcox = cbind(unique(df.res$dataset),
    c(gwas.wilcox$p.value, wood.wilcox$p.value, mgi.wilcox$p.value,
      lui.wilcox$p.value, growth.wilcox$p.value))

write.table(all.wilcox, out.file.wilcox, quote=FALSE,
    sep="\t", col.names=FALSE, row.names=FALSE)

# ----------------------------------------------------------------------------------------
# --- Write LaTeX tables
# ----------------------------------------------------------------------------------------

first.col.names = c("SNP", "Value", "Gene", "p", "adj. p")

stat.sig.sort.clean = cbind(
    paste(stat.sig.sort$chr, stat.sig.sort$start, sep=":"),
    stat.sig.sort[4:ncol(stat.sig.sort)])
names(stat.sig.sort.clean) = c(first.col.names, "")

gwas.sig.sort.clean = cbind(
    paste(gwas.sig.sort$chr, gwas.sig.sort$start, sep=":"),
    gwas.sig.sort[4:ncol(gwas.sig.sort)])
names(gwas.sig.sort.clean) = c(first.col.names, "")

xt.stat  = xtable(stat.sig.sort.clean, digits=c(0,0,3,0,-3,-3,0))
xt.gwas  = xtable(gwas.sig.sort.clean, digits=c(0,0,3,0,-3,-3,0))

print(xt.stat, type="html",  file=gsub(".txt$", ".html", out.file.outliers),
    include.rownames=FALSE)
print(xt.stat, type="latex", file=gsub(".txt$", ".tex",  out.file.outliers),
    include.rownames=FALSE, tabular.environment="longtable",
    floating=FALSE,
    size="\\fontsize{9pt}{10pt}\\selectfont")

print(xt.gwas, type="html",  file=gsub(".txt$", ".html", out.file.gwas),
    include.rownames=FALSE)
print(xt.gwas, type="latex", file=gsub(".txt$", ".tex",  out.file.gwas),
    include.rownames=FALSE, tabular.environment="longtable",
    floating=FALSE,
    size="\\fontsize{9pt}{10pt}\\selectfont")
