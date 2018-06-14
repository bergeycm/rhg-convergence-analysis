#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

library(ggplot2)

# ========================================================================================
# --- Plot stats, separated those that overlap known prior regions of interest
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Bring in results of parse_windowed_stats
# ----------------------------------------------------------------------------------------

stats = read.table(file="results/basic_stats.anno.txt", sep="\t", header=TRUE)

# ----------------------------------------------------------------------------------------
# --- Do plotting
# ----------------------------------------------------------------------------------------

out.dir = "reports/distribution_comparisons/"
dir.create(out.dir, showWarnings=FALSE)

# --- Perry et al 2014

# GWAS
p = ggplot(stats, aes(x=fst.all.WEIGHTED_FST, ..density.., 
        fill=overlap.perry.gwas)) + 
    geom_histogram(alpha=0.5, position="identity") + 
    xlab(expression(Average~F[ST])) + guides(fill=FALSE) + theme_bw()

ggsave(p, file=paste0(out.dir, "perry.gwas.fst.all.pdf"), height=7, width=7)

p = ggplot(stats, aes(x=fst.ns.WEIGHTED_FST, ..density.., 
        fill=overlap.perry.gwas)) + 
    geom_histogram(alpha=0.2, position="identity") + 
    xlab(expression(Average~F[ST]))

ggsave(p, file=paste0(out.dir, "perry.gwas.fst.ns.pdf"), height=7, width=7)

# BayeScan
p = ggplot(stats, aes(x=fst.all.WEIGHTED_FST, ..density.., 
        fill=overlap.perry.bs.batwa.bakiga)) + 
    geom_histogram(alpha=0.5, position="identity") + 
    xlab(expression(Average~F[ST])) + guides(fill=FALSE) + theme_bw()

ggsave(p, file=paste0(out.dir, "perry.bayescan.fst.all.pdf"), height=7, width=7)

p = ggplot(stats, aes(x=fst.ns.WEIGHTED_FST, ..density.., 
        fill=overlap.perry.bs.batwa.bakiga)) + 
    geom_histogram(alpha=0.5, position="identity") + 
    xlab(expression(Average~F[ST])) + guides(fill=FALSE) + theme_bw()

ggsave(p, file=paste0(out.dir, "perry.bayescan.fst.ns.pdf"), height=7, width=7)

# iHS
p = ggplot(stats, aes(x=tajd.ns.agr.TajimaD, ..density.., 
        fill=overlap.perry.ihs.bakiga)) + 
    geom_histogram(alpha=0.5, position="identity") + 
    xlab("Nonsynonymous Tajima's D - Bakiga") + guides(fill=FALSE) + theme_bw()

ggsave(p, file=paste0(out.dir, "perry.ihs-bakiga.tajd.ns.pdf"), height=7, width=7)

p = ggplot(stats, aes(x=tajd.ns.rhg.TajimaD, ..density.., 
        fill=overlap.perry.ihs.batwa)) + 
    geom_histogram(alpha=0.5, position="identity") + 
    xlab("Nonsynonymous Tajima's D - Batwa") + guides(fill=FALSE) + theme_bw()

ggsave(p, file=paste0(out.dir, "perry.ihs-batwa.tajd.ns.pdf"), height=7, width=7)

# --- Jarvis et al 2012 and Lachance et al 2012 data

p = ggplot(stats, aes(x=fst.all.WEIGHTED_FST, ..density.., 
        fill=overlap.jl)) + 
    geom_histogram(alpha=0.5, position="identity") + 
    xlab(expression(Average~F[ST])) + guides(fill=FALSE) + theme_bw()

ggsave(p, file=paste0(out.dir, "jarvis-lachance.gwas.fst.all.pdf"), height=7, width=7)

# ----------------------------------------------------------------------------------------
# --- Plot distribution for candidate gene lists
# ----------------------------------------------------------------------------------------

stats.genes.expand = read.table(file="results/basic_stats.genes.txt")

omim     = read.table("data/Wood_etal_2014/OMIM_height_genes.txt")
stopgain = read.table("results/stopgain.genes.txt")

stats.genes.expand$is.omim     = stats.genes.expand$genes %in% omim$V1
stats.genes.expand$is.stopgain = stats.genes.expand$genes %in% stopgain$V1

p = ggplot(stats.genes.expand, aes(x=fst.ns.WEIGHTED_FST, ..density.., fill=is.omim)) + 
    geom_histogram(alpha=0.5, position="identity", bins=100) + 
    xlab(expression(Average~F[ST])) + guides(fill=FALSE) + theme_bw()

ggsave(p, file=paste0(out.dir, "omim.skeletal.fst.ns.pdf"), height=7, width=7)

p = ggplot(stats.genes.expand, aes(x=fst.syn.WEIGHTED_FST, ..density.., fill=is.omim)) + 
    geom_histogram(alpha=0.5, position="identity", bins=100) + 
    xlab(expression(Average~F[ST])) + guides(fill=FALSE) + theme_bw()

ggsave(p, file=paste0(out.dir, "omim.skeletal.fst.syn.pdf"), height=7, width=7)

p = ggplot(stats.genes.expand, aes(x=fst.all.WEIGHTED_FST, ..density.., fill=is.omim)) + 
    geom_histogram(alpha=0.5, position="identity", bins=100) + 
    xlab(expression(Average~F[ST])) + guides(fill=FALSE) + theme_bw()

ggsave(p, file=paste0(out.dir, "omim.skeletal.fst.all.pdf"), height=7, width=7)
