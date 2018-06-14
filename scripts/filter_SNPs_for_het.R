#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Screen out possible paralogous sites based on heterozygosity
# ----------------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
hwe.in = args[1]    # e.g. "results/HWE_info.eAGR.hwe"

hwe = read.table(hwe.in, header=TRUE)

counts = data.frame(do.call(rbind,
    strsplit(as.character(hwe$OBS.HOM1.HET.HOM2.), '/', fixed=TRUE)))
names(counts) = c("obs.hom1", "obs.het", "obs.hom2")

counts = data.frame(sapply(counts, as.numeric))

counts$N = 2 * (counts$obs.hom1 + counts$obs.het + counts$obs.hom2)
counts$num.ind = counts$N / 2

counts$p = (2 * counts$obs.hom1 + counts$obs.het) / counts$N
counts$q = (2 * counts$obs.hom2 + counts$obs.het) / counts$N

# Expected heterozygosity using simple method
counts$exp.het.simp = 2 * counts$p * counts$q

# Expected heterozygotes using simple method
counts$exp.het.ind.simp = counts$exp.het.simp * counts$num.ind

# Allele count for the two alleles
counts$p.count = counts$p * counts$N
counts$q.count = counts$q * counts$N

# p allele component of sum
p.component = (counts$p.count * (counts$p.count - 1)) / (counts$N * (counts$N - 1))
# q allele component of sum
q.component = (counts$q.count * (counts$q.count - 1)) / (counts$N * (counts$N - 1))

# Expected heterozygosity using Hohenlohe method
counts$exp.het.hohen = 1 - (p.component + q.component)

# Expected heterozygotes using Hohenlohe method
counts$exp.het.ind.hohen = counts$exp.het.hohen * counts$num.ind

# Observed homozygosity
counts$obs.het.prop = counts$obs.het / counts$num.ind

# FIS - simple method
counts$fis.simp = 1 - (counts$obs.het.prop / counts$exp.het.simp)

# FIS - Hohenlohe et al method
counts$fis.hohen = 1 - (counts$obs.het.prop / counts$exp.het.hohen)

filter.fail.idx = which(counts$obs.het.prop > 0.5 |
    counts$fis.simp < -0.5 | counts$fis.hohen < -0.5)

to.remove = hwe[filter.fail.idx,1:2]
to.remove$END = to.remove$POS + 1

# Write BED file of SNPs that fail filter for this population
write.table(to.remove, file=paste0(hwe.in, ".filtered.bed"),
    quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
