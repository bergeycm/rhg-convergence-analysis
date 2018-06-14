#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ========================================================================================
# --- Remove SNPs that did not survive interpolation of genetic map
# ========================================================================================

args = commandArgs(trailingOnly=TRUE)

orig.snp.list      = args[1]
predicted.snp.list = args[2]
in.haps.file       = args[3]

orig    = read.table(orig.snp.list)
predict = read.table(predicted.snp.list)

surviving.lines = which(orig$V1 %in% predict$V4)

in.haps = read.table(in.haps.file)

haps.reduced = in.haps[surviving.lines,]

write.table(haps.reduced, file=paste0(in.haps.file, ".reduced"),
    sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
