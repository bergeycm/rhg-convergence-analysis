#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Compute 95th percentile value of S site distribution of PBS stats
# ----------------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

in.bed     = args[1]  # results/pbs.Batwa.GBR.anno.syn.bed
if (length(args) > 1) {
    percentile = as.numeric(args[2])
} else {
    percentile = 0.95
}

stat.exon.s = read.table(in.bed)

# Get genome-wide Nth percentile of S SNPs (default = 0.95)
overall.s.percentile.cutoff = quantile(stat.exon.s[,4], percentile)

out.file = gsub(".bed$", paste0(".perc", percentile * 100, ".txt"), in.bed)

write.table(overall.s.percentile.cutoff, file=out.file,
    quote=FALSE, row.names=FALSE, col.names=FALSE)
