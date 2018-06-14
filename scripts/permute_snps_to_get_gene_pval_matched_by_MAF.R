#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Permute to get p-value for each gene
# ----------------------------------------------------------------------------------------

library(plyr)
library(parallel)

args = commandArgs(trailingOnly=TRUE)
stat.in = args[1]                          # e.g. results/pbs.Batwa.GBR.anno.bed

print(paste("Processing file ", stat.in, " - ", date()))

stat = read.table(stat.in)

# Remove SNPs that are not associated with genes or with infinite values
stat = stat[stat$V5 != ".",]
stat = stat[abs(stat$V4) != Inf,]

print(paste("After filtration, there are ", nrow(stat), "rows"))

# Duplicate lines for SNPs associated with multiple genes

dups = stat[grepl(",", stat$V5),]

split.line = function (line) {

    genes = strsplit(line$V5, split=',', fixed=TRUE)[[1]]
    full.lines = cbind(line[,1:4], genes)
    return(full.lines)

}

dups.sing = ldply(1:nrow(dups), function (x) {
    if (x %% 10000 == 0) { print(paste("On row", x, "of", nrow(dups))) }
    split.line(dups[x,])
})

dups.sing.uniq = unique(dups.sing)

names(stat) = names(dups.sing.uniq) = c("chr", "start", "end", "val", "gene")

stat.all = rbind(stat[!grepl(",", stat$V5),], dups.sing.uniq)

# ----------------------------------------------------------------------------------------
# --- Read in global MAF info and add column to data.frame
# ----------------------------------------------------------------------------------------

# Use African MAFs for Bayenv for now

if (grepl("Batwa", stat.in) | grepl("Bakiga", stat.in) |
    grepl("twa",   stat.in) | grepl("kiga",   stat.in) |
    grepl("eAGR_eRHG", stat.in) | grepl("bayenv", stat.in)) {

    frq.in = "results/AGRHUM_EASTERN_100x267251.frq"
} else {
    frq.in = "results/GreatApe_Indiafinal_Austosome.sm.recode.frq"
}

frq = read.table(frq.in, sep="\t", skip=1)
names(frq) = c("chr", "start", "N_ALLELES", "N_CHR", "allele1.frq", "allele2.frq")

frq$allele1.frq = gsub(".*:", "", frq$allele1.frq)
frq$allele2.frq = gsub(".*:", "", frq$allele2.frq)

frq$maf = apply(frq[5:6], 1, min)

stat.all = merge(stat.all, frq[c("chr", "start", "maf")])

stat.all$maf = as.numeric(stat.all$maf)

# Get average by gene (after removing SNPs not in MAF input file)
avg.by.gene = aggregate(stat.all$val, list(gene=stat.all$gene), mean)

names(avg.by.gene)[2] = "actual.avg"

# ----------------------------------------------------------------------------------------
# --- Do randomization
# ----------------------------------------------------------------------------------------

perm.and.avg = function() {

    maf.bin.size = 0.01

    bin.starts = seq(from=0, to=0.5, by=maf.bin.size)
    bin.starts = bin.starts[-c(length(bin.starts))]

    # Keep ones with MAF = 0.5
    bin.starts = c(bin.starts, 0.5)

    stat.all.split = split(stat.all, cut(stat.all$maf, bin.starts, include.lowest=TRUE))

    # Jumble each subset of the data and reassemble
    stat.rand = do.call(rbind, lapply(stat.all.split, function (stat.subset) {
        stat.subset.rand = transform(stat.subset, gene = sample(gene))
    }))

    rand.avg.by.gene = aggregate(stat.rand$val, list(gene=stat.rand$gene), mean)

    rand.avg.by.gene[order(rand.avg.by.gene$gene),]
}

num.trials = 100000

#   perm.avgs.by.gene = do.call(cbind, mclapply(1:num.trials, function(x) {
#       if (x %% 100 == 0) { print(paste("On trial", x, "of", num.trials, " - ", date())) }
#       perm.and.avg()$x
#   }, mc.cores=19))

cl = makeCluster(detectCores() - 1)
clusterExport(cl, list("perm.and.avg", "stat.all", "num.trials"))
perm.avgs.by.gene = do.call(cbind, parLapply(cl = cl, 1:num.trials, function(x) {
   if (x %% 100 == 0) { print(paste("On trial", x, "of", num.trials, " - ", date())) }
   perm.and.avg()$x
}))

save(list=c("perm.avgs.by.gene"), file=paste0(stat.in, "_parallel_test_mafCor.Rdata"))

avg.by.gene$p = 1 - (rowSums(avg.by.gene$actual.avg >= perm.avgs.by.gene) / num.trials)

p.val.out = gsub("\\.\\.", ".", gsub("anno(.*).bed", "\\1.pvals.txt", stat.in))
p.val.out = gsub("results", "results_mafCor", p.val.out)

write.table(avg.by.gene[,c(1,3)], file=p.val.out,
    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
