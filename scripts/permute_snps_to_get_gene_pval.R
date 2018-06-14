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

# Get average by gene
avg.by.gene = aggregate(stat.all$val, list(gene=stat.all$gene), mean)

names(avg.by.gene)[2] = "actual.avg"

# Do randomization

perm.and.avg = function() {

    # Randomize gene column (randomize SNPs and their associated genes,
    # keeping total number of SNPs per gene constant)
    stat.rand = transform(stat.all, gene = sample(gene))

    rand.avg.by.gene = aggregate(stat.rand$val, list(gene=stat.rand$gene), mean)
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

save(list=c("perm.avgs.by.gene"), file=paste0(stat.in, "_parallel_test.Rdata"))

avg.by.gene$p = 1 - (rowSums(avg.by.gene$actual.avg >= perm.avgs.by.gene) / num.trials)

p.val.out = gsub("\\.\\.", ".", gsub("anno(.*).bed", "\\1.pvals.txt", stat.in))

write.table(avg.by.gene[,c(1,3)], file=p.val.out,
    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
