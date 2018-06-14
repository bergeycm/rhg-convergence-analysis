#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Permute to get joint p-value for each gene
# ----------------------------------------------------------------------------------------

library(plyr)
library(parallel)

args = commandArgs(trailingOnly=TRUE)
stat.in.1 = args[1]                        # e.g. results/pbs.Batwa.GBR.anno.bed
stat.in.2 = args[2]                        # e.g. results/pbs.Anda.LWK.anno.bed
out.file  = args[3]                        # e.g. results/pbs.Batwa.Anda.jointp.txt

read.stat = function (stat.in) {

    stat = read.table(stat.in)

    # Remove SNPs that are not associated with genes or with infinite values
    stat = stat[stat$V5 != ".",]
    stat = stat[abs(stat$V4) != Inf,]

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

    stat.all
}

stat.all.1 = read.stat(stat.in.1)
stat.all.2 = read.stat(stat.in.2)

# Get average by gene
avg.by.gene.1 = aggregate(stat.all.1$val, list(gene=stat.all.1$gene), mean)
avg.by.gene.2 = aggregate(stat.all.2$val, list(gene=stat.all.2$gene), mean)

names(avg.by.gene.1)[2] = names(avg.by.gene.2)[2] = "actual.avg"

# Merge two stats data.frames together
avg.by.gene = merge(avg.by.gene.1, avg.by.gene.2, by = "gene", all = FALSE)

# Do randomization

perm.and.avg = function() {

    # Randomize gene column (randomize SNPs and their associated genes,
    # keeping total number of SNPs per gene constant)
    stat.rand.1 = transform(stat.all.1, gene = sample(gene))
    stat.rand.2 = transform(stat.all.2, gene = sample(gene))

    rand.avg.by.gene.1 = aggregate(stat.rand.1$val, list(gene=stat.rand.1$gene), mean)
    rand.avg.by.gene.2 = aggregate(stat.rand.2$val, list(gene=stat.rand.2$gene), mean)

    rand.avg.by.gene = merge(avg.by.gene,
        merge(rand.avg.by.gene.1, rand.avg.by.gene.2,
            by = "gene", all = FALSE), by="gene", all.x=TRUE)
}

num.trials = 10000

cl = makeCluster(detectCores() - 1)
clusterExport(cl, list("perm.and.avg", "stat.all.1", "stat.all.2",
                       "num.trials", "avg.by.gene"))

more.extreme.than.both.by.gene = do.call(cbind, parLapply(cl = cl, 1:num.trials,
    function(x) {
        if (x %% 100 == 0) { print(paste("On trial", x, "of", num.trials, " - ", date())) }
        this.perm = perm.and.avg()
        return((this.perm$x.x < this.perm$actual.avg.x) &
               (this.perm$x.y < this.perm$actual.avg.y))
    }
))

avg.by.gene$p = 1 - (rowSums(more.extreme.than.both.by.gene) / num.trials)

write.table(avg.by.gene, file= out.file,
    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
