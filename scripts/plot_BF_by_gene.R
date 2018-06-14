#!/usr/bin/env Rscript

library(ggplot2)
library(ggrepel)

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Plot BF index values by gene
# ----------------------------------------------------------------------------------------

stats.simple = read.table("results/bayenv_BF_pass.pvals.txt")
names(stats.simple) = c("gene", "BF")

# ----------------------------------------------------------------------------------------

# --- Get gene positions and add to stats
#refgene = read.table("refGene/refGene.sort.simple.gene.gtf")
refgene = read.table("refGene/refGene.sort.simple.justGenes.gtf")
# Get exon minimum and maximum
refgene$exon.min = with(refgene, pmin(V4, V5))
refgene$exon.max = with(refgene, pmax(V4, V5))
# Get gene minimums and maximums and center
gene.mins = aggregate(refgene$exon.min, by=list(refgene$V9), FUN=min)
gene.maxs = aggregate(refgene$exon.max, by=list(refgene$V9), FUN=max)
gene.info = merge(gene.mins, gene.maxs, by="Group.1")
names(gene.info) = c("gene", "min", "max")
gene.info$center = (gene.info$min + gene.info$max) / 2
# Add chromosome
gene.info = merge(gene.info, refgene[,c(1,9)], by.x="gene", by.y="V9")[,c(1,5,2:4)]
names(gene.info)[2] = "chr"

stats.simple[which(stats.simple$BF == 0),]$BF = 1e-5
stats.simple$neg.log.p = -log(stats.simple$BF, base=10)

stats = merge(stats.simple, gene.info)

# ----------------------------------------------------------------------------------------

stats = stats[stats$gene != ".",]
stats = stats[which(! is.na(stats$BF)),]
stats = stats[which(! is.na(stats$center)),]

stats = stats[abs(stats$BF) != Inf,]

chrs = paste("chr", 1:22, sep="")

just.pos = unique(stats[,c("chr", "center")])

just.chrs = gsub("chr", "", just.pos$chr)
just.chrs.t = table(just.chrs)[order(as.numeric(names(table(just.chrs))))]

all.snp.pos = do.call(rbind, lapply(1:22, function (x) {
    chr = paste0("chr", x)
    this.chr.pos = sort(unique(stats[stats$chr == chr,]$center))

    # Shift chromosomes to cumulative index
    if (x == 1) {
        cum.len = 0
    } else {
        cum.len = sum(just.chrs.t[1:(x-1)])
    }
    snp.pos = cbind(x, this.chr.pos, cum.len + (1:length(this.chr.pos)))
    snp.pos
}))

all.snp.pos = data.frame(all.snp.pos)
names(all.snp.pos) = c("chr", "center", "index")
all.snp.pos$chr = paste0("chr", all.snp.pos$chr)

stats.idx = merge(stats, all.snp.pos, by=c("chr", "center"))
stats.idx.srt = stats.idx[order(stats.idx$neg.log.p, decreasing=TRUE),]

stats.idx.srt$top.gene = "NO"

stats.idx.srt[which(stats.idx.srt$neg.log.p > 3),]$top.gene = "HIGH"
stats.idx.srt[which(stats.idx.srt$neg.log.p == 5),]$top.gene = "TOP"

# Grab top and high SNPs plus a subset (size 10k) of the rest
stats.idx.subset = rbind(stats.idx.srt[stats.idx.srt$top.gene != "NO",],
                         stats.idx.srt[sample(which(stats.idx.srt$top.gene == "NO"), 10000),])

stats.idx.subset = unique(stats.idx.subset)

top.genes = stats.idx.subset[stats.idx.subset$top.gene == "TOP",]

p = ggplot(stats.idx.subset, aes(index, neg.log.p,
        col=as.numeric(gsub("chr", "", chr)) %% 2)) +
    geom_hline(yintercept=0) +
    geom_point() +
    geom_label_repel(data=top.genes,
        aes(index, neg.log.p, label=gene, fill=stat),
        col='white', fill='orange', size=2,
        segment.color = '#cccccc',
        nudge_y=1, direction='y', max.iter=10, force=10) +
    xlab("SNP") +
    ylab("Bayenv BF index [-log(BF index)]") +
    theme_bw() +
    guides(fill=FALSE, color=FALSE) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black")) +
    ylim(c(0,8))

ggsave(p, file="reports/BF_by_gene.pdf", height=7, width=12)
