#!/usr/bin/env Rscript

library(ggplot2)
library(ggrepel)

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Plot convergent PBS values by SNP
# ----------------------------------------------------------------------------------------

pbs.conv = read.table("results/joint_SNP_pvalues.txt")
stat.ids = c("pbs.batwa", "pbs.anda")
stats = pbs.conv[pbs.conv$stat1 %in% stat.ids & pbs.conv$stat2 %in% stat.ids,]

# ----------------------------------------------------------------------------------------

stats = stats[stats$genes.x != ".",]

stats$chr   = gsub(":.*", "", stats$snp.info)
stats$start = as.numeric(gsub(".*:", "", stats$snp.info))

chrs = paste("chr", 1:22, sep="")

just.pos = unique(stats[,11:12])

just.chrs = gsub("chr", "", just.pos$chr)
just.chrs.t = table(just.chrs)[order(as.numeric(names(table(just.chrs))))]

all.snp.pos = do.call(rbind, lapply(1:22, function (x) {
    chr = paste0("chr", x)
    this.chr.pos = sort(unique(stats[stats$chr == chr,]$start))

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
names(all.snp.pos) = c("chr", "start", "index")
all.snp.pos$chr = paste0("chr", all.snp.pos$chr)

stats.idx = merge(stats, all.snp.pos, by=c("chr", "start"))
stats.idx.srt = stats.idx[order(stats.idx$joint.p.edington, decreasing=FALSE),]

stats.idx.srt$top.gene = "NO"

stats.idx.srt[1:100,]$top.gene = "HIGH"
stats.idx.srt[1:10,]$top.gene = "TOP"

# Grab top and high SNPs plus a subset (size 10k) of the rest
stats.idx.subset = rbind(stats.idx.srt[stats.idx.srt$top.gene != "NO",],
                         stats.idx.srt[sample(which(stats.idx.srt$top.gene == "NO"), 10000),])

stats.idx.subset$genes = gsub(",.*", "", stats.idx.subset$genes)

top.snps = stats.idx.subset[stats.idx.subset$top.gene == "TOP",]

p = ggplot(stats.idx.subset, aes(index, -log(joint.p.edington, base=10),
        col=as.numeric(gsub("chr", "", chr)) %% 2)) +
    geom_hline(yintercept=0) +
    geom_point() +
    geom_label_repel(data=top.snps,
        aes(index, -log(joint.p.edington, base=10), label=genes, fill=stat, angle=90),
        col='white', fill='orange', size=2,
        segment.color = '#cccccc') +
    xlab("SNP") +
    ylab("Convergent PBS (-log10)") +
    theme_bw() +
    guides(fill=FALSE, color=FALSE) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))

ggsave(p, file="reports/convergent_PBS_p_by_SNP.pdf", height=7, width=12)
