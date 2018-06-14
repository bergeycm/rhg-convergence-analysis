#!/usr/bin/env Rscript

library(ggplot2)
library(ggrepel)

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Plot PBS values by SNP
# ----------------------------------------------------------------------------------------

# --- Batwa

pbs.batwa = read.table("results/pbs.Batwa.GBR.anno.bed")
names(pbs.batwa) = c("chr", "start", "end", "pbs", "genes")
pbs.batwa$genes = gsub(",.*", "", pbs.batwa$genes)
pbs.batwa$stat = "pbs.batwa"

# --- Bakiga

pbs.bakiga = read.table("results/pbs.Bakiga.GBR.anno.bed")
names(pbs.bakiga) = c("chr", "start", "end", "pbs", "genes")
pbs.bakiga$genes = gsub(",.*", "", pbs.bakiga$genes)
pbs.bakiga$stat = "pbs.bakiga"

# --- Anda

pbs.anda = read.table("results/pbs.Anda.LWK.anno.bed")
names(pbs.anda) = c("chr", "start", "end", "pbs", "genes")
pbs.anda$genes = gsub(",.*", "", pbs.anda$genes)
pbs.anda$stat = "pbs.anda"

# --- BR1

pbs.br1 = read.table("results/pbs.BR1.LWK.anno.bed")
names(pbs.br1) = c("chr", "start", "end", "pbs", "genes")
pbs.br1$genes = gsub(",.*", "", pbs.br1$genes)
pbs.br1$stat = "pbs.br1"

# ----------------------------------------------------------------------------------------

stats = rbind(rbind(rbind(pbs.batwa, pbs.bakiga), pbs.anda), pbs.br1)

# While here, write table of SNP counts by gene
gene.SNP.counts = table(stats$genes,stats$stat)
write.table(data.frame(gene.SNP.counts), file="results/gene_SNP_counts.txt")

save.image("gene.stats.Rdata")
load("gene.stats.Rdata")

# refgene = read.table("refGene/refGene.sort.simple.genes.gtf")

stats = stats[stats$genes != ".",]

stats = stats[abs(stats$pbs) != Inf,]

chrs = paste("chr", 1:22, sep="")

just.pos = unique(stats[,1:2])

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
stats.idx.srt = stats.idx[order(stats.idx$pbs, decreasing=TRUE),]

stats.idx.srt$top.gene = "NO"

stats.idx.srt[head(which(stats.idx.srt$stat == "pbs.anda"),   n=100),]$top.gene = "HIGH"
stats.idx.srt[head(which(stats.idx.srt$stat == "pbs.bakiga"), n=100),]$top.gene = "HIGH"
stats.idx.srt[head(which(stats.idx.srt$stat == "pbs.batwa"),  n=100),]$top.gene = "HIGH"
stats.idx.srt[head(which(stats.idx.srt$stat == "pbs.br1"),    n=100),]$top.gene = "HIGH"

stats.idx.srt[head(which(stats.idx.srt$stat == "pbs.anda"),   n=5),]$top.gene = "TOP"
stats.idx.srt[head(which(stats.idx.srt$stat == "pbs.bakiga"), n=5),]$top.gene = "TOP"
stats.idx.srt[head(which(stats.idx.srt$stat == "pbs.batwa"),  n=5),]$top.gene = "TOP"
stats.idx.srt[head(which(stats.idx.srt$stat == "pbs.br1"),    n=5),]$top.gene = "TOP"

# Grab top and high SNPs plus a subset (size 10k) of the rest
stats.idx.subset = rbind(stats.idx.srt[stats.idx.srt$top.gene != "NO",],
                         stats.idx.srt[sample(which(stats.idx.srt$top.gene == "NO"), 10000),])

stats.idx.subset$genes = gsub(",.*", "", stats.idx.subset$genes)

stats.idx.subset$stat = factor(stats.idx.subset$stat,
    levels=c("pbs.batwa", "pbs.bakiga", "pbs.anda", "pbs.br1"),
    labels=c("Batwa", "Bakiga", "Andamanese", "Brahmin"))

top.snps = stats.idx.subset[stats.idx.subset$top.gene == "TOP",]

p = ggplot(stats.idx.subset, aes(index, pbs,
        col=as.numeric(gsub("chr", "", chr)) %% 2)) +
    geom_hline(yintercept=0) +
    geom_point() +
    geom_label_repel(data=top.snps,
        aes(index, pbs, label=genes, fill=stat, angle=90),
        col='white', fill='orange', size=2,
        nudge_y = ifelse(top.snps$stat %in% c("Batwa", "Bakiga"), 30, 20),
        segment.color = '#cccccc') +
    facet_grid(stat ~ .) +
    xlab("SNP") +
    ylab("PBS") +
    theme_bw() +
    guides(fill=FALSE, color=FALSE) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))

ggsave(p, file="reports/PBS_by_SNP.pdf", height=7, width=12)
