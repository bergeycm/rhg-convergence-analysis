#!/usr/bin/env Rscript

library(ggplot2)
library(ggrepel)

options(stringsAsFactors=FALSE)

args = commandArgs(trailingOnly=TRUE)
cor.type = args[1]

if (is.na(cor.type)) {
    cor.type = ""
}

# ----------------------------------------------------------------------------------------
# --- Plot PBS index values by gene
# ----------------------------------------------------------------------------------------

# --- Batwa

pbs.batwa = read.table(paste0("results", cor.type, "/pbs.Batwa.GBR.pvals.txt"))
names(pbs.batwa) = c("gene", "pbs")
pbs.batwa$stat = "pbs.batwa"

# --- Bakiga

pbs.bakiga = read.table(paste0("results", cor.type, "/pbs.Bakiga.GBR.pvals.txt"))
names(pbs.bakiga) = c("gene", "pbs")
pbs.bakiga$stat = "pbs.bakiga"

# --- Anda

pbs.anda = read.table(paste0("results", cor.type, "/pbs.Anda.LWK.pvals.txt"))
names(pbs.anda) = c("gene", "pbs")
pbs.anda$stat = "pbs.anda"

# --- BR1

pbs.br1 = read.table(paste0("results", cor.type, "/pbs.BR1.LWK.pvals.txt"))
names(pbs.br1) = c("gene", "pbs")
pbs.br1$stat = "pbs.br1"

# ----------------------------------------------------------------------------------------

stats.simple = rbind(rbind(rbind(pbs.batwa, pbs.bakiga), pbs.anda), pbs.br1)

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

stats.simple[which(stats.simple$pbs == 0),]$pbs = 1e-5
stats.simple$neg.log.p = -log(stats.simple$pbs, base=10)

stats = merge(stats.simple, gene.info)

# ----------------------------------------------------------------------------------------

stats = stats[stats$gene != ".",]
stats = stats[which(! is.na(stats$pbs)),]
stats = stats[which(! is.na(stats$center)),]

stats = stats[abs(stats$pbs) != Inf,]

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

stats.idx.subset$stat = factor(stats.idx.subset$stat,
    levels=c("pbs.batwa", "pbs.bakiga", "pbs.anda", "pbs.br1"),
    labels=c("Batwa", "Bakiga", "Andamanese", "Brahmin"))

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
    facet_grid(stat ~ .) +
    xlab("SNP") +
    ylab("PBS index [-log(PBS index)]") +
    theme_bw() +
    guides(fill=FALSE, color=FALSE) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black")) +
    ylim(c(0,8))

ggsave(p, file=paste0("reports", cor.type, "/PBS_by_gene.pdf"), height=7, width=12)

# ========================================================================================
# --- Plot selection stats by gene length or proxies to search for gene length bias
# ========================================================================================

stats.idx.u = unique(stats.idx)

stats.idx.u$stat = factor(stats.idx.u$stat,
    levels=c("pbs.batwa", "pbs.bakiga", "pbs.anda", "pbs.br1"),
    labels=c("Batwa", "Bakiga", "Andamanese", "Brahmin"))

# --- Plot PBS index by gene length

p = ggplot(stats.idx.u, aes(max - min, neg.log.p)) +
    geom_point(alpha=0.1) +
    facet_grid(stat ~ .) +
    theme_bw() +
    ylab("-log(PBS index)") +
    xlab("Gene length") +
    geom_smooth(method='lm')

ggsave(p, file=paste0("reports", cor.type, "/PBS_by_gene_length.pdf"), height=7, width=12)

# Test for relationship via linear regression
lentest.len = stats.idx.u$max - stats.idx.u$min
lentest.neg.log.p = stats.idx.u$neg.log.p

lentest.fit = lm(lentest.neg.log.p ~ lentest.len)
summary(lentest.fit)

# --- Also plot by number of SNPs (-log index)

gene.SNP.counts = read.table("results/gene_SNP_counts.txt")
names(gene.SNP.counts) = c("gene", "stat", "SNP.count")

gene.SNP.counts = gene.SNP.counts[gene.SNP.counts$gene != ".",]

gene.SNP.counts$stat = factor(gene.SNP.counts$stat,
    levels=c("pbs.batwa", "pbs.bakiga", "pbs.anda", "pbs.br1"),
    labels=c("Batwa", "Bakiga", "Andamanese", "Brahmin"))

stats.idx.u.ct = merge(stats.idx.u, gene.SNP.counts, by=c("gene", "stat"))

p = ggplot(stats.idx.u.ct, aes(SNP.count, neg.log.p)) +
    geom_point(alpha=0.1) +
    facet_grid(stat ~ .) +
    theme_bw() +
    ylab("-log(PBS index)") +
    xlab("SNP count") +
    geom_smooth(method='lm')

ggsave(p, file=paste0("reports", cor.type, "/PBS_by_SNP_count_log.pdf"),
    height=7, width=12)


p = ggplot(stats.idx.u.ct, aes(SNP.count, pbs)) +
    geom_point(alpha=0.1) +
    facet_grid(stat ~ .) +
    theme_bw() +
    ylab("PBS index") +
    xlab("SNP count") +
    geom_smooth(method='lm')

ggsave(p, file=paste0("reports", cor.type, "/PBS_by_SNP_count.pdf"),
    height=7, width=12)

# ========================================================================================
# --- Additional plots
# ========================================================================================

# --- Plot PBS index variance by gene length

stats.idx.u$tot.len = stats.idx.u$max - stats.idx.u$min

bin.width = 10000
bin.starts = seq(from = 0, to=max(stats.idx.u$tot.len), by=bin.width)
bin.vars = do.call(rbind, lapply(bin.starts, function (x) {
    this.pbs = stats.idx.u[stats.idx.u$tot.len > x & stats.idx.u$tot.len < x + bin.width,]
    var(this.pbs$neg.log.p)
}))

bin.vars.df = data.frame(start=bin.starts, var=bin.vars)

p = ggplot(bin.vars.df, aes(start, var)) +
    geom_point() +
    ylab("Variance in -log(PBS selection index)") +
    xlab("Total gene length (bp) - binned")

ggsave(p, file=paste0("reports", cor.type, "/PBS_var_by_gene_length.pdf"),
    height=7, width=12)

# --- Also plot variance by length with equal sized bins
# --- (plotting variance against average size of bin)

srt.len = stats.idx.u[order(stats.idx.u$tot.len),]

bin.length = 50
start.idx = seq(from=1, to=nrow(srt.len), by=bin.length)

eqbin.vars = do.call(rbind, lapply(start.idx, function (x) {
    this.pbs = srt.len[x:(x + bin.length),]
    return(list(bin.num=x,
                avg.len=mean(this.pbs$tot.len),
                pbs.var=var(this.pbs$neg.log.p)))
}))

eqbin.vars = data.frame(eqbin.vars)

eqbin.vars$avg.len = as.numeric(eqbin.vars$avg.len)
eqbin.vars$pbs.var = as.numeric(eqbin.vars$pbs.var)

p = ggplot(eqbin.vars, aes(avg.len, pbs.var)) +
    geom_point() +
    ylab("Variance in -log(PBS selection index)") +
    xlab(paste0("Mean gene length (bp) for bin of ", bin.length, " SNPs"))

ggsave(p, file=paste0("reports", cor.type, "/PBS_var_by_gene_length_SNP_bins.pdf"),
    height=7, width=12)

# --- Plot standard deviation of the p-values by SNP count

std.devs = sapply(1:20, function (x) {

    these.SNPs = stats.idx.u.ct[stats.idx.u.ct$SNP.count == x,]
    #pop.ct = table(these.SNPs$stat)
    #avg.pop.ct = mean(pop.ct)

    #std.err = sd(these.SNPs$pbs, na.rm=TRUE) /
    #    sqrt(length(these.SNPs$pbs[!is.na(these.SNPs$pbs)]))

    if (length(these.SNPs$pbs) >= 20) {
        std.dev = sapply(1:100, function (x) {
            sd(sample(these.SNPs$pbs, size=20, replace=TRUE))
        })
    } else {
        std.dev = NA
    }
    return(quantile(std.dev, c(0.05, 0.5, 0.95), na.rm=TRUE))
})

sd.df = data.frame(cbind(1:20, t(std.devs)))
names(sd.df) = c("SNP.count", "quant5", "median", "quant95")

p = ggplot(sd.df, aes(SNP.count, median)) +
    geom_segment(aes(y=quant5, yend=quant95, x=SNP.count, xend=SNP.count),
        col='darkgrey') +
    geom_point() +
    theme_bw() +
    xlab("SNPs per gene") +
    ylab("Standard deviation of PBS index\n(Equal subsamples)")

ggsave(p, file=paste0("reports", cor.type, "/PBS_stddev_by_SNP_count.pdf"),
    height=7, width=7)

# --- Plot mean of the p-values by SNP count

means = sapply(1:20, function (x) {

    these.SNPs = stats.idx.u.ct[stats.idx.u.ct$SNP.count == x,]
    #pop.ct = table(these.SNPs$stat)
    #avg.pop.ct = mean(pop.ct)

    #std.err = sd(these.SNPs$pbs, na.rm=TRUE) /
    #    sqrt(length(these.SNPs$pbs[!is.na(these.SNPs$pbs)]))

    if (length(these.SNPs$pbs) >= 20) {
        this.mean = sapply(1:100, function (x) {
            mean(sample(these.SNPs$pbs, size=20, replace=TRUE))
        })
    } else {
        this.mean = NA
    }
    return(quantile(this.mean, c(0.05, 0.5, 0.95), na.rm=TRUE))
})

means.df = data.frame(cbind(1:20, t(means)))
names(means.df) = c("SNP.count", "quant5", "median", "quant95")

p = ggplot(means.df, aes(SNP.count, median)) +
    geom_segment(aes(y=quant5, yend=quant95, x=SNP.count, xend=SNP.count),
        col='darkgrey') +
    geom_point() +
    theme_bw() +
    xlab("SNPs per gene") +
    ylab("Mean of PBS index\n(Equal subsamples)")

ggsave(p, file=paste0("reports", cor.type, "/PBS_mean_by_SNP_count.pdf"),
    height=7, width=7)

# --- Plot standard error of the p-values by SNP count

std.errs = sapply(1:20, function (x) {

    these.SNPs = stats.idx.u.ct[stats.idx.u.ct$SNP.count == x,]
    pop.ct = table(these.SNPs$stat)
    avg.pop.ct = mean(pop.ct)

    std.err = sd(these.SNPs$pbs, na.rm=TRUE) #/
       # sqrt(length(these.SNPs$pbs[!is.na(these.SNPs$pbs)]))
})

p = ggplot(data.frame(std.errs = std.errs), aes(1:20, std.errs)) +
    geom_point() +
    theme_bw() +
    ylab("Standard Error of PBS index") +
    xlab("SNPs per gene")

ggsave(p, file=paste0("reports", cor.type, "/PBS_stderr_by_SNP_count.pdf"),
    height=7, width=7)

# --- Do boxplot of selection index

stats.idx.u.ct.sm = stats.idx.u.ct[stats.idx.u.ct$SNP.count >= 1 &
                                   stats.idx.u.ct$SNP.count < 30,]

ct.tbl = table(stats.idx.u.ct.sm$stat, stats.idx.u.ct.sm$SNP.count)

library(reshape2)
ct.tbl.long = melt(ct.tbl, id.vars=c("subject", "sex"))
names(ct.tbl.long) = c("stat", "SNP.count", "total.genes")

stats.idx.u.ct.sm = merge(stats.idx.u.ct.sm, ct.tbl.long)

p = ggplot(stats.idx.u.ct.sm, aes(SNP.count, pbs, group=SNP.count, fill=total.genes)) +
    geom_boxplot(outlier.color="black", outlier.shape=1,
        outlier.size=0.5, outlier.alpha=0.25) +
    scale_fill_gradient("Gene count", low = "white", high = "black") +
    facet_grid(stat ~ .) +
    theme_bw() +
    ylab("PBS index") +
    xlab("SNPs per gene") +
    geom_smooth(method='lm', se=TRUE, aes(group=1), lty=2, col='red')

ggsave(p, file=paste0("reports", cor.type, "/PBS_boxplots_by_SNP_count.pdf"),
    height=7, width=12)

# --- Plot significant gene count by SNP length

p = ggplot(stats.idx.u.ct, aes(SNP.count, fill=pbs < 0.05)) +
    geom_histogram(binwidth=1) +
    xlab("Gene SNP count") +
    ylab("Total Genes") +
    theme_bw() +
    xlim(c(0,25))

ggsave(p, file=paste0("reports", cor.type, "/PBS_sig_gene_SNP_count_histogram.pdf"))

# --- Plot proportion of genes significant by SNP count and population

sig.prop = stats.idx.u.ct[,c("stat", "pbs", "SNP.count")]
sig.prop$sig = sig.prop$pbs < 0.05

sig.prop.pop = aggregate(sig.prop$sig, by=list(sig.prop$stat, sig.prop$SNP.count),
    FUN=mean)
names(sig.prop.pop) = c("pop", "SNP.ct", "prop.sig")

sig.prop.pop = sig.prop.pop[sig.prop.pop$SNP.ct >= 1 & sig.prop.pop$SNP.ct < 21,]

p = ggplot(sig.prop.pop, aes(SNP.ct, prop.sig)) +
    geom_point() +
    facet_grid(pop ~ .) +
    xlab("Gene SNP count") +
    ylab("Proportion significant (p<0.05)") +
    theme_bw() +
    #xlim(c(1,25)) +
    geom_hline(yintercept=0.05, col='red', lty=2)

ggsave(p, file=paste0("reports", cor.type, "/PBS_prop_sig_gene_SNP_count_by_pop.pdf"))

# --- Plot proportion of genes significant by SNP count

sig.tbl = table(stats.idx.u.ct$pbs < 0.05, stats.idx.u.ct$SNP.count)
prop.sig = sig.tbl["TRUE",] / colSums(sig.tbl)
prop.sig.df = data.frame(SNP.ct = as.numeric(names(prop.sig)),
                         this.prop.sig=as.numeric(prop.sig))
p = ggplot(prop.sig.df[2:26,], aes(SNP.ct, this.prop.sig)) +
    geom_point() +
    xlab("Gene SNP count") +
    ylab("Proportion significant (p<0.05)") +
    theme_bw() +
    xlim(c(1,25)) +
    geom_hline(yintercept=0.05, col='red', lty=2)

ggsave(p, file=paste0("reports", cor.type, "/PBS_prop_sig_gene_SNP_count.pdf"))
