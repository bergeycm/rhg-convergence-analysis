#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

library(ggplot2)

# ========================================================================================
# --- Plot stats in various ways
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Bring in results of parse_windowed_stats
# ----------------------------------------------------------------------------------------

stats = read.table(file="results/basic_stats.anno.txt", sep="\t", header=TRUE)

stats.ok = stats[stats$exon_filter_pass,]

# ----------------------------------------------------------------------------------------
# --- Plot Fst across genome, showing top 5%
# ----------------------------------------------------------------------------------------

stats.ok$chr = factor(gsub("chr", "", stats.ok$chr), levels=1:22)

fst.99 = quantile(stats.ok$fst.all.WEIGHTED_FST, 0.99, na.rm=TRUE)

stats.ok$fst.top1 = stats.ok$fst.all.WEIGHTED_FST > fst.99

p = ggplot(stats.ok, aes(win_start, fst.all.WEIGHTED_FST, col=fst.top1)) + 
    geom_point(alpha=0.2) + facet_grid(chr ~ .) +
    ylab(expression(Average~F[ST])) + xlab("Position along chromosome") +
    theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# ----------------------------------------------------------------------------------------
# --- Plot pi ratio across genome
# ----------------------------------------------------------------------------------------

# Compare pi between populations

p = ggplot(stats, aes(pi.all.agr.PI, pi.all.rhg.PI)) + 
    geom_point(alpha=0.3) + geom_abline(slope=1, col='red', lty=3)

# ----------------------------------------------------------------------------------------
# --- Plot Fst by number of variants in window
# ----------------------------------------------------------------------------------------

p = ggplot(stats, aes(fst.all.N_VARIANTS, fst.all.WEIGHTED_FST)) + 
    geom_point(alpha=0.2) + 
    xlab("Number of variants in window") + 
    ylab(expression(Average~F[ST])) + 
    geom_hline(yintercept=mean(stats$fst.all.WEIGHTED_FST, na.rm=TRUE), col='red', lty=2)

# ----------------------------------------------------------------------------------------
# --- Plot SNP class counts
# ----------------------------------------------------------------------------------------

snp.counts = read.table("reports/SNP_class_counts.txt")
names(snp.counts) = c("num", "class")
snp.counts = snp.counts[order(snp.counts$num, decreasing=TRUE),]
snp.counts = snp.counts[grepl("3b", snp.counts$class) == FALSE,]
snp.counts$class = factor(snp.counts$class, levels=snp.counts$class)

p = ggplot(snp.counts, aes(x=class, y=num, fill=class)) + 
    geom_bar(stat="identity",position="dodge") + coord_flip() + 
    ylab("SNP class") + xlab("Count")

# ----------------------------------------------------------------------------------------
# --- Plot PolyPhen prediction counts
# ----------------------------------------------------------------------------------------

polyphen.counts = read.table("reports/Polyphen_prediction_counts.txt")
polyphen.counts[polyphen.counts$V2 == "B",]$V2 = "benign"
polyphen.counts[polyphen.counts$V2 == "P",]$V2 = "possibly damaging"
polyphen.counts[polyphen.counts$V2 == "D",]$V2 = "probably damaging"
names(polyphen.counts) = c("num", "class")

p = ggplot(polyphen.counts, aes(x=class, y=num, fill=class)) + 
    geom_bar(stat="identity",position="dodge") + coord_flip() + 
    ylab("PolyPhen prediction") + xlab("Count")
