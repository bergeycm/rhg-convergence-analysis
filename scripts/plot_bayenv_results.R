#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Plot Bayenv results and annotate with gene info
# ----------------------------------------------------------------------------------------

library(ggplot2)

bf = read.table("results/bayenv_BF_full.txt")
names(bf) = c("chr", "start", "label", "BF")

miss = read.table("results/all.pops.merged.clean.overly.missing.uniq.txt")
names(miss) = c("chr", "start")

miss$too.missing = TRUE
bf.m = merge(bf, miss, sort=FALSE, all.x=TRUE)
bf.m[is.na(bf.m$too.missing),]$too.missing = FALSE

bf.pass = subset(bf.m, too.missing == FALSE)

p = ggplot(bf.pass, aes(start, log10(BF), label=label)) +
    geom_point() +
    geom_text(data=subset(bf, BF > 10),
            aes(start, log10(BF), label=label),
            check_overlap = TRUE, angle = 45, vjust = 0, nudge_y = 0.05, size=3) +
    xlab("SNP Index") + ylab("log10(Bayes Factor)") +
    geom_hline(yintercept=log10(10), lty=3) +
    geom_hline(yintercept=log10(30), lty=3, col='red') +
    theme_bw() +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(p, file="reports/bayenv_BFs.pdf", height=8, width=12)

high = bf.pass[bf.pass$BF > 30,]
high[order(high$BF, decreasing=TRUE),]

write.table(bf.pass, file="results/bayenv_BF_pass.txt",
    sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Sort
bf.pass$chr.num = as.numeric(gsub("chr", "", bf.pass$chr))
bf.pass = bf.pass[order(bf.pass$chr.num, bf.pass$start),]

bf.bed = data.frame(bf.pass[,1:2], bf.pass[,2]+1, bf.pass[,4])
write.table(bf.bed, file="results/bayenv_BF_pass.bed",
    sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# From https://www.r-bloggers.com/what-does-a-bayes-factor-feel-like/:

#       BF   Interpretation
# --------   -----------------------------
#     >100   Extreme evidence for H1
# 30 – 100   Very strong evidence for H1
# 10 –  30   Strong evidence for H1
#  3 –  10   Moderate evidence for H1
#  1 –   3   Anecdotal evidence for H1
#        1   No evidence
