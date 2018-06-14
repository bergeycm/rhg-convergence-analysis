#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Plot genotype frequencies by population for a given SNP
# ----------------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

chr      = args[1]
pos      = args[2]
vcf      = args[3]      # data/all.pops.merged.clean.recode.vcf
pop.list = args[4]      # data/all.pop.list.txt, column 1 is ID, column 2 is population

pops = read.table(pop.list)

vcf.inds = read.csv(pipe(paste("grep 'CHROM'", vcf, "| tr '\t' ','")), header=FALSE)
vcf.inds = data.frame(t(vcf.inds[1,-c(1:9)]))
names(vcf.inds) = "IID"

vcf.pops = merge(vcf.inds, pops, by.x="IID", by.y="V1", sort=FALSE, all.x=TRUE)
vcf.pops = vcf.pops[order(match(vcf.pops$IID, vcf.inds$IID)),]

awk.cmd = paste0("awk '{ if ($1 == \"chr", chr, "\" && $2 == ", pos, ") print $0 }' ", vcf)
this.snp.line = read.csv(pipe(awk.cmd), header=FALSE, sep="\t")

names(this.snp.line) = c("chr", "pos", "rs", "REF", "PASS", "VAL", "PASS", "", "COLS",
                         paste0("ind", 1:(ncol(this.snp.line) - 9)))

this.snp.line[,10:ncol(this.snp.line)] = gsub(":.*", "",
                                              this.snp.line[,10:ncol(this.snp.line)])

vcf.pops.geno = cbind(vcf.pops, t(this.snp.line[,10:ncol(this.snp.line)]))
names(vcf.pops.geno) = c("ind", "pop", "geno")

# Remove individuals not in the populations file (e.g. Birhor that slipped through)
vcf.pops.geno = vcf.pops.geno[!is.na(vcf.pops.geno$pop),]

out.file = paste0(gsub('.vcf', '', gsub('.*/', 'results/genotype_freqs-', vcf)),
    "-chr", chr, "_", pos, ".txt")

# Make table of genotyping results
sink(out.file)
    table(vcf.pops.geno[,2:3])
sink()

# Plot genotypes
plot.title = paste0("chr", chr, ":", pos, " - ", this.snp.line$rs)
p = ggplot(vcf.pops.geno, aes(geno, ..count.., fill=pop)) +
    geom_bar(position="dodge") +
    ggtitle(plot.title) +
    xlab("Genotype")

ggsave(p, filename = gsub(".txt$", ".pdf", out.file), height=7, width=7)
