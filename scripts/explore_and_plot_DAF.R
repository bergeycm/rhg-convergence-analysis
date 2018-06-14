#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

library(ggplot2)

# ========================================================================================
# --- Compute and plot derived allele frequency
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Bring in SNP's ancestral allele info
# ----------------------------------------------------------------------------------------

ancest  = read.table("data/SNPAncestralAllele.bcp")
alleles = read.table("data/Allele.bcp", sep="\t", fill=TRUE)

vcf = read.table("data/AGRHUM_EASTERN_100x267251.vcf")
vcf$rsid = gsub("rs", "", vcf$V3)

tmp  = merge(vcf, ancest,  by.x="rsid", by.y="V1")
tmp2 = merge(tmp, alleles, by.x="V2.y", by.y="V1")

ref.info = cbind(tmp2[,3:7], tmp2$V2)
names(ref.info) = c("chr", "pos", "rs.id", "ref", "alt", "aa")

write.table(ref.info, file="results/ref.info.txt", quote=FALSE, sep="\t")

# ----------------------------------------------------------------------------------------
# --- Bring in population-specific frequency info
# ----------------------------------------------------------------------------------------

agr.frq = read.table("results/AGRHUM_EASTERN_100x267251.AGR.frq",
    sep="\t", row.names=NULL)
rhg.frq = read.table("results/AGRHUM_EASTERN_100x267251.RHG.frq",
    sep="\t", row.names=NULL)

names(agr.frq) = names(rhg.frq) = c("chr", "pos", "n_alleles", "n_chr", "freq1", "freq2")

names(agr.frq)[5:6] = paste0("AGR", names(agr.frq)[5:6])
names(rhg.frq)[5:6] = paste0("RHG", names(rhg.frq)[5:6])

frq = merge(agr.frq, rhg.frq)

frq.info = merge(frq, ref.info, by=c("chr", "pos"))

# Why are there duplicate rows?
frq.info.u = unique(frq.info)

frq.info.u$AGRfreq1 = gsub(".*:", "", frq.info.u$AGRfreq1)
frq.info.u$AGRfreq2 = gsub(".*:", "", frq.info.u$AGRfreq2)
frq.info.u$RHGfreq1 = gsub(".*:", "", frq.info.u$RHGfreq1)
frq.info.u$RHGfreq2 = gsub(".*:", "", frq.info.u$RHGfreq2)

# ----------------------------------------------------------------------------------------
# --- Bring in PBS info and compute DAF
# ----------------------------------------------------------------------------------------

pbs.batwa.N = read.table("results/pbs.Batwa.GBR.anno.nonsyn.bed")
pbs.batwa.S = read.table("results/pbs.Batwa.GBR.anno.syn.bed")

pbs.batwa.N$type = "N"
pbs.batwa.S$type = "S"

pbs.batwa = rbind(pbs.batwa.N, pbs.batwa.S)

names(pbs.batwa) = c("chr", "pos", "end", "pbs", "genes", "type")

frq.info.u.NS = merge(frq.info.u, pbs.batwa, all.x=TRUE)

ref.eq.aa = which(frq.info.u.NS$ref == frq.info.u.NS$aa)
alt.eq.aa = which(frq.info.u.NS$alt == frq.info.u.NS$aa)

frq.info.u.NS$daf.rhg = NA
frq.info.u.NS[ref.eq.aa,]$daf.rhg = frq.info.u.NS[ref.eq.aa,]$RHGfreq2
frq.info.u.NS[alt.eq.aa,]$daf.rhg = frq.info.u.NS[alt.eq.aa,]$RHGfreq1

frq.info.u.NS$daf.agr = NA
frq.info.u.NS[ref.eq.aa,]$daf.agr = frq.info.u.NS[ref.eq.aa,]$AGRfreq2
frq.info.u.NS[alt.eq.aa,]$daf.agr = frq.info.u.NS[alt.eq.aa,]$AGRfreq1

# ----------------------------------------------------------------------------------------
# --- Plot DAF distribution and PBS by DAF
# ----------------------------------------------------------------------------------------

# RHG DAF density
p = ggplot(frq.info.u.NS[which(!is.na(frq.info.u.NS$type)),], aes(daf.rhg, col=type)) +
    geom_density() + xlab("Batwa DAF")

ggsave("reports/DAF.batwa.density.pdf", height=5, width=7)

# AGR DAF density
p = ggplot(frq.info.u.NS[which(!is.na(frq.info.u.NS$type)),], aes(daf.agr, col=type)) +
    geom_density() + xlab("Bakiga DAF")

ggsave("reports/DAF.bakiga.density.pdf", height=5, width=7)

# ----------------------------------------------------------------------------------------

# RHG PBS by DAF
p = ggplot(frq.info.u.NS[which(!is.na(frq.info.u.NS$type)),],
        aes(daf.rhg, pbs, col=type)) +
    geom_point(alpha=0.25) + facet_grid(type ~ .) +
    xlab("Batwa DAF") + ylab("Batwa PBS")

ggsave("reports/DAF.PBS.batwa.pdf", height=7, width=7)

# AGR PBS by DAF
p = ggplot(frq.info.u.NS[which(!is.na(frq.info.u.NS$type)),],
        aes(daf.agr, pbs, col=type)) +
    geom_point(alpha=0.25) + facet_grid(type ~ .) +
    xlab("Bakiga DAF") + ylab("Bakiga PBS")

ggsave("reports/DAF.PBS.bakiga.pdf", height=7, width=7)

# ----------------------------------------------------------------------------------------

# RHG DAF by AGR DAF
p = ggplot(frq.info.u.NS[which(!is.na(frq.info.u.NS$type)),],
        aes(daf.agr, daf.rhg, col=type)) +
    geom_point(alpha=0.25) + facet_grid(type ~ .) +
    xlab("Bakiga DAF") + ylab("Batwa DAF")

ggsave("reports/DAF.batwa.bakiga.pdf", height=7, width=7)
