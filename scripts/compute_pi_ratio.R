#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

library(ggplot2)

# ----------------------------------------------------------------------------------------
# --- Compute pi ratio
# ----------------------------------------------------------------------------------------

pi.files = list.files(path="results", pattern="pi_.*.sites.pi")

pops = gsub("pi_(.*)\\.sites.pi", "\\1", pi.files)

pi.data = lapply(pops, function(x) {
    this.pi = read.table(paste0("results/pi_", x, ".sites.pi"), header=TRUE)
    this.pi$pop = x
    return(this.pi)
})
names(pi.data) = pops

# --- Make graph that shows relationship between transformed value and its input pi values

a = seq(0, 1, 0.01)
sim = data.frame(pi.sub = sample(a, 1000000, replace=TRUE),
    pi.all = sample(a, 1000000, replace=TRUE))
sim = unique(sim)
sim$pi.ratio = 4 * (sim$pi.all - sim$pi.sub) / (1 + sim$pi.all + sim$pi.sub) ^ 2
p = ggplot(sim, aes(pi.all, pi.sub, col=pi.ratio)) +
    geom_point() +
    scale_colour_gradientn(colors = rainbow(7),
        guide = guide_colorbar(title=expression("Scaled"~pi))) +
    xlab(expression(pi~"in all three populations combined")) +
    ylab(expression(pi~"in single population")) +
    theme_bw()

ggsave(p, file="reports/transformed_pi_values.pdf", height=7, width=7)

# --- Compute pi ratio for each population

lapply(pops, function (x) {

    if (grepl("_all_pops", x)) {
        return()
    }

    if (x %in% c("Batwa", "Bakiga")) {
        all.pop.name = "GBR_all_pops"
    } else {
        all.pop.name = "LWK_all_pops"
    }

    both.pops = merge(pi.data[[x]], pi.data[[all.pop.name]],
        by=c("CHROM", "POS"), suffixes=c(".sub",".all"))
    both.pops$pi.ratio = 4 * (both.pops$PI.all - both.pops$PI.sub) /
        (1 + both.pops$PI.all + both.pops$PI.sub) ^ 2

    both.pops$end = both.pops$POS + 1

    # Sort
    both.pops$chr.num = gsub("chr", "", both.pops$CHROM)
    both.pops = both.pops[order(both.pops$chr.num, both.pops$POS),]

    snp.list.bed = both.pops[,c(1,2,8,7)]

    write.table(snp.list.bed, file=paste0("results/pi_ratio.", x, ".bed"),
        quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
})
