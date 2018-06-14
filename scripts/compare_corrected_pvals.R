#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Compare corrected PBS selection index values to originals
# ----------------------------------------------------------------------------------------

library(ggplot2)

stat.prefixes = c("bayenv_BF_pass",
    "pbs.Batwa.GBR", "pbs.Bakiga.GBR",
    "pbs.Anda.LWK", "pbs.BR1.LWK")

all.stats = lapply(stat.prefixes, function(stat.prefix) {

    stat.in = paste0("results/", stat.prefix, ".pvals.txt")
    stat.sizeCor = gsub("results", "results_sizeCor", stat.in)
    stat.mafCor = gsub("results", "results_mafCor", stat.in)

    orig = read.table(stat.in)
    cor.size = read.table(stat.sizeCor)
    cor.maf = read.table(stat.mafCor)

    names(orig) = names(cor.size) =  names(cor.maf) = c("gene", "pval")

    cor.size$cor.type = "size"
    cor.maf$cor.type = "maf"

    cor.both = rbind(cor.size, cor.maf)

    this.stat = merge(orig, cor.both,
            by="gene", all=TRUE, suffixes=c(".orig", ".cor"))

    names(this.stat)[names(this.stat) == "pval"] = "pval.maf"

    this.stat
})

names(all.stats) = stat.prefixes

all.stats.named = do.call(rbind,
    lapply(grep("pbs", names(all.stats), value=TRUE), function (stat.prefix) {
        all.stats[[stat.prefix]]$prefix = stat.prefix
        all.stats[[stat.prefix]]
    })
)

all.stats.named$cor.type[all.stats.named$cor.type == "maf"]  = "MAF-based correction"
all.stats.named$cor.type[all.stats.named$cor.type == "size"] = "size-based correction"

all.stats.named$prefix = gsub("\\..*", "", gsub("pbs\\.", "", all.stats.named$prefix))
all.stats.named$prefix[all.stats.named$prefix == "Anda"] = "Andamanese"
all.stats.named$prefix[all.stats.named$prefix == "BR1"]  = "Brahmin"

p = ggplot(all.stats.named, aes(log(pval.orig, base=10), log(pval.cor, base=10),
        col=abs(pval.orig - pval.cor) / pval.orig)) +
    geom_point(alpha=0.5) +
    geom_abline(slope=1, intercept=0, lty=2) +
    facet_grid(cor.type ~ prefix) +
    theme_bw() +
    scale_color_gradient(low = "grey", high = "red", guide=FALSE) +
    xlab(expression(paste("Original"~log[10], "p"))) +
    ylab(expression(paste("Corrected"~log[10], "p"))) +
    coord_fixed()

ggsave(p, file="reports/corrected_pval_comparison.pdf")

all.r.sq = do.call(rbind, lapply(names(all.stats), function (stat.prefix) {

    cor.subset = all.stats[[stat.prefix]]

    cor.subset$cor.type = factor(cor.subset$cor.type)

    subset.size = cor.subset[cor.subset$cor.type == "size",]
    subset.maf  = cor.subset[cor.subset$cor.type == "maf",]

    r.sq.size = cor(subset.size$pval.orig, subset.size$pval.cor) ^ 2
    r.sq.maf  = cor(subset.maf$pval.orig,  subset.maf$pval.cor)  ^ 2

    c(stat.prefix =stat.prefix, r.sq.size=r.sq.size, r.sq.maf=r.sq.maf)

}))
