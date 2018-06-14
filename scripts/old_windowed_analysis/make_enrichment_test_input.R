w#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

dir.create("results/gene_lists")

# ----------------------------------------------------------------------------------------
# --- Make input files for enrichment tests
# ----------------------------------------------------------------------------------------

stats = read.table("results/basic_stats.anno.txt", header=TRUE)

genes = read.table("results/windows_100kb.genes.bed")
names(genes) = c(names(stats)[1:3], "genes")

genes$win_start = genes$win_start + 1

stats.genes = merge(stats, genes)
# Remove windows with no genes or that fail exon_filter_pass filter
stats.genes = stats.genes[stats.genes$genes != '.' & stats.genes$exon_filter_pass,]

# Expand genes (and duplicate other columns) so that each gene is on its own line

stats.genes.expand = do.call(rbind, lapply(1:nrow(stats.genes), function (i) {
    x = stats.genes[i,]
    these.genes = strsplit(x$genes, ',')[[1]]
    # Duplicate stuff in other columns, once per gene
    other.cols = x[-c(ncol(x))][rep(1,length(these.genes)),]
    other.cols$genes = these.genes
    return(other.cols)
}))

# Write file of stats with genes too
write.table(stats.genes.expand, file="results/basic_stats.genes.txt", 
    sep="\t", row.names=FALSE, quote=FALSE)

# ----------------------------------------------------------------------------------------

stats.of.interest = c("fst.all.WEIGHTED_FST",
                      "fst.syn.WEIGHTED_FST",
                      "fst.ns.WEIGHTED_FST",
                      "pi.all.agr.PI",
                      "pi.all.rhg.PI",
                      "pi.syn.agr.PI",
                      "pi.syn.rhg.PI",
                      "pi.ns.agr.PI",
                      "pi.ns.rhg.PI",
                      "pi.all.ratio",
                      "pi.syn.ratio",
                      "pi.ns.ratio",
                      "tajd.all.agr.TajimaD",
                      "tajd.all.rhg.TajimaD",
                      "tajd.syn.agr.TajimaD",
                      "tajd.syn.rhg.TajimaD",
                      "tajd.ns.agr.TajimaD",
                      "tajd.ns.rhg.TajimaD")

gene.lists = lapply(stats.of.interest, function(x) {
    this.gene.list = cbind(stats.genes.expand$genes, stats.genes.expand[x])
    # Remove NAs
    this.gene.list = this.gene.list[!is.na(this.gene.list[,2]),]
    out.file = paste0("results/gene_lists/gene_list_for_enrichment_", x, ".txt")
    write.table(this.gene.list, file=out.file, 
        sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
    return(this.gene.list)
})

names(gene.lists) = stats.of.interest

# ----------------------------------------------------------------------------------------

# --- Get outlier windows

overrep.gene.lists = lapply(stats.of.interest, function(x) {
    this.gene.list = gene.lists[[x]]
    this.gene.list.sort = this.gene.list[order(this.gene.list[,2]),]
    quants = quantile(this.gene.list.sort[,2], c(0.01, 0.99))

    top = this.gene.list.sort[this.gene.list.sort[,2] > quants[2],]
    bot = this.gene.list.sort[this.gene.list.sort[,2] < quants[1],]

    out.file.top = paste0("results/gene_lists/gene_list_for_overrep.top0.01.", x, ".txt")
    write.table(top, file=out.file.top, 
        sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)

    out.file.bot = gsub("top", "bot", out.file.top)
    write.table(bot, file=out.file.bot, 
        sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
})

names(overrep.gene.lists) = stats.of.interest
