#!/usr/bin/Rscript

library(ggplot2)
library(RColorBrewer)

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly = TRUE)
adm.prefix = args[1]
fam.file   = args[2]

# ========================================================================================
# === Analyze and plot ADMIXTURE results
# ========================================================================================

adm.inds = read.table(fam.file)$V1

# Ugandans
agr.ids = read.table("data/eAGR_IDs.txt")
rhg.ids = read.table("data/eRHG_IDs.txt")

gbr.ids = read.table("data/1000genomes/subset_GBR.txt")
lwk.ids = read.table("data/1000genomes/subset_LWK.txt")

ind.info = data.frame(ind=adm.inds, pop="")

if (grepl("AGRHUM", fam.file)) {
    ind.info[ind.info$ind %in% agr.ids$V1,]$pop = "Bakiga (AGR)"
    ind.info[ind.info$ind %in% rhg.ids$V1,]$pop = "Batwa (RHG)"
    ind.info[ind.info$ind %in% gbr.ids$V1,]$pop = "GBR 1000G (AGR)"

    ind.info$pop = factor(ind.info$pop,
        levels = c("Batwa (RHG)", "Bakiga (AGR)", "GBR 1000G (AGR)"))

} else {
    ind.info[grepl("BIR", ind.info$ind),]$pop = "Birhor (AGR)"
    ind.info[grepl("BR1", ind.info$ind),]$pop = "Brahmin (AGR)"
    ind.info[grepl("ONG", ind.info$ind),]$pop = "Onge (RHG)"
    ind.info[grepl("JAR", ind.info$ind),]$pop = "Jarawa (RHG)"
    ind.info[ind.info$ind %in% lwk.ids$V1,]$pop = "LWK 1000G (AGR)"

    ind.info$pop = factor(ind.info$pop,
        levels = c("Onge (RHG)", "Jarawa (RHG)", "Brahmin (AGR)", "LWK 1000G (AGR)"))
}

ind.info = ind.info[!is.na(ind.info$pop),]

# ----------------------------------------------------------------------------------------

make.adm.plot = function (k) {

    this.adm.file = paste0(adm.prefix, k, ".Q")

    adm = read.table(this.adm.file)
    names(adm) = paste0("ADM_", 1:as.numeric(k))
    adm = cbind(adm.inds, adm)

    adm = merge(adm, ind.info, by.x="adm.inds", by.y="ind")

    num.ind = nrow(adm)

    order.list = append(
                    append(
                        list(adm$pop),
                        lapply(1:10, function (x) { adm[[paste0("ADM_", x)]] })),
                    list(adm$adm.inds)
        )

    # Re-order to match now that rearranging done
    ordered.pops = levels(factor(adm$pop))

    # Remove NULL for missing components (greater than current K)
    order.list = order.list[!sapply(order.list, is.null)]

    # Nice trick to "convert" list to ellipsis from http://stackoverflow.com/a/18545884
    new.order = do.call(function(...) order(...), order.list)

    # Figure out where to draw lines between populations
    cats = as.character(adm[new.order,]$pop)
    last.for.pop = which(cats != c(cats[-1], NA))

    pal = brewer.pal(as.numeric(k), "Set3")

    pdf.out = paste0(adm.prefix, k, ".pdf")

    pdf(file=pdf.out, width=15, height=5)

        par(mar=c(8,4,4,2))

        mp = barplot(t(as.matrix(adm[new.order, grep("ADM",names(adm))])),
            col=pal,
            xlab="", ylab="",
            border=NA, xaxt='n', yaxt='n',
            main=paste0("ADMIXTURE (k=", k, ")"))

        abline(v=mp[last.for.pop] + 0.5 * mp[1], lwd=3)

        cat.midpoints = c(last.for.pop, nrow(adm)) - 0.5 * table(adm$pop)
        axis(1, at=mp[cat.midpoints], labels=ordered.pops, las=1, cex.axis=1)

    dev.off()
}

lapply(2:6, make.adm.plot)
