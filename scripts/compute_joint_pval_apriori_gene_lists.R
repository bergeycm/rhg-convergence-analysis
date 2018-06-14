#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Compute joint p-values for a priori gene lists
# ----------------------------------------------------------------------------------------

library(metap)

overrep.batwa = read.table("reports/pbs.Batwa.GBR.anno.chisq.overrep.flt.txt",
                           header=TRUE, sep="\t")
overrep.anda  = read.table("reports/pbs.Anda.LWK.anno.chisq.overrep.flt.txt",
                           header=TRUE, sep="\t")

overrep.both.p = cbind(overrep.batwa[,1:2],overrep.batwa$pvalue, overrep.anda$pvalue,
    unlist(lapply(
        apply(cbind(overrep.batwa$pvalue, overrep.anda$pvalue), 1, sump),
        function(x) {x$p})
    ),
    unlist(lapply(
        apply(cbind(overrep.batwa$pvalue, overrep.anda$pvalue), 1, sumlog),
        function(x) {x$p})
    ))

names(overrep.both.p) = c("dataset", "p.val.cutoffs", "batwa.p", "anda.p",
    "edg.p", "fish.p")

sig = overrep.both.p[overrep.both.p$fish.p < 0.05 & overrep.both.p$p.val.cutoffs == 0.05,]

# ----------------------------------------------------------------------------------------

enrich.batwa = read.table("reports/pbs.Batwa.GBR.anno.proptest.flt.txt",
                           header=TRUE, sep=" ")
enrich.anda  = read.table("reports/pbs.Anda.LWK.anno.proptest.flt.txt",
                           header=TRUE, sep=" ")

# Set to 1 if NA to exclude from analysis, effectively
if (sum(is.na(enrich.batwa$prop.pvals)) > 0) {
    enrich.batwa[is.na(enrich.batwa$prop.pvals),]$prop.pvals = 1
}
if (sum(is.na(enrich.anda$prop.pvals)) > 0) {
    enrich.anda[is.na(enrich.anda$prop.pvals),]$prop.pvals = 1
}
enrich.both.p = cbind(enrich.batwa[,1:2],enrich.batwa$prop.pvals, enrich.anda$prop.pvals,
    unlist(lapply(
        apply(cbind(enrich.batwa$prop.pvals, enrich.anda$prop.pvals), 1, sump),
        function(x) {x$p})
    ),
    unlist(lapply(
        apply(cbind(enrich.batwa$prop.pvals, enrich.anda$prop.pvals), 1, sumlog),
        function(x) {x$p})
    ))

names(enrich.both.p) = c("dataset", "p.val.cutoffs", "batwa.p", "anda.p",
    "edg.p", "fish.p")

sig = enrich.both.p[enrich.both.p$fish.p < 0.05 & enrich.both.p$p.val.cutoffs == 0.05,]
