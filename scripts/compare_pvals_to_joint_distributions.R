#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Test the observed p-values for Batwa PBS and for Anda PBS
# --- separately against their random counterparts
# --- to figure out how often the observe values are more extreme than both
# ----------------------------------------------------------------------------------------

# For this GO term, how often is random noise more extreme than both the
# Anda and Batwa values individually?

args = commandArgs(trailingOnly=TRUE)
pop1      = args[1]   # e.g. Batwa.GBR
pop2      = args[2]   # e.g. Anda.LWK
ont       = args[3]   # e.g. BP
test.type = args[4]   # "overrep" or "enrich"
cor.type  = args[5]   # Empty string or "_sizeCor" or "_mafCor"

if (is.na(cor.type)) {
    cor.type = ""
}

in.file.pop1 = paste0("results", cor.type, "/pbs.", pop1, ".pvals.", ont,
    ".", test.type, ".results.txt")
in.file.pop2 = paste0("results", cor.type, "/pbs.", pop2, ".pvals.", ont,
    ".", test.type, ".results.txt")

# --- Read in actual p-values

obs.pop1 = read.table(in.file.pop1, header=TRUE, sep="\t")
obs.pop2 = read.table(in.file.pop2, header=TRUE, sep="\t")

# --- Read in bootstrapped p-values

# With gene-PBS permutation
bs.pop1.orig = read.table(paste0("reports", cor.type, "/pbs.", pop1,
    ".pvals.", ont, ".bootstrapped.GO.p-values.", test.type, ".txt"))
bs.pop2.orig = read.table(paste0("reports", cor.type, "/pbs.", pop2,
    ".pvals.", ont, ".bootstrapped.GO.p-values.", test.type, ".txt"))

# Same, but with GO-gene permutation
bs.pop1.gg = read.table(paste0("reports", cor.type, "/pbs.", pop1,
    ".pvals.", ont, ".bootstrapped-geneGO.GO.p-values.", test.type, ".txt"))
bs.pop2.gg = read.table(paste0("reports", cor.type, "/pbs.", pop2,
    ".pvals.", ont, ".bootstrapped-geneGO.GO.p-values.", test.type, ".txt"))

bs.names = c("GO.ID", "Term", paste0("bs.iter.", 1:(ncol(bs.pop1.orig) - 2)))

names(bs.pop1.orig) = names(bs.pop2.orig) =
    names(bs.pop1.gg) = names(bs.pop2.gg) = bs.names

# --- For each GO, see how often both pops are more extreme than
# --- randomly generated values

all.gos = sort(unique(c(obs.pop1$GO.ID, obs.pop2$GO.ID)))

compute.empirical.p = function (go, bs.pop1, bs.pop2) {

    # Return NA if this GO not assessed in one population
    if ((! go %in% bs.pop1$GO.ID) || (! go %in% bs.pop2$GO.ID)) {
        return (NA)
    }

    # Distribution of p-values
    pop1.p.dist = bs.pop1[bs.pop1$GO.ID == go, -c(1:2)]
    pop2.p.dist = bs.pop2[bs.pop2$GO.ID == go, -c(1:2)]

    # Actual p-values
    p.col.name = "classicFisher"
    if (test.type == "enrich") {
        p.col.name = "classicKS"
    }

    p.pop1 = obs.pop1[obs.pop1$GO.ID == go,][[p.col.name]]
    p.pop2 = obs.pop2[obs.pop2$GO.ID == go,][[p.col.name]]

    # # Boolean string - how many random values are more extreme than both by chance?
    # both.bool = (pop1.p.dist <= p.pop1) & (pop2.p.dist <= p.pop2)

    # Boolean string:
    #   # 1 - how often our values are both lower than the random iterations
    #   both.bool = !((p.pop1 <= pop1.p.dist) & (p.pop2 <= pop2.p.dist))

    # how often BOTH of the random iterations are lower than our observed values
    #   both.bool = (pop1.p.dist <= p.pop1) & (pop2.p.dist <= p.pop2)

    # how often:
    #     (1-random.batwa)*(1-random.andaman) >= (1-obs.batwa)*(1-obs.andaman)
    both.bool = ((1 - pop1.p.dist) * (1 - pop2.p.dist)) >=
                ((1 -      p.pop1) * (1 -      p.pop2))

    # Compute empirical p-value
    p = sum(both.bool) / length(both.bool)
}

all.p.vals.orig = lapply(all.gos, function (go) {
    compute.empirical.p (go, bs.pop1 = bs.pop1.orig, bs.pop2 = bs.pop2.orig)
})

all.p.vals.gg   = lapply(all.gos, function (go) {
    compute.empirical.p (go, bs.pop1 = bs.pop1.gg,   bs.pop2 = bs.pop2.gg)
})

go.p.vals.orig = data.frame(GO.ID = all.gos, joint.bs.p = unlist(all.p.vals.orig))
go.p.vals.gg   = data.frame(GO.ID = all.gos, joint.bs.p = unlist(all.p.vals.gg))

# Add this column to original data.frames.
obs.pop1.w.go = merge(merge(obs.pop1, go.p.vals.orig, by="GO.ID", sort=FALSE),
                      go.p.vals.gg, by="GO.ID", sort=FALSE,
                      suffix=c(".genePBS", ".geneGO"))
obs.pop2.w.go = merge(merge(obs.pop2, go.p.vals.orig, by="GO.ID", sort=FALSE),
                      go.p.vals.gg, by="GO.ID", sort=FALSE,
                      suffix=c(".genePBS", ".geneGO"))

write.table(obs.pop1.w.go, file=gsub(".txt", ".wGOiterp.txt", in.file.pop1))
write.table(obs.pop2.w.go, file=gsub(".txt", ".wGOiterp.txt", in.file.pop2))

# Make merged data.frame of both
obs.both.w.go = merge(obs.pop1.w.go, obs.pop2.w.go,
    by=c("GO.ID", "Term", "joint.bs.p.genePBS", "joint.bs.p.geneGO"),
    suffix=paste0(".", c(pop1, pop2)))
obs.both.w.go = obs.both.w.go[order(obs.both.w.go$joint.bs.p.genePBS),]

combined.out.file = paste0("results", cor.type, "/convergent_permuted_pval_",
    pop1, "-", pop2, ".", ont, ".", test.type, ".txt")
write.table(obs.both.w.go, file=combined.out.file, sep="\t", row.names=FALSE)
