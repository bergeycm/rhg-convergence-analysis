#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ========================================================================================
# --- Functions for bootstrapping input to GO to generate empirical p-value distributions
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Do bootstrapped overrep test iteration - permute PBS-gene relationship
# ----------------------------------------------------------------------------------------

do.overrep.bs.iter = function () {

    # Randomize gene list
    gene.names = names(all.genes)
    all.genes.rand = sample(all.genes)
    names(all.genes.rand) = gene.names

    GOdata.bs = new("topGOdata",
        description = description,
        ontology = ont,
        allGenes = all.genes.rand,
        geneSel = top.interest.genes,
        annot = annFUN.org,
        ID = "alias",
        mapping = "org.Hs.eg",
        nodeSize = 50)

    result.fisher.classic.bs  = runTest(GOdata.bs, "classic",  "fisher")

    # Get all p-values
    allGO.bs = usedGO(object = GOdata.bs)

    allRes.overrep.full.bs = GenTable(GOdata.bs,
        classicFisher = result.fisher.classic.bs,
        orderBy = "classicFisher",
        topNodes = length(allGO.bs))

    allRes.overrep.full.bs = allRes.overrep.full.bs[order(allRes.overrep.full.bs$GO.ID),]

    return(allRes.overrep.full.bs[,c(1,2,6)])
}

bs.to.get.GO.pval.dist.overrep = function (num.iter) {

    cl = makeCluster(detectCores() - 1)
    clusterExport(cl, list("do.overrep.bs.iter", "all.genes", "num.iter", "description",
                           "ont", "top.interest.genes", "annFUN.org",
                           "GOBPTerm", "GOMFTerm", "GOCCTerm",
                           "runTest", "is.pval", "cutoff", "usedGO",
                           "GenTable"), envir=environment())

    go.res.bp = parLapply(cl = cl, 1:num.iter, function (x) { do.overrep.bs.iter() })

    bs.p = do.call(cbind, lapply(1:length(go.res.bp), function (x) {
        as.numeric(go.res.bp[[x]][[3]])
    }))

    bs.p = data.frame(cbind(go.res.bp[[1]][[1]], go.res.bp[[1]][[2]], bs.p))
}

# ----------------------------------------------------------------------------------------
# --- Do bootstrapped enrich test iteration - permute PBS-gene relationship
# ----------------------------------------------------------------------------------------

do.enrich.bs.iter = function () {

    # Randomize gene list
    gene.names = names(all.genes)
    all.genes.rand = sample(all.genes)
    names(all.genes.rand) = gene.names

    GOdata.bs = new("topGOdata",
        description = description,
        ontology = ont,
        allGenes = all.genes.rand,
        geneSel = top.interest.genes,
        annot = annFUN.org,
        ID = "alias",
        mapping = "org.Hs.eg",
        nodeSize = 50)

    result.ks.classic.bs  = runTest(GOdata.bs, "classic",  "ks")

    # Get all p-values
    allGO.bs = usedGO(object = GOdata.bs)

    allRes.enrich.full.bs = GenTable(GOdata.bs,
        classicKS = result.ks.classic.bs,
        orderBy = "classicKS",
        topNodes = length(allGO.bs))

    allRes.enrich.full.bs = allRes.enrich.full.bs[order(allRes.enrich.full.bs$GO.ID),]

    return(allRes.enrich.full.bs[,c(1,2,6)])
}

bs.to.get.GO.pval.dist.enrich = function (num.iter) {

    cl = makeCluster(detectCores() - 1)
    clusterExport(cl, list("do.enrich.bs.iter", "all.genes", "num.iter", "description",
                           "ont", "top.interest.genes", "annFUN.org",
                           "GOBPTerm", "GOMFTerm", "GOCCTerm",
                           "runTest", "is.pval", "cutoff", "usedGO",
                           "GenTable"), envir=environment())

    go.res.bp = parLapply(cl = cl, 1:num.iter, function (x) { do.enrich.bs.iter() })

    bs.p = do.call(cbind, lapply(1:length(go.res.bp), function (x) {
        as.numeric(go.res.bp[[x]][[3]])
    }))

    bs.p = data.frame(cbind(go.res.bp[[1]][[1]], go.res.bp[[1]][[2]], bs.p))
}

# ----------------------------------------------------------------------------------------
# --- Do bootstrapped overrep test iteration - permute gene-GO relationship
# ----------------------------------------------------------------------------------------

do.overrep.bs.iter.GOgene = function () {

    # Randomize GO-to-gene mapping
    org.Hs.eg.db.gene2go.rand = org.Hs.eg.db.gene2go
    names(org.Hs.eg.db.gene2go.rand) = sample(names(org.Hs.eg.db.gene2go.rand))

    GOdata.bs = new("topGOdata",
        description = description,
        ontology = ont,
        allGenes = all.genes,
        geneSel = top.interest.genes,
        annot = annFUN.gene2GO,
        gene2GO = org.Hs.eg.db.gene2go.rand,
        nodeSize = 50)

    result.fisher.classic.bs  = runTest(GOdata.bs, "classic",  "fisher")

    # Get all p-values
    allGO.bs = usedGO(object = GOdata.bs)

    allRes.overrep.full.bs = GenTable(GOdata.bs,
        classicFisher = result.fisher.classic.bs,
        orderBy = "classicFisher",
        topNodes = length(allGO.bs))

    allRes.overrep.full.bs = allRes.overrep.full.bs[order(allRes.overrep.full.bs$GO.ID),]

    return(allRes.overrep.full.bs[,c(1,2,6)])
}

bs.to.get.GO.pval.dist.overrep.GOgene = function (num.iter) {

    org.Hs.eg.db.gene2go.full = select(org.Hs.eg.db, keys=names(all.genes),
        columns=c("SYMBOL","GO"), keytype="SYMBOL")

    org.Hs.eg.db.gene2go.list = split(org.Hs.eg.db.gene2go.full,
                                      org.Hs.eg.db.gene2go.full$SYMBOL)

    org.Hs.eg.db.gene2go = lapply(org.Hs.eg.db.gene2go.list, function (gene.gos) {
        unique(gene.gos[gene.gos$ONTOLOGY == ont,]$GO)
    })

    cl = makeCluster(detectCores() - 1)
    clusterExport(cl, list("do.overrep.bs.iter.GOgene", "all.genes",
                           "num.iter", "description",
                           "ont", "top.interest.genes",
                           "annFUN.gene2GO", "org.Hs.eg.db.gene2go",
                           "GOBPTerm", "GOMFTerm", "GOCCTerm",
                           "runTest", "is.pval", "cutoff", "usedGO",
                           "GenTable"), envir=environment())

    go.res.bp = parLapply(cl = cl, 1:num.iter,
        function (x) { do.overrep.bs.iter.GOgene() })

    bs.p = do.call(cbind, lapply(1:length(go.res.bp), function (x) {
        as.numeric(go.res.bp[[x]][[3]])
    }))

    bs.p = data.frame(cbind(go.res.bp[[1]][[1]], go.res.bp[[1]][[2]], bs.p))
}

# ----------------------------------------------------------------------------------------
# --- Do bootstrapped enrich test iteration - permute gene-GO relationship
# ----------------------------------------------------------------------------------------

do.enrich.bs.iter.GOgene = function () {

    # Randomize GO-to-gene mapping
    org.Hs.eg.db.gene2go.rand = org.Hs.eg.db.gene2go
    names(org.Hs.eg.db.gene2go.rand) = sample(names(org.Hs.eg.db.gene2go.rand))

    GOdata.bs = new("topGOdata",
        description = description,
        ontology = ont,
        allGenes = all.genes,
        geneSel = top.interest.genes,
        annot = annFUN.gene2GO,
        gene2GO = org.Hs.eg.db.gene2go.rand,
        nodeSize = 50)

    result.ks.classic.bs  = runTest(GOdata.bs, "classic",  "ks")

    # Get all p-values
    allGO.bs = usedGO(object = GOdata.bs)

    allRes.enrich.full.bs = GenTable(GOdata.bs,
        classicKS = result.ks.classic.bs,
        orderBy = "classicKS",
        topNodes = length(allGO.bs))

    allRes.enrich.full.bs = allRes.enrich.full.bs[order(allRes.enrich.full.bs$GO.ID),]

    return(allRes.enrich.full.bs[,c(1,2,6)])
}

bs.to.get.GO.pval.dist.enrich.GOgene = function (num.iter) {

    org.Hs.eg.db.gene2go.full = select(org.Hs.eg.db, keys=names(all.genes),
        columns=c("SYMBOL","GO"), keytype="SYMBOL")

    org.Hs.eg.db.gene2go.list = split(org.Hs.eg.db.gene2go.full,
                                      org.Hs.eg.db.gene2go.full$SYMBOL)

    org.Hs.eg.db.gene2go = lapply(org.Hs.eg.db.gene2go.list, function (gene.gos) {
        unique(gene.gos[gene.gos$ONTOLOGY == ont,]$GO)
    })

    cl = makeCluster(detectCores() - 1)
    clusterExport(cl, list("do.enrich.bs.iter.GOgene", "all.genes",
                           "num.iter", "description",
                           "ont", "top.interest.genes",
                           "annFUN.gene2GO", "org.Hs.eg.db.gene2go",
                           "GOBPTerm", "GOMFTerm", "GOCCTerm",
                           "runTest", "is.pval", "cutoff", "usedGO",
                           "GenTable"), envir=environment())

    go.res.bp = parLapply(cl = cl, 1:num.iter,
        function (x) { do.enrich.bs.iter.GOgene() })

    bs.p = do.call(cbind, lapply(1:length(go.res.bp), function (x) {
        as.numeric(go.res.bp[[x]][[3]])
    }))

    bs.p = data.frame(cbind(go.res.bp[[1]][[1]], go.res.bp[[1]][[2]], bs.p))
}
