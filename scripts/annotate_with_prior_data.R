#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ========================================================================================
# --- Match windows with previously ID'd regions/genes/SNPs of interest
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Bring in results of parse_windowed_stats
# ----------------------------------------------------------------------------------------

stats = read.table(file="results/basic_stats.txt", sep="\t", header=TRUE)

# ----------------------------------------------------------------------------------------
# --- Bring in outside info on regions of interest
# ----------------------------------------------------------------------------------------

# --- Perry et al 2014 data
perry = list()

# GWAS associated regions
perry[['gwas']] = read.table(paste0("data/Perry_etal_2014/GWAS_assoc_regions/",
    "Perry_etal_2014-S2-batwa_assoc.hg19.bed"))

# BayeScan
perry[['bayescan']] = list()
perry[['bayescan']][['baka-vs-nzebi-bzime']] = read.table(paste0("data/Perry_etal_2014/",
        "BayeScan_outlier_regions/",
        "Perry_etal_2014-S3-baka-vs-nzebi-bzime_BayeScan_outliers.hg19.bed"), 
    sep="\t", fill=TRUE)
perry[['bayescan']][['batwa-vs-bakiga']] = read.table(paste0("data/Perry_etal_2014/",
        "BayeScan_outlier_regions/",
        "Perry_etal_2014-S3-batwa-vs-bakiga_BayeScan_outliers.hg19.bed"), 
    sep="\t", fill=TRUE)

# iHS
perry[['ihs']] = list()
pops = c('baka', 'bakiga', 'batwa', 'nzebi-nzime')
for (pop in pops) {
    perry[['ihs']][[pop]] = read.table(paste0("data/Perry_etal_2014/iHS_outlier_regions/",
        "Perry_etal_2014-S3-", pop, "_iHS_outliers.hg19.bed"),
        sep="\t", fill=TRUE, row.names=NULL)
}

# --- Jarvis et al 2012 and Lachance et al 2012 data
jl = read.table("data/Jarvis_Lachance-pygmy_assoc.bed", header=TRUE)

# --- Wood et al 2014 data
# OMIM genes

# ----------------------------------------------------------------------------------------
# --- Figure out which of our windows overlap these prior ROIs
# ----------------------------------------------------------------------------------------

find.overlap = function(chr, start, end) {
    matching.stat.row.ids = which(chr == stats$chr & 
        end >= stats$win_start & start <= stats$win_end)
    return(matching.stat.row.ids)
}

# --- Perry et al 2014 data

perry.overlap.ids = list()

perry.overlap.ids[['gwas']] = unlist(sapply(1:nrow(perry[['gwas']]), 
    function(x) { 
        find.overlap(perry[['gwas']][x,1], perry[['gwas']][x,2], perry[['gwas']][x,3])
    }))

perry.overlap.ids[['bayescan']] = list()

for (pair in c("baka-vs-nzebi-bzime", "batwa-vs-bakiga")) {
    
    perry.overlap.ids[['bayescan']][[pair]] = unlist(sapply(
        1:nrow(perry[['bayescan']][[pair]]), 
        function(x) { 
            find.overlap(perry[['bayescan']][[pair]][x,1], 
                perry[['bayescan']][[pair]][x,2], 
                perry[['bayescan']][[pair]][x,3])
        }
    ))
}

perry.overlap.ids[['ihs']] = list()

for (pop in pops) {
    perry.overlap.ids[['ihs']][[pop]] = unlist(sapply(
        1:nrow(perry[['ihs']][[pop]]), 
        function(x) { 
            find.overlap(perry[['ihs']][[pop]][x,1], 
                perry[['ihs']][[pop]][x,2], 
                perry[['ihs']][[pop]][x,3])
        }
    ))
}

# --- Jarvis et al 2012 and Lachance et al 2012 data

jl.overlap.ids = unlist(sapply(1:nrow(jl), 
    function(x) { find.overlap(jl[x,1], jl[x,2], jl[x,3]) }))

# ----------------------------------------------------------------------------------------
# --- Add this overlap info to data.frame
# ----------------------------------------------------------------------------------------

# --- Perry et al 2014 data

# GWAS
stats$overlap.perry.gwas = FALSE
stats[perry.overlap.ids[['gwas']],]$overlap.perry.gwas = TRUE

# BayeScan
stats$overlap.perry.bs.baka.nzebibzime = FALSE
matches = perry.overlap.ids[['bayescan']][['baka-vs-nzebi-bzime']]
stats[matches,]$overlap.perry.bs.baka.nzebibzime = TRUE
rm("matches")

stats$overlap.perry.bs.batwa.bakiga = FALSE
matches = perry.overlap.ids[['bayescan']][['batwa-vs-bakiga']]
stats[matches,]$overlap.perry.bs.batwa.bakiga = TRUE
rm("matches")

# iHS
for (pop in pops) {
    col.name = paste0('overlap.perry.ihs.', pop)
    stats[col.name] = FALSE
    matches = perry.overlap.ids[['ihs']][[pop]]
    stats[matches,][col.name] = TRUE
    rm(list=c("matches", "col.name"))
}

# --- Jarvis et al 2012 and Lachance et al 2012 data

stats$overlap.jl = FALSE
stats[jl.overlap.ids,]$overlap.jl = TRUE

# ----------------------------------------------------------------------------------------
# --- Write file of windowed stats with annotations
# ----------------------------------------------------------------------------------------

write.table(stats, file="results/basic_stats.anno.txt", 
    sep="\t", row.names=FALSE, quote=FALSE)
