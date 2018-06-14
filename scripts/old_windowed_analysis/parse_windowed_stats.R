#!/usr/bin/env Rscript

# ========================================================================================
# --- Gather together all stats
# ========================================================================================

classes = c('all', 'syn', 'ns')
groups = c('agr', 'rhg')

stats.l = list()

# ----------------------------------------------------------------------------------------
# --- Number of exonic bases
# ----------------------------------------------------------------------------------------

stats.l[['cov']] = read.table("results/windows_100kb.exon_cov.bed", header=TRUE)

# ----------------------------------------------------------------------------------------
# --- Fst
# ----------------------------------------------------------------------------------------

# --- Per-SNP

stats.l[['fst']] = list()

stats.l[['fst']][['all']] = read.table("results/eAGR_eRHG.weir.fst", header=TRUE)

# fst.all = stats.l[['fst']][['all']]
# fst.all[fst.all$WEIR_AND_COCKERHAM_FST > 0.5,]

stats.l[['fst']][['syn']] = read.table("results/eAGR_eRHG.syn.weir.fst",
    header=TRUE)
stats.l[['fst']][['ns']]  = read.table("results/eAGR_eRHG.nonsyn.weir.fst", 
    header=TRUE)

# --- Windowed

stats.l[['fst.win']] = list()

stats.l[['fst.win']][['all']] = read.table("results/eAGR_eRHG.windowed.weir.fst", 
    header=TRUE)

stats.l[['fst.win']][['syn']] = read.table("results/eAGR_eRHG.syn.windowed.weir.fst",
    header=TRUE)
stats.l[['fst.win']][['ns']]  = read.table("results/eAGR_eRHG.nonsyn.windowed.weir.fst", 
    header=TRUE)

# ----------------------------------------------------------------------------------------
# --- Nucleotide diversity (windowed)
# ----------------------------------------------------------------------------------------

stats.l[['pi']] = list()

stats.l[['pi']][['all']] = list()
stats.l[['pi']][['all']][['agr']] = read.table("results/eAGR.windowed.pi", header=TRUE)
stats.l[['pi']][['all']][['rhg']] = read.table("results/eRHG.windowed.pi", header=TRUE)

stats.l[['pi']][['syn']] = list()
stats.l[['pi']][['syn']][['agr']] = read.table("results/eAGR.syn.windowed.pi", 
    header=TRUE)
stats.l[['pi']][['syn']][['rhg']] = read.table("results/eRHG.syn.windowed.pi", 
    header=TRUE)

stats.l[['pi']][['ns']] = list()
stats.l[['pi']][['ns']][['agr']] = read.table("results/eAGR.nonsyn.windowed.pi", 
    header=TRUE)
stats.l[['pi']][['ns']][['rhg']] = read.table("results/eRHG.nonsyn.windowed.pi", 
    header=TRUE)

# ----------------------------------------------------------------------------------------
# --- Tajima's D
# ----------------------------------------------------------------------------------------

stats.l[['tajd']] = list()

stats.l[['tajd']][['all']] = list()
stats.l[['tajd']][['all']][['agr']] = read.table("results/eAGR.Tajima.D", header=TRUE)
stats.l[['tajd']][['all']][['rhg']] = read.table("results/eRHG.Tajima.D", header=TRUE)

stats.l[['tajd']][['syn']] = list()
stats.l[['tajd']][['syn']][['agr']] = read.table("results/eAGR.syn.Tajima.D", header=TRUE)
stats.l[['tajd']][['syn']][['rhg']] = read.table("results/eRHG.syn.Tajima.D", header=TRUE)

stats.l[['tajd']][['ns']] = list()
stats.l[['tajd']][['ns']][['agr']] = read.table("results/eAGR.nonsyn.Tajima.D", 
    header=TRUE)
stats.l[['tajd']][['ns']][['rhg']] = read.table("results/eRHG.nonsyn.Tajima.D", 
    header=TRUE)

# ========================================================================================
# --- Bring all stats together in data.frame
# ========================================================================================

# --- Exonic coverage
stats = stats.l[['cov']]
stats$win_start = stats$win_start + 1

stats$exon_filter_pass = TRUE
stats[stats$num_exonic_bases < 500,]$exon_filter_pass = FALSE

# --- Fst

for (x in classes) {
    names(stats.l[['fst.win']][[x]])[4:6] = paste0("fst.", x, ".", 
        names(stats.l[['fst.win']][[x]])[4:6])
}

merge.fst = function(i = 1) {
    if (i < length(stats.l[['fst.win']])) {
        return(merge(stats.l[['fst.win']][[i]], merge.fst(i + 1),
            sort=FALSE, all=TRUE))
    } else {
        return(stats.l[['fst.win']][[i]])
    }
}

fst.merge = merge.fst()

names(fst.merge)[1:3] = names(stats)[1:3]

stats = merge(stats, fst.merge, sort=FALSE, all=TRUE)

# --- pi

for (x in classes) {
    for (grp in groups) {
        names(stats.l[['pi']][[x]][[grp]])[4:5] = paste0("pi.", x, ".", grp, ".",
            names(stats.l[['pi']][[x]][[grp]])[4:5])
    }
}

merge.stat = function(stat, i = 1) {
    both.pops = merge(stats.l[[stat]][[i]][['agr']], stats.l[[stat]][[i]][['rhg']],
        sort=FALSE, all=TRUE)
    if (i < length(stats.l[[stat]])) {
        return(merge(both.pops, merge.stat(stat, i + 1),
            sort=FALSE, all=TRUE))
    } else {
        return(both.pops)
    }
}

pi.merge = merge.stat(stat='pi')

names(pi.merge)[1:3] = names(stats)[1:3]

stats = merge(stats, pi.merge, sort=FALSE, all=TRUE)

# --- Compute and add pi ratio

# Fix SNPs that have NA for pi (set to zero)
for (x in classes) {
    for (grp in groups) {
        col.index = which(names(stats) == paste0('pi.', x, ".", grp, ".PI"))
        if (sum(is.na(stats[col.index])) > 0) {
            stats[is.na(stats[col.index]),][col.index] = 0
        }
    }
}

# Compute ratio of pi between populations (store in AGR data.frame)
for (x in classes) {

    this.agr.pi = stats[[paste0('pi.', x, ".agr.PI")]]
    this.rhg.pi = stats[[paste0('pi.', x, ".rhg.PI")]]

    # Add 0.00005 to both to avoid dividing by zero
    this.agr.pi = this.agr.pi + 0.00005
    this.rhg.pi = this.rhg.pi + 0.00005
    
    this.ratio = this.agr.pi / this.rhg.pi
    
    stats[[paste0('pi.', x, ".ratio")]] = this.ratio
}

# --- Tajima's D

for (x in classes) {
    for (grp in groups) {
        names(stats.l[['tajd']][[x]][[grp]])[3:4] = paste0("tajd.", x, ".", grp, ".",
            names(stats.l[['tajd']][[x]][[grp]])[3:4])
        all.starts = stats.l[['tajd']][[x]][[grp]]$BIN_START
        stats.l[['tajd']][[x]][[grp]]$BIN_START = all.starts + 1
        stats.l[['tajd']][[x]][[grp]]$BIN_END = all.starts + 100000
    }
}

tajd.merge = merge.stat(stat='tajd')

names(tajd.merge)[1:3] = names(stats)[1:3]

stats = merge(stats, tajd.merge, sort=FALSE, all=TRUE)

# ========================================================================================
# --- Compute averages for all stats, filtering for number of exonic bases
# ========================================================================================

# Remove windows with too few exonic bases
stats.flt = stats[stats$exon_filter_pass,]

# Average all the things
means = colMeans(stats.flt[4:ncol(stats.flt)], na.rm = TRUE)

# Count number of passing windows
win.ct = nrow(stats.flt)

# Count number of passing windows that have no variants
no.var.win.ct = nrow(stats.flt[which(stats.flt$pi.all.agr.PI == 0 & 
                     stats.flt$pi.all.rhg.PI == 0),])

# --- Write out info

sink(file='reports/basic_stat_averages.txt')

paste("Overall mean pi - AGR -", 
    sprintf("%.6f", 100 * means[['pi.all.agr.PI']]), "%")
paste("Overall mean pi - RHG -", 
    sprintf("%.6f", 100 * means[['pi.all.rhg.PI']]), "%")

paste("Mean pi at synonymous sites - AGR -", 
    sprintf("%.6f", 100 * means[['pi.syn.agr.PI']]), "%")
paste("Mean pi at synonymous sites - RHG -", 
    sprintf("%.6f", 100 * means[['pi.syn.rhg.PI']]), "%")

paste("Mean pi at nonsynonymous sites - AGR -", 
    sprintf("%.6f", 100 * means[['pi.ns.agr.PI']]), "%")
paste("Mean pi at nonsynonymous sites - RHG -", 
    sprintf("%.6f", 100 * means[['pi.ns.rhg.PI']]), "%")

paste("---")

paste("Overall mean Tajima's D - AGR -", 
    sprintf("%.6f", means[['tajd.all.agr.TajimaD']]))
paste("Overall mean Tajima's D - RHG -", 
    sprintf("%.6f", means[['tajd.all.rhg.TajimaD']]))

paste("Mean Tajima's D at synonymous sites - AGR -", 
    sprintf("%.6f", means[['tajd.syn.agr.TajimaD']]))
paste("Mean Tajima's D at synonymous sites - RHG -", 
    sprintf("%.6f", means[['tajd.syn.rhg.TajimaD']]))

paste("Mean Tajima's D at nonsynonymous sites - AGR -", 
    sprintf("%.6f", means[['tajd.ns.agr.TajimaD']]))
paste("Mean Tajima's D at nonsynonymous sites - RHG -", 
    sprintf("%.6f", means[['tajd.ns.rhg.TajimaD']]))

paste("---")

paste("Number of passing windows without polymorphism -", 
    sprintf("%.2f", 100 * no.var.win.ct / win.ct), "%")

sink()

# --- Write all merged stats to file

write.table(stats, file="results/basic_stats.txt", sep="\t", row.names=FALSE, quote=FALSE)
