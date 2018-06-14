#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Compute PBS
# ----------------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
pop.out = args[1]

if (pop.out == "GBR") {
    pop.of.interest = "Batwa"
    pops.sister     = c("Bakiga")
} else {
    pop.of.interest = "Anda"
    pops.sister     = c("BIR", "BR1")
}

for (pop.sister in pops.sister) {

    pops = sort(c(pop.of.interest, pop.sister, pop.out))

    pop.pairs = combn(pops, 2, simplify=FALSE)

    names(pop.pairs) = do.call(c, lapply(pop.pairs, function(x) {
        paste(x, collapse="_")
    } ))

    fst = list()

    fst = lapply(names(pop.pairs), function (x) {
        this.fst = read.table(paste0("results/PBS_", x, ".weir.fst"), header=TRUE)
        this.fst$pop1 = pop.pairs[x][[1]][1]
        this.fst$pop2 = pop.pairs[x][[1]][2]
        return(this.fst)
    })

    names(fst) = names(pop.pairs)

    # Merge all FST dataframes together
    merge.cols = names(fst[[1]])[1:2]
    tmp = merge(fst[[ names(pop.pairs)[1] ]], fst[[ names(pop.pairs)[2] ]], by=merge.cols)
    all.fst = merge(tmp, fst[[ names(pop.pairs)[3] ]], by=merge.cols)

    names(all.fst)[3:ncol(all.fst)] = apply(expand.grid(c("FST", "pop1", "pop2"),
        names(pop.pairs)), 1, paste, collapse=".")

    # Compute time since divergence using "classical transformation by Cavalli-Sforza"
    # See Yi et al 2010

    all.T = data.frame(sapply(names(pop.pairs), function(x) {
        this.T = -1 * log(1 - all.fst[[paste0("FST.", x)]])
    }))

    names(all.T) = paste0("T.", names(all.T))

    all.fst = cbind(all.fst, all.T)

    na.fst.rows = unique(sort(do.call(c, apply(all.T, 2, function(x) { which(is.na(x)) }))))
    all.fst.noNA = all.fst[-na.fst.rows,]

    # Sum of branches containing focus population, minus branch that lacks focus pop,
    # divided by two
    str.sister.of.interest = paste(sort(c(pop.sister,      pop.of.interest)), collapse="_")
    str.sister.out         = paste(sort(c(pop.sister,      pop.out)),         collapse="_")
    str.of.interest.out    = paste(sort(c(pop.of.interest, pop.out)),         collapse="_")

    str.sister.of.interest = paste0("T.", str.sister.of.interest)
    str.sister.out         = paste0("T.", str.sister.out        )
    str.of.interest.out    = paste0("T.", str.of.interest.out   )

    all.fst.noNA$pbs.of.interest =
        (all.fst.noNA[[str.sister.of.interest]] + all.fst.noNA[[str.of.interest.out]] -
        all.fst.noNA[[str.sister.out]]) / 2

    all.fst.noNA$pbs.sister =
        (all.fst.noNA[[str.sister.of.interest]] + all.fst.noNA[[str.sister.out]] -
        all.fst.noNA[[str.of.interest.out]]) / 2

    # Write all Fst, T, and PBS to file
    all.fst.noNA.bed = cbind(all.fst.noNA[,1:2],
                             all.fst.noNA$POS + 1,
                             all.fst.noNA[,3:ncol(all.fst.noNA)])
    write.table(all.fst.noNA.bed, file=paste0("results/", pop.out, ".FST_for_PBS.txt"),
        sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

    # Sort in correct order
    all.fst.noNA$chr.num = gsub("chr", "", all.fst.noNA$CHROM)
    all.fst.sort = all.fst.noNA[order(as.numeric(all.fst.noNA$chr.num), all.fst.noNA$POS),]

    pbs.bed.of.interest =
        cbind(all.fst.sort[,1:2], all.fst.sort$POS + 1, all.fst.sort$pbs.of.interest)
    pbs.bed.sister =
        cbind(all.fst.sort[,1:2], all.fst.sort$POS + 1, all.fst.sort$pbs.sister)

    write.table(pbs.bed.of.interest,
        file=paste0("results/pbs.", pop.of.interest, ".", pop.out, ".bed"),
        row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    write.table(pbs.bed.sister,
        file=paste0("results/pbs.", pop.sister,      ".", pop.out, ".bed"),
        row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

    # ----------------------------------------------------------------------------------------
    # --- Plot PBS
    # ----------------------------------------------------------------------------------------

    rg = read.table("refGene/refGene.sort.simple.justGenes.gtf")

    names(rg)[c(1,4,5,9)] = c("chr", "start", "end", "gene")

    find.overlap = function(chr, pos) {
        matching.genes = rg[which(chr == rg$chr &
            pos >= rg$start & pos <= rg$end),]$gene
        matching.genes = unique(as.character(matching.genes))
        return(paste0(matching.genes, collapse=','))
    }

    # Just annotate the outliers
    all.fst.noNA$is.outlier = FALSE
    outliers =
        all.fst.noNA$pbs.of.interest > quantile(all.fst.noNA$pbs.of.interest, probs=0.99) |
        all.fst.noNA$pbs.sister      > quantile(all.fst.noNA$pbs.sister,      probs=0.99)
    all.fst.noNA[outliers,]$is.outlier = TRUE

    pbs.overlap = do.call(rbind, lapply(1:nrow(all.fst.noNA),
        function(x) {
            if (all.fst.noNA[x,]$is.outlier) {
                find.overlap(all.fst.noNA[x,1], all.fst.noNA[x,2])
            } else {
                NA
            }
        }
    ))

    all.fst.noNA$gene = pbs.overlap

    library(ggplot2)

    # --- Thin dataset
    all.fst.noNA.noInf = all.fst.noNA[
                           (all.fst.noNA$pbs.of.interest %in% c(Inf, -Inf) == FALSE) &
                           (all.fst.noNA$pbs.sister      %in% c(Inf, -Inf) == FALSE),
                         ]
    samp = sample(1:nrow(all.fst.noNA.noInf), 50000)
    all.fst.thin = all.fst.noNA.noInf[samp,]

    # --- Plot PBS of target population vs. PBS of sister population

    p = ggplot(all.fst.thin, aes(pbs.sister, pbs.of.interest)) +
        geom_point(alpha=0.1) +
        geom_text(data=subset(all.fst.noNA.noInf,
            pbs.sister != Inf & pbs.of.interest != Inf &
            (pbs.of.interest > quantile(pbs.of.interest, 0.999) |
            pbs.sister > quantile(pbs.sister, 0.999))),
            aes(pbs.sister, pbs.of.interest, label=gene),
                check_overlap = FALSE, size=2, col='red', vjust=-1) +
        theme_bw() + xlab(paste("PBS", pop.of.interest)) + ylab(paste("PBS", pop.sister))

    ggsave(paste0("reports/PBS_outliers.", pop.out, ".pdf"), plot=p, height=8, width=8)

    # --- Plot each PBS against Fst

    p = ggplot(all.fst.thin,
            aes_string(
                paste0("FST.", gsub("T.", "", str.sister.of.interest)),
                "pbs.of.interest"
            )) +
        geom_point(alpha=0.1) +
        geom_text(data=subset(all.fst.noNA.noInf,
            pbs.sister != Inf & pbs.of.interest != Inf &
            (pbs.of.interest > quantile(pbs.of.interest, 0.999) |
            pbs.sister > quantile(pbs.sister, 0.999))),
            aes_string(
                paste0("FST.", gsub("T.", "", str.sister.of.interest)),
                "pbs.of.interest",
                label="gene"
            ),
            check_overlap = FALSE, size=2, col='red', vjust=-1) +
        theme_bw() +
        xlim(0,1) +
        xlab(bquote("F"[ST] ~ .(pop.sister) ~ "vs." ~ .(pop.of.interest))) +
        ylab(paste("PBS", pop.of.interest))

    ggsave(paste0("reports/PBS_", pop.sister, "_vs_Fst.pdf"), plot=p, height=8, width=8)

    p = ggplot(all.fst.thin,
            aes_string(
                paste0("FST.", gsub("T.", "", str.sister.of.interest)),
                "pbs.sister"
            )) +
        geom_point(alpha=0.1) +
        geom_text(data=subset(all.fst.noNA.noInf,
            pbs.sister != Inf & pbs.of.interest != Inf &
            (pbs.of.interest > quantile(pbs.of.interest, 0.999) |
            pbs.sister > quantile(pbs.sister, 0.999))),
            aes_string(
                paste0("FST.", gsub("T.", "", str.sister.of.interest)),
                "pbs.sister",
                label="gene"
            ),
            check_overlap = FALSE, size=2, col='red', vjust=-1) +
        theme_bw() +
        xlim(0,1) +
        xlab(bquote("F"[ST] ~ .(pop.sister) ~ "vs." ~ .(pop.of.interest))) +
        ylab(paste("PBS", pop.sister))

    ggsave(paste0("reports/PBS_", pop.of.interest, "_vs_Fst.pdf"), plot=p,
        height=8, width=8)
}
