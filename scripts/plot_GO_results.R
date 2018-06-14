#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Plot GO results
# ----------------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
in.file = args[1]   # E.g. "results/pbs.Batwa.GBR.pvals.BP.results.txt"
out.pdf = args[2]   # E.g. "reports/pbs.Batwa.GBR.pvals.BP.GO.overrep.pdf"

go.res = read.table(in.file, header=TRUE, sep="\t", quote="~")

p = ggplot(go.res, aes(Expected, Significant)) +
    geom_abline(slope=1, lty=2, color="darkgrey") +
    geom_point(color="lightgrey") +
    geom_point(data=subset(go.res, classicFisher < 0.1),
        aes(Expected, Significant, color=classicFisher)) +
    geom_segment(data=subset(go.res, go.res$classicFisher == min(go.res$classicFisher)),
        aes(x=Expected, y=Significant,
            xend=(Expected + 10), yend=(Significant + 30))) +
    geom_label(data=subset(go.res, go.res$classicFisher == min(go.res$classicFisher)),
        aes(Expected, Significant, label=GO.ID, fill=classicFisher),
        col='white', size=3, nudge_y=30, nudge_x=10) +
    xlab("Expected number of genes") +
    ylab("Observed number of genes") +
    theme_bw() +
    theme(legend.position="none",
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    xlim(c(0,100)) + ylim(c(0,100))

ggsave(p, file=out.pdf, height=5, width=5)
