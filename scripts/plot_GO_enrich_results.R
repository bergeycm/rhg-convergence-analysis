#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Plot Anda and Batwa GO enrichment results
# ----------------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)

library(ggplot2)
library(ggrepel)

library(GO.db)

args = commandArgs(trailingOnly=TRUE)
in.pop1     = args[1]   # E.g. "results/pbs.Batwa.GBR.pvals.BP.enrich.results.txt"
in.pop2     = args[2]   # E.g. "results/pbs.Anda.LWK.pvals.BP.enrich.results.txt"
joint.go    = args[3]   # E.g. "results/convergent_permuted_pval_Batwa.GBR-Anda.LWK.BP.enrich.txt"
out.pdf     = args[4]   # E.g. "reports/pbs.Batwa-Anda.pvals.BP.GO.enrich.pdf"

# BP:
# in.pop1     = "results/pbs.Batwa.GBR.pvals.BP.enrich.results.txt"
# in.pop2     = "results/pbs.Anda.LWK.pvals.BP.enrich.results.txt"
# joint.go    = "results/convergent_permuted_pval_Batwa.GBR-Anda.LWK.BP.enrich.txt"
# out.pdf     = "reports/pbs.Batwa-Anda.pvals.BP.GO.enrich.pdf"

# MF:
# in.pop1     = "results/pbs.Batwa.GBR.pvals.MF.enrich.results.txt"
# in.pop2     = "results/pbs.Anda.LWK.pvals.MF.enrich.results.txt"
# joint.go    = "results/convergent_permuted_pval_Batwa.GBR-Anda.LWK.MF.enrich.txt"
# out.pdf     = "reports/pbs.Batwa-Anda.pvals.MF.GO.enrich.pdf"

# BP and MF (RHG):
# in.pop1     = "results/pbs.Batwa.GBR.pvals.BPMF.enrich.results.txt"
# in.pop2     = "results/pbs.Anda.LWK.pvals.BPMF.enrich.results.txt"
# joint.go    = "results/convergent_permuted_pval_Batwa.GBR-Anda.LWK.BPMF.enrich.txt"
# out.pdf     = "reports/pbs.Batwa-Anda.pvals.BPMF.GO.enrich.pdf"

# BP and MF (AGR):
# in.pop1     = "results/pbs.Bakiga.GBR.pvals.BPMF.enrich.results.txt"
# in.pop2     = "results/pbs.BR1.LWK.pvals.BPMF.enrich.results.txt"
# joint.go    = "results/convergent_permuted_pval_Bakiga.GBR-BR1.LWK.BPMF.enrich.txt"
# out.pdf     = "reports/pbs.Bakiga-BR1.pvals.BPMF.GO.enrich.pdf"

# ========================================================================================

pal.orange = c("#EC7E29", "#FFB57C", "#FC9C53", "#C45C0D", "#9B4300")
pal.blue   = c("#26579B", "#658BBF", "#406BA6", "#134181", "#093066")
pal.green  = c("#1FB539", "#67D379", "#3FC156", "#0A9622", "#007715")

# ========================================================================================

go.res.pop1 = read.table(in.pop1, header=TRUE, sep="\t", quote="~")
go.res.pop2 = read.table(in.pop2,  header=TRUE, sep="\t", quote="~")

go.res = merge(go.res.pop1, go.res.pop2, by="GO.ID", suffixes=c(".Pop1", ".Pop2"))

go.res$log.p.Pop1 = -log(go.res$classicKS.Pop1, base=10)
go.res$log.p.Pop2 = -log(go.res$classicKS.Pop2, base=10)

# Compute log values for SEs
go.res$log.p.Pop1.SElow = -log(go.res$classicKS.Pop1 - go.res$std.err.Pop1, base=10)
go.res$log.p.Pop1.SEup  = -log(go.res$classicKS.Pop1 + go.res$std.err.Pop1, base=10)
go.res$log.p.Pop2.SElow = -log(go.res$classicKS.Pop2 - go.res$std.err.Pop2, base=10)
go.res$log.p.Pop2.SEup  = -log(go.res$classicKS.Pop2 + go.res$std.err.Pop2, base=10)

joint = read.table(joint.go, header=TRUE)
joint = joint[,c("GO.ID", "joint.bs.p.genePBS")]
go.res.jt = merge(go.res, joint, by="GO.ID", suffixes = c("", ".Joint"), all.x=TRUE)

# Remove the biggest categories
go.res.jt = go.res.jt[!(go.res.jt$Term.Pop1 %in%
    c("biological_process", "molecular_function", "cellular_component")),]

# awk 'BEGIN {FS="\t"}{ if ($7 < 0.01) print $0}' \
#     results/pbs.Batwa.GBR.pvals.BP.enrich.results.txt

go.res.jt$label = ""
interesting.go = c()
if (grepl("Batwa", out.pdf)) {
    if (grepl("pvals\\..*BP.*\\.GO", out.pdf)) {
        # Both high
        go.res.jt[go.res.jt$GO.ID == "GO:0035265",]$label = "Organ growth"
        go.res.jt[go.res.jt$GO.ID == "GO:1903828",]$label = "Negative regulation of\ncellular protein localization"
        go.res.jt[go.res.jt$GO.ID == "GO:0048738",]$label = "Cardiac muscle\ntissue development"
        # (Not significant:)
        # Same as below: go.res.jt[go.res.jt$GO.ID == "GO:0016202",]$label = "Regulation of striated\nmuscle tissue development"
        #go.res.jt[go.res.jt$GO.ID == "GO:1901861",]$label = "Regulation of\nmuscle tissue development"
        go.res.jt[go.res.jt$GO.ID == "GO:0035265",]$label = "Organ growth"
        #go.res.jt[go.res.jt$GO.ID == "GO:0048634",]$label = "Regulation of\nmuscle organ development"
        # Batwa high
        go.res.jt[go.res.jt$GO.ID == "GO:0003231",]$label = "Cardiac ventricle\ndevelopment"
        # Anda high
        go.res.jt[go.res.jt$GO.ID == "GO:0035051",]$label = "Cardiocyte differentiation"
    }

    if (grepl("pvals\\..*MF.*\\.GO", out.pdf)) {
        # Both high
        go.res.jt[go.res.jt$GO.ID == "GO:0019838",]$label = "Growth factor\nbinding"
        # Batwa high
        # (None)
        # Anda high
        #go.res.jt[go.res.jt$GO.ID == "GO:0005085",]$label = "Guanyl-nucleotide\nexchange factor\nactivity"
    }

} else if (grepl("Bakiga", out.pdf)) {
    if (grepl("pvals\\..*BP.*\\.GO", out.pdf)) {
        # Both high
        go.res.jt[go.res.jt$GO.ID == "GO:0002521",]$label = "Leukocyte\ndifferentiation"
        go.res.jt[go.res.jt$GO.ID == "GO:0046777",]$label = "Protein\nautophosphorylation"
        # Bakiga high
        go.res.jt[go.res.jt$GO.ID == "GO:0006665",]$label = "Sphingolipid\nmetabolic process"
        go.res.jt[go.res.jt$GO.ID == "GO:0031532",]$label = "Actin cytoskeleton\nreorganization"
        go.res.jt[go.res.jt$GO.ID == "GO:0006952",]$label = "Defense response"
        go.res.jt[go.res.jt$GO.ID == "GO:0006909",]$label = "Phagocytosis"
        # BR1 high
        go.res.jt[go.res.jt$GO.ID == "GO:0044242",]$label = "Cellular lipid\ncatabolic process"
        go.res.jt[go.res.jt$GO.ID == "GO:0001654",]$label = "Eye development"
    }

    if (grepl("pvals\\..*MF.*\\.GO", out.pdf)) {
        # Both high
        #go.res.jt[go.res.jt$GO.ID == "GO:0008236",]$label = "Serine-type\npeptidase activity"
        # Bakiga high
        go.res.jt[go.res.jt$GO.ID == "GO:0003779",]$label = "Actin binding"
        # BR1 high
        go.res.jt[go.res.jt$GO.ID == "GO:0031406",]$label = "Carboxylic acid binding"
        go.res.jt[go.res.jt$GO.ID == "GO:0019903",]$label = "Protein phosphatase binding"
        go.res.jt[go.res.jt$GO.ID == "GO:0017171",]$label = "Serine hydrolase activity"
    }
}

interesting.go = go.res.jt[go.res.jt$label != "",]$GO.ID

go.res.jt$to.label = FALSE
go.res.jt[go.res.jt$GO.ID %in% interesting.go,]$to.label = TRUE

# ----------------------------------------------------------------------------------------

# --- Find growth genes

growth.go = NA
growth.go.of.interest = c()

if (grepl("BP", out.pdf)) {

    # 'growth'
    growth.go = "GO:0040007"

    growth.go.of.interest = c(growth.go.of.interest,
        growth.go, GOBPOFFSPRING[[growth.go]])

}

if (grepl("MF", out.pdf)) {

    # 'growth factor binding'
    growth.go.1 = "GO:0019838"
    # 'growth factor receptor binding'
    growth.go.2 = "GO:0070851"

    growth.go.descendants = c(GOMFOFFSPRING[[growth.go.1]], GOMFOFFSPRING[[growth.go.2]])

    # Add 'growth hormone receptor activity' and 'growth factor activity'
    growth.go.of.interest = c(growth.go.of.interest,
        unique(c(growth.go.1, growth.go.2, growth.go.descendants,
            "GO:0004903", "GO:0008083")))
}

go.res.jt$growth.go = FALSE
if (sum(go.res.jt$GO.ID %in% growth.go.of.interest) > 0) {
    go.res.jt[go.res.jt$GO.ID %in% growth.go.of.interest,]$growth.go = TRUE
}

# ----------------------------------------------------------------------------------------

# Bring in labels which have been manually edited
label.data = read.table("data/Batwa-Andamanese_convergence_figures_label_positions.txt",
    sep="\t", header=TRUE)
label.data = label.data[grepl("enrich", label.data$dataset),]

go.res.jt = merge(go.res.jt, label.data[,c(1,5:6)], all=TRUE)

# ----------------------------------------------------------------------------------------

p = ggplot(subset(go.res.jt, Annotated.Pop1 > 10 & Annotated.Pop2 > 10),
        aes(log.p.Pop1, log.p.Pop2, size=Annotated.Pop1)) +
    geom_vline(xintercept=0, lty=2, color="darkgrey") +
    geom_hline(yintercept=0, lty=2, color="darkgrey")

# Add segments for labels
p = p + geom_segment(data=subset(go.res.jt, to.label == TRUE),
    aes(x    = log.p.Pop1,   y    = log.p.Pop2,
        xend = label.loc.x,  yend = label.loc.y),
    color="black", lwd=0.3, lty=1)

# Add special segment for outlier
if (out.pdf == "reports/pbs.Bakiga-BR1.pvals.BPMF.GO.enrich.pdf") {
    df = data.frame(x1 = 0.5, x2 = 0.3909686, y1 = 3.3, y2 = 3.5)

    p = p + geom_segment(data=df,
        aes(x    = x1,  y    = y1,
            xend = x2,  yend = y2),
        color="black", lwd=0.3, lty=1,
        arrow = arrow(length = unit(0.01, "npc")))

} else if (out.pdf == "reports_sizeCor/pbs.Bakiga-BR1.pvals.BPMF.GO.enrich.pdf") {
    df = data.frame(x1 = 0.5, x2 = 0.3978113, y1 = 3.3, y2 = 3.5)

    p = p + geom_segment(data=df,
        aes(x    = x1,  y    = y1,
            xend = x2,  yend = y2),
        color="black", lwd=0.3, lty=1,
        arrow = arrow(length = unit(0.01, "npc")))
}

data.to.add = subset(go.res.jt,
    (classicKS.Pop1 < 0.01 | classicKS.Pop2 < 0.01))
if (nrow(data.to.add) > 0) {
    p = p + geom_segment(data=data.to.add,
        aes(x    = log.p.Pop1.SElow, y    = log.p.Pop2,
            xend = log.p.Pop1.SEup,  yend = log.p.Pop2),
        color=pal.blue[1], lwd=0.1, lty=1)
}

data.to.add = subset(go.res.jt,
    (classicKS.Pop1 < 0.01 | classicKS.Pop2 < 0.01))

if (nrow(data.to.add) > 0) {
    p = p + geom_segment(data=data.to.add,
        aes(x    = log.p.Pop1, y    = log.p.Pop2.SElow,
            xend = log.p.Pop1, yend = log.p.Pop2.SEup),
        color=pal.blue[1], lwd=0.1, lty=1)
}

if (nrow(subset(go.res.jt, joint.bs.p.genePBS <= 0.001)) > 0) {
    p = p + geom_segment(data=subset(go.res.jt, joint.bs.p.genePBS <= 0.001 &
            log.p.Pop1 > -4),
        aes(x    = log.p.Pop1.SElow, y    = log.p.Pop2,
            xend = log.p.Pop1.SEup,  yend = log.p.Pop2),
        color=pal.orange[1], lwd=0.1, lty=1) +
        geom_segment(data=subset(go.res.jt, joint.bs.p.genePBS <= 0.001 &
                log.p.Pop2 > -4),
            aes(x    = log.p.Pop1, y    = log.p.Pop2.SElow,
                xend = log.p.Pop1, yend = log.p.Pop2.SEup),
            color=pal.orange[1], lwd=0.1, lty=1)
}

p = p + geom_point(shape=16, color="darkgrey", fill="grey", alpha=0.5)

if (nrow(subset(go.res.jt,
            classicKS.Pop1 < 0.01 | classicKS.Pop2 < 0.01)) > 0) {
    p = p + geom_point(data=subset(go.res.jt,
            classicKS.Pop1 < 0.01 | classicKS.Pop2 < 0.01),
        aes(log.p.Pop1, log.p.Pop2,
            fill="fillblue", color="colblue"), pch=21)
}
if (nrow(subset(go.res.jt, joint.bs.p.genePBS <= 0.001)) > 0) {
    p = p + geom_point(data=subset(go.res.jt, joint.bs.p.genePBS <= 0.001),
        aes(log.p.Pop1, log.p.Pop2,
            fill="fillorange", color="colorange"), pch=21)
}
p = p + geom_point(data=subset(go.res.jt, growth.go == TRUE),
        aes(log.p.Pop1, log.p.Pop2,
            shape="1",
            color='colblack'), size=2) +
    # geom_label_repel(data=subset(go.res.jt, to.label == TRUE),
    #         aes(log.p.Pop1, log.p.Pop2, label=label),
    #     col='black', fill='white', size=4, force=100,
    #     segment.color = 'black',
    #     maxiter=1000,
    #     nudge_x = ifelse(subset(go.res.jt, to.label == TRUE)$log.p.Pop1 > 1,
    #         0.5, 0),
    #     nudge_y = ifelse(subset(go.res.jt, to.label == TRUE)$log.p.Pop2 > 1,
    #         0.5, 0)) +
    theme_bw() +
    theme(legend.justification=c(1,1),
          legend.position=c(0.995,0.99),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())

if (grepl("Batwa", in.pop1)) {
    pop1.name = "Batwa"
    pop2.name = "Andamanese"
} else if (grepl("Bakiga", in.pop1)) {
    pop1.name = "Bakiga"
    pop2.name = "Brahmin"
}

p = p +
    xlab(paste0("Enrichment p-value (-log10)\n", pop1.name)) +
    ylab(paste0("Enrichment p-value (-log10)\n", pop2.name))

p = p +
    scale_x_continuous(limits=c(0,3.5), expand = c(0, 0)) +
    scale_y_continuous(limits=c(0,3.5), expand = c(0, 0))

p = p + coord_fixed() +
    scale_fill_manual(    name="Significant Shifts",
                          values = c("fillblue"=pal.blue[2],
                                     "fillorange"=pal.orange[2],
                                     "fillblack"="black"),
                          breaks=c("fillblue", "fillorange"),
                          labels=c("Population-specific (p < 0.01)",
                                   "Convergent (p <= 0.001)")) +
    guides(fill = guide_legend(override.aes = list(size=5, color=NA))) +
    scale_color_manual(   guide=FALSE,
                          name="",
                          values = c("colblue"=pal.blue[3],
                                     "colorange"=pal.orange[5],
                                     "colblack"="black")) +
    scale_shape_manual(   guide=FALSE,
                          values=c(1,124,21,95)) +
    scale_size_continuous(name="Annotated Genes",
                          range = c(1, 10),
                          breaks=c(50,100,500,1000)) +
    guides(size = guide_legend(override.aes = list(shape=16,
                                                   color="darkgrey",
                                                   fill="grey",
                                                   alpha=0.5)))


p = p + geom_label(data=subset(go.res.jt, to.label == TRUE),
            aes(label.loc.x, label.loc.y, label=label),
        col='black', fill='white', size=4)

ggsave(p, file=out.pdf, height=8, width=8, useDingbats=FALSE)

# Adjust label placement
# tmp = subset(go.res.jt, to.label == TRUE)[,c("GO.ID", "Term.Pop1", "log.p.Pop1", "log.p.Pop2")]
# write.table(tmp, file="tmp.txt", sep="\t", quote=FALSE, row.names=FALSE)
# This just returns points as fraction of plot area:
# built = ggplot_gtable(ggplot_build(p))
# built$grobs[[6]]$children$geom_label_repel.labelrepeltree.16$data
