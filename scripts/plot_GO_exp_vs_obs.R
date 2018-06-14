#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Plot Anda and Batwa GO against one another
# ----------------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)

library(ggplot2)
library(ggrepel)

library(GO.db)

args = commandArgs(trailingOnly=TRUE)
in.pop1     = args[1]   # E.g. "results/pbs.Batwa.GBR.pvals.BP.overrep.results.txt"
in.pop2     = args[2]   # E.g. "results/pbs.Anda.LWK.pvals.BP.overrep.results.txt"
joint.go    = args[3]   # E.g. "results/convergent_permuted_pval_Batwa.GBR-Anda.LWK.BP.overrep.txt"
out.pdf     = args[4]   # E.g. "reports/pbs.Batwa-Anda.pvals.BP.GO.overrep.pdf"

# BP:
# in.pop1     = "results/pbs.Batwa.GBR.pvals.BP.overrep.results.txt"
# in.pop2     = "results/pbs.Anda.LWK.pvals.BP.overrep.results.txt"
# joint.go    = "results/convergent_permuted_pval_Batwa.GBR-Anda.LWK.BP.overrep.txt"
# out.pdf     = "reports/pbs.Batwa-Anda.pvals.BP.GO.overrep.pdf"

# MF:
# in.pop1     = "results/pbs.Batwa.GBR.pvals.MF.overrep.results.txt"
# in.pop2     = "results/pbs.Anda.LWK.pvals.MF.overrep.results.txt"
# joint.go    = "results/convergent_permuted_pval_Batwa.GBR-Anda.LWK.MF.overrep.txt"
# out.pdf     = "reports/pbs.Batwa-Anda.pvals.MF.GO.overrep.pdf"

# BP and MF (RHG):
# in.pop1     = "results/pbs.Batwa.GBR.pvals.BPMF.overrep.results.txt"
# in.pop2     = "results/pbs.Anda.LWK.pvals.BPMF.overrep.results.txt"
# joint.go    = "results/convergent_permuted_pval_Batwa.GBR-Anda.LWK.BPMF.overrep.txt"
# out.pdf     = "reports/pbs.Batwa-Anda.pvals.BPMF.GO.overrep.pdf"

# BP and MF (AGR):
# in.pop1     = "results/pbs.Bakiga.GBR.pvals.BPMF.overrep.results.txt"
# in.pop2     = "results/pbs.BR1.LWK.pvals.BPMF.overrep.results.txt"
# joint.go    = "results/convergent_permuted_pval_Bakiga.GBR-BR1.LWK.BPMF.overrep.txt"
# out.pdf     = "reports/pbs.Bakiga-BR1.pvals.BPMF.GO.overrep.pdf"

# ========================================================================================

pal.orange = c("#EC7E29", "#FFB57C", "#FC9C53", "#C45C0D", "#9B4300")
pal.blue   = c("#26579B", "#658BBF", "#406BA6", "#134181", "#093066")
pal.green  = c("#1FB539", "#67D379", "#3FC156", "#0A9622", "#007715")

# ========================================================================================

go.res.pop1 = read.table(in.pop1, header=TRUE, sep="\t", quote="~")
go.res.pop2 = read.table(in.pop2, header=TRUE, sep="\t", quote="~")

go.res = merge(go.res.pop1, go.res.pop2, by="GO.ID", suffixes=c(".Pop1", ".Pop2"))

go.res$log.ratio.Pop1 = log(go.res$Significant.Pop1 / go.res$Expected.Pop1, base=2)
go.res$log.ratio.Pop2 = log(go.res$Significant.Pop2 / go.res$Expected.Pop2, base=2)

# Compute log values for CIs
if (grepl("NvsS", out.pdf) == FALSE) {
    go.res$log.ratio.Pop1.CI95low = log(go.res$CI.95.low.Pop1, base=2)
    go.res$log.ratio.Pop1.CI95up  = log(go.res$CI.95.up.Pop1,  base=2)
    go.res$log.ratio.Pop2.CI95low = log(go.res$CI.95.low.Pop2, base=2)
    go.res$log.ratio.Pop2.CI95up  = log(go.res$CI.95.up.Pop2,  base=2)
}

joint = read.table(joint.go, header=TRUE)
joint = joint[,c("GO.ID", "joint.bs.p.genePBS")]
go.res.jt = merge(go.res, joint, by="GO.ID", suffixes = c("", ".Joint"), all.x=TRUE)

# Remove the biggest categories
go.res.jt = go.res.jt[!(go.res.jt$Term.Pop1 %in%
    c("biological_process", "molecular_function", "cellular_component")),]

go.res.jt$label = ""

if (grepl("Batwa", out.pdf)) {
    if (grepl("pvals\\..*BP.*\\.GO", out.pdf)) {
        # Batwa high
        go.res.jt[go.res.jt$GO.ID == "GO:0007517",]$label = "Muscle organ development"
        go.res.jt[go.res.jt$GO.ID == "GO:0045926",]$label = "Negative regulation\nof growth"
        # Anda high
        go.res.jt[go.res.jt$GO.ID == "GO:0006302",]$label = "Double-strand break repair"
        #go.res.jt[go.res.jt$GO.ID == "GO:0070085",]$label = "Glycosylation"
        # (Others exist that are more significant)
        go.res.jt[go.res.jt$GO.ID == "GO:0045596",]$label = "Negative regulation of\ncell differentiation"
        # Both high
        #go.res.jt[go.res.jt$GO.ID == "GO:0048705",]$label = "Skeletal system morphogenesis"
        #go.res.jt[go.res.jt$GO.ID == "GO:1901617",]$label = "Organic hydroxy compound\nbiosynthetic process"
        #go.res.jt[go.res.jt$GO.ID == "GO:0007034",]$label = "Vacuolar transport"
        go.res.jt[go.res.jt$GO.ID == "GO:0035108",]$label = "Limb morphogenesis"
        go.res.jt[go.res.jt$GO.ID == "GO:0030326",]$label = "Embryonic limb\nmorphogenesis"
    }

    if (grepl("pvals\\..*MF.*\\.GO", out.pdf)) {
        # Batwa high
        #go.res.jt[go.res.jt$GO.ID == "GO:0003723",]$label = "RNA binding"
        #go.res.jt[go.res.jt$GO.ID == "GO:0043167",]$label = "Ion binding"
        # Anda high
        go.res.jt[go.res.jt$GO.ID == "GO:0043130",]$label = "Ubiquitin binding"
        #go.res.jt[go.res.jt$GO.ID == "GO:0008233",]$label = "Peptidase activity"
        # Both high
        # (None)
    }
} else if (grepl("Bakiga", out.pdf)) {
    if (grepl("pvals\\..*BP.*\\.GO", out.pdf)) {
        # Bakiga high
        go.res.jt[go.res.jt$GO.ID == "GO:0002283",]$label = "Neutrophil activation\ninvolved in immune response"
        #go.res.jt[go.res.jt$GO.ID == "GO:0002446",]$label = "Neutrophil mediated\nimmunity"
        #go.res.jt[go.res.jt$GO.ID == "GO:0042119",]$label = "Neutrophil activation"
        #go.res.jt[go.res.jt$GO.ID == "GO:0036230",]$label = "Granulocyte activation"
        #go.res.jt[go.res.jt$GO.ID == "GO:0002446",]$label = "Neutrophil mediated immunity"
        #go.res.jt[go.res.jt$GO.ID == "GO:0090407",]$label = "organophosphate biosynthetic process"
        #go.res.jt[go.res.jt$GO.ID == "GO:0043299",]$label = "leukocyte degranulation"
        #go.res.jt[go.res.jt$GO.ID == "GO:0002275",]$label = "myeloid cell activation involved in immunity"
        #go.res.jt[go.res.jt$GO.ID == "GO:0002444",]$label = "myeloid leukocyte mediated immunity"
        #go.res.jt[go.res.jt$GO.ID == "GO:0045055",]$label = "regulated exocytosis"
        # BR1 high
        go.res.jt[go.res.jt$GO.ID == "GO:0046777",]$label = "Protein autophosphorylation"
        #go.res.jt[go.res.jt$GO.ID == "GO:0006260",]$label = "DNA replication"
        #go.res.jt[go.res.jt$GO.ID == "GO:0031323",]$label = "Regulation of cellular metabolic process"
        #go.res.jt[go.res.jt$GO.ID == "GO:0050790",]$label = "regulation of catalytic activity"
        #go.res.jt[go.res.jt$GO.ID == "GO:0007169",]$label = "transmembrane receptor protein tyrosine ..."
        #go.res.jt[go.res.jt$GO.ID == "GO:0002682",]$label = "regulation of immune system process"
        #go.res.jt[go.res.jt$GO.ID == "GO:0002761",]$label = "regulation of myeloid leukocyte differen..."
        #go.res.jt[go.res.jt$GO.ID == "GO:0060249",]$label = "anatomical structure homeostasis"
        #go.res.jt[go.res.jt$GO.ID == "GO:0014066",]$label = "regulation of phosphatidylinositol 3-kin..."
        # Both high
        #go.res.jt[go.res.jt$GO.ID == "GO:0002761",]$label = "Regulation of\nmyeloid leukocyte differentiation"
    }

    if (grepl("pvals\\..*MF.*\\.GO", out.pdf)) {
        # Bakiga high
        #go.res.jt[go.res.jt$GO.ID == "GO:0004866",]$label = "Endopeptidase\ninhibitor activity"
        #go.res.jt[go.res.jt$GO.ID == "GO:0032403",]$label = "Protein complex binding"
        # BR1 high
        go.res.jt[go.res.jt$GO.ID == "GO:0031406",]$label = "Carboxylic acid binding"
        #go.res.jt[go.res.jt$GO.ID == "GO:0008083",]$label = "Growth factor activity"
        # Both high
        #go.res.jt[go.res.jt$GO.ID == "GO:0070851",]$label = "Growth factor receptor binding"
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

go.res.jt$on.axis = 21
go.res.jt[go.res.jt$log.ratio.Pop2 == -Inf,]$on.axis = 124
go.res.jt[go.res.jt$log.ratio.Pop1 == -Inf,]$on.axis = 95
go.res.jt$on.axis = factor(go.res.jt$on.axis, levels = c(21,124,95,1))

# ----------------------------------------------------------------------------------------

# Bring in labels which have been manually edited
label.data = read.table("data/Batwa-Andamanese_convergence_figures_label_positions.txt",
    sep="\t", header=TRUE)
label.data = label.data[grepl("overrep", label.data$dataset),]

go.res.jt = merge(go.res.jt, label.data[,c(1,5:6)], all=TRUE)

# ----------------------------------------------------------------------------------------

p = ggplot(subset(go.res.jt, Annotated.Pop1 > 10 & Annotated.Pop2 > 10),
        aes(log.ratio.Pop1, log.ratio.Pop2, size=Annotated.Pop1)) +
    geom_vline(xintercept=0, lty=2, color="darkgrey") +
    geom_hline(yintercept=0, lty=2, color="darkgrey")

# Add segments for labels
p = p + geom_segment(data=subset(go.res.jt, to.label == TRUE),
    aes(x    = log.ratio.Pop1,  y    = log.ratio.Pop2,
        xend = label.loc.x,     yend = label.loc.y),
    color="black", lwd=0.3, lty=1)

if (grepl("NvsS", out.pdf) == FALSE) {

    data.to.add = subset(go.res.jt,
        (classicFisher.Pop1 < 0.01 | classicFisher.Pop2 < 0.01) &
            log.ratio.Pop1 > -4)
    if (nrow(data.to.add) > 0) {
        p = p + geom_segment(data=data.to.add,
            aes(x    = log.ratio.Pop1.CI95low, y    = log.ratio.Pop2,
                xend = log.ratio.Pop1.CI95up,  yend = log.ratio.Pop2),
            color=pal.blue[1], lwd=0.1, lty=1)
    }

    data.to.add = subset(go.res.jt,
        (classicFisher.Pop1 < 0.01 | classicFisher.Pop2 < 0.01) &
            log.ratio.Pop2 > -4)

    if (nrow(data.to.add) > 0) {
        p = p + geom_segment(data=data.to.add,
            aes(x    = log.ratio.Pop1, y    = log.ratio.Pop2.CI95low,
                xend = log.ratio.Pop1, yend = log.ratio.Pop2.CI95up),
            color=pal.blue[1], lwd=0.1, lty=1)
    }

    if (nrow(subset(go.res.jt, joint.bs.p.genePBS <= 0.001)) > 0) {
        p = p + geom_segment(data=subset(go.res.jt, joint.bs.p.genePBS <= 0.001 &
                log.ratio.Pop1 > -4),
            aes(x    = log.ratio.Pop1.CI95low, y    = log.ratio.Pop2,
                xend = log.ratio.Pop1.CI95up,  yend = log.ratio.Pop2),
            color=pal.orange[1], lwd=0.1, lty=1) +
            geom_segment(data=subset(go.res.jt, joint.bs.p.genePBS <= 0.001 &
                    log.ratio.Pop2 > -4),
                aes(x    = log.ratio.Pop1, y    = log.ratio.Pop2.CI95low,
                    xend = log.ratio.Pop1, yend = log.ratio.Pop2.CI95up),
                color=pal.orange[1], lwd=0.1, lty=1)
    }
}

p = p + geom_point(aes(shape=on.axis), color="darkgrey", fill="grey", alpha=0.5)

if (nrow(subset(go.res.jt,
            classicFisher.Pop1 < 0.01 | classicFisher.Pop2 < 0.01)) > 0) {
    p = p + geom_point(data=subset(go.res.jt,
            classicFisher.Pop1 < 0.01 | classicFisher.Pop2 < 0.01),
        aes(log.ratio.Pop1, log.ratio.Pop2,
            fill="fillblue", color="colblue"), pch=21)
}
if (nrow(subset(go.res.jt, joint.bs.p.genePBS <= 0.001)) > 0) {
    p = p + geom_point(data=subset(go.res.jt, joint.bs.p.genePBS <= 0.001),
        aes(log.ratio.Pop1, log.ratio.Pop2,
            fill="fillorange", color="colorange"), pch=21)
} else {
    # Add fake point off graph so it shows up in legend
    fake.df = go.res.jt[1,]
    fake.df[1,]$log.ratio.Pop1 = 10
    fake.df[1,]$log.ratio.Pop2 = 10
    p = p + geom_point(data=fake.df,
        aes(log.ratio.Pop1, log.ratio.Pop2,
            fill="fillorange", color="colorange"), pch=21)
}
p = p + geom_point(data=subset(go.res.jt, growth.go == TRUE),
        aes(log.ratio.Pop1, log.ratio.Pop2,
            shape="1",
            color='colblack'), size=2) +
    # geom_label_repel(data=subset(go.res.jt, to.label == TRUE),
    #         aes(log.ratio.Pop1, log.ratio.Pop2, label=label),
    #     col='black', fill='white', size=4, force=250,
    #     segment.color = 'black',
    #     nudge_x = ifelse(subset(go.res.jt, to.label == TRUE)$log.ratio.Pop1 > 0,
    #         1, 0),
    #     nudge_y = ifelse(subset(go.res.jt, to.label == TRUE)$log.ratio.Pop2 > 0,
    #         1, 0)) +
    theme_bw() +
    theme(legend.justification=c(0,0),
          legend.position=c(0.01,0.01),
          legend.background=element_rect(fill = NA),
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

if (grepl("NvsS", out.pdf) == FALSE) {
    p = p +
        xlab(paste0("Significant Genes : Expected Genes (log2)\n", pop1.name)) +
        ylab(paste0("Significant Genes : Expected Genes (log2)\n", pop2.name))
} else {
    p = p +
        xlab(paste0("Significant Nonsyn. SNPs : Expected Nonsyn. SNPs (log2)\n",
            pop1.name)) +
        ylab(paste0("Significant Nonsyn. SNPs : Expected Nonsyn. SNPs (log2)\n",
            pop2.name))
}

# Extra space for N-vs-S labels
if (grepl("NvsS", out.pdf) == FALSE) {
    p = p + xlim(c(-3.5,3.5)) + ylim(c(-3.5,3.5))
} else {
    p = p + xlim(c(-4,5)) + ylim(c(-4,5))
}

p = p + coord_fixed() +
    scale_fill_manual(    name="Significant\nEnrichment",
                          values = c("fillblue"=pal.blue[2],
                                     "fillorange"=pal.orange[2],
                                     "fillblack"="black"),
                          breaks=c("fillblue", "fillorange"),
                          labels=c("Population-specific (p < 0.01)",
                                   "Convergent (p <= 0.001)"),
                          drop=FALSE) +
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
# tmp = subset(go.res.jt, to.label == TRUE)[,c("GO.ID", "Term.Pop1", "log.ratio.Pop1", "log.ratio.Pop2")]
# write.table(tmp, file="tmp.txt", sep="\t", quote=FALSE, row.names=FALSE)
# built = ggplot_gtable(ggplot_build(p))
# built$grobs[[6]]$children$geom_label_repel.labelrepeltree.16$data
