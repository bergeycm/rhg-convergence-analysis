#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Plot PBS schematic
# ----------------------------------------------------------------------------------------

library("ape")

# ----------------------------------------------------------------------------------------
# --- Get branch lengths for outliers and average
# ----------------------------------------------------------------------------------------

pbs.gbr = read.table("results/GBR.FST_for_PBS.txt.anno.bed")
pbs.lwk = read.table("results/LWK.FST_for_PBS.txt.anno.bed")

# Bakiga, Batwa, GBR
names(pbs.gbr) = c("chr", "start", "end",
                   "Fst_focal_sister", "focal_sister_pop1", "focal_sister_pop2",
                   "Fst_sister_out", "sister_out_pop1", "sister_out_pop2",
                   "Fst_focal_out", "focal_out_pop1", "focal_out_pop2",
                   "T_sister_focal", "T_sister_out", "T_focal_out",
                   "PBS_focal", "PBS_sister", "genes")

# Alphabetical order is focal, sister, out, so column order is a bit different
# Anda, BR1, LWK
names(pbs.lwk) = c("chr", "start", "end",
                   "Fst_focal_sister", "focal_sister_pop1", "focal_sister_pop2",
                   "Fst_focal_out", "focal_out_pop1", "focal_out_pop2",
                   "Fst_sister_out", "sister_out_pop1", "sister_out_pop2",
                   "T_sister_focal", "T_focal_out", "T_sister_out",
                   "PBS_focal", "PBS_sister", "genes")

# Compute PBS for outgroup

pbs.gbr$PBS_out = (pbs.gbr$T_focal_out + pbs.gbr$T_sister_out -
                   pbs.gbr$T_sister_focal) / 2
pbs.lwk$PBS_out = (pbs.lwk$T_focal_out + pbs.lwk$T_sister_out -
                   pbs.lwk$T_sister_focal) / 2

pbs.lwk = pbs.lwk[abs(pbs.lwk$PBS_focal) + abs(pbs.lwk$PBS_sister) + abs(pbs.lwk$PBS_out) != Inf,]

# --- Compute average PBS

# GBR
avg.T.gbr.focal_sister = -log(1 - mean(pbs.gbr$Fst_focal_sister))
avg.T.gbr.sister_out   = -log(1 - mean(pbs.gbr$Fst_sister_out))
avg.T.gbr.focal_out    = -log(1 - mean(pbs.gbr$Fst_focal_out))

avg.PBS.gbr.focal = (avg.T.gbr.focal_sister + avg.T.gbr.focal_out - avg.T.gbr.sister_out) / 2
avg.PBS.gbr.sister = (avg.T.gbr.focal_sister + avg.T.gbr.sister_out - avg.T.gbr.focal_out) / 2
avg.PBS.gbr.out = (avg.T.gbr.focal_out + avg.T.gbr.sister_out - avg.T.gbr.focal_sister) / 2

if (avg.PBS.gbr.focal < 0)  avg.PBS.gbr.focal  = 0
if (avg.PBS.gbr.sister < 0) avg.PBS.gbr.sister = 0
if (avg.PBS.gbr.out < 0)    avg.PBS.gbr.out    = 0

# LWK
avg.T.lwk.focal_sister = -log(1 - mean(pbs.lwk$Fst_focal_sister))
avg.T.lwk.sister_out   = -log(1 - mean(pbs.lwk$Fst_sister_out))
avg.T.lwk.focal_out    = -log(1 - mean(pbs.lwk$Fst_focal_out))

avg.PBS.lwk.focal = (avg.T.lwk.focal_sister + avg.T.lwk.focal_out - avg.T.lwk.sister_out) / 2
avg.PBS.lwk.sister = (avg.T.lwk.focal_sister + avg.T.lwk.sister_out - avg.T.lwk.focal_out) / 2
avg.PBS.lwk.out = (avg.T.lwk.focal_out + avg.T.lwk.sister_out - avg.T.lwk.focal_sister) / 2

if (avg.PBS.lwk.focal < 0)  avg.PBS.lwk.focal  = 0
if (avg.PBS.lwk.sister < 0) avg.PBS.lwk.sister = 0
if (avg.PBS.lwk.out < 0)    avg.PBS.lwk.out    = 0

# --- Grab outliers to show

out.PBS.gbr = pbs.gbr[order(pbs.gbr$PBS_focal, decreasing=TRUE),][1:2,c(16:17,19)]
out.PBS.gbr$PBS_sister[out.PBS.gbr$PBS_sister < 0] = 0

out.PBS.lwk = pbs.lwk[order(pbs.lwk$PBS_focal, decreasing=TRUE),][1:2,c(16:17,19)]
out.PBS.lwk$PBS_sister[out.PBS.lwk$PBS_sister < 0] = 0

# ----------------------------------------------------------------------------------------
# --- Make trees for PBS values
# ----------------------------------------------------------------------------------------

pbs.afr.avg.tre  = read.tree(text = paste0("((Batwa:", avg.PBS.gbr.focal, ",",
                                           "Bakiga:",  avg.PBS.gbr.sister, "):0,",
                                           "British:", avg.PBS.gbr.out, ");"))
pbs.afr.out1.tre = read.tree(text = paste0("((Batwa:", out.PBS.gbr$PBS_focal[1], ",",
                                           "Bakiga:",  out.PBS.gbr$PBS_sister[1], "):0,",
                                           "British:", out.PBS.gbr$PBS_out[1], ");"))
pbs.afr.out2.tre = read.tree(text = paste0("((Batwa:", out.PBS.gbr$PBS_focal[2], ",",
                                           "Bakiga:",  out.PBS.gbr$PBS_sister[2], "):0,",
                                           "British:", out.PBS.gbr$PBS_out[2], ");"))

pbs.asn.avg.tre  = read.tree(text = paste0("((Andamanese:", avg.PBS.lwk.focal, ",",
                                           "Brahmin:",      avg.PBS.lwk.sister, "):0,",
                                           "Kenyan:",       avg.PBS.lwk.out, ");"))
pbs.asn.out1.tre = read.tree(text = paste0("((Andamanese:", out.PBS.lwk$PBS_focal[1], ",",
                                           "Brahmin:",      out.PBS.lwk$PBS_sister[1], "):0,",
                                           "Kenyan:",       out.PBS.lwk$PBS_out[1], ");"))
pbs.asn.out2.tre = read.tree(text = paste0("((Andamanese:", out.PBS.lwk$PBS_focal[2], ",",
                                           "Brahmin:",      out.PBS.lwk$PBS_sister[2], "):0,",
                                           "Kenyan:",       out.PBS.lwk$PBS_out[2], ");"))

# ----------------------------------------------------------------------------------------
# --- Plot PBS schematic
# ----------------------------------------------------------------------------------------

pdf("reports/pbs_schematic.pdf", height=6, width=8)

par(mfrow=c(2,3))
par(mar=c(1,1,1,1))

# African Average
plot(pbs.afr.avg.tre, type="unrooted", label.offset=0.01,
    rotate.tree=180, no.margin=FALSE, cex=1, edge.width=4,
    main="Mean PBS")
# African outlier 1
plot(pbs.afr.out1.tre, type="unrooted", label.offset=0.01,
    rotate.tree=180, no.margin=FALSE, cex=1,
    main=paste("Outlier in", gsub(",.*", "",
        pbs.gbr[order(pbs.gbr$PBS_focal, decreasing=TRUE),][1,]$genes)))
# African outlier 2
plot(pbs.afr.out2.tre, type="unrooted", label.offset=0.01,
    rotate.tree=180, no.margin=FALSE, cex=1,
    main=paste("Outlier in", gsub(",.*", "",
        pbs.gbr[order(pbs.gbr$PBS_focal, decreasing=TRUE),][2,]$genes)))

# Asian Average
plot(pbs.asn.avg.tre, type="unrooted", label.offset=0.01,
    rotate.tree=180, no.margin=FALSE, cex=1, edge.width=4,
    main="Mean PBS")
# Asian outlier 1
plot(pbs.asn.out1.tre, type="unrooted", label.offset=0.01,
    rotate.tree=180, no.margin=FALSE, cex=1,
    main=paste("Outlier in", gsub(",.*", "",
        pbs.lwk[order(pbs.lwk$PBS_focal, decreasing=TRUE),][1,]$genes)))
# Asian outlier 2
plot(pbs.asn.out2.tre, type="unrooted", label.offset=0.01,
    rotate.tree=180, no.margin=FALSE, cex=1,
    main=paste("Outlier in", gsub(",.*", "",
        pbs.lwk[order(pbs.lwk$PBS_focal, decreasing=TRUE),][2,]$genes)))

dev.off()
