#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Plot map
# ----------------------------------------------------------------------------------------

library(globe)

pdf("reports/sampling_map.pdf")

globeearth(eye=list(lat=20, lon=54), top=place("northpole"))

globepoints(loc=list(lat=-1.0806, lon=29.6614), col='red', bg='black', pch=21)
globepoints(loc=list(lat=11.7401, lon=92.6586), col='red', bg='black', pch=21)
globepoints(loc=list(lat=27,      lon=81),      col='red', bg='black', pch=21)

globepoints(loc=list(lat=-1.0806, lon=29.6614), col='black', cex=2)
globepoints(loc=list(lat=11.7401, lon=92.6586), col='black', cex=2)
globepoints(loc=list(lat=27,      lon=81),      col='black', cex=2)

dev.off()
