#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Find indivudals with less admixture for iHS analyses
# ----------------------------------------------------------------------------------------

# --- Bring in info from Perry et al 2014 to get admixture proportions
ind.info.file.2014 = paste0("data/Perry_etal_2014/",
    "Perry_etal_2014-individual_info-height_map.csv")

info.2014 = read.csv(ind.info.file.2014, comment.char="#")

# --- Bring in exome info file to know who was sequenced
ind.info.exomes = read.table("data/East_AGR_RHG.txt", sep="\t", header=TRUE)

# --- Merge
orig.case = ind.info.exomes$Indiv
ind.info.exomes$Indiv = tolower(ind.info.exomes$Indiv)

ind.info = merge(ind.info.exomes, info.2014, by.x="Indiv", by.y="Sample")

# --- Make pop. files, excluding individuals who are Batwa
# --- but have nearly "pure" Bakiga ancestry
info.batwa  = ind.info[which(ind.info$Admixture.Proportion.Batwa >= 0.9),]
info.bakiga = ind.info[which(ind.info$Admixture.Proportion.Batwa <= 0.1 &
                             ind.info$Comment != "Bakiga_in_PCA"),]

write.table(info.batwa$Indiv,  file="data/eRHG_IDs.iHS.txt",
    quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(info.bakiga$Indiv, file="data/eAGR_IDs.iHS.txt",
    quote=FALSE, row.names=FALSE, col.names=FALSE)

# Also write lists containing IDs used in VCF file of exome variants
write.table(info.batwa$ID,  file="data/eRHG_IDs.iHS.VCF.txt",
    quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(info.bakiga$ID, file="data/eAGR_IDs.iHS.VCF.txt",
    quote=FALSE, row.names=FALSE, col.names=FALSE)
