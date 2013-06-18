#!/usr/bin/Rscript

# ---------- SET UP OPTIONS ----------
options(width=60)
library(GOSemSim)
library(org.Sc.sgd.db)
library(GO.db)
setwd('/Users/kdoroschak/Documents/DREU/GOSemSim')

# ---------- OPEN FILES ----------
biogrid = read.delim("biogrid_ppi_nom-3.ppi", header=FALSE)
numPPIs = nrow(biogrid)
biogrid <- data.frame(lapply(biogrid, as.character), stringsAsFactors=FALSE)
biogrid["GOSemSimScore"] <- NA


# ---------- COMPUTE SIMILARITY SCORES ----------
matrixcount = 0
processedcount = 0
for (i in 1:numPPIs) {
  score = data.frame(geneSim=0.0)
  gene1 = biogrid[i,1]
  gene2 = biogrid[i,2]
  biogrid[i,3] = 0.0
  try((score = geneSim(gene1, gene2, ont = "MF", organism = "yeast",
          measure = "Wang", drop = "NULL", combine = "avg")), silent = TRUE)
  if (!is.matrix(score$geneSim)) {
    #print(score$geneSim)
    biogrid[i,3] = score$geneSim
  }
  else {
    matrixcount = matrixcount + 1
  }
  processedcount = processedcount + 1
}
print(matrixcount)
print(processedcount)
print(numPPIs)

# ---------- WRITE TO FILE ---------
write.table(biogrid, file = "biogrid_ppi_nom-3.gosemsim.prob", 
            sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
