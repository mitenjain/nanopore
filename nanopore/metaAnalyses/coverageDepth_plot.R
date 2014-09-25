#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
library(lattice)

inFile <- args[1]
outFile <- args[2]
depthFile = read.delim(inFile, sep="\t")
pdf(args[2])
par(mfrow=c(1,1))
plot(unlist(depthFile[2]), unlist(depthFile[3]), main="Coverage across reference", xlab="Position ac
ross reference", ylab="Coverage", pch = 20, cex = 1.0, col = "red")
dev.off()
