#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

library(stats)
library(lattice)

f <- args[1]
out <- args[2]
inf <- args[3]

png(paste(out, "levelplot.png", sep="_"), width=10000, height=10000)

data <- read.table(f, row.names = 1, header = T)

levelplot(as.matrix(data), scales=list(x=list(rot=45, cex=0.9), 
	y=list(rot=45, cex=0.9)), col.regions=colorRampPalette(c("white","red"))(256), 
	main=out, xlab="Read Kmer", ylab="Reference Kmer")

#dev.off()

#png(paste(out, "heatmap.png", sep="_"), width=5000, height=5000)

#title <- paste(inf, "clustered")

#heatmap(as.matrix(data), main = title)

#dev.off()
