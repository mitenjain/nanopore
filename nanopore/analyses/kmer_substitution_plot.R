#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

library(stats)
library(lattice)

f <- args[1]
out <- args[2]
inf <- args[3]


#first make levelplot of all data
png(paste(out, "levelplot.png", sep="_"), width=10000, height=10000, type="cairo")

data <- read.table(f, row.names = 1, header = T)

levelplot(as.matrix(data), scales=list(x=list(rot=45, cex=0.9), 
	y=list(rot=45, cex=0.9)), col.regions=colorRampPalette(c("white","red"))(256), 
	main=out, xlab="Read Kmer", ylab="Reference Kmer")

dev.off()

#second make clustered heatmap of all data; error handling for too sparse matrix

png(paste(out, "heatmap.png", sep="_"), width=5000, height=5000, type="cairo")

title <- paste(inf, "clustered")
r = tryCatch({
    heatmap(as.matrix(data), main = title)
    dev.off()
    }, warning = function(w) {
        dev.off()
        file.remove(paste(out, "heatmap.png", sep="_"))
    }, error = function(e) {
        dev.off()
        file.remove(paste(out, "heatmap.png", sep="_"))
    })



