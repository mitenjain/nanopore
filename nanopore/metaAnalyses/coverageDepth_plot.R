#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
library(lattice)

inFile <- args[1]
outFile <- args[2]
depthFile = read.delim(inFile, sep="\t")
pdf(args[2])

par(mfrow <- c(1,1))
plot(smooth.spline(unlist(depthFile[2]), unlist(depthFile[3])), main="Coverage across reference", xlab="Position across reference", ylab="Coverage", type = "l", col = "red")

hist(unlist(depthFile[3]), main="Coverage histogram", xlab="Coverage", ylab="Frequency", col = "blue")

plot(density(unlist(depthFile[3])), main="Coverage density", xlab="Coverage", ylab="Frequency", col = "green")

cov_change <- diff(unlist(depthFile[3]))
plot(cov_change, main="Coverage derivative across reference", xlab="Position across reference", ylab="Coverage Change", type = "l", col = "blue")

dev.off()
