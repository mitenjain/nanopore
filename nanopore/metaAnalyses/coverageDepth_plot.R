#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
library(lattice)

inFile <- args[1]
outFile_1 <- args[2]
outFile_2 <- args[3]

depthFile = read.delim(inFile, sep="\t")

outliers <- length(boxplot.stats(c(unlist(depthFile[3], dpois(unlist(depthFile[3]), mean(unlist(depthFile[3]), "pois", "mle")$estimate))))$out)
proportion <- round((outliers / length(unlist(depthFile[3]))) * 100, 2)

pdf(args[2])
par(mfrow <- c(1,1))

plot(density(unlist(depthFile[3])), main=paste("Coverage density\n Number of outliers = ", outliers, ", Genome proportion = ", proportion, "%"), xlab="Coverage", ylab="Frequency", col = "green")
hist(unlist(depthFile[3]), main=paste("Coverage histogram\n Number of outliers = ", outliers, ", Genome proportion = ", proportion, "%"), xlab="Coverage", ylab="Frequency", col = "blue")

dev.off()

pdf(args[3])
par(mfrow <- c(1,1))

plot(smooth.spline(unlist(depthFile[2]), unlist(depthFile[3])), main="Coverage across reference", xlab="Position across reference", ylab="Coverage", type = "l", col = "red")
cov_change <- diff(unlist(depthFile[3]))
plot(cov_change, main="Coverage derivative across reference", xlab="Position across reference", ylab="Coverage Change", type = "l", col = "blue")

dev.off()
