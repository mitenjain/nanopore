#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
library(lattice)

data <- read.table(args[1])
outf <- args[2]
cols <- as.numeric(args[3])
rows <- as.numeric(args[4])
#set it to variant calling algorithms x proportion held out plots
par(mfrow=c(cols, rows))

fprs <- data[seq(1, length(data), 2)]
tprs <- data[seq(2, length(data), 2)]

for (i in 1:dim(fprs)[1]) {
    xyplot(tprs[i,]~fpr[i,], xlab="False positive rate", ylab="True positive rate", main=paste("variantCaller:\n", fprs[i,][1], "\nProportion held out: ", fprs[i,][2], sep=" "), type="l", xlim=c(0,1), ylim=c(0,1))

}

dev.off()