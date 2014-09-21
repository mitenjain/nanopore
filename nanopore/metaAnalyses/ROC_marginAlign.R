#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
library(lattice)

data <- read.table(args[1], row.names=1)
outf <- args[2]
cols <- as.numeric(args[3])
rows <- as.numeric(args[4])
#set it to variant calling algorithms x proportion held out plots
par(mfrow=c(cols, rows))

for (i in 1:cols*rows) {
    if (i %% 2) {
        tpr <- data[i,]
    }
    else if (i %% 3) {
        fpr <- data[i,]
    }
    else {
        xyplot(tpr~fpr, xlab="False positive rate", ylab="True positive rate", main=rowname(data[i,]))
    }
}

dev.off()