#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
library(lattice)

data <- read.table(args[1], row.names=1)
outf <- args[2]
cols <- as.numeric(args[3])
rows <- as.numeric(args[4])
#set it to variant calling algorithms x proportion held out plots
par(mfrow=c(cols, rows))

fprs <- data[seq(1, length(data), 3)]
tprs <- data[seq(2, length(data), 3)]

for (i in 1:dim(fprs)[1]) {
    xyplot(tprs[i,]~fpr[i,], xlab="False positive rate", ylab="True positive rate", main=rowname(fprs[i,]), type="l", xlim=c(0,1), ylim=c(0,1))

}

dev.off()