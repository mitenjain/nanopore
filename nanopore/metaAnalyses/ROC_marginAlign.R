#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
library(lattice)

data <- read.table(args[1])
outf <- pdf(args[2])
algorithms <- as.numeric(args[3])

#plot every algorithm on one page, one page per held out amount
#each plot contains ROC curves for each coverage tested
if (algorithms %% 2 == 0) {
    par(mfrow=c(ceiling(sqrt(algorithms)), floor(sqrt(algorithms))))
    plots <- ceiling(sqrt(algorithms)) * floor(sqrt(algorithms))
} else {
    par(mfrow=c(ceiling(sqrt(algorithms)), ceiling(sqrt(algorithms))))
    plots <- ceiling(sqrt(algorithms)) * ceiling(sqrt(algorithms))
}


for (page in seq(algorithms * algorithms, dim(data)[1], algorithms * algorithms)) {
    for (pos in seq(algorithms, page, algorithms)) {
        count <- 0
        tprs <- data[seq(pos-10, pos+1, 2),]
        fprs <- data[seq(pos-11, pos, 2),]
        coverages <- tprs[,4]
        algorithm <- tprs[,2][1]
        held_out <- tprs[,3][1]
        tprs <- tprs[,-(1:4)]
        fprs <- fprs[,-(1:4)]
        matplot(t(fprs), t(tprs), type="l", col=c(1:length(coverages)), xlim=c(0,1), ylim=c(0,1), main=paste("VariantCaller:\n", algorithm, "\nProportionHeldOut: ", held_out, sep=""), cex.main=0.5, cex.axis=0.5, xlab="False Positive Rate", ylab="True Positive Rate")
        legend("topright", legend=coverages, col=c(1:length(coverages)), cex=0.35, pch="-", title="Coverage")
        count <- count + 1
    }
    if (count < plots) {
        while (count < plots) {
            plot.new()
            count <- count + 1
        }
    }
}
dev.off()