#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)


if (file.info(args[1])$size != 0) {
    cols <- max(count.fields(args[1], sep=","))
    dist <- read.table(args[1], fill=T, sep=",", row.names=1, col.names=paste("V",seq_len(cols)))
    dist <- dist[order(rownames(dist)),]

    if (dim(dist)[2] > 2) {
        pdf(args[2])
        name <- args[3]

        plot(density(as.numeric(dist[1,]), na.rm=T, adjust=0.7), xlab="Identity", xlim=c(0,1), main=paste(args[3], "Identity Distribution", sep=" "))
    }
}