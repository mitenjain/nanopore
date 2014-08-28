#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)


f <- args[1]
out <- args[2]

dist <- read.table(args[1])

if (dim(data[1]) > 1) {
    r <- topo.colors(length(rownames(dist)))
    n <- 0
    plot(dist[1,], type="n", main="Convergence of Likelihoods", xlab="Trials", ylab="Running Likelihoods")
    for (i in 1:length(rownames(dist))) {
        n <- n + 1
        lines(dist[i,], col =r[n], lty=n%%2)

    dev.off()
    }
}