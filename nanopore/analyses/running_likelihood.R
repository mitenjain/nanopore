#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)


f <- args[1]
out <- args[2]

dist <- read.table(f)

if (dim(dist)[1] > 1) {

    pdf(out)

    r <- topo.colors(length(rownames(dist)))
    n <- 0
    plot(x=seq(1,dim(dist)[2]),y=as.vector(dist[1,]/max(dist[1,])), type="n", main="Convergence of Likelihoods", xlab="Trials", ylab="Running Likelihoods")
    for (i in 1:length(rownames(dist))) {
        n <- n + 1
        lines(x=seq(1,dim(dist)[2]),y=as.vector(dist[i,]/max(dist[i,])), col =r[n], lty=n%%2)
    }
    dev.off()
}