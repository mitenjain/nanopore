#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)


f <- args[1]
out <- args[2]
tryCatch({
    dist <- read.table(f)

    if (dim(dist)[1] > 1) {

        m <- max(dist)

        pdf(out)

        r <- topo.colors(length(rownames(dist)))
        n <- 0
        plot(x=seq(1,dim(dist)[2]),y=as.vector(dist[1,]/m), type="n", main="Convergence of Likelihoods", xlab="Trials", ylab="Running Log Likelihoods Ratio")
        for (i in 1:length(rownames(dist))) {
            n <- n + 1
            lines(x=seq(1,dim(dist)[2]),y=as.vector(dist[i,]/m), col =r[n])
        }

        n <- 0
        plot(x=seq(1,dim(dist)[2]),y=log(as.vector(dist[1,]/m)), type="n", main="Convergence of Likelihoods", xlab="Trials", ylab="Running Log-Log Likelihoods Ratio", ylim=c(0,0.01))
        for (i in 1:length(rownames(dist))) {
            n <- n + 1
            lines(x=seq(1,dim(dist)[2]),y=log(as.vector(dist[i,]/m)), col =r[n])
        }

        dev.off()
    }
}, error = function(err) { 
})