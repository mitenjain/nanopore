#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

library(lattice)

f <- args[1]
out <- args[2]
inf <- args[3]



myPanel <- function(x, y, z, ...) {
    panel.levelplot(x, y, z, ...)
    panel.text(x, y, paste(100 * round(exp(-z),4), "%", sep=""))
}

d <- read.table(f, header = T, row.names = 1)

if ( dim(d)[1] > 0 && sum(d) > 0) {

    pdf(out)

    p <- levelplot(as.matrix(-log(d)), main=inf, xlab="Reference bases", ylab="Read bases", panel = myPanel, col.regions=colorRampPalette(c("white","red"))(256))

    print(p)

    dev.off()

}
