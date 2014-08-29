#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

library(lattice)

f <- args[1]
out <- args[2]


myPanel <- function(x, y, z, ...) {
    panel.levelplot(x, y, z, ...)
    panel.text(x, y, round(z,3))
}

d <- read.table(f, header = T)

if ( dim(d)[1] > 0 && sum(d) > 0) {

    pdf(out)

    p <- levelplot(t(as.matrix(d[1,])), ylab="Emission Probability", xlab="Nucleotide", main="Insertion Emission Probabilities", panel = myPanel, col.regions=colorRampPalette(c("white","red"))(256))

    print(p)

    p2 <- levelplot(t(as.matrix(d[2,])), ylab="Emission Probability", xlab="Nucleotide", main="Deletion Emission Probabilities", panel = myPanel, col.regions=colorRampPalette(c("white","red"))(256))

    print(p2)

    dev.off()

}
