#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

library(lattice)

f <- args[1]
out <- args[2]
inf <- args[3]

pdf(out)

myPanel <- function(x, y, z, ...) {
    panel.levelplot(x, y, z, ...)
    panel.text(x, y, 100 * round(exp(-z),4))
}

d <- read.table(f, header = T, row.names = 1)

if ( sum(rowSums(d)) > 0 ) {

	levelplot(as.matrix(-log(d)), main=inf, xlab="Read bases", ylab="Reference bases", panel = myPanel, col.regions=colorRampPalette(c("white","red"))(256))
	#use the version below if you want two colors
	#levelplot(as.matrix(-log(d)), main=inf, xlab="Read bases", ylab="Reference bases", panel = myPanel)

	dev.off()

}