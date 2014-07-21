#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

data <- read.table(args[1], sep="\t", header=T)


if ( dim(data)[1] > 0 ) {

	pdf(args[2])


	par(mfrow=c(3,2))

	hist(100 * data$alignmentIdentity, main="Alignment Identity", xlab="% Identity")
	hist(100 * data$alignmentCoverage, main="Alignment Coverage", xlab="% Coverage")
	hist(100 * data$readIdentity, main="Read Identity", xlab="% Identity")
	hist(100 * data$readCoverage, main="Read Coverage", xlab="% Coverage")
	hist(100 * data$referenceIdentity, main="Reference Identity", xlab="% Identity")
	hist(100 * data$referenceCoverage, main="Reference Coverage", xlab="% Coverage")


	dev.off()
}
