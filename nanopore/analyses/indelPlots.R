#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

cols <- max(count.fields(args[1], sep="\t"))
indels <- read.table(args[1], fill=T, sep="\t", row.names=1, col.names=paste("V",seq_len(cols)))


if (dim(indels)[2] > 1) {
	
	pdf(args[2])
	par(mfrow=c(2,1))

	hist(as.numeric(indels["readInsertionLengths",]), main="Read Insertion\nLength Distribution",xlab="Insertion Length")
	hist(as.numeric(indels["readDeletionLengths",]), main="Read Deletion\nLength Distribution",xlab="Deletion Length")

	plot(x=as.numeric(indels["distributionReadSequenceLengths",]),y=as.numeric(indels["distributionNumberReadInsertions",]), main="Insertions vs. Read Length", xlab="Read Length", ylab="Read Insertions", pch=19, col="blue")
	plot(x=as.numeric(indels["distributionReadSequenceLengths",]),y=as.numeric(indels["distributionNumberReadDeletions",]), main="Deletions vs. Read Length", xlab="Read Length", ylab="Read Deletions", pch=19, col="blue")

	plot(density(as.numeric(indels["distributionMedianReadInsertionLengths",])), main="Distribution of Insertion Lengths", xlab="Median Insertion Lengths")
	plot(density(as.numeric(indels["distributionMedianReadDeletionLengths",])), main="Distribution of Deletion Lengths", xlab="Median Deletion Lengths")
	
	dev.off()

}