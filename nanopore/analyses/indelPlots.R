#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

indels <- read.table(args[1], fill=T, sep="\t", header=T, na.strings="None", comment.char="")


if (dim(indels)[1] > 2) {

	pdf(args[2])
	par(mfrow=c(2,1))

	hist(as.numeric(indels$readInsertionLengths), main="Read Insertion\nLength Distribution",xlab="Insertion Length", breaks="FD")
	hist(as.numeric(indels$readDeletionLengths), main="Read Deletion\nLength Distribution",xlab="Deletion Length", breaks="FD")

	hist(as.numeric(indels$readInsertionLengths), main="Read Insertion\nLength Distribution",xlab="Insertion Length", breaks="FD", xlim=c(0,10))
	hist(as.numeric(indels$readDeletionLengths), main="Read Deletion\nLength Distribution",xlab="Deletion Length", breaks="FD", xlim=c(0,10))

	plot(x=as.numeric(indels$ReadSequenceLengths),y=as.numeric(indels$NumberReadInsertions), main="Insertions vs. Read Length", xlab="Read Length", ylab="Read Insertions", pch=19, col="blue")
	plot(x=as.numeric(indels$ReadSequenceLengths),y=as.numeric(indels$NumberReadDeletions), main="Deletions vs. Read Length", xlab="Read Length", ylab="Read Deletions", pch=19, col="blue")

	if (! "NaN" %in% indels$MedianReadInsertionLengths && ! "NaN" %in% indels$MedianReadDeletionLengths) {
		plot(density(as.numeric(indels$MedianReadInsertionLengths), na.rm=T), main="Distribution of Median Insertion Lengths", xlab="Median Insertion Lengths")
		plot(density(as.numeric(indels$MedianReadDeletionLengths), na.rm=T), main="Distribution of Median Deletion Lengths", xlab="Median Deletion Lengths")
	}

	dev.off()

}
