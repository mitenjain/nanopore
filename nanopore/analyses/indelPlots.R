#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

indels <- read.table(args[1], fill=T, sep="\t", header=T, na.strings="None", comment.char="")

g <- function(p, x) {p*(1-p)^(x-1)} #geometric parameterized to start at 1

if (dim(indels)[1] > 2) {

	pdf(args[2])
	par(mfrow=c(2,1))

	if ( ! is.null(indels$readInsertionLengths) && ! is.null(indels$readDeletionLengths) && length(indels$readInsertionLengths[!is.na(indels$readInsertionLengths)]) > 1 && length(indels$readDeletionLengths[!is.na(indels$readDeletionLengths)]) ) {

		insertionLengths <- indels$readInsertionLengths[!is.na(indels$readInsertionLengths)]
		deletionLengths <- indels$readDeletionLengths[!is.na(indels$readDeletionLengths)]

		hist(insertionLengths, main="Read Insertion\nLength Distribution",xlab="Insertion Length", breaks="FD", xlim=c(1,10))
		phat <- 1 / mean(insertionLengths)
		curve(length(insertionLengths)*g(phat, x), xlim=c(1,10), add=T)

		hist(deletionLengths, main="Read Deletion\nLength Distribution",xlab="Deletion Length", breaks="FD", xlim=c(1,10))
		phat <- 1 / mean(deletionLengths)
		curve(length(insertionLengths)*g(phat, x), xlim=c(1,10), add=T)

	}

	plot(x=as.numeric(indels$ReadSequenceLengths),y=as.numeric(indels$NumberReadInsertions), main="Insertions vs. Read Length", xlab="Read Length", ylab="Read Insertions", pch=19, col="blue")
	plot(x=as.numeric(indels$ReadSequenceLengths),y=as.numeric(indels$NumberReadDeletions), main="Deletions vs. Read Length", xlab="Read Length", ylab="Read Deletions", pch=19, col="blue")

	if (! "NaN" %in% indels$MedianReadInsertionLengths && ! "NaN" %in% indels$MedianReadDeletionLengths) {
		if (! length(indels$MedianReadDeletionLengths[!is.na(indels$MedianReadDeletionLengths)]) > 1 && ! length(indels$MedianReadInsertionLengths[!is.na(indels$MedianReadInsertionLengths)]) > 1 ) {
			barplot( c(indels$MedianReadDeletionLengths[!is.na(indels$MedianReadDeletionLengths)], indels$MedianReadInsertionLengths[!is.na(indels$MedianReadInsertionLengths)]), main="Median Insertion Lengths", names.arg=c("Deletions", "Insertions"))
		}
		else {
			plot(density(as.numeric(indels$MedianReadInsertionLengths), na.rm=T), main="Distribution of Median Insertion Lengths", xlab="Median Insertion Lengths")
			plot(density(as.numeric(indels$MedianReadDeletionLengths), na.rm=T), main="Distribution of Median Deletion Lengths", xlab="Median Deletion Lengths")
		}
	}

	dev.off()

}