#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

summary <- read.csv(args[1], header=T, row.names=1)
name <- paste(args[3], args[4], sep="\n")

if (dim(summary)[1] >= 1) {

	r <- rainbow(length(rownames(summary)))

	pdf(paste(args[2], "pdf", sep="."))

	par(mfrow=c(2,2))

	#barplot of unmapped read counts per mapper
	q <- barplot(height=summary$UnmappedReadCount, xaxt="n", col=r, main=paste(name,"Unmapped Reads",sep="\n"), ylab="Counts", cex.main=0.9)
        text(cex=0.25, x=q-.25,y=-1.25, rownames(summary), xpd=T, srt=45)

	#scatterplot of median read coverage vs avg posterior probability
	plot(100 * summary$MedianReadCoverage, summary$AveragePosteriorMatchProbability, ylab="Average Match Probability", xlab="Median Read Coverage", main=name, col=r, pch=19, xlim=c(0,100), cex.main=0.9)
	text(100 * summary$MedianReadCoverage, summary$AveragePosteriorMatchProbability, cex=0.9, pos=2, labels=rownames(summary))

	#scatterplot of median read coverage vs median identity
	plot(100 * summary$MedianReadCoverage, 100 * summary$MedianIdentity, ylab="Median Identity", xlab="Median Read Coverage", main=name, col=r, pch=19, xlim=c(0,100), cex.main=0.9)
	text(100 * summary$MedianReadCoverage, 100 * summary$MedianIdentity, cex=0.9, pos=2, labels=rownames(summary))

	#scatterplot of median read coverage vs median # of indels per base
	plot(100 * summary$MedianReadCoverage, summary$MedianDeletionsPerReadBase + summary$MedianInsertionsPerReadBase, ylab="Median Indels Per Read Base", xlab="Median Read Coverage", main=name, col=r, pch=19, xlim=c(0,100), cex.main=0.9)
	text(100 * summary$MedianReadCoverage, summary$MedianDeletionsPerReadBase + summary$MedianInsertionsPerReadBase, cex=0.9, pos=2, labels=rownames(summary))

	dev.off()

}