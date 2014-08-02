#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

summary <- read.csv(args[1], header=T)
dist <- read.table(args[2], fill=T, sep=",")

if (dim(dist)[2] > 2) {
	
	pdf(args[3])
        #par(mfrow=c(2,3))
	#this is all so hacky - first we make a color rainbow
	n <- 1
	r <- rainbow(length(rownames(dist)))
	#then we find the biggest y value
	m <- 0
	for (i in rownames(dist)) {
		x <- density(as.numeric(dist[i,]), na.rm=T)
		m <- max(m, max(x$y))
	}
	#now we iterate over all of the mappers and plot based on the biggest y value and varying colors
	plot(density(as.numeric(dist[1,]), na.rm=T), xlab="Read Coverage", xlim=c(0,1), main="Coverage by Mapper", ylim=c(0,m), col=r[n])
	for (i in rownames(dist)[-1]) {
		n <- n + 1
		lines(density(as.numeric(dist[i,]), na.rm=T), col = r[n])
	}

	legend(x="topright", y="topright", col=r, legend=rownames(dist), lty=1)

	#barplot of unmapped read counts per mapper
	barplot(height=summary$UnmappedReadCount, names.arg=summary$Mapper, col=r, main="Unmapped Reads by Mapper", ylab="Counts")

	#scatterplot of median read coverage vs avg posterior probability
	plot(100 * summary$MedianReadCoverage, summary$AveragePosteriorMatchProbability, ylab="Average Match Probability", xlab="Median Read Coverage", main="Match Probability vs. Read Coverage", col="blue", pch=19, xlim=c(0,100))
	text(100 * summary$MedianReadCoverage, summary$AveragePosteriorMatchProbability, cex=0.7, pos=3, labels=summary$Mapper)

	#scatterplot of median read coverage vs median identity
	plot(100 * summary$MedianReadCoverage, summary$MedianIdentity, ylab="Median Identity", xlab="Median Read Coverage", main="Median Identity vs. Read Coverage", col="blue", pch=19, xlim=c(0,100))
	text(100 * summary$MedianReadCoverage, summary$MedianIdentity, cex=0.7, pos=3, labels=summary$Mapper)

	#scatterplot of median read coverage vs median # of indels per base
	plot(100 * summary$MedianReadCoverage, summary$MedianDeletionsPerReadBase + summary$MedianInsertionsPerReadBase, ylab="Median Indels Per Read Base", xlab="Median Read Coverage", main="Median Indels vs. Read Coverage", col="blue", pch=19, xlim=c(0,100))
	text(100 * summary$MedianReadCoverage, summary$MedianDeletionsPerReadBase + summary$MedianInsertionsPerReadBase, cex=0.7, pos=3, labels=summary$Mapper)

	dev.off()

}
