#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

dist <- read.table(args[1], fill=T, sep=",", row.names=1)

if (dim(dist)[2] > 2) {
	
	pdf(args[2])
	#this is all so hacky - first we make a color rainbow
	n <- 1
	r <- rainbow(length(rownames(dist)))
	#then we find the biggest y value
	m <- 0
	for (i in 1:length(rownames(dist))) {
		x <- density(as.numeric(dist[i,]), na.rm=T)
		m <- max(m, max(x$y))
	}
	#now we iterate over all of the mappers and plot based on the biggest y value and varying colors
	plot(density(as.numeric(dist[1,]), na.rm=T), xlab="Read Coverage", xlim=c(0,1), main="Coverage by Mapper", ylim=c(0,m), col=r[n])
	for (i in 2:length(rownames(dist))) {
		n <- n + 1
		lines(density(as.numeric(dist[i,]), na.rm=T), col = r[n])
	}

	legend(x="topleft", y="topleft", col=r, legend=rownames(dist), lty=1)

	dev.off()

}
