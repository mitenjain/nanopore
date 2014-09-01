#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

if (file.info(args[1])$size != 0) {
	cols <- max(count.fields(args[1], sep=","))
	dist <- read.table(args[1], fill=T, sep=",", row.names=1, col.names=paste("V",seq_len(cols)))
	dist <- dist[order(rownames(dist)),]


	if (dim(dist)[2] > 2) {
		
		tmp <- dist
		dist <- vector()
		#throw out rows with single-value columns
		for (i in 1:length(rownames(tmp))) {
			r <- tmp[i,]
			if (length(r[!is.na(r)]) > 1) {
				dist <- rbind(dist, r)
			}
		}

		pdf(args[2])
		#this is all so hacky - first we make a color topo.colors
		n <- 1
		r <- topo.colors(length(rownames(dist)))
		#then we find the biggest y value
		m <- 0
		for (i in 1:length(rownames(dist))) {
			x <- density(as.numeric(dist[i,]), na.rm=T)
			m <- max(m, max(x$y))
		}
		#now we iterate over all of the mappers and plot based on the biggest y value and varying colors
		plot(density(as.numeric(dist[1,]), na.rm=T, adjust=0.7), xlab="Identity", xlim=c(0,1), main=paste(args[2],"Identity by Mapper", sep="\n"), ylim=c(0,m), col=r[n])
		if (dim(dist)[1] > 1) {
			for (i in 2:length(rownames(dist))) {
				n <- n + 1
				lines(density(as.numeric(dist[i,]), na.rm=T, adjust=0.7), col = r[n])
			}
		}
		legend(x="top", col=r, legend=rownames(dist), cex=0.7, lty=1)

		dev.off()

	}
}