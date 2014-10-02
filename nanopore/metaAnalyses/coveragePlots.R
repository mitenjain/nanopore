#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

if (file.info(args[1])$size != 0) {
    cols <- max(count.fields(args[1], sep=","))
    dist <- read.table(args[1], fill=T, sep=",", row.names=1, col.names=paste("V",seq_len(cols)))
    dist <- dist[order(rownames(dist)),]
}
if (! is.null(dim(dist)) && dim(dist)[2] > 2 && dim(dist)[1] > 1) {
        
    tmp <- dist
    dist <- vector()
    #throw out rows with single-value columns
    for (i in 1:length(rownames(tmp))) {
        r <- tmp[i,]
        if (length(r[!is.na(r)]) > 1) {
            dist <- rbind(dist, r)
        }
    }

    pdf(args[3])

    hists <- list()
    xmax <- 0
    ymax <- 0
    b <- 0
    for (i in 1:length(rownames(dist))){
        b <- max(b, nclass.FD(dist[i,][!is.na(dist[i,])]))
    }
    for (i in 1:length(rownames(dist))){
        hists[[i]] <- hist(dist[i,][!is.na(dist[i,])], plot=F, breaks=b)
        xmax <- max(xmax, hists[[i]]$mids)
        ymax <- max(ymax, hists[[i]]$counts)
    }
    colmap <- expand.grid(1:8, 15:(15+ceiling(length(hists)/8)))
    #whatever
    plot(1, xlim=c(0,xmax), ylim=c(0,ymax), main=paste(args[2],"Identity by Mapper", sep="\n"), xlab="Identity", type="n", ylab="Frequency")


    q <- vector()
    for (i in 1:length(hists)) {
        x <- hists[[i]]
        points(x$mids, x$counts, col=colmap[i,1], pch=colmap[i,2], cex=0.4)
        lines(x$mids, x$counts, col=colmap[i,1], pch=colmap[i,2])
        q[i] <- paste(formatC(round(100*mean(dist[i,][!is.na(dist[i,])]),1),format="f",digits=1), sep= "")
    }
    col <- rep(c(rep(1:8, length.out=length(hists)),rep("NA",times=length(hists))), times=ceiling(length(hists)/8))
    legend(x="top", col=col, pch=colmap[,2], legend=cbind(rownames(dist),q), ncol=2, title="\t\tAverage % Identity", cex=0.55,bty="n")
    
    #hacky way to determine the number of replicates seen
    num <- 1
    for (i in 1:length(rownames(dist))){
        num <- as.numeric(strsplit(rownames(dist)[i], "[.]")[[1]][2])
        if (! is.na(num)) {
            times <- max(times, num)
        }
    }

    #this plot plots the replicates as the same color with different pch and one (average) line
    #build a new color map where we have one color for every aligner group
    colmap <- vector()
    for (i in 1:8){
        colmap <- c(colmap, rep(i, times=num))
    }
    pchmap <- rep(15:(15+num-1), times=ceiling(length(hists)/num))

    #build these bins to force each replicate to have the same midpoints
    bins <- seq(0:(b-1))/b
    bins <- bins[bins <= xmax]


    #force each replicate to have the same midpoints
    hists2 <- list()
    for (i in 1:length(rownames(dist))){
        hists2[[i]] <- hist(dist[i,][!is.na(dist[i,])], plot=F, breaks=bins)
        xmax <- max(xmax, hists2[[i]]$mids)
        ymax <- max(ymax, hists2[[i]]$counts)
    }    


    plot(1, xlim=c(0,xmax), ylim=c(0,ymax), main=paste(args[2],"Identity by Mapper", sep="\n"), xlab="Identity", type="n", ylab="Frequency")
    
    for (i in seq(num, length(hists2), num)) {
        x <- hists2[[i]]
        y <- hists2[[i-1]]
        z <- hists2[[i-2]]
        points(z$mids, z$counts, col=colmap[i-2], pch=pchmap[i-2], cex=0.4)
        points(y$mids, y$counts, col=colmap[i-1], pch=pchmap[i-1], cex=0.4)
        points(x$mids, x$counts, col=colmap[i], pch=pchmap[i], cex=0.4)
        lines(z$mids, z$counts, col=colmap[i-2])
        lines(y$mids, y$counts, col=colmap[i-1])
        lines(x$mids, x$counts, col=colmap[i])
    }

    avgs <- vector()
    names <- vector()
    for (i in seq(num, length(hists2), num)){
        tmp <- mean(c(dist[i,][!is.na(dist[i,])], dist[i-1,][!is.na(dist[i-1,])], dist[i-2,][!is.na(dist[i-2,])]))
        tmp <- paste(formatC(round(100*tmp,1),format="f",digits=1), sep= "")
        avgs <- c(avgs, tmp)
        names <- c(names, rownames(dist)[i-2])
    }

    coltmp <- rep(c(1:length(names), rep("NA", times=length(names))), times=10)
    legend("top", col=coltmp, pch="-", cex=0.55, legend=cbind(names, avgs), ncol=2, title="\t\tAverage % Identity", bty="n")

    dev.off()
}