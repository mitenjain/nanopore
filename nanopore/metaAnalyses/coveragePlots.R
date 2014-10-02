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
    colmap <- expand.grid(1:length(hists), 15:(15+ceiling(length(hists)/8)))
    #whatever
    plot(1, xlim=c(0,xmax), ylim=c(0,ymax), main=paste(args[2],"Identity by Mapper", sep="\n"), xlab="Identity", type="n", ylab="Frequency")


    q <- vector()
    for (i in 1:length(hists)) {
        x <- hists[[i]]
        points(x$mids, x$counts, col=colmap[i,1], pch=colmap[i,2], cex=0.3)
        lines(x$mids, x$counts, col=colmap[i,1], pch=colmap[i,2])
        q[i] <- paste(formatC(round(100*mean(dist[i,][!is.na(dist[i,])]),1),format="f",digits=1), "%\t\tn=", as.character(length(dist[i,][!is.na(dist[i,])])) , sep= "")
    }

    legend(x="top", col=c(colmap[,1], rep("NA", times=length(hists))), pch=colmap[,2], legend=cbind(rownames(dist),q), ncol=2, title="\t\t\tAverage % Identity", cex=0.55,bty="n")
    


    dev.off()
}