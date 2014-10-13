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
        tmp <- as.numeric(strsplit(rownames(dist)[i], "[.]")[[1]][2])
        if (! is.na(tmp)) {
            num <- max(tmp, num)
        }
    }
    if (num > 1) {
        #this plot plots the replicates as the same color with different pch and one (average) line
        #build a new color map where we have one color for every aligner group
        colmap <- vector()
        for (i in 1:8){
            colmap <- c(colmap, rep(i, times=num))
        }
        pchmap <- rep(15:(15+num-1), times=ceiling(length(hists)/num))

        plot(1, xlim=c(0,xmax), ylim=c(0,ymax), main=paste(args[2],"Identity by Mapper", sep="\n"), xlab="Identity", type="n", ylab="Frequency")
        if (length(hists) > num){
            for (i in seq(num, length(hists), num)) {
                for (j in 0:(num-1)) {
                    x <- hists[[i-j]]
                points(x$mids, x$counts, col=colmap[i-j], pch=pchmap[i-j], cex=0.4)
                lines(x$mids, x$counts, col=colmap[i-j])
                }
            }   
            
            avgs <- vector()
            names <- vector()
            for (i in seq(num, length(hists), num)){
                tmp <- mean(c(dist[i,][!is.na(dist[i,])], dist[i-1,][!is.na(dist[i-1,])], dist[i-2,][!is.na(dist[i-2,])]))
                tmp <- paste(formatC(round(100*tmp,1),format="f",digits=1), sep= "")
                avgs <- c(avgs, tmp)
                names <- c(names, rownames(dist)[i-2])
            }

            coltmp <- rep(c(1:length(names), rep("NA", times=length(names))), times=10)
            legend("top", col=coltmp, pch=pchmap, cex=0.55, legend=cbind(names, avgs), ncol=2, title="\t\tAverage % Identity", bty="n")
        
            #now we plot them added together
            hists2 <- list()
            xmax2 <- 0
            ymax2 <- 0
            for (i in seq(1, length(rownames(dist)), num)){
                t <- ceiling(i/num)
                tmp <- vector()
                for (q in 0:(num-1)) {
                    tmp <- c(tmp,dist[i+q,][!is.na(dist[i+q,])])
                }
                hists2[[t]] <- hist(tmp, plot=F, breaks=b)
                xmax2 <- max(xmax2, hists2[[t]]$mids)
                ymax2 <- max(ymax2, hists2[[t]]$counts)
            }
            plot(1, xlim=c(0,xmax2), ylim=c(0,ymax2), main=paste(args[2],"Identity by Mapper", sep="\n"), xlab="Identity", type="n", ylab="Frequency")
            for (i in 1:length(hists2)) {
                x <- hists2[[i]]
                points(x$mids, x$counts, col=i, cex=0.4, pch=16)
                lines(x$mids, x$counts, col=i)
            }
            legend("top", col=1:8, pch=c(rep("-", times=length(hists2)), rep(NA, times=length(hists2))), cex=0.55, legend=cbind(names, avgs), ncol=2, title="\t\tAverage % Identity", bty="n")
        
        }   
    }
    dev.off()
}