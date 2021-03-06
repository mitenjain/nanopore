#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

#thanks SO: http://stackoverflow.com/questions/6602881/text-file-to-list-in-r
raw <- strsplit(readLines(args[1]), "[[:space:]]+")
data <- lapply(raw, tail, n = -1)
names(data) <- lapply(raw, head, n = 1)
data <- lapply(data, as.numeric)

library(lattice)

for (i in 1:length(data)) {
    if ( sum(data[[i]], na.rm=T) == 0 ) {
        quit()
    }
}

if ( length(data$MappedReadLengths) > 1 && length(data$UnmappedReadLengths) > 1) {
    #open a pdf
    pdf(args[2])
    par(mfrow=c(2,2))
    #plot the mapped/unmapped read lengths
    #future - stack them; remove very long unmapped reads?
    lengths <- c(data$MappedReadLengths, data$UnmappedReadLengths)
    lengths.sort <- lengths[order(lengths)]
    #remove the 5% longest reads so they don't skew the graph
    xmax <- lengths.sort[round(0.95*length(lengths.sort))]
    b <- max(nclass.FD(data$MappedReadLengths), nclass.FD(data$UnmappedReadLengths))
    m <- hist(data$MappedReadLengths, breaks = b, plot=F)
    u <- hist(data$UnmappedReadLengths, breaks = b, plot=F)
    plot(m, main=paste("Mapped Read Length Distribution\nn = ", length(data$MappedReadLengths), sep=" "), xlab="Read Length", xlim=c(0,xmax), cex.main=0.8)
    plot(u, main=paste("Unmapped Read Length Distribution\nn = ", length(data$UnmappedReadLengths), sep=" "), xlab="Read Length", xlim=c(0,xmax), cex.main=0.8)
    #plot read coverage distribution
    hist(data$ReadCoverage, breaks="FD", main="Read Coverage Distribution", xlab="Read Coverage", cex.main=0.8)
    hist(data$ReadIdentity, breaks="FD", main="Read Identity Distribution", xlab="Read Identity", cex.main=0.8)
    
    #put new plots on next page
    par(mfrow=c(1,1))
    
    #stacked histogram of mapped/unmapped length distribution
    ymax <- max(u$counts, m$counts)
    plot(m, col=rgb(1,0,0,0.5), xlim=c(0,xmax), ylim=c(0,ymax), main="Mapped and Unmapped Read Length Distributions", xlab="Read Length")
    plot(u, col=rgb(0,0,1,0.5), add=T, xlim=c(0,xmax), ylim=c(0,ymax), main="", xlab="", ylab="")
    legend("topleft", pch=15, legend=c(paste("Mapped n =", length(data$MappedReadLengths)), paste("Unmapped n =", length(data$UnmappedReadLengths))), col=c(rgb(1,0,0),rgb(0,0,1)))
    
    p1 <- xyplot(data$ReadIdentity~data$MappedReadLengths, main=list("Read Identity vs.\nRead Length", cex=0.95), xlab="Read Length", ylab="Read Identity", grid=T, panel=panel.smoothScatter)
    p2 <- xyplot((data$InsertionsPerBase+data$DeletionsPerBase)~data$MappedReadLengths, main=list("Indels Per Aligned Base vs.\nRead Length", cex=0.95), xlab="Read Length", ylab="Indels Per Base", grid=T, panel=panel.smoothScatter)
    p3 <- xyplot(data$MismatchesPerReadBase~data$MappedReadLengths, main=list("Mismatches Per Aligned Base vs.\nRead Length", cex=0.95), ylab="Mismatches Per Aligned Base", xlab="Read Length", grid=T, panel=panel.smoothScatter)
    
    #print the graphs
    print(p1, position=c(0, 0.5, 0.5, 1), more=T)
    print(p2, position=c(0.5, 0.5, 1, 1), more=T)
    print(p3, position=c(0, 0, 0.5, 0.5))

    #now time to do more graphs
    #filter out any read identities that are outside of 3 sd's from mean for fit lines
    #but still plot the points
    m <- mean(data$ReadIdentity)
    s <- sd(data$ReadIdentity)
    inliers <- which(data$ReadIdentity >= m - 3 * s & data$ReadIdentity <= m + 3 * s)
    
    p1 <- xyplot(data$MismatchesPerReadBase~(data$InsertionsPerBase+data$DeletionsPerBase), 
        main=list("Mismatches Per Aligned Base vs.\nIndels Per Aligned Base", cex=0.95), 
        xlab="Indels Per Aligned Base", ylab="Mismatches Per Aligned Base", grid=T, panel=function(...){
        panel.smoothScatter(...)
        panel.lmline(..., lwd=1.2)
        },
        key=simpleKey(paste("R^2 =", as.character(round(summary.lm(
        lm(data$MismatchesPerReadBase~(data$InsertionsPerBase+data$DeletionsPerBase)))$adj.r.squared,3)), sep=" "),
        points=F, corner=c(0,1), cex=0.9))
    p2 <- xyplot((data$InsertionsPerBase+data$DeletionsPerBase)~data$ReadIdentity, 
        main=list("Indels Per Aligned Base vs.\nRead Identity", cex=0.95), xlab="Read Identity", 
        ylab="Indels Per Base", grid=T, panel=function(...){
        panel.smoothScatter(...)
        panel.abline(lwd=1.2, lm((data$InsertionsPerBase+data$DeletionsPerBase)[inliers]~data$ReadIdentity[inliers]))
        },
        key=simpleKey(paste("R^2 =", as.character(round(summary.lm(
        lm((data$InsertionsPerBase+data$DeletionsPerBase)~data$ReadIdentity))$adj.r.squared,3)), sep=" "),
        points=F, corner=c(1,1), cex=0.9))
    p3 <- xyplot(data$MismatchesPerReadBase~data$ReadIdentity, main=list("Mismatches Per Aligned Base vs.\nRead Identity",cex=0.95), 
        ylab="Mismatches Per Aligned Base", xlab="Read Identity", grid=T, panel=function(...){
        panel.smoothScatter(...)
        panel.abline(lwd=1.2, lm(data$MismatchesPerReadBase[data$ReadIdentity>0.6]~data$ReadIdentity[data$ReadIdentity>0.6]))
        },
        key=simpleKey(paste("R^2 =", as.character(round(summary.lm(
        lm((data$InsertionsPerBase+data$DeletionsPerBase)~data$ReadIdentity))$adj.r.squared,3)), sep=" "),
        points=F, corner=c(0,0), cex=0.9))
    p4 <- xyplot(data$InsertionsPerBase~data$DeletionsPerBase, main=list("Insertions Per Aligned Base vs.\nDeletions Per Aligned Base",cex=0.95), 
        xlab="Deletions Per Aligned Base", ylab="Insertions Per Aligned Base", grid=T, panel=function(...){
            panel.smoothScatter(...)
            panel.lmline(..., lwd=1.2)
        },
        key=simpleKey(paste("R^2 =", as.character(round(summary.lm(
        lm(data$InsertionsPerBase~data$DeletionsPerBase))$adj.r.squared,3)), sep=" "),
        points=F, corner=c(0,1), cex=0.9))

    
    #print the graphs
    print(p1, position=c(0, 0.5, 0.5, 1), more=T)
    print(p2, position=c(0.5, 0.5, 1, 1), more=T)
    print(p3, position=c(0, 0, 0.5, 0.5), more=T)
    print(p4, position=c(0.5, 0, 1, 0.5))
    
    
    
    dev.off()
}