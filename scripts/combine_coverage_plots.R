#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

#thanks SO: http://stackoverflow.com/questions/6602881/text-file-to-list-in-r
raw1 <- strsplit(readLines(args[1]), "[[:space:]]+")
data1 <- lapply(raw1, tail, n = -1)
names(data1) <- lapply(raw1, head, n = 1)
data1 <- lapply(data1, as.numeric)

raw2 <- strsplit(readLines(args[2]), "[[:space:]]+")
data2 <- lapply(raw2, tail, n = -1)
names(data2) <- lapply(raw2, head, n = 1)
data2 <- lapply(data2, as.numeric)

raw3 <- strsplit(readLines(args[3]), "[[:space:]]+")
data3 <- lapply(raw3, tail, n = -1)
names(data3) <- lapply(raw3, head, n = 1)
data3 <- lapply(data3, as.numeric)

pdf(args[4])


library(lattice)

data <- list()
data$MappedReadLengths <- c(data1$MappedReadLengths,data2$MappedReadLengths,data3$MappedReadLengths)
data$UnmappedReadLengths <- c(data1$UnmappedReadLengths,data2$UnmappedReadLengths,data3$UnmappedReadLengths)
data$MismatchesPerReadBase <- c(data1$MismatchesPerReadBase, data2$MismatchesPerReadBase, data3$MismatchesPerReadBase)
data$ReadIdentity <- c(data1$ReadIdentity, data2$ReadIdentity, data3$ReadIdentity)
data$InsertionsPerBase <- c(data1$InsertionsPerBase, data2$InsertionsPerBase, data3$InsertionsPerBase)
data$DeletionsPerBase <- c(data1$DeletionsPerBase, data2$DeletionsPerBase, data3$DeletionsPerBase)

   #stacked histogram of mapped/unmapped length distribution
    #b <- max(nclass.FD(data$MappedReadLengths), nclass.FD(data$UnmappedReadLengths))
    b <- max(c(data$MappedReadLengths, data$UnmappedReadLengths))
    m <- hist(data$MappedReadLengths, breaks = seq(1,b+100,100), plot=F)
    u <- hist(data$UnmappedReadLengths, breaks = seq(1,b+100,100), plot=F)
    lengths <- c(data$MappedReadLengths, data$UnmappedReadLengths)
    lengths.sort <- lengths[order(lengths)]
    xmax <- lengths.sort[round(0.99*length(lengths.sort))]+1000
    ymax <- max(u$counts, m$counts)
    plot(m, col=rgb(1,0,0,0.5), xlim=c(0,xmax), ylim=c(0,ymax), main="Mapped and Unmapped Read Length Distributions", xlab="Read Length")
    plot(u, col=rgb(0,0,1,0.5), add=T, xlim=c(0,xmax), ylim=c(0,ymax), main="", xlab="", ylab="")
    legend("top", pch=15, legend=c(paste("Mapped n =", length(data$MappedReadLengths)), paste("Unmapped n =", length(data$UnmappedReadLengths))), col=c(rgb(1,0,0),rgb(0,0,1)))
    
    p1 <- xyplot(data$ReadIdentity~data$MappedReadLengths, main=list("Read Identity vs.\nRead Length", cex=0.95),
        xlab="Read Length", ylab="Read Identity", grid=T, panel=function(...){
        panel.smoothScatter(...)})
        panel.lmline(..., lwd=1.2)
        },
        key=simpleKey(paste("R^2 =", as.character(round(summary.lm(
        lm(data$MismatchesPerReadBase~(data$InsertionsPerBase+data$DeletionsPerBase)))$adj.r.squared,3)), sep=" "),
        points=F, corner=c(0,1), cex=0.9))
    p2 <- xyplot((data$InsertionsPerBase+data$DeletionsPerBase)~data$MappedReadLengths, main=list("Indels Per Aligned Base vs.\nRead Length",
        cex=0.95), xlab="Read Length", ylab="Indels Per Base", grid=T, panel=function(...){
        panel.smoothScatter(...)})
        panel.lmline(..., lwd=1.2)
        },
        key=simpleKey(paste("R^2 =", as.character(round(summary.lm(
        lm(data$MismatchesPerReadBase~(data$InsertionsPerBase+data$DeletionsPerBase)))$adj.r.squared,3)), sep=" "),
        points=F, corner=c(0,1), cex=0.9))
    p3 <- xyplot(data$MismatchesPerReadBase~data$MappedReadLengths, main=list("Mismatches Per Aligned Base vs.\nRead Length",
        cex=0.95), ylab="Mismatches Per Aligned Base", xlab="Read Length", grid=T, panel=function(...){
        panel.smoothScatter(...)})
        panel.lmline(..., lwd=1.2)
        },
        key=simpleKey(paste("R^2 =", as.character(round(summary.lm(
        lm(data$MismatchesPerReadBase~(data$InsertionsPerBase+data$DeletionsPerBase)))$adj.r.squared,3)), sep=" "),
        points=F, corner=c(0,1), cex=0.9))
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
        panel.smoothScatter(...)})
        #panel.lmline(..., lwd=1.2)
        #},
        #key=simpleKey(paste("R^2 =", as.character(round(summary.lm(
        #lm(data$MismatchesPerReadBase~(data$InsertionsPerBase+data$DeletionsPerBase)))$adj.r.squared,3)), sep=" "),
        #points=F, corner=c(0,1), cex=0.9))
    p2 <- xyplot((data$InsertionsPerBase+data$DeletionsPerBase)~data$ReadIdentity, 
        main=list("Indels Per Aligned Base vs.\nRead Identity", cex=0.95), xlab="Read Identity", 
        ylab="Indels Per Base", grid=T, panel=function(...){
        panel.smoothScatter(...)})
        #panel.abline(lwd=1.2, lm((data$InsertionsPerBase+data$DeletionsPerBase)[inliers]~data$ReadIdentity[inliers]))
        #},
        #key=simpleKey(paste("R^2 =", as.character(round(summary.lm(
        #lm((data$InsertionsPerBase+data$DeletionsPerBase)~data$ReadIdentity))$adj.r.squared,3)), sep=" "),
        #points=F, corner=c(1,1), cex=0.9))
    p3 <- xyplot(data$MismatchesPerReadBase~data$ReadIdentity, main=list("Mismatches Per Aligned Base vs.\nRead Identity",cex=0.95), 
        ylab="Mismatches Per Aligned Base", xlab="Read Identity", grid=T, panel=function(...){
        panel.smoothScatter(...)})
        #panel.abline(lwd=1.2, lm(data$MismatchesPerReadBase[data$ReadIdentity>0.6]~data$ReadIdentity[data$ReadIdentity>0.6]))
        #},
        #key=simpleKey(paste("R^2 =", as.character(round(summary.lm(
        #lm((data$InsertionsPerBase+data$DeletionsPerBase)~data$ReadIdentity))$adj.r.squared,3)), sep=" "),
        #points=F, corner=c(0,0), cex=0.9))
    p4 <- xyplot(data$InsertionsPerBase~data$DeletionsPerBase, main=list("Insertions Per Aligned Base vs.\nDeletions Per Aligned Base",cex=0.95), 
        xlab="Deletions Per Aligned Base", ylab="Insertions Per Aligned Base", grid=T, panel=function(...){
            panel.smoothScatter(...)})
        #    panel.lmline(..., lwd=1.2)
        #},
        #key=simpleKey(paste("R^2 =", as.character(round(summary.lm(
        #lm(data$InsertionsPerBase~data$DeletionsPerBase))$adj.r.squared,3)), sep=" "),
        #points=F, corner=c(0,1), cex=0.9))

    
    #print the graphs
    print(p1, position=c(0, 0.5, 0.5, 1), more=T)
    print(p2, position=c(0.5, 0.5, 1, 1), more=T)
    print(p3, position=c(0, 0, 0.5, 0.5), more=T)
    print(p4, position=c(0.5, 0, 1, 0.5))
    
    
dev.off()