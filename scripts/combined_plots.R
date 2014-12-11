#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

#thanks SO: http://stackoverflow.com/questions/6602881/text-file-to-list-in-r
raw1 <- strsplit(readLines(args[1]), "[[:space:]]+")
data1 <- lapply(raw1, tail, n = -1)
names(data1) <- lapply(raw1, head, n = 1)
data1 <- lapply(data1, as.numeric)

pdf(args[2])


library(lattice)

data <- list()
data$MappedReadLengths <- data1$length
data$MismatchesPerReadBase <- data1$mismatches
data$ReadIdentity <- data1$identity
data$DeletionsPerBase <- data1$deletions
data$InsertionsPerBase <- data1$insertions

    #filter out any read identities that are outside of 3 sd's from mean for fit lines
    #but still plot the points
    id.inliers <- which(data$ReadIdentity >= mean(data$ReadIdentity) - 2 * sd(data$ReadIdentity) & data$ReadIdentity <= mean(data$ReadIdentity) + 2 * sd(data$ReadIdentity))
    len.inliers <- which(data$MappedReadLengths >= mean(data$MappedReadLengths) - 2 * sd(data$MappedReadLengths) & data$MappedReadLengths <= mean(data$MappedReadLengths) + 2 * sd(data$MappedReadLengths))
    mismatches.inliers <- which(data$MismatchesPerReadBase >= mean(data$MismatchesPerReadBase) - 2 * sd(data$MismatchesPerReadBase) & data$MismatchesPerReadBase <= mean(data$MismatchesPerReadBase) + 2 * sd(data$MismatchesPerReadBase))
    deletions.inliers <- which(data$DeletionsPerBase >= mean(data$DeletionsPerBase) - 2 * sd(data$DeletionsPerBase) & data$DeletionsPerBase <= mean(data$DeletionsPerBase) + 2 * sd(data$DeletionsPerBase))
    insertions.inliers <- which(data$InsertionsPerBase >= mean(data$InsertionsPerBase) - 2 * sd(data$InsertionsPerBase) & data$InsertionsPerBase <= mean(data$InsertionsPerBase) + 2 * sd(data$InsertionsPerBase))
    inliers <- intersect(insertions.inliers, intersect(deletions.inliers, intersect(mismatches.inliers, intersect(id.inliers, len.inliers))))

    p1 <- xyplot(data$ReadIdentity~data$MappedReadLengths, main=list("Read Identity vs.\nRead Length", cex=0.95), xlab="Read Length", ylab="Read Identity", grid=T, 
        panel=function(...){
        panel.smoothScatter(...)
        panel.abline(lwd=1.2, lm(data$ReadIdentity[inliers]~data$MappedReadLengths[inliers]))
        },
        key=simpleKey(paste("R^2 =", as.character(round(summary.lm(
        lm(data$ReadIdentity[inliers]~data$MappedReadLengths[inliers]))$adj.r.squared,3)), sep=" "),
        points=F, corner=c(1,1), cex=0.9))
    p2 <- xyplot((data$InsertionsPerBase+data$DeletionsPerBase)~data$MappedReadLengths, main=list("Indels Per Aligned Base vs.\nRead Length", cex=0.95), xlab="Read Length", ylab="Indels Per Base", grid=T,
        panel=function(...){
        panel.smoothScatter(...)
        panel.abline(lwd=1.2, lm((data$InsertionsPerBase+data$DeletionsPerBase)[inliers]~data$MappedReadLengths[inliers]))
        },
        key=simpleKey(paste("R^2 =", as.character(round(summary.lm(
        lm((data$InsertionsPerBase+data$DeletionsPerBase)[inliers]~data$MappedReadLengths[inliers]))$adj.r.squared,3)), sep=" "),
        points=F, corner=c(1,1), cex=0.9))
    p3 <- xyplot(data$MismatchesPerReadBase~data$MappedReadLengths, main=list("Mismatches Per Aligned Base vs.\nRead Length", cex=0.95), ylab="Mismatches Per Aligned Base", xlab="Read Length", grid=T,
        panel=function(...){
        panel.smoothScatter(...)
        panel.abline(lwd=1.2, lm(data$MismatchesPerReadBase[inliers]~data$MappedReadLengths[inliers]))
        },
        key=simpleKey(paste("R^2 =", as.character(round(summary.lm(
        lm(data$MismatchesPerReadBase[inliers]~data$MappedReadLengths[inliers]))$adj.r.squared,3)), sep=" "),
        points=F, corner=c(1,1), cex=0.9))
    
    #print the graphs
    print(p1, position=c(0, 0.5, 0.5, 1), more=T)
    print(p2, position=c(0.5, 0.5, 1, 1), more=T)
    print(p3, position=c(0, 0, 0.5, 0.5))

    #now time to do more graphs

    
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
