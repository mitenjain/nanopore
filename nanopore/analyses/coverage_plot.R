#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

#thanks SO: http://stackoverflow.com/questions/6602881/text-file-to-list-in-r
raw <- strsplit(readLines(args[1]), "[[:space:]]+")
data <- lapply(raw, tail, n = -1)
names(data) <- lapply(raw, head, n = 1)
data <- lapply(data, as.numeric)

library(lattice)

if ( length(data$MappedReadLengths) > 1 && length(data$UnmappedReadLengths) > 1) {
    #open a pdf
    pdf(args[2])
    par(mfrow=c(2,2))
    #plot the mapped/unmapped read lengths
    #future - stack them; remove very long unmapped reads?
    lengths <- cbind(data$MappedReadLengths, data$UnmappedReadLengths)
    lengths.sort <- lengths[order(lengths)]
    #remove the 5% longest reads so they don't skew the graph
    xmax <- lengths.sort[round(0.95*length(lengths.sort))]
    m <- hist(data$MappedReadLengths, breaks = "FD")   
    u <- hist(data$UnmappedReadLengths, breaks = "FD")
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
    legend("topleft", pch=15, legend=c("Mapped", "Unmapped"), col=c(rgb(1,0,0),rgb(0,0,1)))
    

    p1 <- xyplot(data$ReadCoverage~data$MatchIdentity, main="Read Coverage vs. Match Identity", ylab="Read Coverage", xlab="Match Identity", grid=T, panel=panel.smoothScatter)
    p2 <- xyplot(data$ReadIdentity~data$MappedReadLengths, main="Read Identity vs. Read Length", xlab="Read Length", ylab="Read Identity", grid=T, panel=panel.smoothScatter)
    p3 <- xyplot((data$InsertionsPerBase+data$DeletionsPerBase)~data$MappedReadLengths, main="Indels Per Base vs. Read Length", xlab="Read Length", ylab="Indels Per Base", grid=T, panel=panel.smoothScatter)
    p4 <- xyplot(data$ReadIdentity~(data$InsertionsPerBase+data$DeletionsPerBase), main="Read Identity vs. Indels Per Base", ylab="Read Identity", xlab="Indels Per Base", grid=T, panel=panel.smoothScatter)
    
    #print the graphs
    print(p1, position=c(0, 0.5, 0.5, 1), more=T)
    print(p2, position=c(0.5, 0.5, 1, 1), more=T)
    print(p3, position=c(0, 0, 0.5, 0.5), more=T)
    print(p4, position=c(0.5, 0, 1, 0.5))
    
    dev.off()
}
