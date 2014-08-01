#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

#thanks SO: http://stackoverflow.com/questions/6602881/text-file-to-list-in-r
raw <- strsplit(readLines(args[1]), "[[:space:]]+")
data <- lapply(raw, tail, n = -1)
names(data) <- lapply(raw, head, n = 1)
data <- lapply(data, as.numeric)

library(lattice)

if ( length(data$MappedReadLengths) > 1 ) {
    #open a pdf
    pdf(args[2])
    par(mfrow=c(2,2))
    #plot the mapped/unmapped read lengths
    #future - stack them; remove very long unmapped reads?
    hist(data$MappedReadLengths, breaks = "FD", main="Mapped Read Length Distribution", xlab="Read Length")
    hist(data$UnmappedReadLengths, breaks = "FD", main="Unmapped Read Length Distribution", xlab="Read Length")
    #plot read coverage distribution
    hist(data$ReadCoverage, breaks="FD", main="Read Coverage Distribution", xlab="Read Coverage")
    #put new plots on next page
    par(mfrow=c(1,1))
    p1 <- xyplot(data$ReadCoverage~data$ReadIdentity, main="Read Coverage vs. Read Identity", ylab="Read Coverage", xlab="Read Identity", grid=T, panel=panel.smoothScatter)
    p2 <- xyplot(data$ReadCoverage~data$MappedReadLengths, main="Read Coverage vs. Read Length", ylab="Read Coverage", xlab="Read Length", grid=T, panel=panel.smoothScatter)
    p3 <- xyplot(data$MappedReadLengths~data$ReadIdentity, main="Read Length vs. Read Identity", ylab="Read Length", xlab="Read Identity", grid=T, panel=panel.smoothScatter)
    p4 <- xyplot(data$ReadCoverage~(data$InsertionsPerBase+data$DeletionsPerBase), main="Coverage vs. Indels Per Base", ylab="Read Coverage", xlab="Indels Per Base", grid=T, panel=panel.smoothScatter)

    #print the graphs
    print(p1, position=c(0, 0.5, 0.5, 1), more=T)
    print(p2, position=c(0.5, 0.5, 1, 1), more=T)
    print(p3, position=c(0, 0, 0.5, 0.5), more=T)
    print(p4, position=c(0.5, 0, 1, 0.5))

    dev.off()
}