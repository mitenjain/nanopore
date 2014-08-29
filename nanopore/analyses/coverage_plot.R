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
    xmax <- lengths.sort[round(0.95*length(lengths.sort))]
    hist(data$MappedReadLengths, breaks = "FD", main=paste("Mapped Read Length Distribution n= ", length(data$MappedReadLengths), sep=" "), xlab="Read Length", xlim=c(0,xmax), cex.main=0.8)
    hist(data$UnmappedReadLengths, breaks = "FD", main=paste("Unmapped Read Length Distribution n= ", length(data$UnmappedReadLengths), sep=" "), xlab="Read Length", xlim=c(0,xmax), cex.main=0.8)
    #plot read coverage distribution
    hist(data$ReadCoverage, breaks="FD", main="Read Coverage Distribution", xlab="Read Coverage")

    #plot relative density of mapped to unmapped reads
    #first, find max x value to expect
    xmax <- max(data$UnmappedReadLengths,data$MappedReadLengths)
    #generate a combined data set of all reads
    tot <- rbind(c(data$UnmappedReadLengths,data$MappedReadLengths))
    #find density of mapped and combined
    totDens <- density(tot, from=0, to=xmax, n=1000, adjust=0.5)
    mapDens <- density(data$MappedReadLengths, from=0, to=xmax, n=1000, adjust=0.5)
    #find porportionality difference in # of members in each group
    mapProp <- length(data$MappedReadLengths)/length(tot)
    plot(mapDens$x, mapProp*mapDens$y/totDens$y, type="l", xlab="Read Length", ylab="Relative Density", main="Relative Density\nMapped To Unmapped Read Length", ylim=c(0,2))

    #put new plots on next page
    par(mfrow=c(1,1))
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
