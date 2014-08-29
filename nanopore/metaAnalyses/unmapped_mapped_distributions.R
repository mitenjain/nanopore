#!/usr/bin/env Rscript
#
args <- commandArgs(trailingOnly = T)
mapped <- read.table(args[1])
unmapped <- read.table(args[2])
analysis <- args[4]

if (length(unmapped[,1]) >= 1) {
    pdf(args[3])
    
    #remove outliers so graph is prettier
    repoutliers <- function(x){ med=median(x); mad=mad(x); x[x>med+3*mad | x<med-3*mad]=NA; return(x)}
    mapped <- apply(mapped, 2, repoutliers)
    unmapped <- apply(unmapped, 2, repoutliers)

    mapped.density <- density(mapped[,1], na.rm=T)
    unmapped.density <- density(unmapped[,1], na.rm=T)
    xmax <- max(unmapped.density$x, mapped.density$x)
    ymax <- max(unmapped.density$y, mapped.density$y)
    plot(mapped.density, xlim=c(0,xmax), ylim=c(0,ymax), xlab="Read Length", main=paste(analysis, "Mappable and Unmappable Read Lengths", sep=" "), col="blue")
    lines(unmapped.density, col="red")
    legend(x="topright", col=c("blue","red"), legend=c("Mapped","Unmapped"), pch="-")

    lengths <- c(length(unmapped[,1]), length(mapped[,1]))
    barplot(lengths, names.arg=c("Unmapped","Mapped"), col=c("red","blue"), main=paste(analysis, " # Of Mapped and Unmapped Reads"))

    dev.off()
}