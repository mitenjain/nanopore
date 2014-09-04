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

    m <- hist(mapped, breaks = "FD", plot=F)
    u <- hist(unmapped, breaks = "FD", plot=F)
    combined <- c(mapped, unmapped)
    combined.sort <- combined[order(combined)]
    xmax <- combined.sort[round(0.95*length(combined.sort))]
    ymax <- max(m$counts, u$counts)
    plot(m, col=rgb(1,0,0,0.5), xlim=c(0,xmax), ylim=c(0,ymax), main="Mapped and Unmapped Read Length Distributions", xlab="Read Length")
    plot(u, col=rgb(0,0,1,0.5), add=T, xlim=c(0,xmax), ylim=c(0,ymax), main="", xlab="", ylab="")
    legend("topleft", pch=15, legend=c("Mapped", "Unmapped"), col=c(rgb(1,0,0),rgb(0,0,1)))

    lengths <- c(length(unmapped[,1]), length(mapped[,1]))
    barplot(lengths, names.arg=c("Unmapped","Mapped"), col=c("red","blue"), main=paste(analysis, " # Of Mapped and Unmapped Reads"))


    dev.off()
}