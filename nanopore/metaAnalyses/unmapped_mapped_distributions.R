#!/usr/bin/env Rscript
#
args <- commandArgs(trailingOnly = T)
unmapped <- read.table(args[1])
mapped <- read.table(args[2])
analysis <- args[4]

if (length(unmapped[,1]) >= 1 && length(mapped[,1]) >= 1) {
    pdf(args[3])

    m <- hist(mapped[,1], breaks = "FD", plot=F)
    u <- hist(unmapped[,1], breaks = "FD", plot=F)
    combined <- c(mapped[,1], unmapped[,1])
    combined <- combined[!is.na(combined)]
    combined.sort <- combined[order(combined)]
    xmax <- combined.sort[round(0.98*length(combined.sort))]
    ymax <- max(m$counts, u$counts)
    plot(m, col=rgb(1,0,0,0.5), xlim=c(0,xmax), ylim=c(0,ymax), main="Mapped and Unmapped Read Length Distributions", xlab="Read Length")
    plot(u, col=rgb(0,0,1,0.5), add=T, xlim=c(0,xmax), ylim=c(0,ymax), main="", xlab="", ylab="")
    legend("topright", pch=15, legend=c(paste("Mapped n =", dim(mapped)[1]), paste("Unmapped n =", dim(unmapped)[1])), col=c(rgb(1,0,0),rgb(0,0,1)))

    #stacked histogram without removing outliers
    ymax <- max(u$counts, m$counts)
    xmax <- combined.sort[length(combined.sort)]
    plot(m, col=rgb(1,0,0,0.5), xlim=c(0,xmax), ylim=c(0,ymax), main="Mapped and Unmapped Read Length Distributions", xlab="Read Length")
    plot(u, col=rgb(0,0,1,0.5), add=T, xlim=c(0,xmax), ylim=c(0,ymax), main="", xlab="", ylab="")
    legend("topright", pch=15, legend=c(paste("Mapped n =", dim(mapped)[1]), paste("Unmapped n =", dim(unmapped)[1])), col=c(rgb(1,0,0),rgb(0,0,1)))


    lengths <- c(length(unmapped[,1]), length(mapped[,1]))
    q <- barplot(lengths, names.arg=c("Unmapped","Mapped"), col=c("red","blue"), main=paste(analysis, " # Of Mapped and Unmapped Reads"))
    text(y=lengths+0.02, x=q, labels=paste(as.character(round(lengths/sum(lengths)*100),3), "%"), xpd=T)

    dev.off()
}