#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

data <- read.table(args[1], row.names=1, header=T)

if (dim(data)[1] > 1) {
    
    pdf(args[2])

    sorted <- data[order(data$ReadCount, decreasing=T),]
    q<- barplot(as.matrix(t(sorted)), main="Sorted Channel Mappability", xlab="Channel", ylab="Read Counts", legend.text=T, xaxt="n")
    text(cex=0.2, x=q-.25,y=-1.25, rownames(sorted), xpd=T, srt=45)

    barplot(t(data)[2,]/(t(data)[1,]+t(data)[2,])*100, main="% Mappable Reads Per Channel", xlab="Channel", ylab="% Mappable")

    par(mfrow=(c(1,2)))

    barplot(t(data)[1,], main="Total Read Counts", xlab="Channel", ylab="Read Counts")
    barplot(t(data)[2,], main="Mappable Read Counts", xlab="Channel", ylab="Read Counts")

    dev.off()

}