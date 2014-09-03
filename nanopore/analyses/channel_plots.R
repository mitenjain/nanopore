#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

data <- read.table(args[1], row.names=1, header=T)

if (dim(data)[1] > 1) {
    
    sorted <- data[order(data$ReadCount, decreasing=T),]
    sorted <- t(sorted[sorted$ReadCount > 0,])

    png(args[3], height=3000, width=3000, type="cairo")

    q<- barplot(sorted, main=paste("Sorted Channel Mappability", paste("# Reporting = ", length(sorted[1,]), sep=""), sep="\n"), xlab="Channel", ylab="Read Counts", legend.text=T, xaxt="n", col=c("blue","red"), args.legend=c(cex=3), cex.names=3)
    text(cex=0.5, x=q-.25, y=-1.25, colnames(sorted), xpd=T, srt=65)

    dev.off()

    pdf(args[2])

    sorted.percent <- sorted["MappableReadCount",]/sorted["ReadCount",]
    sorted.percent <- sorted.percent[order(sorted.percent, decreasing=T)]
    sorted.percent <- sorted.percent[sorted.percent > 0]
    sorted.percent <- sorted.percent[!is.na(sorted.percent)]

    q<- barplot(sorted.percent, main="Sorted Channel Percent Mappability", xlab="Channel", ylab="Read Counts", xaxt="n")
    text(cex=0.27, x=q-.25,y=-0.005, names(sorted.percent), xpd=T, srt=45)


    #do linear regression
    reg <- lm(sorted["MappableReadCount",]~sorted["ReadCount",])
    #plot scatterplot
    plot(sorted["MappableReadCount",]~sorted["ReadCount",], pch=20, col="blue", xlab="Total Read Count", ylab="Mappable Read Count", main="Mappable vs Total Reads\nReporting Channels Only")
    #add regression line
    abline(reg)
    #add R2
    legend("topleft",legend=c(paste("R^2 = ", round(summary.lm(reg)$adj.r.squared,4))))

    barplot(t(data)[2,]/(t(data)[1,]+t(data)[2,])*100, main="% Mappable Reads Per Channel", xlab="Channel", ylab="% Mappable")

    par(mfrow=(c(1,2)))

    barplot(t(data)[1,], main="Total Read Counts", xlab="Channel", ylab="Read Counts")
    barplot(t(data)[2,], main="Mappable Read Counts", xlab="Channel", ylab="Read Counts")

    dev.off()

}
