#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

summary <- read.csv(args[1], header=T, row.names=1)
name <- args[2]
summary <- summary[order(rownames(summary)),]



if (dim(summary)[1] >= 1) {


    num <- 1
    names <- vector()
    for (i in 1:length(rownames(summary))){
        tmp <- as.numeric(strsplit(rownames(summary)[i], "[.]")[[1]][2])
        n <- as.character(strsplit(rownames(summary)[i], "[.]")[[1]][1])
        if (! is.na(tmp)) {
            num <- max(tmp, num)
        }
        else {
            names <- c(names, n)
        }
    }
    colmap <- vector()
    for (i in 1:8){
        colmap <- c(colmap, rep(i, times=num))
    }

    pchmap <- rep(15:(15+num-1), times=ceiling(length(hists)/num))

    pdf(args[3])


    ################################################
    ####Looking at insertion/deletion/match rates
    ################################################
    
    plot(summary$AvgInsertionsPerReadBase, summary$AvgDeletionsPerReadBase, xlab="Avg. Insertions Per Aligned Read Base", ylab="Avg. Deletions Per Aligned Read Base", main=name, col=colmap, pch=pchmap, xlim=c(0,0.2), ylim=c(0,0.2), cex.main=0.9)
    legend("topright", cex=0.75, legend=names, pch="*", col=1:8, pt.cex=1.2)
    
    #scatterplot of Avg identity coverage vs Avg # of indels per base
    plot(100 * summary$AvgIdentity, summary$AvgDeletionsPerReadBase + summary$AvgInsertionsPerReadBase, ylab="Avg. Indels Per Aligned Read Base", xlab="Avg. Match Identity", main=name, col=colmap, pch=pchmap, xlim=c(0,100), cex.main=0.9)
    legend("topright", cex=0.75, legend=names, pch="*", col=1:8)

    q <- barplot(height=summary$AvgInsertionsPerReadBase, xaxt="n", col=colmap, main=paste(name,"Avg. Insertions Per Aligned Base",sep="\n"), ylab="Avg. Insertions Per Aligned Base", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)
    
    q <- barplot(height=summary$AvgDeletionsPerReadBase, xaxt="n", col=colmap, main=paste(name,"Avg. Deletions Per Aligned Base",sep="\n"), ylab="Avg. Deletions Per Aligned Base", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)
    
    q <- barplot(height=summary$AvgIdentity, xaxt="n", col=colmap, main=paste(name,"Identity of Aligned Bases",sep="\n"), ylab="Proportion of matched aligned bases.", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)
    
    ################################################
    ####Looking at identity - demonstrating chaining and EM are useful improvements
    ################################################
    
    q <- barplot(height=summary$AvgIdentity, xaxt="n", col=colmap, main=paste(name,"Identity (inc. gaps)",sep="\n"), ylab="Identity (inc. gaps)", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)
    
    #scatterplot of Avg read coverage vs Avg match identity
    plot(100 * summary$AvgReadCoverage, 100 * summary$AvgIdentity, ylab="Avg. Match Identity", xlab="Avg. Read Coverage", main=name, col=colmap, pch=pchmap, xlim=c(0,100), ylim=c(0,100), cex.main=0.9)
    legend("topright", cex=0.75, legend=names, pch="*", col=1:8)

    ################################################
    ####Looking at unmapped reads - demonstrating last is doing a great job
    ################################################

    #barplot of unmapped read counts per mapper
    q <- barplot(height=summary$NumberOfMappedReads/summary$NumberOfReads, xaxt="n", col=colmap, main=paste(name,"Prop. Mapped Reads",sep="\n"), ylab="Prop. Reads Mapped", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)

    ################################################
    ####Combining mapped reads + identity to get measure of pipeline productivity.
    ################################################
    
    #scatterplot of Avg. identity vs Proportion reads mapped.
    plot(100 * summary$AvgIdentity, 100 * summary$NumberOfMappedReads/summary$NumberOfReads, xlab="Avg. Identity", ylab="Prop. Reads Mapped", main=name, col=colmap, pch=20, xlim=c(0,100), ylim=c(0,100), cex.main=0.9)
    text(100 * summary$AvgIdentity, 100 * summary$NumberOfMappedReads/summary$NumberOfReads, cex=0.65, pos=2, labels=rownames(summary))
    
    q <- barplot(height=(summary$NumberOfMappedReads/summary$NumberOfReads*summary$AvgIdentity), xaxt="n", col=colmap, main=paste(name,"Prop. Of Sequenced Bases Aligned with Identity",sep="\n"), ylab="Prop. Of Sequenced Bases Aligned with Identity", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)



    dev.off()

}
