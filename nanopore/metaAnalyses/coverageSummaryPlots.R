#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

summary <- read.csv(args[1], header=T, row.names=1)
name <- args[2]
summary <- summary[order(rownames(summary)),]

if (dim(summary)[1] >= 1) {

    r <- topo.colors(length(rownames(summary)))

    pdf(args[3])


    ################################################
    ####Looking at insertion/deletion/match rates
    ################################################
    
    plot(summary$AvgInsertionsPerReadBase, summary$AvgDeletionsPerReadBase, xlab="Avg. Insertions Per Aligned Base", ylab="Avg. Deletions Per Aligned Base", main=name, col=r, pch=c(19,18,15,17), xlim=c(0,0.2), ylim=c(0,0.2), cex.main=0.9)
    text(summary$AvgInsertionsPerReadBase, summary$AvgDeletionsPerReadBase, cex=0.75, pos=2, labels=rownames(summary))
    #legend(x="top", legend=rownames(summary), col=r, pch=c(19,18,15,17), cex=0.75)
    
    #scatterplot of Avg identity coverage vs Avg # of indels per base
    plot(100 * summary$AvgMismatchesPerReadBase, summary$AvgDeletionsPerReadBase + summary$AvgInsertionsPerReadBase, ylab="Avg. Indels Per Aligned Base", xlab="Avg. Mismatches Per Aligned Base", main=name, col=r, pch=20, xlim=c(0,100), cex.main=0.9)
    text(100 * summary$AvgMismatchesPerReadBase, summary$AvgDeletionsPerReadBase + summary$AvgInsertionsPerReadBase, cex=0.75, pos=2, labels=rownames(summary))

    q <- barplot(height=summary$AvgInsertionsPerReadBase, xaxt="n", col=r, main=paste(name,"Avg. Insertions Per Aligned Base",sep="\n"), ylab="Avg. Insertions Per Aligned Base", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)
    
    q <- barplot(height=summary$AvgDeletionsPerReadBase, xaxt="n", col=r, main=paste(name,"Avg. Deletions Per Aligned Base",sep="\n"), ylab="Avg. Deletions Per Aligned Base", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)
    
    q <- barplot(height=summary$AvgMismatchesPerReadBase, xaxt="n", col=r, main=paste(name,"Avg. Mismatches Per Aligned Base",sep="\n"), ylab="Avg. Mismatches Per Aligned Base", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)
    
    ################################################
    ####Looking at identity - demonstrating chaining and EM are useful improvements
    ################################################
    
    q <- barplot(height=summary$AvgIdentity, xaxt="n", col=r, main=paste(name,"Identity",sep="\n"), ylab="Identity", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)
    
    #scatterplot of Avg read coverage vs Avg match identity
    plot(100 * summary$AvgReadCoverage, 100 * summary$AvgMismatchesPerReadBase, ylab="Avg. Mismatches Per Aligned Base", xlab="Avg. Read Coverage", main=name, col=r, pch=20, xlim=c(0,100), ylim=c(0,100), cex.main=0.9)
    text(100 * summary$AvgReadCoverage, 100 * summary$AvgMismatchesPerReadBase, cex=0.75, pos=2, labels=rownames(summary))

    ################################################
    ####Looking at unmapped reads - demonstrating last is doing a great job
    ################################################

    #barplot of unmapped read counts per mapper
    q <- barplot(height=summary$NumberOfMappedReads/summary$NumberOfReads, xaxt="n", col=r, main=paste(name,"Prop. Mapped Reads",sep="\n"), ylab="Prop. Reads Mapped", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)

    ################################################
    ####Combining mapped reads + identity to get measure of pipeline productivity.
    ################################################
    
    #scatterplot of Avg. identity vs Proportion reads mapped.
    plot(100 * summary$AvgIdentity, 100 * summary$NumberOfMappedReads/summary$NumberOfReads, xlab="Avg. Identity", ylab="Prop. Reads Mapped", main=name, col=r, pch=20, xlim=c(0,100), ylim=c(0,100), cex.main=0.9)
    text(100 * summary$AvgIdentity, 100 * summary$NumberOfMappedReads/summary$NumberOfReads, cex=0.75, pos=2, labels=rownames(summary))
    
    q <- barplot(height=(summary$NumberOfMappedReads/summary$NumberOfReads*summary$AvgIdentity), xaxt="n", col=r, main=paste(name,"Prop. Of Sequenced Bases Aligned with Identity",sep="\n"), ylab="Prop. Of Sequenced Bases Aligned with Identity", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)



    dev.off()

}
