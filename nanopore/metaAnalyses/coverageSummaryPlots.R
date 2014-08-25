#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

summary <- read.csv(args[1], header=T, row.names=1)
name <- paste(args[3], args[4], sep="\n")
summary <- summary[order(rownames(summary)),]

if (dim(summary)[1] >= 1) {

	r <- topo.colors(length(rownames(summary)))

	pdf(paste(args[2], "pdf", sep="."))

	#par(mfrow=c(2,2))

	#scatterplot of Avg read coverage vs avg posterior probability
	#plot(100 * summary$AvgReadCoverage, summary$AveragePosteriorMatchProbability, ylab="Avg. Match Probability", xlab="Avg. Read Coverage", main=name, col=r, pch=20, xlim=c(0,100), cex.main=0.9)
	#text(100 * summary$AvgReadCoverage, summary$AveragePosteriorMatchProbability, cex=0.75, pos=2, labels=rownames(summary))
    
    ################################################
    ####Looking at insertion/deletion/match rates
    ################################################
    
	plot(summary$AvgInsertionsPerReadBase, summary$AvgDeletionsPerReadBase, ylab="Avg. Insertions Per Aligned Read Base", xlab="Avg. Deletions Per Aligned Read Base", main=name, col=r, pch=c(19,18,15,17), xlim=c(0,0.5), ylim=c(0,0.5), cex.main=0.9)
    text(summary$AvgInsertionsPerReadBase, summary$AvgDeletionsPerReadBase, cex=0.75, pos=2, labels=rownames(summary))
	#legend(x="top", legend=rownames(summary), col=r, pch=c(19,18,15,17), cex=0.75)
    
    #scatterplot of Avg identity coverage vs Avg # of indels per base
	plot(100 * summary$AvgMatchIdentity, summary$AvgDeletionsPerReadBase + summary$AvgInsertionsPerReadBase, ylab="Avg. Indels Per Aligned Read Base", xlab="Avg. Match Identity", main=name, col=r, pch=20, xlim=c(0,100), cex.main=0.9)
    text(100 * summary$AvgMatchIdentity, summary$AvgDeletionsPerReadBase + summary$AvgInsertionsPerReadBase, cex=0.75, pos=2, labels=rownames(summary))

    q <- barplot(height=summary$AvgInsertionsPerReadBase, xaxt="n", col=r, main=paste(name,"Avg. Insertions Per Read Base",sep="\n"), ylab="Avg. Insertions Per Read Base", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)
    
    q <- barplot(height=summary$AvgDeletionsPerReadBase, xaxt="n", col=r, main=paste(name,"Avg. Deletions Per Read Base",sep="\n"), ylab="Avg. Deletions Per Read Base", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)
    
    q <- barplot(height=summary$AvgMatchIdentity, xaxt="n", col=r, main=paste(name,"Identity of Aligned Bases",sep="\n"), ylab="Proportion of matched aligned bases.", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)
    
    ################################################
    ####Looking at identity - demonstrating chaining and EM are useful improvements
    ################################################
    
    q <- barplot(height=summary$AvgIdentity, xaxt="n", col=r, main=paste(name,"Identity (inc. gaps)",sep="\n"), ylab="Identity (inc. gaps)", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)
    
    #scatterplot of Avg read coverage vs Avg match identity
	plot(100 * summary$AvgReadCoverage, 100 * summary$AvgMatchIdentity, ylab="Avg. Match Identity", xlab="Avg. Read Coverage", main=name, col=r, pch=20, xlim=c(0,100), ylim=c(0,100), cex.main=0.9)
	text(100 * summary$AvgReadCoverage, 100 * summary$AvgMatchIdentity, cex=0.75, pos=2, labels=rownames(summary))

    ################################################
    ####Looking at unmapped reads - demonstrating last is doing a great job
    ################################################

    #barplot of unmapped read counts per mapper
	q <- barplot(height=summary$MappedReadCount/(summary$UnmappedReadCount+summary$MappedReadCount), xaxt="n", col=r, main=paste(name,"Prop. Mapped Reads",sep="\n"), ylab="Prop. Reads Mapped", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)

    ################################################
    ####Combining mapped reads + identity to get measure of pipeline productivity.
    ################################################
    
    #scatterplot of Avg. identity vs Proportion reads mapped.
	plot(100 * summary$AvgIdentity, 100 * summary$MappedReadCount/(summary$UnmappedReadCount+summary$MappedReadCount), xlab="Avg. Identity", ylab="Prop. Reads Mapped", main=name, col=r, pch=20, xlim=c(0,100), ylim=c(0,100), cex.main=0.9)
	text(100 * summary$AvgIdentity, 100 * summary$MappedReadCount/(summary$UnmappedReadCount+summary$MappedReadCount), cex=0.75, pos=2, labels=rownames(summary))
    
    q <- barplot(height=((summary$MappedReadCount/(summary$UnmappedReadCount+summary$MappedReadCount))*summary$AvgIdentity), xaxt="n", col=r, main=paste(name,"Prop. Of Sequenced Bases Aligned with Identity",sep="\n"), ylab="Prop. Of Sequenced Bases Aligned with Identity", cex.main=0.9)
    text(cex=0.8, x=q, y=-0.01, rownames(summary), xpd=T, srt=90)
    
    #legend(x="top", legend=rownames(summary), col=r, pch=c(19,18,15,17), cex=0.75)


    #legend(x="top", legend=rownames(summary), col=r, pch=c(19,18,15,17), cex=0.75)
	#text(100 * summary$AvgIdentityCoverage, summary$AvgDeletionsPerReadBase + summary$AvgInsertionsPerReadBase, cex=0.75, pos=2, labels=rownames(summary))

	#scatterplot of Avg read coverage vs avg posterior probability
	#plot(100 * summary$AvgReadCoverage, summary$AveragePosteriorMatchProbability, ylab="Avg. Match Probability", xlab="Avg. Read Coverage", main=name, col=r, pch=c(19,18,15,17), xlim=c(0,100), cex.main=0.9)
	#legend(x="top", legend=rownames(summary), col=r, pch=c(19,18,15,17), cex=0.75)

	#scatterplot of Avg read coverage vs Avg identity
	#plot(100 * summary$AvgReadCoverage, 100 * summary$AvgIdentity, ylab="Avg. Identity", xlab="Avg. Read Coverage", main=name, col=r, pch=c(19,18,15,17), xlim=c(0,100), cex.main=0.9)
	#legend(x="top", legend=rownames(summary), col=r, pch=c(19,18,15,17), cex=0.75)




	dev.off()

}
