#!/usr/bin/env Rscript
#
args <- commandArgs(trailingOnly = T)
blast <- args[1]
unmapped <- args[2]

if (length(nums[,1]) >= 1) {
    pdf(args[3])
    q <- barplot(c(blast, unmapped)), names.arg=c("Blast Hits", "Still Unmapped"), col=c("red","blue"), main="Percent of Reads"))
    text(y=c(blast,unmapped)+0.015, x=q, labels=paste(as.character(round(nums*100, 3)), "%"), xpd=T)
    dev.off()
}
