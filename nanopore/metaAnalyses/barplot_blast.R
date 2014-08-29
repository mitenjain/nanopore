#!/usr/bin/env Rscript
#
args <- commandArgs(trailingOnly = T)
nums <- read.table(args[1])

if (length(nums[,1]) >= 1) {
    pdf(args[2])
    q <- barplot(as.matrix(nums[1,]), names.arg=c("Blast Hits", "Originally Mapped","Still Unmapped"), col=c("red","blue","purple"), main=paste(args[3], " # Of Mapped and Unmapped Reads"))
    text(y=nums+0.01, x=q, labels=paste(as.character(round(nums*100, 4)), "%"), xpd=T)
    dev.off()
}
