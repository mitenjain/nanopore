#!/usr/bin/env Rscript
#
args <- commandArgs(trailingOnly = T)
blast <- args[1]
unmapped <- args[2]
mapped <- args[3]

if (length(nums[,1]) >= 1) {
    pdf(args[3])
    q <- barplot(c(blast, unmapped, mapped)), names.arg=c("Blast Hits", "Still Unmapped", "Originally Mapped"), col=c("red","blue"), main="Percent of Reads"))
    text(y=c(blast, unmapped, mapped)+0.015, x=q, labels=paste(as.character(round(c(blast, unmapped, mapped)*100, 3)), "%"), xpd=T)
    dev.off()
}
