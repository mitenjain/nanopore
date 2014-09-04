#!/usr/bin/env Rscript
#
args <- commandArgs(trailingOnly = T)
blast <- round(as.numeric(args[1]), 3)
unmapped <- round(as.numeric(args[2]), 3)
mapped <- round(as.numeric(args[3]), 3)
readType <- args[4]

pdf(args[5])
q <- barplot(c(blast, unmapped, mapped), names.arg=c("Blast Hits", "Still Unmapped", "Originally Mapped"), col=c("red","blue","green"), main=paste(readType, "Percent of Reads", sep=" "))
text(y=c(blast, unmapped, mapped)+0.02, x=q, labels=paste(as.character(round(c(blast, unmapped, mapped)*100, 3)), "%"), xpd=T)
dev.off()
