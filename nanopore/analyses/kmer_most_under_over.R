#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

if (file.info(args[1])$size != 0) {
    kmers <- read.table(args[1])
    colnames(kmers) <- c("Kmer", "refCount", "RefFreq", "ReadCount", "ReadFreq", "LogRatio", "RefFreq", "LogRatioAgain")

    top <- head(kmers[order(kmers$LogRatio),], n=10L)
    bot <- head(kmers[order(kmers$LogRatio, decreasing=T),], n=10L)

    write.table(top[,c(-7,-8)], file=args[2], sep="\t")
    write.table(bot[,c(-7,-8)], file=args[3], sep="\t")
}