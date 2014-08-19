#!/usr/bin/env Rscript
#
args <- commandArgs(trailingOnly = T)
nums <- read.table(args[1])

if (length(nums[,1]) >= 1) {
    pdf(args[2])
    barplot(nums, names.arg=c("Blast Hits", "Mapped","Unapped"), col=c("red","blue","purple"), main=paste(args[3], " # Of Mapped and Unmapped Reads"))
    dev.off()
}