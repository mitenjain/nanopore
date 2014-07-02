#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

stats <- read.table(args[1], row.names=1)
#implement stats later
insert <- scan(args[2])
delete <- scan(args[3])
out <- args[4]
png(out)

#plotrix contains a function to add tables to plots
#library(plotrix)

#set it for two plots and extra room for the table
#op <- par(oma=c(0,0,0,7), mfrow=c(1,2), xpd=T)
par(mfrow=c(1,2))

hist(insert, main="Insertions (to reference)", xlab="Size")
hist(delete, main="Deletions (from reference)", xlab="Size")

dev.off()
