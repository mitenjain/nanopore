#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

stats <- read.table(args[1], row.names=1)
#implement stats later
insert <- scan(args[2])
delete <- scan(args[3])
out <- args[4]


pdf(out)


par(mfrow=c(1,2))


hist(insert, main="Size of insertions (to reference)", xlab="Size (bp)")
hist(delete, main="Size of deletions (from reference)", xlab="Size (bp)")


dev.off()
