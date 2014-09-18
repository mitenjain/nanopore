#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

data <- read.table(args[1], row.names=1, header=T)
outf <- args[2]
outsig <- args[3]

library(stats)

trials <- sum(data$readCount)

test <- function(x, p, n){binom.test(x, n, p, alternative="two.sided", conf.level=0.99)$p.value}

pvals <- mapply(test, data$readCount, data$refFraction, MoreArgs=list(n=trials))

adjusted <- p.adjust(pvals, method="bonferroni")

finished <- cbind(data, adjusted)

colnames(finished) <- c(colnames(finished)[1:dim(finished)[2]-1], "p_value")

write.table(finished, outf)

significant <- finished[finished$p_value <= 0.05,]

top <- head(significant[order(significant$foldChange),], n=20L)
bot <- head(significant[order(significant$foldChange, decreasing=T),], n=20L)

write.table(rbind(top,bot), args[3])