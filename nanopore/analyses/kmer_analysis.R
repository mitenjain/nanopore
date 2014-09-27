#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

data <- read.table(args[1], row.names=1, header=T)
outf <- args[2]
outsig <- args[3]
library(stats)


num_trials <- 5000
#trial size is a max of 10,000
if (sum(data$refCount/50) > 10000) {
    trial_size <- 10000
} else {
    trial_size <- sum(data$refCount)/50
}


trial_fn <- function(d, t) {
    matrix(table(factor(sample(1024, t, prob=d, replace=T), levels=1:1024)))[,1]
}

ref <- replicate(num_trials, trial_fn(d=data$refFraction, t=trial_size))

read <- replicate(num_trials, trial_fn(d=data$readFraction, t=trial_size))

p_values <- rep(0, 1024)

for (i in 1:1024) {
    p_values[i] <- ks.test(ref[i,], read[i,])$p.value
}

adjusted_p_value <- p.adjust(p_values, method="bonferroni")
#combine original data frame with two new vectors
finished <- cbind(data, p_values, adjusted_p_value)
#write full dataset
write.table(finished, outf)
#find significant hits
significant <- finished[finished$adjusted_p_value <= 0.05,]
#sort the significant hits by fold change
ordered <- significant[order(significant$logFoldChange),]
#report the top 20 and bottom 20 significant hits
top <- head(ordered, n=20L)
bot <- apply(tail(ordered, n=20L), 2, rev)

write.table(rbind(top,bot), outsig)