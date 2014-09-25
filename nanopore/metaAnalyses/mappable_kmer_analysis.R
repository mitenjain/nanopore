#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

data <- read.table(args[1], row.names=1, header=T)
outf <- args[2]
outsig <- args[3]
library(stats)


#10,000 trials
num_trials <- 5000
#we want each trial to be around 1/50th of the number of kmers seen in the mappable
trial_size <- round(sum(data$mappableCount/50))


counts <- rep(0, 1024)
#samples from the read population
trial_fn <- function(data) {
    tmp <- matrix(table(factor(sample(1024, trial_size, prob=data$mappableFraction, replace=T), levels=1:1024)))[,1]
    round(tmp/sum(tmp),4) == round(data$unmappableFraction,4)
}

for (i in 1:num_trials){
    counts <- counts + trial_fn(data)
}

p_values <- rep(0, 1024)
#loop over each kmer, count the number of successful trials, then run a binomial test on that
for (kmer in 1:1024) {
    p_values[kmer] <- binom.test(counts[kmer], num_trials, data$unmappableFraction[kmer], alternative="two.sided", conf.level=0.95)$p.value
}
#do a bonferroni correction on the pvalues
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
bot <- tail(ordered, n=20L)

write.table(rbind(top,bot), outsig)
