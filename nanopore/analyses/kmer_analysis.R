#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

data <- read.table(args[1], row.names=1, header=T)
outf <- args[2]
outsig <- args[3]
library(stats)

#turn the count of reads into a vector representing the number of times each kmer is seen
#i.e. there will be 100 1's if AAAAA was seen 100 times in the read set
#this lets us sample without replacement
counts <- vector()
for (i in 1:1024) {
    counts <- c(counts, rep(i, times=data[i,]$readCount))
}

#10,000 trials
num_trials <- 10000
#we want each trial to be around 1/25th of the number of kmers seen in the reads
trial_size <- max(round(length(counts)/25), 10000)


#samples from the read population
trial_fn <- function(counts) {
   replicate(num_trials, sample(counts, size=trial_size, replace=F), simplify=F)
}
#runs binomial exact test
test <- function(x, p, n){
    binom.test(x, n, p, alternative="two.sided", conf.level=0.95)$p.value
}
#makes a count table out of the trials, counting non-present things also
tableize <- function(x) {
    tmp <- matrix(table(factor(x, levels=1:1024)))[,1]
    tmp/sum(tmp)
}
count_success <- function(x, real) {
    std <- sd(x)
    #std/10 just to remove float rounding issues and stuff
    length(x[x <= real+std/10 && x >= real-std/10])
}
#generate a trial dataset by sampling from the counts vector without replacement
trials <- trial_fn(counts)
#count the number of times each kmer was found in the trial dataset
trial_table <- sapply(trials, tableize)

p_values <- rep(0, 1024)

for (kmer in 1:1024) {
    s <- count_success(trial_table[kmer,], real=data[kmer,]$refFraction)
    p_values[kmer] <- binom.test(x=s, n=num_trials, p=data[kmer,]$refFraction)$p.value
}

adjusted_p_value <- p.adjust(p_values, method="bonferroni")

finished <- cbind(data, p_values, adjusted_p_value)

write.table(finished, outf)

significant <- finished[finished$adjusted_p_value <= 0.05,]

ordered <- significant[order(significant$logFoldChange),]

top <- head(ordered, n=20L)
bot <- tail(ordered, n=20L)

write.table(rbind(top,bot), outsig)
