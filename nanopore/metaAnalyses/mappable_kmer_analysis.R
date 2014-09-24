#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

data <- read.table(args[1], row.names=1, header=T)
outf <- args[2]
outsig <- args[3]
library(stats)

mySample <- function(values, size, nElementsPerValue){
  nElementsPerValue <- as.integer(nElementsPerValue)
  if(sum(nElementsPerValue) < size)
    stop("Total number of elements per value is lower than the sample size")
  if(length(values) != length(nElementsPerValue))
    stop("nElementsPerValue must have the same length of values")
  if(any(nElementsPerValue < 0))
    stop("nElementsPerValue cannot contain a negative numbers")

  # remove values having zero elements inside
  nElementsPerValue <- nElementsPerValue[which(nElementsPerValue > 0)]
  values <- values[which(nElementsPerValue > 0)]

  # pre-allocate the result vector
  res <- rep.int(0.0,size)
  for(i in 1:size){
    idx <- sample(1:length(values),size=1,replace=F,prob=nElementsPerValue)
    res[i] <- values[idx]
    # remove sampled value from nElementsPerValue
    nElementsPerValue[idx] <- nElementsPerValue[idx] - 1
    # if zero elements remove also from values
    if(nElementsPerValue[idx] == 0){
      values <- values[-idx]
      nElementsPerValue <- nElementsPerValue[-idx]
    }
  }
  return(res)
}

#10,000 trials
num_trials <- 10000
#we want each trial to be around 1/50th of the number of kmers seen in the reads
trial_size <- round(length(counts)/50)


#samples from the read population
trial_fn <- function(data, num_trials) {
    replicate(num_trials, mySample(1:1024, trial_size, data$readCount), simplify=F)
}
#runs binomial exact test
test <- function(x, p, n){
    binom.test(x, n, p, alternative="two.sided", conf.level=0.95)$p.value
}
#makes a count table out of the trian replicates 
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
trials <- trial_fn(data, num_trials)
#count the number of times each kmer was found in the trial dataset
trial_table <- sapply(trials, tableize)

p_values <- rep(0, 1024)

for (kmer in 1:1024) {
    s <- count_success(trial_table[kmer,], real=data[kmer,]$unmappableFraction)
    p_values[kmer] <- binom.test(x=s, n=num_trials, p=data[kmer,]$unmappableFraction)$p.value
}

adjusted_p_value <- p.adjust(p_values, method="bonferroni")

finished <- cbind(data, p_values, adjusted_p_value)

write.table(finished, outf)

significant <- finished[finished$adjusted_p_value <= 0.05,]

ordered <- significant[order(significant$logFoldChange),]

top <- head(ordered, n=20L)
bot <- tail(ordered, n=20L)

write.table(rbind(top,bot), outsig)
