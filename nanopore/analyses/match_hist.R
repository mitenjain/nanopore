#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

data <- read.table(args[1], sep="\t")

pdf(args[2])

hist(t(data), main = "Average Posterior Match Probability", xlab="Probability")

dev.off()
