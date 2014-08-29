#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

data <- read.table(args[1], sep="\t")

if ( dim(data)[1] > 0 ) {

    pdf(args[2])

    hist(t(data), main = "Average Posterior Match Probability", xlab="Probability")

    dev.off()

}