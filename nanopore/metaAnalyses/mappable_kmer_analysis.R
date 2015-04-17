#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

data <- read.table(args[1], row.names=1, header=T)
outf <- args[2]
outsig <- args[3]
outplot <- args[4]
library(stats)
library(lattice)


num_trials <- 1000
trial_size <- 5000


trial_fn <- function(d, t) {
    matrix(table(factor(sample(1024, t, prob=d, replace=T), levels=1:1024)))[,1]
}
if (sum(data$mappableFraction) != 0 & sum(data$unmappableFraction) != 0) {
    mappable <- replicate(num_trials, trial_fn(d=data$mappableFraction, t=trial_size))

    unmappable <- replicate(num_trials, trial_fn(d=data$unmappableFraction, t=trial_size))

    p_values <- rep(0, 1024)

    for (i in 1:1024) {
        p_values[i] <- ks.test(mappable[i,], unmappable[i,])$p.value
    }

    adjusted_p_value <- p.adjust(p_values, method="bonferroni")
    #combine original data frame with two new vectors
    finished <- cbind(data, p_values, adjusted_p_value)
    #write full dataset
    write.table(finished, outf)
    #find significant hits
    significant <- finished[finished$adjusted_p_value <= 0.05,]
    pdf(outplot)
    xyplot(finished$adjusted_p_value~finished$logFoldChange, xlab="Log Fold Change", ylab="Adjusted P Value", main="Unmappable vs. Mappable Volcano Plot")
    dev.off()
    #sort the significant hits by fold change
    ordered <- significant[order(significant$logFoldChange),]
    #report the top 20 and bottom 20 significant hits
    top <- head(ordered, n=20L)
    bot <- apply(tail(ordered, n=20L), 2, rev)

    write.table(rbind(top,bot), outsig)
}