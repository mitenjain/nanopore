#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
library(lattice)

inFile <- args[1]
pdf(args[2])

depthFile = read.delim(inFile, sep="\t", header=F)

outliers <- length(boxplot.stats(c(unlist(unlist(depthFile[3]), dpois(unlist(depthFile[3]), mean(unlist(depthFile[3]))))))$out)
proportion <- round((outliers / length(unlist(depthFile[3]))) * 100, 2)

par(mfrow <- c(1,1))
obs_cov <- unlist(depthFile[3])
avg <- mean(unlist(depthFile[3]))
dev <- sd(unlist(depthFile[3]))

samples <- length(unlist(depthFile[3]))

xmax <- max(obs_cov) + 1000

x <- seq(0,max(obs_cov))
y <- dnorm(x,mean=avg, sd=dev)

hist(obs_cov, freq=F, main=paste("Coverage distribution\n Number of under-represented sites = ", outliers, "(", proportion, "%)"), xlab="Coverage", ylab="Density", xlim=c(0, xmax), col = "green")
lines(x,y, lwd=1, col="red")
legend("topleft",legend=c("Expected","Observed"),col=c("red","green"),pch="-")

dev.off()

pdf(args[3])
par(mfrow <- c(1,1))

plot(smooth.spline(unlist(depthFile[2]), unlist(depthFile[3])), main="Coverage across reference", xlab="Position across reference", ylab="Coverage", type = "l", col = "red")
cov_change <- diff(unlist(depthFile[3]))

dev.off()

