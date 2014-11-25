#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
library(lattice)

inFile <- args[1]
pdf(args[2])

library(boot)
depthFile = read.delim(inFile, sep="\t", header=F)

obs_cov <- unlist(depthFile[3])
avg <- mean(unlist(depthFile[3]))
dev <- sd(unlist(depthFile[3]))

samples <- length(obs_cov)

outlier_pos <- boxplot.stats(c(unlist(obs_cov, dpois(obs_cov, avg))))$out

outliers <- length(outlier_pos)
proportion <- round((outliers / length(unlist(depthFile[3]))) * 100, 2)
xmax <- max(obs_cov) + 1000

par(mfrow <- c(1,1))
hist(obs_cov, freq=F, main=paste("Coverage distribution\n Number of under-represented sites = ", outliers, "(", proportion, "%)"), xlab="Coverage", ylab="Density", xlim=c(0, xmax), col = "green", breaks="FD")
library(MASS)
d <- fitdistr(obs_cov, "Weibull")
lines(1:max(obs_cov) + 1000,dweibull(1:max(obs_cov) + 1000,shape=d$estimate[1],scale=d$estimate[2]), col="red")

legend("topleft",legend=c("Expected","Observed"),col=c("red","green"),pch="-")

dev.off()

pdf(args[3])
par(mfrow <- c(1,1))
plot(spline(unlist(depthFile[2]), obs_cov/max(obs_cov)), main="Coverage across reference", xlab="Position across reference", ylab="Normalized Coverage", type = "l", col = "red")

dev.off()

write(paste(names(outlier_pos), outlier_pos), args[4])
