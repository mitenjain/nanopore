#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
library(lattice)

data <- read.table(args[1])

#find # of algorithms, coverages and proportions heldout used
algorithms <- length(unique(data[,2]))
coverage <- length(unique(data[,4]))
heldout <- length(unique(data[,3]))

#set it to 4 plots per page because we are assuming 4 heldouts
#if there are not 4 heldouts then this won't work properly

#loop over every algorithm block, which is heldout*coverage*2
#this corresponds to each row of plots
for (i in seq(1, algorithms*heldout*coverage*2, heldout*coverage*2)) {
    #open a pdf for each algorithm, put into the folder for this readtype/mapper combination
    pdf(paste(args[2], data[,2][i], args[3], sep=""))
    par(mfrow=c(2,2))
    #loop over every fpr/tpr block, which happens every coverage*2
    #this corresponds to one plot in a row
    for (j in seq(i, i+heldout*coverage*2-heldout-coverage, coverage*2)) {
        #pull out these fprs and tprs (which alternate down the whole file)
        fprs <- data[seq(j+2, j+coverage*2-1, 2),]
        tprs <- data[seq(j+3, j+coverage*2, 2),]
        #find the coverages for this plot *should always be the same*
        coverages <- tprs[,4]
        #find this trials algorithm *should be the same for each row*
        algorithm <- tprs[,2][1]
        #find the proportion heldout for this plot, turn into a text percentage
        held_out <- paste(formatC(100*round(tprs[,3][1],6), format="f", digits=2), "%", sep="")
        #remove the non-data columns in preparation for plotting
        tprs <- tprs[,-(1:4)]
        fprs <- fprs[,-(1:4)]
        #plot and draw legend
        matplot(t(fprs), t(tprs), type="l", col=c(1:length(coverages)), main=paste("VariantCaller:\n", algorithm, "\nProportionHeldOut: ", held_out, sep=""), cex.main=0.8, cex.axis=0.7, xlab="False Positive Rate", ylab="True Positive Rate")
        legend("topleft", legend=coverages, col=c(1:length(coverages)), cex=0.8, pch="-", title="Coverage")
    }
    dev.off()
}


