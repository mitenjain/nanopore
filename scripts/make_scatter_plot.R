
args <- commandArgs(trailingOnly = T)

data <- read.table(args[1],row.names=1,header=T)

pdf(args[2])

col1 <- rgb(1,0,0,0.7)
col2 <- rgb(0,1,0,0.7)
col3 <- rgb(0,0,1,0.7)
col4 <- rgb(0,0,0,0.7)
cols <- c(rep(col1,3),rep(col2,3),rep(col3,3),rep(col4,3))

plot(data$AvgInsert+data$AvgDelete,data$avgMismatch,pch=c(15,16,17),col=cols,cex=1.5,ylab="Average Mismatch Rate",xlab="Average Indel Rate",main="Mismatch vs. Indel")
legend("topright",legend=rownames(data),pch=c(15,16,17),col=cols,cex=0.7)

plot(data$AvgInsert,data$AvgDelete,pch=c(15,16,17),col=cols,cex=1.5,ylab="Avg Deleletions Per Aligned Read Base",xlab="Average Insertions Per Aligned Read Base",main="Insertions vs. Deletions")
legend("bottomright",legend=rownames(data),pch=c(15,16,17),col=cols,cex=0.7)

dev.off()