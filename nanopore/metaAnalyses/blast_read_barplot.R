blast <- read.table(args[4])
lengths <- c(length(unmapped[,1]), length(mapped[,1]), length(blast[,1]))
barplot(lengths, names.arg=c("Unmapped","Mapped","BLAST"), col=c("red", "blue", "yellow"), main="# Of Mapped and Unmapped Reads\nAnd BLAST hits")
