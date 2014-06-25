#!/usr/bin/env Rscript
if(!require("optparse")) { install.packages("optparse", repos = "http://cran.r-project.org") }
suppressPackageStartupMessages(library("optparse"))

# Last modified: 14-VI-2012

option_list <- list(
                 make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
                         help="Print extra output [default]"),
                 make_option(c("-q", "--quietly"), action="store_false",
                         dest="verbose", help="Print little output"),
                 make_option(c("--data1"), type="character", action="store",
                         help="REQUIRED. First file with counts.",default=NULL,
                         metavar="file_counts"),
		 make_option(c("--data2"), type="character", action="store",
                         help="Optional. Second file with counts.",default=NULL,
                         metavar="file_counts"),
		 make_option(c("--name1"), type="character", action="store",
                         help="Optional. Name for the first sample",default="Sample 1",
                         metavar="sample_name"),
		 make_option(c("--name2"), type="character", action="store",
                         help="Optional. Name for the second sample.",default="Sample 2",
                         metavar="sample_name"),
                 make_option(c("--info"), type="character", action="store",
                         help="Optional. Info file.",default=NULL, metavar="file_info"),
		 make_option(c("--groups"), type="character", action="store",
                         help="Optional. File with groups for the --info file." , default=NULL, metavar="file_groups"),
		 make_option(c("-k", "--threshold"), type="integer", action="store",
                         help="Optional. Threshold for the number of counts.",default=5,
                         metavar="counts_threshold"),
		 make_option(c("-o", "--dirOut"), type="character", action="store",
                         help="Optional. Output folder.",default="./",
                         metavar="folder_output"),
		make_option("--homesrc", type="character", action="store",
                         help="DEVELOPMENTAL. NOT TO BE USED IN THIS VERSION", default="./",
                         metavar="home_src folder")
                 )

opt <- parse_args(OptionParser(option_list=option_list))

#HOMESRC <- "/home/fdgarcia/qualimap/sonia"
HOMESRC <- opt$homesrc
source(file.path(HOMESRC, "qualimapRfunctions.r"))

#datos1 <- "counts1.txt"
datos1 <- opt$data1
if(is.null(datos1)){
	stop("--data1 is a REQUIRED argument")
}
datos2 <- opt$data2  # with the same features than datos1 and in the same order


if(!file.exists(opt$dirOut)){
	dir.create(opt$dirOut, recursive=TRUE)
}

# cutoff for the number of counts to consider a biological feature as detected
k <- opt$threshold 


cm <- 1.5   # cex.main
cl <- 1.2   # cex.lab
ca <- 1     # cex.axis
cc <- 1.2   # cex


nom1 <- opt$name1
nom2 <- opt$name2




if (!is.null(opt$info)){
#file including features ID in "datos1" and biotypes or other classification
  infobio <- opt$info # NO header
  if (!is.null(opt$groups)){
    agru.biotype <- TRUE
		# two columns file: 1) all biotypes in infobio   2) group for biotypes
    file.agru <- opt$groups # NO header
  }else{
    file.agru <- NULL
    agru.biotype <- FALSE
  }
}




###########################################################################


#### DATA

fichero1 <- read.delim(datos1, header = FALSE, sep = "\t", as.is = TRUE)

misdatos1 <- fichero1[,2]
names(misdatos1) <- fichero1[,1]

if (is.null(datos2)) {

    misdatos2 <- NULL

} else {

  fichero2 <- read.delim(datos2, header = FALSE, sep = "\t", as.is = TRUE)

  misdatos2 <- fichero2[,2]
  names(misdatos2) <- fichero2[,1]
}



#### BIOTYPES
if (!is.null(opt$info)) {

  infofile <- read.delim(infobio, header = FALSE, sep = "\t", as.is = TRUE)
  myinfo <- infofile[,2]
  names(myinfo) <- infofile[,1]


  ## Completing files: data & info

  genes1 <- as.character(names(misdatos1))
  genes2 <- as.character(names(misdatos2))
  genesI <- as.character(names(myinfo))

  if (length(genes2) > 0) {
    allgenes <- union(genes1, genes2)
  } else { allgenes <- genes1 }

  
  # Completing myinfo
  unknown <- setdiff(allgenes, genesI) # genes in data without biotype 

  if(length(unknown) > 0) {
    unk <- rep("unknown", length(unknown))
    names(unk) <- unknown
    myinfo <- c(myinfo, unk)   # adding unknown genes (without biotype)    
    write.table(unknown, file = file.path(opt$dirOut, "nogroup.txt"),
                sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE) # file for unknown genes
  }

  genesII <- names(myinfo)


  # Completing my data
  notdet1 <- setdiff(genesII, genes1) # genes1 not in complete info -> counts=0
  if (length(genes2) > 0) {
    notdet2 <- setdiff(genesII, genes2)
    # genes2 not in complete info -> counts=0
  } else { notdet2 <- NULL }
  
  if(length(notdet1) > 0) {
    ceros1 <- rep(0, length(notdet1))
    names(ceros1) <- notdet1
    misdatos1 <- c(misdatos1, ceros1)   # adding genes with 0s to data
    nom <- paste(nom1, "nocounts.txt", sep = ".")
    write.table(notdet1, file = file.path(opt$dirOut, nom),
                sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE) # file for added genes
  }

  if(length(notdet2) > 0) {
    ceros2 <- rep(0, length(notdet2))
    names(ceros2) <- notdet2
    misdatos2 <- c(misdatos2, ceros2)   # adding genes with 0s to data
    nom <- paste(nom2, "nocounts.txt", sep = ".")
    write.table(notdet2, file = file.path(opt$dirOut, nom),
                sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE) # file for added genes
  }


  # OK files
  ok1 <- intersect(genes1, genesII)   # genes1 with counts and biotype
  
  nom <- paste(nom1, "OK.txt", sep = ".")
  write.table(ok1, file = file.path(opt$dirOut, nom),
              sep = "\t", row.names = FALSE,
              col.names = FALSE, quote = FALSE) # file for added genes

  if (!is.null(notdet2)) {

    ok2 <- intersect(genes2, genesII)   # genes2 with counts and biotype

    nom <- paste(nom2, "OK.txt", sep = ".")
    write.table(ok2, file = file.path(opt$dirOut, nom),
                sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE) # file for added genes
  }

 

  ## Sorting files according to genes (all with the same genes)
  misdatos1 <- misdatos1[genesII]

  if(!is.null(misdatos2)) {
    misdatos2 <- misdatos2[genesII]
  }
  

  

  ## Preparing biotype files
  
  if (agru.biotype) { # If you want to group biotypes

    biogroups <- read.delim(file.agru, header = FALSE, sep = "\t", as.is = TRUE)
    biogrupos <- unique(biogroups)
  
  # Creating a list for groupes of biotypes
    bio.agru <- vector("list", length = length(table(biogrupos[,2])))
    names(bio.agru) <- names(table(biogrupos[,2]))

    for (i in 1:length(bio.agru)) {

      bio.agru[[i]] <- biogrupos[is.element(biogrupos[,2],
                                            names(bio.agru)[i]),1]
    }

    biotipus <- bio.agru

    bio.agru2 <- NULL

    for (i in 1:length(bio.agru)) {
      bio.agru2 <- c(bio.agru2, rep(names(bio.agru)[i], length(bio.agru[[i]])))
    }
    names(bio.agru2) <- unlist(bio.agru)

  } else {  # No grouping

    biotipus <- as.list(names(table(myinfo)))
    names(biotipus) <- names(table(myinfo))

    bio.agru2 <- names(table(myinfo))
    names(bio.agru2) <- bio.agru2

  }





###########################################################################
###########################################################################


####  FEATURES DETECTION per BIOTYPE

  if (agru.biotype) {

    mybioagru <- bio.agru2[myinfo]

  } else { mybioagru <- myinfo }

  detect <- misdatos1 > k

  genome <- 100*table(bio.agru2[myinfo])/sum(table(bio.agru2[myinfo]))

  ordre <- order(genome, decreasing = TRUE)

  perdet1 <- genome*table(mybioagru, detect)[,2]/
             table(bio.agru2[myinfo])[rownames(table(mybioagru, detect))]
  perdet2 <- 100*table(mybioagru, detect)[,2]/sum(table(mybioagru, detect)[,2])

  ceros <- rep(0, length(genome))

  biotable <- as.matrix(rbind(genome[ordre], perdet1[ordre],
                              perdet2[ordre], ceros))
  rownames(biotable) <- c("genome", "detectionVSgenome", "detectionVSsample",
                        "ceros")

  if (is.null(misdatos2)) {

    # Saving data in a txt file
    mybiotable <- cbind(colnames(biotable), t(biotable[-4,]))
    colnames(mybiotable)[1] <- "biotype"
    
    write.table(mybiotable,
                file.path(opt$dirOut, "DetectionPerGroup.txt"),
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

    # Plot (1 sample)
    png(filename = file.path(opt$dirOut, "DetectionPerGroup.png"),
        width = 3*480*2, height = 3*480, units = "px", pointsize = 3*12)

    par(mar = c(11, 4, 2, 2), bg = "#e6e6e6")

    barplot(biotable[c(1,3),], main = "Detection per group",
          xlab = "", ylab = "%features", axis.lty = 1, legend = FALSE,
          beside = TRUE, col = c("grey", 2), las = 2,
          ylim = c(0, 100), border = c("grey", 2))

    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
       col = "white")

    barplot(biotable[c(1,3),], main = "Detection per group",
          xlab = "", ylab = "%features", axis.lty = 1, legend = FALSE,
          beside = TRUE, col = c("grey", 2), las = 2, add = TRUE,
          ylim = c(0, 100), border = c("grey", 2))

    barplot(biotable[c(2,4),], main = "Detection per group",
          xlab = "", ylab = "%features", axis.lty = 1, legend = FALSE,
          beside = TRUE, col = c(2, 1), las = 2, density = 30,
          ylim = c(0, 100), border = 2, add = TRUE)

    legend(x = "top", bty = "n", horiz = TRUE,
         fill = c("grey", 2, 2), density = c(NA,30,NA),
         border = c("grey", 2, 2),
         legend = c("% in genome", "detected", "% in sample"))

    garbage <- dev.off()

  } else {

    detect2 <- misdatos2>k

    perdet3 <- genome*table(mybioagru, detect2)[names(genome),2]/
             table(bio.agru2[myinfo])[names(genome)]
    perdet4 <- 100*table(mybioagru, detect2)[,2]/
             sum(table(mybioagru, detect2)[,2])

    biotable2 <- as.matrix(rbind(genome[ordre], perdet3[ordre],
                               perdet4[ordre], ceros))
    rownames(biotable2) <- c("genome", "detectionVSgenome", "detectionVSsample",
                          "ceros")


    # Saving data in a txt file
    mybiotable <- cbind(colnames(biotable), t(biotable[-4,]))
    colnames(mybiotable)[1] <- "biotype"
    mybiotable2 <- t(biotable2[-c(1,4),])
    mybiotable2 <- mybiotable2[rownames(mybiotable),]

    mybiotable12 <- cbind(mybiotable, mybiotable2)
    colnames(mybiotable12)[-c(1,2)] <- apply(cbind(colnames(mybiotable12)[-c(1,2)],
                                                   rep(1:2, each = 2)), 1, paste, collapse = "_")
    
    write.table(mybiotable12,
                file.path(opt$dirOut, "DetectionPerGroup.txt"),
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

    
    # Plot (2 samples)
    png(filename = file.path(opt$dirOut, "DetectionPerGroup.png"),
      width = 3*480*2, height = 3*480*2, units = "px", pointsize = 3*12)

    par(mar = c(11, 4, 2, 2), mfrow = c(2,1), bg = "#e6e6e6")

    # Datos1
    barplot(biotable[c(1,3),],
          main = paste("Detection per group:", nom1),
          xlab = "", ylab = "%features", axis.lty = 1, legend = FALSE,
          beside = TRUE, col = c("grey", 2), las = 2,
          ylim = c(0, 100), border = c("grey", 2))

    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
       col = "white")

    barplot(biotable[c(1,3),],
          main = paste("Detection per group:", nom1),
          xlab = "", ylab = "%features", axis.lty = 1, legend = FALSE,
          beside = TRUE, col = c("grey", 2), las = 2, add = TRUE,
          ylim = c(0, 100), border = c("grey", 2))

    barplot(biotable[c(2,4),],
          main = paste("Detection per group:", nom1),
          xlab = "", ylab = "%features", axis.lty = 1, legend = FALSE,
          beside = TRUE, col = c(2, 1), las = 2, density = 30,
          ylim = c(0, 100), border = 2, add = TRUE)

    legend(x = "top", bty = "n", horiz = TRUE,
         fill = c("grey", 2, 2), density = c(NA,30,NA),
         border = c("grey", 2, 2),
         legend = c("% in genome", "detected", "% in sample"))


    # Datos2
    barplot(biotable2[c(1,3),],
          main = paste("Detection per group:", nom2),
          xlab = "", ylab = "%features", axis.lty = 1, legend = FALSE,
          beside = TRUE, col = c("grey", 4), las = 2,
          ylim = c(0, 100), border = c("grey", 4))

    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
       col = "white")

    barplot(biotable2[c(1,3),],
          main = paste("Detection per group:", nom2),
          xlab = "", ylab = "%features", axis.lty = 1, legend = FALSE,
          beside = TRUE, col = c("grey", 4), las = 2, add = TRUE,
          ylim = c(0, 100), border = c("grey", 4))

    barplot(biotable2[c(2,4),],
          main = paste("Detection per group:", nom2),
          xlab = "", ylab = "%features", axis.lty = 1, legend = FALSE,
          beside = TRUE, col = c(4, 1), las = 2, density = 30,
          ylim = c(0, 100), border = 4, add = TRUE)

    legend(x = "top", bty = "n", horiz = TRUE,
         fill = c("grey", 4, 4), density = c(NA,30,NA),
         border = c("grey", 4, 4),
         legend = c("% in genome", "detected", "% in sample"))

    garbage <- dev.off()

  }
  
}





###########################################################################
###########################################################################


#### SATURATION GLOBAL PLOT

png(filename = paste(opt$dirOut, "GlobalSaturation.png", sep=""),
    width = 3*480, height = 3*480, units = "px", pointsize = 3*12)

par(mar = c(5,6,4,6), bg = "#e6e6e6") 

satur.plot2(datos1 = misdatos1, datos2 = misdatos2,
            myoutput = file.path(opt$dirOut, "GlobalSaturation.txt"),
            tit = "Global Saturation Plot", k = k,
            #yleftlim = NULL, yrightlim = NULL,
            cex.main = cm, cex.lab = cl, cex.axis = ca, cex = cc,
            legend = c(nom1, nom2))

garbage <- dev.off()








############################################################################



#### SATURATION PLOTS per BIOTYPE

if (!is.null(opt$info)) {
  
  satbio.mine <- saturbio.dat(datos1 = misdatos1, datos2 = misdatos2, k = k,
                              infobio = myinfo, biotypes = biotipus)

  mytable1 <- satbio.mine$depth1
  mytable2 <- satbio.mine$depth2

  for (i in 1:length(biotipus)) {
    
    mytable1 <- cbind(mytable1, satbio.mine$cond1[[i]], satbio.mine$newdet1[[i]])    
    mytable2 <- cbind(mytable2, satbio.mine$cond2[[i]], satbio.mine$newdet2[[i]])    

    png(paste(opt$dirOut, names(biotipus)[i],".png", sep = ""),
        width = 3*480, height = 3*480, pointsize = 3*12)

    par(mar = c(5,6,4,6), bg = "#e6e6e6")
    
    titulo <- paste("Number of features with reads >", k)

    saturbio.plot(depth1 = satbio.mine$depth1, depth2 = satbio.mine$depth2,
                  sat1 = satbio.mine$cond1[[i]], sat2 = satbio.mine$cond2[[i]],
                  newdet1 = satbio.mine$newdet1[[i]],
                  newdet2 = satbio.mine$newdet2[[i]],
                  bionum = satbio.mine$bionum[i],
                  main = paste(names(biotipus)[i],
                    " (", satbio.mine$bionum[i], ")", sep = ""),
                  legend = c(nom1, nom2),
                  yleftlim = c(0, satbio.mine$bionum[i]),
                  ylabL = titulo,
                  ylabR = "New detections per million reads",
                  cex = cc, cex.main = cm, cex.lab = cl)
    
    garbage <- dev.off()
    
    }
  
  colnames(mytable1) <- c("depth",
                          apply(cbind(rep(c("detec","newdetec"), length(satbio.mine$cond1)),
                                      rep(names(satbio.mine$cond1), each = 2)),
                                1, paste, collapse = "_"))
  
  write.table(mytable1, file.path(opt$dirOut, "SaturationPerGroup_1.txt"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  if(!is.null(mytable2)) {
    
    colnames(mytable2) <- c("depth",
                            apply(cbind(rep(c("detec","newdetec"), length(satbio.mine$cond2)),
                                        rep(names(satbio.mine$cond2), each = 2)),
                                  1, paste, collapse = "_"))

    write.table(mytable2, file.path(opt$dirOut, "SaturationPerGroup_2.txt"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
    
    }





############################################################################
############################################################################



#### COUNTS DISTRIBUTION -- per BIOTYPE
  
  countdist.mine <- countsbio.dat(misdatos1, misdatos2, k = k, nbox = 1,
                                  infobio = myinfo, biotypes = biotipus)
  
  countlist <- countdist.mine$cond1[[1]]

  for (i in 2:length(countdist.mine$cond1)) {
    countlist <- c(countlist, countdist.mine$cond1[[i]])
    }
  names(countlist) <- names(countdist.mine$cond1)

  ysup <- max(sapply(countlist, quantile, probs = 0.75, na.rm = TRUE),
              na.rm = TRUE)*1.3

  if (is.null(misdatos2)) {

    png(paste(opt$dirOut, "counts_boxplot.png", sep=""),
        width = 3*480, height = 3*480, pointsize = 3*12)

    par(mar = c(11,4,6,2), bg = "#e6e6e6")
    
    } else {
      
      countlist2 <- countdist.mine$cond2[[1]]
      
      for (i in 2:length(countdist.mine$cond2)) {
        countlist2 <- c(countlist2, countdist.mine$cond2[[i]])
        }
      names(countlist2) <- names(countdist.mine$cond2)
      
      ysup2 <- max(sapply(countlist2, quantile, probs = 0.75, na.rm = TRUE),
                 na.rm = TRUE)*1.3
      
      ysup <- max(ysup, ysup2, na.rm = TRUE)
      
      png(paste(opt$dirOut, "counts_boxplot.png", sep=""),
          width = 3*480*2, height = 3*480, pointsize = 3*12)
      
      par(mar = c(11,4,6,2), mfrow = c(1,2), bg = "#e6e6e6")
      
      }
  
  boxplot(countlist[ordre], las = 2,
          main = paste("Counts distribution per group:", nom1),
          ylim = c(k, ysup), col = 2)

  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")

  boxplot(countlist[ordre], las = 2,
          main = paste("Counts distribution per group:", nom1),
          ylim = c(k, ysup), col = 2, add = TRUE,
          ylab = "Read counts of detected features")

  cuantos <- sapply(countlist, length)
  cuantos1 <- which(cuantos == 1)
  sumant <- sapply(countlist, sum, na.rm = TRUE)
  sumant0 <- which(sumant == 0)
  cuantosNA <- intersect(cuantos1, sumant0)
  cuantos[cuantosNA] <- 0
  mtext(cuantos[ordre], 3, at = 1:length(countlist), cex = 1, las = 2)


  # txt file
  write.table(unlist(countlist), file.path(opt$dirOut, "Counts_boxplot_1.txt"),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = TRUE)
  
  
  
  if (!is.null(misdatos2)) {
    
    write.table(unlist(countlist2), file.path(opt$dirOut, "Counts_boxplot_2.txt"),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = TRUE)

    boxplot(countlist2[ordre], las = 2,
            main = paste("Counts distribution per group:", nom2),
            ylim = c(k, ysup), col = 4)

    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
         col = "white")

    boxplot(countlist2[ordre], las = 2,
            main = paste("Counts distribution per group:", nom2),
            ylim = c(k, ysup), col = 4, add = TRUE,
            ylab = "Read counts of detected features")

    cuantos <- sapply(countlist2, length)
    cuantos1 <- which(cuantos == 1)
    sumant <- sapply(countlist2, sum, na.rm = TRUE)
    sumant0 <- which(sumant == 0)
    cuantosNA <- intersect(cuantos1, sumant0)
    cuantos[cuantosNA] <- 0
    mtext(cuantos[ordre], 3, at = 1:length(countlist2), cex = 1, las = 2)
    
    }
  
  garbage <- dev.off()






############################################################################


  
#### COUNTS DISTRIBUTION -- per BIOTYPE + sequencing depth
  
  countbio.mine <- countsbio.dat(misdatos1, misdatos2, k = k,
                                 infobio = myinfo, biotypes = biotipus)

  mytable1 <- c("depth", "--", countbio.mine$depth1)
  mytable2 <- c("depth", "--", countbio.mine$depth2)
  
  mystats <- c("lower_whisker", "lower_quartile", "median", "upper_quartile", "upper_whisker")

  for (i in 1:length(biotipus)) {
    
    mytable1 <- rbind(mytable1,
                      cbind(mystats, names(countbio.mine$cond1)[i],
                            sapply(countbio.mine$cond1[[i]],
                                   function(x) { boxplot(x, plot = FALSE)$stats })))
    mytable1 <- rbind(mytable1, c("number", names(countbio.mine$cond1)[i],
                            sapply(countbio.mine$cond1[[i]],
                                   function(x) { boxplot(x, plot = FALSE)$n })))

    mytable2 <- rbind(mytable2,
                      cbind(mystats, names(countbio.mine$cond2)[i],
                            sapply(countbio.mine$cond2[[i]],
                                   function(x) { boxplot(x, plot = FALSE)$stats })))
    mytable2 <- rbind(mytable2, c("number", names(countbio.mine$cond2)[i],
                            sapply(countbio.mine$cond2[[i]],
                                   function(x) { boxplot(x, plot = FALSE)$n })))
    

    png(paste(opt$dirOut, names(biotipus)[i],"_boxplot.png", sep = ""),
        width = 3*480, height = 3*480, pointsize = 3*12)

    countbio.plot(depth1 = countbio.mine$depth1,
                  depth2 = countbio.mine$depth2,
                  sat1 = countbio.mine$cond1[[i]],
                  sat2 = countbio.mine$cond2[[i]],
                  bionum = countbio.mine$bionum[[i]], las = 1, cex = 1.2,
                  legend = c(nom1, nom2),
                  main = names(biotipus)[i],
                  ylab = "Read counts of detected features",
                  xlab = "Sequencing Depth (million reads)")

    

    garbage <- dev.off()
    }
  
  colnames(mytable1) <- c("stat", "group", paste("depth", 1:length(countbio.mine$depth1), sep = ""))

  write.table(mytable1, file.path(opt$dirOut, "CountsPerGroup_boxplot_1.txt"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  
  if(!is.null(misdatos2)) {
    
    colnames(mytable2) <- c("stat", "group", paste("depth", 1:length(countbio.mine$depth2), sep = ""))

    write.table(mytable2, file.path(opt$dirOut, "CountsPerGroup_boxplot_2.txt"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
    
    }
  
  }






############################################################################


#### CORRELATION PLOT


if (!is.null(misdatos2)) {
  
  png(paste(opt$dirOut, "correlation_plot.png", sep = ""),
      width = 3*480, height = 3*480, pointsize = 3*12)
  
  cor.plot.2D(misdatos1, misdatos2, noplot = 0.001, log.scale = TRUE, 
              xlab = paste("log2(", nom1, "+1)", sep = ""),
              ylab = paste("log2(", nom2, "+1)", sep = ""))
  
  garbage <- dev.off()
  
  }












