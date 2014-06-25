####################################################################
###################      Qualimap R Functions   ####################
####################################################################

# By Sonia Tarazona, Fernando Garcia-Alcalde and Konstantin Okonechnikov
#



# This function creates input data based on the description file

load.counts.data <- function(input.desc) {
    
    input.data <- read.table(input.desc, sep= "\t", stringsAsFactors=FALSE)
    
    sample.names <- character()
    conditions <- character()
    counts <- data.frame()
    gene.names <- NULL
    for (i in 1:nrow(input.data)) {
        sample.names <- c(sample.names, input.data[i,1])
        conditions <- c(conditions, input.data[i,2])
        path <- input.data[i,3]
        data.column <- as.numeric(input.data[i,4])
        sample.counts <- read.table(path, sep = "\t")
        if (is.null(gene.names)) {
            gene.names <- as.character(sample.counts[,1])
            head(gene.names)
            counts <- sample.counts[data.column]
        } else {
            counts <- cbind(counts, sample.counts[,data.column])
        }
    }
    head(counts)
    head(gene.names)
    rownames(counts) <- gene.names
    colnames(counts) <- sample.names
    
    attr(counts, "factors") <- as.factor(conditions)

    counts
    
    
}


#****************************************************************************#


## Function to intersect multiple sets

int.mult <- function(lista, todos = NULL) {
  
  if(is.null(todos)) {
    todos <- unlist(lista)
  }

  comunes <- todos

  for(i in 1:length(lista)) {
    comunes <- intersect(comunes, lista[[i]])
  }

  comunes
}




#****************************************************************************#




## To count number of non-zero elements

noceros <- function (x, num = TRUE, k = 0) {
  
  nn <- length(which(x > k))
  
  if (num) {
    nn
    
  } else {
    if(nn > 0) { which(x > k) } else { NULL }
  }
}




#***************************************************************************#





## Saturation Plot (GLOBAL)


satur.plot <- function (datos1, datos2, ylim = NULL, k = 0, tit = "Saturation",
                        cex.main = cex.main, cex.lab = cex.lab,
                        cex.axis = cex.axis, cex = cex,
                        legend = c(deparse(substitute(datos1)),
                                   deparse(substitute(datos2)))) {


# For datos1
  n1 <- ncol(as.matrix(datos1))

  
  if (n1 > 1) {

    muestra1 <- as.list(1:n1)
    for (i in 2:n1) {
      
      combi1 <- combn(n1, i, simplify = FALSE)

      if( length(combi1) > 20 ) {
        sub20 <- sample(1:length(combi1), size = 20, replace = FALSE)
        combi1 <- combi1[sub20]
      }
      
      muestra1 <- append(muestra1, combi1)
    }

    varias1 <- vector("list", length = length(muestra1))
    names(varias1) <- sapply(muestra1, function(x) {
      paste("C1.", x, collapse = "", sep = "")})

    for (i in 1:length(muestra1)) {
      varias1[[i]] <- apply(as.matrix(datos1[,muestra1[[i]]]), 1, sum)
    }

    satura1 <- data.frame("muestra" = names(varias1),
                          "seq.depth" = sapply(varias1, sum),
                          "noceros" = sapply(varias1, noceros, k = k))
  }


 
  if (n1 == 1) {
    
    total1 <- sum(datos1)
    satura1 <- NULL
    
    for (i in 1:9) {     # 10%, 20%, ..., 90% reads (apart 100% is calculated)     
      muestra1 <- rmultinom(10, size = round(total1*i/10,0), prob = datos1)
      detec1 <- mean(apply(muestra1, 2, noceros, k = k))
      satura1 <- rbind(satura1, c(round(total1*i/10,0), detec1))
    }
    satura1 <- rbind(satura1, c(total1, noceros(datos1, k = k)))
    colnames(satura1) <- c("seq.depth", "noceros")
    satura1 <- as.data.frame(satura1)
  }




# For datos2

  if (!is.null(datos2)) {

    n2 <- ncol(as.matrix(datos2))

    if (n2 > 1) {

      if (n1 == n2) {
        muestra2 <- muestra1

      } else {

        muestra2 <- as.list(1:n2)

        for (i in 2:n2) {

          combi2 <- combn(n2, i, simplify = FALSE)

          if (length(combi2) > 20) {
            sub20 <- sample(1:length(combi2), size = 20, replace = FALSE)
            combi2 <- combi2[sub20]
          }

          muestra2 <- append(muestra2, combi2)

        }
      }

      varias2 <- vector("list", length = length(muestra2))
      names(varias2) <- sapply(muestra2, function(x) {
        paste("C2.", x, collapse = "", sep = "")})

      for (i in 1:length(muestra2)) {
        varias2[[i]] <- apply(as.matrix(datos2[,muestra2[[i]]]), 1, sum)
      }

      satura2 <- data.frame("muestra" = names(varias2),
                            "seq.depth" = sapply(varias2, sum),
                            "noceros" = sapply(varias2, noceros, k = k))
    }


    if (n2 == 1) {

      total2 <- sum(datos2)
      satura2 <- NULL

      for (i in 1:9) {   # 10%, 20%, ..., 90% reads (apart 100% is calculated)

        muestra2 <- rmultinom(10, size = round(total2*i/10,0), prob = datos2)
        detec2 <- mean(apply(muestra2, 2, noceros, k = k))
        satura2 <- rbind(satura2, c(round(total2*i/10,0), detec2))
      }

      satura2 <- rbind(satura2, c(total2, noceros(datos2, k = k)))
      colnames(satura2) <- c("seq.depth", "noceros")
      satura2 <- as.data.frame(satura2)
    }


    if (is.null(ylim)) {
      ylim <- c(0, NROW(datos1))
    }

    SS1 <- range(satura1$seq.depth/10^6)
    SS2 <- range(satura2$seq.depth/10^6)

    xM <- max(SS1[2], SS2[2])
    xm <- min(SS1[1], SS2[1])

    par(xpd=TRUE, mar=par()$mar+c(0,0,2,0), bg = "#e6e6e6")

    plot(satura1$seq.depth/10^6, satura1$noceros, pch = 16, col = 2,
         ylim = ylim, xlim = c(xm, xM), main = tit, type = "o",
         xlab = "Sequencing Depth (million reads)",
         ylab = paste("Number of features with reads >", k), cex = 1.5, las = 1,
         cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis)

    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
         col = "white")

    points(satura1$seq.depth/10^6, satura1$noceros, pch = 16, col = 2, 
	   type = "o", cex = 1.5)

    points(satura2$seq.depth/10^6, satura2$noceros, pch = 16, col = 4,
           type = "o", cex = 1.5)

    legend(x = (xm+xM)/2, y = ylim[2]+0.15*diff(ylim),
           legend = legend, text.col = c(2,4), bty = "n", xjust = 0.5,
           lty = 1, lwd = 2, col = c(2, 4), cex = cex, horiz = TRUE)

    par(mar=c(5, 4, 4, 2) + 0.1)
    

  } else {

    if (is.null(ylim)) {
      ylim <- c(0, NROW(datos1))
    }

    xlim <- range(satura1$seq.depth/10^6)

    par(bg = "#e6e6e6")

    plot(satura1$seq.depth/10^6, satura1$noceros, pch = 16, col = 2,
         ylim = ylim, xlim = xlim, main = tit, type = "o",
         xlab = "Sequencing Depth (million reads)",
         ylab = paste("Number of features with reads >", k), cex = 1.5, las = 1,
         cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis)

    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
         col = "white")

    points(satura1$seq.depth/10^6, satura1$noceros, pch = 16, col = 2,
	   type = "o", cex = 1.5)
    
  }

}





#***************************************************************************#







###############################################################################
###############################################################################


## Data for saturation plot with different biotypes


saturbio.dat <- function (datos1, datos2, k = 0, infobio, biotypes, nsim = 5) {

  # infobio = vector containing biotype for each gene in "datos"
  # biotypes = list containing groups of biotypes to be studied

  biog <- lapply(biotypes, function(x) { which(is.element(infobio, x)) })

  satura1 <- satura2 <- vector("list", length = length(biotypes))
  names(satura1) <- names(satura2) <- names(biotypes)
  
  newdet1 <- newdet2 <- vector("list", length = length(biotypes))
  names(newdet1) <- names(newdet2) <- names(biotypes)


  # datos1
  n1 <- NCOL(datos1)


  if (n1 > 1) {  # when replicates are available in condition 1

    varias1 <- vector("list", length = n1) # counts for each sequencing depth
    names(varias1) <- paste(1:n1, "rep", sep = "")

    for(i in 1:n1) {
      
      muestra1 <- combn(n1, i, simplify = FALSE)

      if (length(muestra1) > 20) {
        sub20 <- sample(1:length(muestra1), size = 20, replace = FALSE)
        muestra1 <- muestra1[sub20]
      }

      sumrepli <- NULL

      for (com in 1:length(muestra1)) {
        sumrepli <- cbind(sumrepli, rowSums(as.matrix(datos1[,muestra1[[com]]])))
      }

      varias1[[i]] <- as.matrix(sumrepli)
      
    }
  }


  if (n1 == 1) {   # simulating replicates for datos1

    varias1 <- vector("list", length = nsim) # counts for each sequencing depth
    names(varias1) <- paste("dp", 1:nsim, sep = "")
    
    total1 <- sum(datos1)    

    for (i in 1:(nsim-1)) {    # total1*i/nbox reads (100% is calculated apart)

      muestra1 <- rmultinom(10, size = round(total1*i/nsim,0), prob = datos1)
      
      varias1[[i]] <- muestra1
    }

    varias1[[nsim]] <- as.matrix(datos1)
  }

  seq.depth1 <- sapply(varias1, function(x) { mean(colSums(x)) })

  # computing saturation for each biotype (datos1)
  for (j in 1:length(satura1))  {

    conbio1 <- lapply(varias1, function(x) { as.matrix(x[biog[[j]],]) })

    satura1[[j]] <- sapply(conbio1,
                           function(x) { mean(apply(x, 2, noceros, k = k)) })
  }
  

  ## computing detection increasing per million reads: condition 1
  for (j in 1:length(newdet1))  {

    puntos1 <- data.frame("x" = seq.depth1, "y" = satura1[[j]])

    pendi <- NULL

    for(i in 2:nrow(puntos1)) {

      nuevo <- max(0, (puntos1$y[i]-puntos1$y[i-1])/
                   (puntos1$x[i]-puntos1$x[i-1]))
      
      pendi <- c(pendi, nuevo)                 
    }

    pendimil1 <- c(NA, pendi)*1000000

    newdet1[[j]] <- pendimil1

  }


  
  
  #### datos2
  if (!is.null(datos2)) {

    n2 <- NCOL(datos2)

    if (n2 > 1) {   # when replicates are available in condition 2

      varias2 <- vector("list", length = n2) # counts for each sequencing depth
      names(varias2) <- paste(1:n2, "rep", sep = "")

      for(i in 1:n2) {

        muestra2 <- combn(n2, i, simplify = FALSE)

        if (length(muestra2) > 20) {
          sub20 <- sample(1:length(muestra2), size = 20, replace = FALSE)
          muestra2 <- muestra2[sub20]
        }

        sumrepli <- NULL

        for (com in 1:length(muestra2)) {
          sumrepli <- cbind(sumrepli,
                            rowSums(as.matrix(datos2[,muestra2[[com]]])))
        }

        varias2[[i]] <- as.matrix(sumrepli)

      }
    }


    if (n2 == 1) {  # replicates have to be simulated

      varias2 <- vector("list", length = nsim) # counts for each seq depth
      names(varias2) <- paste("dp", 1:nsim, sep = "")

      total2 <- sum(datos2)

      for (i in 1:(nsim-1)) {  # total1*i/nbox reads (100% is calculated apart)

        muestra2 <- rmultinom(10, size = round(total2*i/nsim,0), prob = datos2)

        varias2[[i]] <- muestra2
      }

      varias2[[nsim]] <- as.matrix(datos2)
    }

    seq.depth2 <- sapply(varias2, function(x) { mean(colSums(x)) })

    # computing saturation for each biotype (datos2)
    for (j in 1:length(satura2))  {

      conbio2 <- lapply(varias2, function(x) { as.matrix(x[biog[[j]]]) })

      satura2[[j]] <- sapply(conbio2,
                             function(x) { mean(apply(x, 2, noceros, k = k)) })
    }


    # computing detection increasing per million reads: condition 2
    for (j in 1:length(newdet2))  {

      puntos2 <- data.frame("x" = seq.depth2, "y" = satura2[[j]])

      pendi <- NULL

      for(i in 2:nrow(puntos2)) {

        nuevo <- max(0, (puntos2$y[i]-puntos2$y[i-1])/
                     (puntos2$x[i]-puntos2$x[i-1]))
        
        pendi <- c(pendi, nuevo)
      }

      pendimil2 <- c(NA, pendi)*1000000

      newdet2[[j]] <- pendimil2
    }

    satura <- list("cond1" = satura1, "cond2" = satura2,
                   "bionum" = sapply(biog, length),
                   "depth1" = seq.depth1, "depth2" = seq.depth2,
                   "newdet1" = newdet1, "newdet2" = newdet2)
  

  } else {

    satura <- list("cond1" = satura1, "cond2" = satura2,
                   "bionum" = sapply(biog, length),
                   "depth1" = seq.depth1, "depth2" = NULL,
                   "newdet1" = newdet1, "newdet2" = NULL)
  }
  

  ### Results

  satura

}





#*****************************************************************************#






## Saturation plot with different biotypes

saturbio.plot <- function(depth1, depth2 = NULL, sat1, sat2 = NULL, bionum,
                          newdet1 = NULL, newdet2 = NULL,
                          xlim = NULL, yleftlim = NULL, yrightlim = NULL,
                          main = "Saturation per biotype", lwdL = 2, lwdR = 20,
                          legend = c("sample1", "sample2"), 
                          ylabL = "Number of detected features",
                          ylabR = "New detections per million reads",
                          cex.main = 1, cex.lab = 1, cex.axis = 1, cex = 1) {


  if (bionum > 0) {

    if (!is.null(depth2))  {
  
  # yleftlim for plot and plot.y2
      if (is.null(yleftlim)) {
        yleftlim <- c(min(c(sat1, sat2)), max(c(sat1, sat2)))
      }

  # xlim for plot
      if (is.null(xlim)) {
        SS1 <- range(depth1/10^6)
        SS2 <- range(depth2/10^6)

        xM <- max(SS1[2], SS2[2])
        xm <- min(SS1[1], SS2[1])

        xlim <- c(xm, xM)
      }

      percen1 <- round(100*max(sat1)/bionum, 1)
      percen2 <- round(100*max(sat2)/bionum, 1)


      if(is.null(newdet1)) {

    # PLOT for detections
	#par(bg = "#e6e6e6")

        plot(depth1/10^6, sat1, pch = 16, col = 2, ylim = yleftlim, lwd = lwdL,
             xlim = xlim, main = main, type = "o",
             xlab = "Sequencing Depth (million reads)",
             ylab = ylabL, cex.main = cex.main, cex.lab = cex.lab, las = 1,
             cex.axis = cex.axis, cex = 1.5)

	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "white")

	points(depth1/10^6, sat1, pch = 16, col = 2, lwd = lwdL,
               type = "o", cex = 1.5)

        points(depth2/10^6, sat2, pch = 16, col = 4, lwd = lwdL,
               type = "o", cex = 1.5)

        legend("top",
               legend = paste(legend, ": ",
                        c(percen1, percen2), "% detected", sep = ""),
               text.col = c(2,4), bty = "n", lty = 1, lwd = lwdL,
               col = c(2, 4), cex = 0.9)

      } else {

    # yrightlim for plot.y2
        if (is.null(yrightlim)) {
          maxi <- max(na.omit(c(newdet1, newdet2)))
          if (maxi < 10) {
            yrightlim <- c(0, max(10,maxi))
          } else {
            yrightlim <- c(0, maxi*1.1)
          }
        }

    # PLOT with 2 axis
        plot.y2(x = depth1/10^6, yright = newdet1, yleft = sat1,
                type = c("h", "o"), lwd = c(lwdR, lwdL),
                xlab = "Sequencing depth (million reads)", xlim = xlim,
                yrightlim = yrightlim, yleftlim = yleftlim, 
                yylab = c(ylabR, ylabL), pch = c(1,19), col = c("pink",2),
                main = main, x2 = depth2/10^6, yright2 = newdet2, yleft2 = sat2,
                col2 = c("lightblue1",4), cex.main = cex.main,
                cex.lab = cex.lab, cex.axis = cex.axis, cex2 = c(1.5,1.5))
  
       legend.text <- c(paste(legend[1], "(left)"),
                        paste(legend[2], "(left)"),
                        paste(legend[1], "(right)"),
                        paste(legend[2], "(right)"),
                        paste(percen1,"% detected"),                        
                        paste(percen2,"% detected"))

        legend.pch <- c(16, 16, 15, 15, 1, 1)

        legend.col <- c(2, 4, "pink", "lightblue1", 0, 0)
        
        legend("top",
               legend = legend.text, pch = legend.pch,
               col = legend.col, ncol = 3, bty = "n", cex = 0.9)
      }
      
    } else {



      
      ### Only for 1 sample

      # yleftlim for plot and plot.y2
      if (is.null(yleftlim)) {
        yleftlim <- range(sat1)
      }

      # xlim for plot
      if (is.null(xlim)) {
        xlim <- range(depth1/10^6)
      }

      percen1 <- round(100*max(sat1)/bionum, 1)

      
      if(is.null(newdet1)) {

    # PLOT for detections
	#par(bg = "#e6e6e6")

        plot(depth1/10^6, sat1, pch = 16, col = 4, ylim = yleftlim, lwd = lwdL,
             xlim = xlim, main = main, type = "o",
             xlab = "Sequencing Depth (million reads)",
             ylab = ylabL, cex.main = cex.main, cex.lab = cex.lab, las = 1,
             cex.axis = cex.axis, cex = 1.5)

	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
             col = "white")

	points(depth1/10^6, sat1, pch = 16, col = 2, lwd = lwdL,
	       type = "o", cex = 1.5)

        legend("top",
               legend = paste(legend, ": ", percen1, "% detected", sep = ""),
               text.col = 2, bty = "n", lty = 1, lwd = lwdL,
               col = 4, cex = 0.9)

      } else {

    # yrightlim for plot.y2
        if (is.null(yrightlim)) {
          maxi <- max(na.omit(newdet1))
          if (maxi < 10) {
            yrightlim <- c(0, max(10,maxi))
          } else {
            yrightlim <- c(0, maxi*1.1)
          }
        }

    # PLOT with 2 axis
        plot.y2(x = depth1/10^6, yright = newdet1, yleft = sat1,
                type = c("h", "o"),
                lwd = c(lwdR, lwdL),
                xlab = "Sequencing depth (million reads)", xlim = xlim,
                yrightlim = yrightlim, yleftlim = yleftlim,
                yylab = c(ylabR, ylabL), pch = c(1,19), 
                main = main, col = c("pink",2), cex.main = cex.main,
                cex.lab = cex.lab, cex.axis = cex.axis, cex2 = c(1.5,1.5))

        legend.text <- c("Left axis", "Right axis",
                         paste(percen1,"%detected"))

        legend.pch = c(16, 15, 1)

        legend.col = c(2, "pink", 0)
        
        legend("top",
               legend = legend.text, pch = legend.pch,
               col = legend.col, ncol = 3, cex = 0.9)
      }
    }

  } else {

    print("Biotype not found in this dataset")
  }
}







#*****************************************************************************#







## Counts for detected genes Plot according to BIOTYPES


countsbio.dat <- function (datos1, datos2, k = 0, infobio, biotypes, nbox = 5)
  {

  # infobio = vector containing biotype for each gene in "datos"
  # biotypes = list containing groups of biotypes to be studied
  # nbox = number of different depths to be plotted (only for simulated data)
 

  biog <- lapply(biotypes, function(x) { which(is.element(infobio, x)) })
  # which genes belong to each biotype

  satura1 <- satura2 <- vector("list", length = length(biotypes))
  names(satura1) <- names(satura2) <- names(biotypes)

  
  ### DATOS 1
  n1 <- NCOL(datos1)
  
  ## replicates available for condition 1
  if (n1 > 1) {

    varias1 <- vector("list", length = n1) # counts for each sequencing depth
    names(varias1) <- paste(1:n1, "rep", sep = "")
    
    for(i in 1:n1) {

      muestra1 <- combn(n1, i, simplify = FALSE)

      if (length(muestra1) > 20) {
        sub20 <- sample(1:length(muestra1), size = 20, replace = FALSE)
        muestra1 <- muestra1[sub20]
      }

      sumrepli <- NULL

      for (com in 1:length(muestra1)) {
        sumrepli <- cbind(sumrepli, rowSums(as.matrix(datos1[,muestra1[[com]]])))
      }

      varias1[[i]] <- rowMeans(sumrepli)
    }
  }


  ## replicates simulated for condition 1
  if (n1 == 1) {  # replicates have to be simulated

    varias1 <- vector("list", length = nbox) # counts for each sequencing depth
    names(varias1) <- paste("dp", 1:nbox, sep = "")
    
    total1 <- sum(datos1)

    if (nbox > 1) {

      for (i in 1:(nbox-1)) { # total1*i/nbox reads (100% is calculated apart)

        muestra1 <- rmultinom(10, size = round(total1*i/nbox,0), prob = datos1)

        varias1[[i]] <- rowMeans(muestra1)
      }
    }

    varias1[[nbox]] <- datos1
  }

  seq.depth1 <- sapply(varias1, sum)  # sequencing depth for each new sample

 
  
  ## selecting detected genes for each biotype (datos1)
  for (j in 1:length(satura1))  {

    satura1[[j]] <- vector("list", length = length(varias1))
    names(satura1[[j]]) <- names(varias1)

    # selecting genes in bioclass j for each sample
    conbio1 <- lapply(varias1, function(x) { x[biog[[j]]] })  

    for (i in 1:length(conbio1)) {
            
      # selecting the genes with counts > k
      noK <- noceros(conbio1[[i]], k = k, num = FALSE)

      if (is.null(noK)) {
        satura1[[j]][[i]] <- NA
      } else {
        satura1[[j]][[i]] <- conbio1[[i]][noK]
      }
    }
  }



  #### DATOS 2

  if (!is.null(datos2)) {

    n2 <- NCOL(datos2)
  
  ## replicates available for condition 2
    if (n2 > 1) {

      varias2 <- vector("list", length = n2)
      names(varias2) <- paste(1:n2, "rep", sep = "")

      for (i in 1:n2) {

        muestra2 <- combn(n2, i, simplify = FALSE)

        if (length(muestra2) > 20) {
          sub20 <- sample(1:length(muestra2), size = 20, replace = FALSE)
          muestra2 <- muestra2[sub20]
        }

        sumrepli2 <- NULL

        for (com in 1:length(muestra2)) {
          sumrepli2 <- cbind(sumrepli2,
                             rowSums(as.matrix(datos2[,muestra2[[com]]])))
        }

        varias2[[i]] <- rowMeans(sumrepli2)
      }
    }    

    
  ## replicates simulated for condition 2
    if (n2 == 1) {  # replicates have to be simulated

      varias2 <- vector("list", length = nbox) # counts for each seq depth
      names(varias2) <- paste("dp", 1:nbox, sep = "")

      total2 <- sum(datos2)

      if (nbox > 1) {

        for (i in 1:(nbox-1)) { # total1*i/nbox reads (100% is calculated apart)

          muestra2 <- rmultinom(10,
                                size = round(total2*i/nbox,0), prob = datos2)

          varias2[[i]] <- rowMeans(muestra2)
        }
      }

      varias2[[nbox]] <- datos2
    }

    seq.depth2 <- sapply(varias2, sum)  # sequencing depth for each new sample

  
  ## selecting detected genes for each biotype (datos2)
    for (j in 1:length(satura2))  {

      satura2[[j]] <- vector("list", length = length(varias2))
      names(satura2[[j]]) <- names(varias2)

    # selecting genes in bioclass j for each sample
      conbio2 <- lapply(varias2, function(x) { x[biog[[j]]] })

      for (i in 1:length(conbio2)) {
            
      # selecting the genes with counts > k
        noK <- noceros(conbio2[[i]], k = k, num = FALSE)

        if (is.null(noK)) {
          satura2[[j]][[i]] <- NA
        } else {
          satura2[[j]][[i]] <- conbio2[[i]][noK]
        }
      }
    }

    ## results
    satura <- list("cond1" = satura1, "cond2" = satura2,
                   "bionum" = sapply(biog, length),
                   "depth1" = seq.depth1, "depth2" = seq.depth2)
  } else {

    ## results
    satura <- list("cond1" = satura1, "cond2" = satura2,
                   "bionum" = sapply(biog, length),
                   "depth1" = seq.depth1, "depth2" = NULL)

  }

  satura

}




#****************************************************************************#





## PLOT: Mean length for detected genes Plot according to BIOTYPES

countbio.plot <- function (depth1, depth2, sat1, sat2, bionum, 
                           legend = c("sample1", "sample2"), main = "",
                           ylab = "# counts of detected features",
                           xlab = "Sequencing Depth", ylim = NULL, las = 0,
                           cex.main = 1, cex.lab = 1, cex.axis = 1, cex = 1)  {

  if (bionum == 0) {

    print("Biotype not found in this dataset")
    
  } else {

    if (!is.null(depth2)) {

    # ylim for plot
      if (is.null(ylim)) {
        
        #maxi <- sapply(c(sat1,sat2), max, na.rm = TRUE)
        q75 <- sapply(c(sat1,sat2), quantile, probs = 0.75, na.rm = TRUE)

        #limis <- q75 + 0.2*(maxi - q75)
        
        #ylim <- c(0, max(limis, na.rm = TRUE))
        ylim <- c(0, suppressWarnings(max(q75, na.rm = TRUE))*1.3)
      }

    # boxplots data

      depth <- c(depth1, depth2)
      depth <- round(depth/10^6, 1)
      d.sort <- sort(depth)
      d.order <- order(depth)

      n1 <- length(depth1)
      n2 <- length(depth2)

      databox <- vector("list", length = n1+n2)
      names(databox) <- d.sort

      colo <- NULL

      for (i in 1:(n1+n2)) {

        if ( d.order[i] > n1 ) {

          dd <- d.order[i] - n1
          databox[[i]] <- sat2[[dd]]
          colo <- c(colo, 4)

        } else {

          dd <- d.order[i]
          databox[[i]] <- sat1[[dd]]
          colo <- c(colo, 2)
        }
      }
    } else {

      #### only DATOS1

          # ylim for plot
      if (is.null(ylim)) {
        
        #maxi <- sapply(sat1, max, na.rm = TRUE)
        q75 <- sapply(sat1, quantile, probs = 0.75, na.rm = TRUE)

        #limis <- q75 + 0.001*(maxi - q75)

        #ylim <- c(0, max(limis, na.rm = TRUE))
        ylim <- c(0, suppressWarnings(max(q75, na.rm = TRUE))*1.3)
      }

    # boxplots data

      depth <- round(depth1/10^6,1)
      
      n1 <- length(depth1)      

      databox <- sat1

      colo <- rep(2, n1)

    }

    numNA <- sapply(databox, function(x) { length(which(is.na(x))) })
    cuantos <- sapply(databox, length)

    if (sum(cuantos-numNA) != 0) {

      # BOXPLOT
      par(bg = "#e6e6e6")

      boxplot(databox, col = colo, ylim = ylim, main = main,
              type = "b", xlab = xlab, ylab = ylab, names = depth,
              cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis)

      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")

      boxplot(databox, col = colo, ylim = ylim, main = main,
              type = "b", xlab = xlab, ylab = ylab, names = depth,
              cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, add = TRUE)


      if (!is.null(depth2)) {

        legend("top",
               legend = legend, text.col = c(2,4), bty = "o",               
               bg = "white", fill = c(2, 4), cex = cex, horiz = TRUE)
      }

      cuantos1 <- which(cuantos == 1)
      sumant <- sapply(databox, sum, na.rm = TRUE)
      sumant0 <- which(sumant == 0)
      cuantosNA <- intersect(cuantos1, sumant0)

      cuantos[cuantosNA] <- 0

      mtext(cuantos, 3, at = 1:length(databox), cex = 0.7*cex, las = las)

    } else {
      mismar <- par()$mar
      par(mar = c(2,2,4,2), bg = "#e6e6e6")
      plot(1,1, col = "white", xlab = "", ylab = "", xaxt = "n", yaxt = "n", main  = main)
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
      text(1,1, "No counts found for this class")
      par(mar = mismar)
      #print("No features detected for this biotype")
    }
  }
}





#*************************************************************************#



#############################################################################
##############   Plot with 2 different Y axis (left and right)   ############
#############################################################################


# By Ajay Shah (taken from [R] Plot 2 time series with different y axes
# (left and right),
# in https://stat.ethz.ch/pipermail/r-help/2004-March/047775.html) 

# Modified by: Sonia Tarazona

### PARAMETERS (default):
# x: data to be drawn on X-axis
# yright: data to be drawn on Y right axis
# yleft: data to be drawn on Y left axis
# yrightlim (range(yright, na.rm = TRUE)): ylim for rigth Y-axis
# yleftlim (range(yleft, na.rm = TRUE)): ylim for left Y-axis
# xlab (NULL): Label for X-axis
# yylab (c("","")): Labels for right and left Y-axis
# pch (c(1,2)): Type of symbol for rigth and left data
# col (c(1,2)): Color for rigth and left data
# linky (TRUE): If TRUE, points are connected by lines.
# smooth (0): Friedman's super smoothing
# lwds (1): Line width for smoothed line
# length (10): Number of tick-marks to be drawn on axis
# ...: Other graphical parameters to be added by user (such as main, font, etc.)
###



plot.y2 <- function(x, yright, yleft, yrightlim = range(yright, na.rm = TRUE),
                    yleftlim = range(yleft, na.rm = TRUE),
                    xlim = range(x, na.rm = TRUE),
                    xlab = NULL, yylab = c("",""), lwd = c(2,2),
                    pch = c(1,2), col = c(1,2), type = c("o","o"),
                    linky = TRUE, smooth = 0, cex2 = c(1,1),
                    lwds = 1, length = 10, ...,
                    x2 = NULL, yright2 = NULL, yleft2 = NULL, col2 = c(3,4))
{
  #par(mar = c(5,2,4+2,2), oma = c(0,3,0,3))
  #par(mar = c(5,4,4+2,4))

  ## Plotting RIGHT axis data

  #par(bg = "#e6e6e6")

  plot(x, yright, ylim = yrightlim, axes = FALSE, ylab = "", xlab = xlab,
       xlim = xlim,
       pch = pch[1], type = type[1], lwd = lwd[1], col = col[1], ...)
  
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
       col = "white")

  points(x, yright, pch = pch[1], type = type[1], lwd = lwd[1],
         col = col[1], cex = cex2[1], ...)
 
  axis(4, pretty(yrightlim, length), col = 1, col.axis = 1, las = 1)

  if (!is.null(yright2)) {
    points(x2, yright2, type = type[1], pch = pch[1], lwd = lwd[1],
           col = col2[1], cex = cex2[1], ...)
  }
  
  #if (linky) lines(x, yright, col = col[1], ...)
  
  if (smooth != 0) lines(supsmu(x, yright, span = smooth), col = col[1],
        lwd = lwds, ...)
  
  if(yylab[1]=="") {
    mtext(deparse(substitute(yright)), side = 4, outer = FALSE, line = 4,
          col = 1, las = 3, ...)
  } else {
    mtext(yylab[1], side = 4, outer = FALSE, line = 4, col = 1, las = 3, ...)
  }
  

  par(new = TRUE)#, xpd = TRUE, mar = c(5,2,4+2,2), oma = c(0,3,0,3))

  ## Plotting LEFT axis data
  
  plot(x, yleft, ylim = yleftlim, axes = FALSE, ylab = "" , xlab = xlab,
       xlim = xlim, cex = cex2[2],
       pch = pch[2], type = type[2], lwd = lwd[2], col = col[2], ...)
  
  box()
  
  axis(2, pretty(yleftlim, length), col = 1, col.axis = 1, las = 1)

  if (!is.null(yleft2)) {
    points(x2, yleft2, type = type[2], pch = pch[2], lwd = lwd[2],
           col = col2[2], cex = cex2[2], ...)
  }
  

  #if (linky) lines(x, yleft, col = col[2], ...)
  
  if (smooth != 0) lines(supsmu(x, yleft, span = smooth), col = col[2],
        lwd=lwds, ...)
  
  if(yylab[2] == "") {
    mtext(deparse(substitute(yleft)), side = 2, outer = FALSE, line = 4,
          col = 1, las = 3,...)
  } else {
    mtext(yylab[2], side = 2, outer = FALSE, line = 4, col = 1, las = 3,...)
  }
  
  
  ## X-axis
  axis(1, at = pretty(xlim, length))
  
   
}






###############################################################################
###############################################################################





## Global Saturation Plot including new detection per million reads


satur.plot2 <- function (datos1, datos2, myoutput, yleftlim = NULL, yrightlim = NULL,
                         k = 0, tit = "Saturation", cex.main = cex.main,
                         cex.lab = cex.lab, cex.axis = cex.axis, cex = cex,
                         legend = c(deparse(substitute(datos1)),
                                    deparse(substitute(datos2)))) {


# For datos1
  n1 <- ncol(as.matrix(datos1))

  
  if (n1 > 1) {

    muestra1 <- as.list(1:n1)
    for (i in 2:n1) {
      
      combi1 <- combn(n1, i, simplify = FALSE)

      if( length(combi1) > 20 ) {
        sub20 <- sample(1:length(combi1), size = 20, replace = FALSE)
        combi1 <- combi1[sub20]
      }
      
      muestra1 <- append(muestra1, combi1)
    }

    varias1 <- vector("list", length = length(muestra1))
    names(varias1) <- sapply(muestra1, function(x) {
      paste("C1.", x, collapse = "", sep = "")})

    for (i in 1:length(muestra1)) {
      varias1[[i]] <- apply(as.matrix(datos1[,muestra1[[i]]]), 1, sum)
    }

    satura1 <- data.frame("muestra" = names(varias1),
                          "seq.depth" = sapply(varias1, sum),
                          "noceros" = sapply(varias1, noceros, k = k))
  }


 
  if (n1 == 1) {
    
    total1 <- sum(datos1)
    satura1 <- NULL
    
    for (i in 1:9) {     # 10%, 20%, ..., 90% reads (apart 100% is calculated)     
      muestra1 <- rmultinom(10, size = round(total1*i/10,0), prob = datos1)
      detec1 <- mean(apply(muestra1, 2, noceros, k = k))
      satura1 <- rbind(satura1, c(round(total1*i/10,0), detec1))
    }
    satura1 <- rbind(satura1, c(total1, noceros(datos1, k = k)))
    colnames(satura1) <- c("seq.depth", "noceros")
    satura1 <- as.data.frame(satura1)
  }


  # new detections (sample 1)

  puntos1 <- data.frame("x" = satura1$seq.depth, "y" = satura1$noceros)

  pendi <- NULL

  for (i in 2:nrow(puntos1)) {

    nuevo <- max(0, (puntos1$y[i]-puntos1$y[i-1])/
                   (puntos1$x[i]-puntos1$x[i-1]))

    pendi <- c(pendi, nuevo)
  }

  newdet1 <- c(NA, pendi)*1000000





# For datos2

  if (!is.null(datos2)) {

    n2 <- ncol(as.matrix(datos2))

    if (n2 > 1) {

      if (n1 == n2) {
        muestra2 <- muestra1

      } else {

        muestra2 <- as.list(1:n2)

        for (i in 2:n2) {

          combi2 <- combn(n2, i, simplify = FALSE)

          if (length(combi2) > 20) {
            sub20 <- sample(1:length(combi2), size = 20, replace = FALSE)
            combi2 <- combi2[sub20]
          }

          muestra2 <- append(muestra2, combi2)

        }
      }

      varias2 <- vector("list", length = length(muestra2))
      names(varias2) <- sapply(muestra2, function(x) {
        paste("C2.", x, collapse = "", sep = "")})

      for (i in 1:length(muestra2)) {
        varias2[[i]] <- apply(as.matrix(datos2[,muestra2[[i]]]), 1, sum)
      }

      satura2 <- data.frame("muestra" = names(varias2),
                            "seq.depth" = sapply(varias2, sum),
                            "noceros" = sapply(varias2, noceros, k = k))
    }


    if (n2 == 1) {

      total2 <- sum(datos2)
      satura2 <- NULL

      for (i in 1:9) {   # 10%, 20%, ..., 90% reads (apart 100% is calculated)

        muestra2 <- rmultinom(10, size = round(total2*i/10,0), prob = datos2)
        detec2 <- mean(apply(muestra2, 2, noceros, k = k))
        satura2 <- rbind(satura2, c(round(total2*i/10,0), detec2))
      }

      satura2 <- rbind(satura2, c(total2, noceros(datos2, k = k)))
      colnames(satura2) <- c("seq.depth", "noceros")
      satura2 <- as.data.frame(satura2)
    }


    # new detections (sample 2)

    puntos2 <- data.frame("x" = satura2$seq.depth, "y" = satura2$noceros)

    pendi <- NULL

    for (i in 2:nrow(puntos2)) {

      nuevo <- max(0, (puntos2$y[i]-puntos2$y[i-1])/
                      (puntos2$x[i]-puntos2$x[i-1]))

      pendi <- c(pendi, nuevo)
    }

    newdet2 <- c(NA, pendi)*1000000




    ## PLOT LIMITS

    # yleftlim for plot.y2
    if (is.null(yleftlim)) { 
      yleftlim <- c(0, NROW(datos1))
    }


    # xlim
    SS1 <- range(satura1$seq.depth/10^6)
    SS2 <- range(satura2$seq.depth/10^6)

    xM <- max(SS1[2], SS2[2])
    xm <- min(SS1[1], SS2[1])


    # yrightlim for plot.y2
    if (is.null(yrightlim)) {
      maxi <- max(na.omit(c(newdet1, newdet2)))
      if (maxi < 10) {
        yrightlim <- c(0, max(10,maxi))
      } else {
        yrightlim <- c(0, maxi*1.1)
      }
    }

    

    ## PLOT for SAMPLES 1 and 2

    #par(xpd=TRUE, mar=par()$mar+c(0,0,2,0))

    plot.y2(x = satura1$seq.depth/10^6, yright = newdet1,
            yleft = satura1$noceros,
            type = c("h", "o"), lwd = c(20,2), pch = c(1, 19), 
            yrightlim = yrightlim, yleftlim = yleftlim, 
            xlim = c(xm, xM), main = tit, 
            xlab = "Sequencing Depth (million reads)", col = c("pink", 2),
            yylab = c("New detections per million reads",
                      paste("Number of features with reads >", k)),
            cex2 = c(1.5,1.5),
            cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis,
            x2 = satura2$seq.depth/10^6, yright2 = newdet2,
            yleft2 = satura2$noceros, col2 = c("lightblue1", 4))

    legend.text <- c(paste(legend[1], "(left axis)"),
                     paste(legend[2], "(left axis)"),
                     paste(legend[1], "(right axis)"),
                     paste(legend[2], "(right axis)"))

    legend.pch <- c(16, 16, 15, 15)

    legend.col <- c(2, 4, "pink", "lightblue1")

    legend("top",
           legend = legend.text, pch = legend.pch,
           col = legend.col, ncol = 2, bty = "n")
 
    par(mar=c(5, 4, 4, 2) + 0.1)


    # txt file for sample 1 and 2
    mytxt <- cbind(satura1$seq.depth, satura1$noceros, newdet1,
                   satura2$seq.depth, satura2$noceros, newdet2)
    colnames(mytxt) <- c("depth1", "detec1", "newdetec1",
                         "depth2", "detec2", "newdetec2")

    write.table(mytxt, file = myoutput, quote = FALSE,
                col.names = TRUE, row.names = FALSE, sep = "\t")


  } else {  ## 1 sample

        ## PLOT LIMITS

    # yleftlim for plot.y2
    if (is.null(yleftlim)) { 
      yleftlim <- c(0, NROW(datos1))
    }


    # xlim
    xlim <- range(satura1$seq.depth/10^6)
 

    # yrightlim for plot.y2
    if (is.null(yrightlim)) {
      maxi <- max(na.omit(newdet1))
      if (maxi < 10) {
        yrightlim <- c(0, max(10,maxi))
      } else {
        yrightlim <- c(0, maxi*1.1)
      }
    }
    

    ## PLOT for SAMPLE 1       
    plot.y2(x = satura1$seq.depth/10^6, yright = newdet1,
            yleft = satura1$noceros,
            type = c("h", "o"), lwd = c(20,2), pch = c(1, 19), 
            yrightlim = yrightlim, yleftlim = yleftlim,
            xlim = xlim, main = tit, 
            xlab = "Sequencing Depth (million reads)", col = c("pink", 2),
            yylab = c("New detections per million reads",
                      paste("Number of features with reads >", k)),
            cex = c(1.5,1.5),
            cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis)

    legend.text <- c("Left axis", "Right axis")

    legend.pch <- c(16, 15)

    legend.col <- c(2, "pink")

    legend("top",
           legend = legend.text, pch = legend.pch,
           col = legend.col, ncol = 2, bty = "n")
    
 
    par(mar=c(5, 4, 4, 2) + 0.1)


    # txt file for sample 1
    mytxt <- cbind(satura1$seq.depth, satura1$noceros, newdet1)                   
    colnames(mytxt) <- c("depth1", "detec1", "newdetec1")

    write.table(mytxt, file = myoutput, quote = FALSE,
                col.names = TRUE, row.names = FALSE, sep = "\t")


  }

}





#***************************************************************************#







## Correlation plot between two samples


cor.plot <- function(datos1, datos2, log.scale = FALSE, ...) {

  if (log.scale) {

    datos1 <- log2(datos1 + 1)
    datos2 <- log2(datos2 + 1)

  }
  
  plot(datos1, datos2,...)

  coefLR <- summary(lm(datos2~datos1))$coefficients[,"Estimate"]

  abline(a = coefLR[1], b = coefLR[2], col = 2, lwd = 2)

  mycor <- cor(datos1, datos2)

  myR2 <- round(100*mycor^2,2)
 
  my.y <- min(datos2)+0.1*diff(range(datos2))

  text(max(datos1), my.y, paste("y =", round(coefLR[1],3), "+", round(coefLR[2],3), "x"), col = 2, adj = 1)
  text(max(datos1), my.y-0.05*diff(range(datos2)), paste("R2 =", myR2, "%"), col = 2, adj = 1)

}





## Function filled.contour with legend

filled.contour <- 
function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
    length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
    levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
    col = color.palette(length(levels) - 1), plot.title, plot.axes, 
    key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
    axes = TRUE, frame.plot = axes, key.border = NA, ...) 
{
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
        stop("increasing 'x' and 'y' values expected")
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
    par(las = las)
    mar <- mar.orig
    mar[4L] <- mar[2L]
    mar[2L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
        yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1L], col = col,
         border=key.border)
    if (missing(key.axes)) {
        if (axes) 
            axis(4)
    }
    else key.axes
    box()
    if (!missing(key.title)) 
        key.title
    mar <- mar.orig
    mar[4L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1L || ncol(z) <= 1L) 
        stop("no proper 'z' matrix specified")
    if (!is.double(z)) 
        storage.mode(z) <- "double"
    
    if (R.Version()$major > 2) {    
        .filled.contour(as.double(x), as.double(y), z, as.double(levels), col = col)
    } else {
        .Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels), 
                                                        col = col))
    }

    if (missing(plot.axes)) {
        if (axes) {
            title(main = "", xlab = "", ylab = "")
            Axis(x, side = 1)
            Axis(y, side = 2)
        }
    }
    else plot.axes
    if (frame.plot) 
        box()
    if (missing(plot.title)) 
        title(...)
    else plot.title
    invisible()
}





## Contour plot for correlation

#library(colorRamps)

cor.plot.2D <- function(datos1, datos2, noplot = 0.01,
                        log.scale = TRUE,...) {

  nozeros <- which(rowSums(cbind(datos1, datos2)) > 0)

  datos1 <- datos1[nozeros]
  datos2 <- datos2[nozeros]

  mycor <- round(cor(datos1, datos2),3)
  
  if (log.scale) {

    datos1 <- log2(datos1 + 1)
    datos2 <- log2(datos2 + 1)

  }

  library(MASS)

  limx <- quantile(datos1, c(noplot, 1-noplot))
  limy <- quantile(datos2, c(noplot, 1-noplot))

  d <- kde2d(datos1, datos2, n=100, lims=c(limx, limy))

  filled.contour(d, color.palette = topo.colors,
                 main = paste("Pearson's correlation coefficient =", mycor),
                 n=100, key.title = title(main="Density"),
                 key.axes = axis(4),...)

}
