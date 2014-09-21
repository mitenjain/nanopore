#!/usr/bin/env Rscript
#
args <- commandArgs(trailingOnly = T)

data <- read.table(args[1], header=T, row.names=1)

library(methods)

##############################################################################
#### R script to
####    Produce Venn Diagrams with 1 to 5 groups
####        an extension on the code from the limma package
####        
#### Written By: Matt Settles
####                Postdoctoral Research Associate
####                Washington State University
####
##############################################################################
####
#### Change Log: 
####    Feb 8, 2008: 
####        formalized code
####    Dec 23, 2008:
####        added mixed type to vennCounts
##############################################################################
####
####    Usage:
####    source("http://bioinfo-mite.crb.wsu.edu/Rcode/Venn.R")
####    can change colors now on 4 way plot
##############################################################################

####################################
## Function ellipse
## Add an elipse to the current plot
"ellipse" <- 
function (center, radius, rotate, 
    segments = 360, add = FALSE, xlab = "", ylab = "", las = par("las"), 
    col = palette()[2], lwd = 2, lty = 1, ...) 
{
    # x' = x cosø + y sinø
    # y' = y cosø - x sinø
    if (!(is.vector(center) && 2 == length(center))) 
        stop("center must be a vector of length 2")
    if (!(is.vector(radius) && 2 == length(radius))) 
        stop("radius must be a vector of length 2")
    
    angles <- (0:segments) * 2 * pi/segments  
    rotate <- rotate*pi/180
    ellipse <- cbind(radius[1] * cos(angles), 
                     radius[2] * sin(angles))
    if(rotate != 0)
        ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate),
                          ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
    ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])  
    if (add) 
        lines(ellipse, col = col, lwd = lwd, lty = lty, ...)
    else    plot(ellipse, type = "l", xlim = c(-4, 4), ylim = c(-4, 4),
            xlab = "", ylab = "", axes = FALSE, col = col, lwd = lwd,
            lty = lty, ...)
}

###################################
## Function vennCounts
## Produce venn table object
"vennCounts" <- 
function (x, include = "both") 
{
    x <- as.matrix(x)
    include <- match.arg(include, c("both", "up", "down","mixed"))
    x <- sign(switch(include, both = abs(x), up = x > 0, down = x < 0, mixed = x))
    nprobes <- nrow(x)
    ncontrasts <- ncol(x)
    names <- colnames(x)
    if (is.null(names)) 
        names <- paste("Group", 1:ncontrasts)
    noutcomes <- 2^ncontrasts
    if ( include == "mixed" ) noutcomes <- 3^ncontrasts
    outcomes <- matrix(0, noutcomes, ncontrasts)
    colnames(outcomes) <- names
    for (j in 1:ncontrasts) {
        if( include == "mixed"){
            outcomes[, j] <- rep(-1:1, times = 3^(j - 1), each = 3^(ncontrasts - j))
        } else {
            outcomes[, j] <- rep(0:1, times = 2^(j - 1), each = 2^(ncontrasts - j))
        }   
    }
    xlist <- list()
    for (i in 1:ncontrasts) {
        if( include == "mixed"){
            xlist[[i]] <- factor(x[, ncontrasts - i + 1], levels = c(-1,0, 1))
        } else {
            xlist[[i]] <- factor(x[, ncontrasts - i + 1], levels = c(0, 1))
        }    
    }
    counts <- as.vector(table(xlist))
    structure(cbind(outcomes, Counts = counts), class = "vennCounts")
}

"vennDiagram" <-
function (object, include = "both", names, mar = rep(1, 4), cex = 1.5, 
    lwd = 2, circle.col, counts.col, show.include, ...) 
{
    if (!is(object, "vennCounts")) {
        if (length(include) > 2) 
            stop("Cannot plot Venn diagram for more than 2 sets of counts")
        if (length(include) == 2) 
            object.2 <- vennCounts(object, include = include[2])
        object <- vennCounts(object, include = include[1])
    }
    else if (length(include == 2)) 
        include <- include[1]
    nsets <- ncol(object) - 1
    if (nsets > 4) 
        stop("Can't plot Venn diagram for more than 4 sets")
    if (missing(names)) 
        names <- colnames(object)[1:nsets]
    counts <- object[, "Counts"]
    if (length(include) == 2) 
    
    if (length(include) == 2 && nsets < 4) 
        counts.2 <- object.2[, "Counts"]
    #Setup colors
     if (missing(circle.col) & nsets == 4)
        circle.col <- c("red","blue","orange","green")
    else 
        circle.col <- par("col")
    if (length(circle.col) < nsets) 
        circle.col <- rep(circle.col, length.out = nsets)
    if (missing(counts.col) & nsets == 4) 
        counts.col <- c("red","blue","orange","green")
     else
          counts.col <- par("col") 
    if (length(counts.col) < length(include) & nsets < 4) 
        counts.col <- rep(counts.col, length.out = length(include))
    else if (length(counts.col) < nsets) 
        counts.col <- rep(counts.col, length.out = nsets)
    
    if (missing(show.include)) 
        show.include <- as.logical(length(include) - 1)
    xcentres <- list(0, c(-1, 1), c(-1, 1, 0), c(-0.2,0.2,-1.05,1.05))[[nsets]]
    ycentres <- list(0, c(0, 0), c(1/sqrt(3), 1/sqrt(3), -2/sqrt(3)),c(.20,.20,-0.35,-0.35))[[nsets]]
    centers <- cbind(xcentres,ycentres)
    r1 <- c(1.5, 1.5, 1.5, 1.5)[nsets]
    r2 <- c(1.5, 1.5, 1.5, 2.7)[nsets]
    radius <- c(r1,r2)
    rotate <- list(0, c(0,0), c(0,0,0), c(-45,45,-45,45))[[nsets]]
    
    xtext <- list(-1.2, c(-1.2, 1.2), c(-1.2, 1.2, 0),c(-3.2,3.2,-3.2,3.2))[[nsets]]
    ytext <- list(1.8, c(1.8, 1.8), c(2.4, 2.4, -3),c(3.2,3.2,-3.2,-3.2))[[nsets]]
    old.par <- par(mar = mar)
    on.exit(par(old.par))
    plot(x = 0, y = 0, type = "n", xlim = c(-4.0, 4.0), ylim = c(-4.0, 
        4.0), xlab = "", ylab = "", axes = FALSE, ...)
    for (circle in 1:nsets) {
        ellipse(centers[circle,],radius,rotate[circle],add=TRUE,
                , lwd = lwd, col = circle.col[circle])
        text(xtext[circle], ytext[circle], names[circle], 
            #pos=ifelse(circle != as.integer(circle/2)*2,4,2),
            offset = 0.5, col = circle.col[circle], cex = cex)
    }
#   rect(-3.5, -3.5, 3.5, 3.5)
 
    switch(nsets, {
        #rect(-3, -2.5, 3, 2.5)
        printing <- function(counts, cex, adj, col, leg) {
            col <- col[1]
            text(2.3, -2.1, counts[1], cex = cex, col = col, 
                adj = adj)
            text(0, 0, counts[2], cex = cex, col = col, adj = adj)
            if (show.include) 
                text(-2.3, -2.1, leg, cex = cex, col = col, adj = adj)
        }
    }, {
        #rect(-3, -2.5, 3, 2.5)
        printing <- function(counts, cex, adj, col, leg) {
            col <- col[1]
            text(2.3, -2.1, counts[1], cex = cex, col = col, 
                adj = adj)
            text(1.5, 0.1, counts[2], cex = cex, col = col, adj = adj)
            text(-1.5, 0.1, counts[3], cex = cex, col = col, 
                adj = adj)
            text(0, 0.1, counts[4], cex = cex, col = col, adj = adj)
            if (show.include) 
                text(-2.3, -2.1, leg, cex = cex, col = col, adj = adj)
        }
    }, {
        #rect(-3, -3.5, 3, 3.3)
        printing <- function(counts, cex, adj, col, leg) {
            col <- col[1]
            text(2.5, -3, counts[1], cex = cex, col = col, adj = adj)
            text(0, -1.7, counts[2], cex = cex, col = col, adj = adj)
            text(1.5, 1, counts[3], cex = cex, col = col, adj = adj)
            text(0.75, -0.35, counts[4], cex = cex, col = col, 
                adj = adj)
            text(-1.5, 1, counts[5], cex = cex, col = col, adj = adj)
            text(-0.75, -0.35, counts[6], cex = cex, col = col, 
                adj = adj)
            text(0, 0.9, counts[7], cex = cex, col = col, adj = adj)
            text(0, 0, counts[8], cex = cex, col = col, adj = adj)
            if (show.include) 
                text(-2.5, -3, leg, cex = cex, col = col, adj = adj)
        }
    }, {
        #rect(-3.5, -3.5, 3.5, 3.5)
        printing <- function(counts, cex, adj, col, leg) {
            text(0, -3,         counts[1], cex = cex, col = "black", adj = adj)
            text(2.5, 0,        counts[2], cex = cex, col = col[4], adj = adj)
             lines(c(2.25,2.75),c(-0.2,-0.2),col=col[4])
             
            text(-2.5, 0,       counts[3], cex = cex, col = col[3], adj = adj)
             lines(c(-2.75,-2.25),c(-0.2,-0.2),col=col[3])
             
            text(0, -2.0,       counts[4], cex = cex, col = "black", adj = adj)
             lines(c(-0.25,0.25),c(-2.2,-2.2),col=col[3])
             lines(c(-0.25,0.25),c(-2.25,-2.25),col=col[4])
            
            text(1.3, 2.1,      counts[5], cex = cex, col = col[2], adj = adj)
             lines(c(1.05,1.55),c(1.9,1.9),col=col[2])
            
            text(1.7, 1.2,      counts[6], cex = cex, col = "black", adj = adj)
             lines(c(1.45,1.95),c(1.0,1.0),col=col[2])
             lines(c(1.45,1.95),c(0.95,0.95),col=col[4])
            
            text(-1.6, -1.1,    counts[7], cex = cex, col = "black", adj = adj)
             lines(c(-1.85,-1.35),c(-1.3,-1.3),col=col[2])
             lines(c(-1.85,-1.35),c(-1.35,-1.35),col=col[3])
             
            text(-0.8, -1.55,   counts[8], cex = cex, col = "black", adj = adj)
             lines(c(-0.55,-1.05),c(-1.75,-1.75),col=col[2])
             lines(c(-0.55,-1.05),c(-1.8,-1.8),col=col[3])
             lines(c(-0.55,-1.05),c(-1.85,-1.85),col=col[4])

            text(-1.3, 2.1,     counts[9], cex = cex, col = col[1], adj = adj)
             lines(c(-1.55,-1.05),c(1.9,1.9),col=col[1])
            
            text(1.6, -1.1,     counts[10], cex = cex, col = "black", adj = adj)
             lines(c(1.85,1.35),c(-1.3,-1.3),col=col[1])
             lines(c(1.85,1.35),c(-1.35,-1.35),col=col[4])

            text(-1.7, 1.2,     counts[11], cex = cex, col = "black", adj = adj)
             lines(c(-1.45,-1.95),c(1.0,1.0),col=col[1])
             lines(c(-1.45,-1.95),c(0.95,0.95),col=col[3])
            
            text(0.8, -1.55,    counts[12], cex = cex, col = "black", adj = adj)
             lines(c(0.55,1.05),c(-1.75,-1.75),col=col[1])
             lines(c(0.55,1.05),c(-1.8,-1.8),col=col[3])
             lines(c(0.55,1.05),c(-1.85,-1.85),col=col[4])
             
            text(0, 1.6,        counts[13], cex = cex, col = "black", adj = adj)
             lines(c(-0.25,0.25),c(1.4,1.4),col=col[1])
             lines(c(-0.25,0.25),c(1.35,1.35),col=col[2])
             
            text(0.9, 0.5,      counts[14], cex = cex, col = "black", adj = adj)
             lines(c(1.15,0.65),c(0.3,0.3),col=col[1])
             lines(c(1.15,0.65),c(0.25,0.25),col=col[2])
             lines(c(1.15,0.65),c(0.2,0.2),col=col[4])
             
            text(-0.9, 0.5,     counts[15], cex = cex, col = "black", adj = adj)
             lines(c(-1.15,-0.65),c(0.3,0.3),col=col[1])
             lines(c(-1.15,-0.65),c(0.25,0.25),col=col[2])
             lines(c(-1.15,-0.65),c(0.2,0.2),col=col[3])
             
            text(0, -0.5,       counts[16], cex = cex, col = "black", adj = adj)
             lines(c(-0.25,0.25),c(-0.7,-0.7),col=col[1])
             lines(c(-0.25,0.25),c(-0.75,-0.75),col=col[2])
             lines(c(-0.25,0.25),c(-0.8,-0.8),col=col[3])
             lines(c(-0.25,0.25),c(-0.85,-0.85),col=col[4])           
      } 
    } )
    adj <- c(0.5, 0.5)
    if (length(include) == 2 & nsets < 4) 
        adj <- c(0.5, 0)
    printing(counts, cex, adj, counts.col, include[1])
    if (length(include) == 2 & nsets < 4) 
        printing(counts.2, cex, c(0.5, 1), counts.col[2], include[2])
    invisible()
}

if (dim(data)[2] > 0) {
    pdf(args[2])
    x <- vennCounts(data[,2:dim(data)[2]])
    vennDiagram(x)
    x[,dim(x)[2]] <- round(100 * x[,dim(x)[2]] / sum(x[,dim(x)[2]]), 2)
    vennDiagram(x)
    dev.off()
}