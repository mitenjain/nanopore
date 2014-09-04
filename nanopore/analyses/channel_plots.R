#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

data <- read.table(args[1], row.names=1, header=T)

library(lattice)

#below is hard coded positions on the nanopore
labels <- c(125, 126, 127, 128, 253, 254, 255, 256, 381, 382, 383, 384, 509, 510, 511, 512, 
           121, 122, 123, 124, 249, 250, 251, 252, 377, 378, 379, 380, 505, 506, 507, 508, 
           117, 118, 119, 120, 245, 246, 247, 248, 373, 374, 375, 376, 501, 502, 503, 504, 
           113, 114, 115, 116, 241, 242, 243, 244, 369, 370, 371, 372, 497, 498, 499, 500, 
           109, 110, 111, 112, 237, 238, 239, 240, 365, 366, 367, 368, 493, 494, 495, 496, 
           105, 106, 107, 108, 233, 234, 235, 236, 361, 362, 363, 364, 489, 490, 491, 492, 
           101, 102, 103, 104, 229, 230, 231, 232, 357, 358, 359, 360, 485, 486, 487, 488, 
           97, 98, 99, 100, 225, 226, 227, 228, 353, 354, 355, 356, 481, 482, 483, 484, 
           93, 94, 95, 96, 221, 222, 223, 224, 349, 350, 351, 352, 477, 478, 479, 480, 
           89, 90, 91, 92, 217, 218, 219, 220, 345, 346, 347, 348, 473, 474, 475, 476, 
           85, 86, 87, 88, 213, 214, 215, 216, 341, 342, 343, 344, 469, 470, 471, 472, 
           81, 82, 83, 84, 209, 210, 211, 212, 337, 338, 339, 340, 465, 466, 467, 468, 
           77, 78, 79, 80, 205, 206, 207, 208, 333, 334, 335, 336, 461, 462, 463, 464, 
           73, 74, 75, 76, 201, 202, 203, 204, 329, 330, 331, 332, 457, 458, 459, 460, 
           69, 70, 71, 72, 197, 198, 199, 200, 325, 326, 327, 328, 453, 454, 455, 456, 
           65, 66, 67, 68, 193, 194, 195, 196, 321, 322, 323, 324, 449, 450, 451, 452, 
           61, 62, 63, 64, 189, 190, 191, 192, 317, 318, 319, 320, 445, 446, 447, 448, 
           57, 58, 59, 60, 185, 186, 187, 188, 313, 314, 315, 316, 441, 442, 443, 444, 
           53, 54, 55, 56, 181, 182, 183, 184, 309, 310, 311, 312, 437, 438, 439, 440, 
           49, 50, 51, 52, 177, 178, 179, 180, 305, 306, 307, 308, 433, 434, 435, 436, 
           45, 46, 47, 48, 173, 174, 175, 176, 301, 302, 303, 304, 429, 430, 431, 432, 
           41, 42, 43, 44, 169, 170, 171, 172, 297, 298, 299, 300, 425, 426, 427, 428, 
           37, 38, 39, 40, 165, 166, 167, 168, 293, 294, 295, 296, 421, 422, 423, 424, 
           33, 34, 35, 36, 161, 162, 163, 164, 289, 290, 291, 292, 417, 418, 419, 420, 
           29, 30, 31, 32, 157, 158, 159, 160, 285, 286, 287, 288, 413, 414, 415, 416, 
           25, 26, 27, 28, 153, 154, 155, 156, 281, 282, 283, 284, 409, 410, 411, 412, 
           21, 22, 23, 24, 149, 150, 151, 152, 277, 278, 279, 280, 405, 406, 407, 408, 
           17, 18, 19, 20, 145, 146, 147, 148, 273, 274, 275, 276, 401, 402, 403, 404, 
           13, 14, 15, 16, 141, 142, 143, 144, 269, 270, 271, 272, 397, 398, 399, 400, 
           9, 10, 11, 12, 137, 138, 139, 140, 265, 266, 267, 268, 393, 394, 395, 396, 
           5, 6, 7, 8, 133, 134, 135, 136, 261, 262, 263, 264, 389, 390, 391, 392, 
           1, 2, 3, 4, 129, 130, 131, 132, 257, 258, 259, 260, 385, 386, 387, 388)


if (dim(data)[1] > 1) {
    
    sorted <- data[order(data$ReadCount, decreasing=T),]
    sorted <- t(sorted[sorted$ReadCount > 0,])

    png(args[3], height=3000, width=3000, type="cairo")

    q<- barplot(sorted, main=paste("Sorted Channel Mappability", paste("# Reporting = ", length(sorted[1,]), sep=""), sep="\n"), xlab="Channel", ylab="Read Counts", legend.text=T, xaxt="n", col=c("blue","red"), args.legend=c(cex=3), cex.names=3)
    text(cex=0.5, x=q-.25, y=-1.25, colnames(sorted), xpd=T, srt=65)

    dev.off()

    pdf(args[2])

    sorted.percent <- sorted["MappableReadCount",]/sorted["ReadCount",]
    sorted.percent <- sorted.percent[order(sorted.percent, decreasing=T)]
    sorted.percent <- sorted.percent[sorted.percent > 0]
    sorted.percent <- sorted.percent[!is.na(sorted.percent)]

    q<- barplot(sorted.percent, main="Sorted Channel Percent Mappability", xlab="Channel", ylab="Read Counts", xaxt="n")
    text(cex=0.27, x=q-.25,y=-0.005, names(sorted.percent), xpd=T, srt=45)


    #do linear regression
    reg <- lm(sorted["MappableReadCount",]~sorted["ReadCount",])
    #plot scatterplot
    plot(sorted["MappableReadCount",]~sorted["ReadCount",], pch=20, col="blue", xlab="Total Read Count", ylab="Mappable Read Count", main="Mappable vs Total Reads\nReporting Channels Only")
    #add regression line
    abline(reg)
    #add R2
    legend("topleft",legend=c(paste("R^2 = ", round(summary.lm(reg)$adj.r.squared,4))))

    barplot(t(data)[2,]/(t(data)[1,]+t(data)[2,])*100, main="% Mappable Reads Per Channel", xlab="Channel", ylab="% Mappable")

    par(mfrow=(c(1,2)))

    barplot(t(data)[1,], main="Total Read Counts", xlab="Channel", ylab="Read Counts")
    barplot(t(data)[2,], main="Mappable Read Counts", xlab="Channel", ylab="Read Counts")

    dev.off()

    png(args[4], height=1000, width=1000, type="cairo")

    is.nan.data.frame <- function(x){
        do.call(cbind, lapply(x, is.nan))}
    #need to find % mapped, but then replace the NaNs with 0s
    mapped <- data$MappableReadCount/data$ReadCount
    mapped[is.nan(mapped)] <- 0

    positions <- labels

    for (i in 1:length(positions) ) {
        positions[match(c(i), positions)] <- mapped[i]
    }
    positions <- matrix(positions, nrow=16)

    rotate <- function(x) t(apply(x, 2, rev))

    p <- levelplot(rotate(positions), main ="MinION Channel Percent Read Mappability Layout", col.regions=colorRampPalette(c("white","red"))(256))

    print(p)

    #do it again except with total # of reads
    reads <- data$ReadCount
    reads[is.nan(reads)] <- 0
    positions <- labels

    for (i in 1:length(positions) ) {
        positions[match(c(i), positions)] <- reads[i]
    }
    positions <- matrix(positions, nrow=16)

    rotate <- function(x) t(apply(x, 2, rev))

    p <- levelplot(rotate(positions), main ="MinION Channel Number of Reads", col.regions=colorRampPalette(c("white","red"))(256))

    print(p)

    dev.off()

}
