library(lattice)
library(latticeExtra)
pdf("combined_plots_channel_mappability.pdf")#, type="cairo")
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

template1 <- read.table("/hive/users/benedict/nanoporeM13/nanopore/output/analysis_template/experiment_DD_575_R7.3_His_M13_11_21_14_R7.X_2D_V1.9_pass_template.fastq_combinedReference.fa_LastParamsRealignEm/analysis_ChannelMappability/channel_mappability.tsv", row.names=1, header=T)
template2 <- read.table("/hive/users/benedict/nanoporeM13/nanopore/output/analysis_template/experiment_DD_575_R7.3_M13_His_11_17_14_R7.X_2D_V1.9_pass_template.fastq_combinedReference.fa_LastParamsRealignEm/analysis_ChannelMappability/channel_mappability.tsv", row.names=1, header=T)
template3 <- read.table("/hive/users/benedict/nanoporeM13/nanopore/output/analysis_template/experiment_MA_286_R7.3_His_M13_11_18_14_R7.X_2D_V1.9_pass_template.fastq_combinedReference.fa_LastParamsRealignEm/analysis_ChannelMappability/channel_mappability.tsv", row.names=1, header=T)
complement1 <- read.table("/hive/users/benedict/nanoporeM13/nanopore/output/analysis_complement/experiment_DD_575_R7.3_His_M13_11_21_14_R7.X_2D_V1.9_pass_complement.fastq_combinedReference.fa_LastParamsRealignEm/analysis_ChannelMappability/channel_mappability.tsv", row.names=1, header=T)
complement2 <- read.table("/hive/users/benedict/nanoporeM13/nanopore/output/analysis_complement/experiment_DD_575_R7.3_M13_His_11_17_14_R7.X_2D_V1.9_pass_complement.fastq_combinedReference.fa_LastParamsRealignEm/analysis_ChannelMappability/channel_mappability.tsv", row.names=1, header=T)
complement3 <- read.table("/hive/users/benedict/nanoporeM13/nanopore/output/analysis_complement/experiment_MA_286_R7.3_His_M13_11_18_14_R7.X_2D_V1.9_pass_complement.fastq_combinedReference.fa_LastParamsRealignEm/analysis_ChannelMappability/channel_mappability.tsv", row.names=1, header=T)
twod1 <- read.table("/hive/users/benedict/nanoporeM13/nanopore/output/analysis_2D/experiment_DD_575_R7.3_His_M13_11_21_14_R7.X_2D_V1.9_pass_2D.fastq_combinedReference.fa_LastParamsRealignEm/analysis_ChannelMappability/channel_mappability.tsv", row.names=1, header=T)
twod2 <- read.table("/hive/users/benedict/nanoporeM13/nanopore/output/analysis_2D/experiment_DD_575_R7.3_M13_His_11_17_14_R7.X_2D_V1.9_pass_2D.fastq_combinedReference.fa_LastParamsRealignEm/analysis_ChannelMappability/channel_mappability.tsv", row.names=1, header=T)
twod3 <- read.table("/hive/users/benedict/nanoporeM13/nanopore/output/analysis_2D/experiment_MA_286_R7.3_His_M13_11_18_14_R7.X_2D_V1.9_pass_2D.fastq_combinedReference.fa_LastParamsRealignEm/analysis_ChannelMappability/channel_mappability.tsv", row.names=1, header=T)

data <- cbind(template1, template2, template3, complement1, complement2, complement3, twod1, twod2, twod3)

#first all read counts
counts <- list()
largest <- 0
for (i in seq(1, dim(data)[2], 2)) {
  positions <- labels
  for (j in 1:length(positions)) {
    positions[match(c(j), positions)] <- data[j,i]
  }
  counts[[ceiling(i/2)]] <- matrix(positions, nrow=16)
  largest <- max(largest, max(positions))
}

b <- seq(1, largest+10, by=round(largest/30))

rotate <- function(x) t(apply(x, 2, rev))


t1 <- levelplot(rotate(counts[[1]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
t2 <- levelplot(rotate(counts[[2]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
t3 <- levelplot(rotate(counts[[3]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
c1 <- levelplot(rotate(counts[[4]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
c2 <- levelplot(rotate(counts[[5]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
c3 <- levelplot(rotate(counts[[6]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
w1 <- levelplot(rotate(counts[[7]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
w2 <- levelplot(rotate(counts[[8]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
w3 <- levelplot(rotate(counts[[9]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)

plots <- c(t1, t2, t3, c1, c2, c3, w1, w2, w3)
print(plots)

#now it is time to do mappable read counts
counts <- list()
largest <- 0
for (i in seq(2, dim(data)[2], 2)) { #this is the only change
  positions <- labels
  for (j in 1:length(positions)) {
    positions[match(c(j), positions)] <- data[j,i]
  }
  counts[[ceiling(i/2)]] <- matrix(positions, nrow=16)
  largest <- max(largest, max(positions))
}

b <- seq(1, largest+10, by=round(largest/30))

rotate <- function(x) t(apply(x, 2, rev))

t1 <- levelplot(rotate(counts[[1]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
t2 <- levelplot(rotate(counts[[2]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
t3 <- levelplot(rotate(counts[[3]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
c1 <- levelplot(rotate(counts[[4]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
c2 <- levelplot(rotate(counts[[5]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
c3 <- levelplot(rotate(counts[[6]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
w1 <- levelplot(rotate(counts[[7]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
w2 <- levelplot(rotate(counts[[8]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
w3 <- levelplot(rotate(counts[[9]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)

plots <- c(t1, t2, t3, c1, c2, c3, w1, w2, w3)
print(plots)

#finally, lets do the percentages

counts <- list()
largest <- 0
for (i in seq(1, dim(data)[2], 2)) {
  positions <- labels
  for (j in 1:length(positions)) {
    vals <- data[j, i]/data[j+1, i+1]
    vals[is.nan(vals)] <- 0
    positions[match(c(j), positions)] <- vals
  }
  counts[[ceiling(i/2)]] <- matrix(positions, nrow=16)
  largest <- max(largest, max(positions))
}

b <- seq(0, 1, 0.05)

rotate <- function(x) t(apply(x, 2, rev))

t1 <- levelplot(rotate(counts[[1]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
t2 <- levelplot(rotate(counts[[2]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
t3 <- levelplot(rotate(counts[[3]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
c1 <- levelplot(rotate(counts[[4]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
c2 <- levelplot(rotate(counts[[5]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
c3 <- levelplot(rotate(counts[[6]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
w1 <- levelplot(rotate(counts[[7]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
w2 <- levelplot(rotate(counts[[8]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)
w3 <- levelplot(rotate(counts[[9]]), col.regions=colorRampPalette(c("white","red"))(256), at=b)

plots <- c(t1, t2, t3, c1, c2, c3, w1, w2, w3)
print(plots)
dev.off()