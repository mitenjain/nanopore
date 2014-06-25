#installing deps

if(!require("optparse")) { 
    install.packages("optparse", repos = "http://cran.r-project.org") 
}

if(!require("XML")) { 
    install.packages("XML", repos = "http://cran.r-project.org")
}
   
if(!require("Repitools")) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("Repitools")
}

if(!require("Rsamtools")) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("Rsamtools")
}

if(!require("rtracklayer")) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("rtracklayer")
}

