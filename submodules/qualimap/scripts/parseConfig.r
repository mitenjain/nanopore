#if(!require("XML",quietly=T)) { install.packages("XML", repos = "http://cran.r-project.org") }

suppressPackageStartupMessages(library("XML",quietly=T))

checkChild <- function(node,name, fileXML="fileConfig"){
		if(is.null(node[[name]])){ stop(paste("Node <", as.character(xmlName(node)), "> does not have a child <", name, "> in ", fileXML,sep="")) }

}

#hasChild <- function(node, name, fileXML="fileConfig") {
#  if (is.null(node[[name]])) {
#    return TRUE;
#  } else {
#    return FALSE;
#  }
#}

parseReplicate <- function(replicate){
	checkChild(replicate, "medips")
	medips = as.character(xmlValue(replicate[["medips"]]))
	input = as.character(xmlValue(replicate[["input"]]))
	name = as.character(xmlValue(replicate[["name"]]))
	d <- data.frame(medips=medips,input=input,name=name)
	d
}

parseSample <- function(sample,ord.sample=0){
	#print(xmlName(sample))
	sample.name <- if(!is.null(xmlAttrs(sample)["name"])) as.character(xmlAttrs(sample)["name"]) else paste("","")
	#print(sample.name)
	checkChild(sample, "replicate")
	parsedMedips <- c()
	parsedInput <- c()
	parsedName <- c()
	i <- 1
	for(curRep in xmlChildren(sample)){
		if(xmlName(curRep) == "replicate"){
			d <- parseReplicate(curRep)
			#print(as.character(d$medips))
			#print(as.character(d$input))
			parsedMedips <- c(parsedMedips,as.character(d$medips))
			parsedInput <- c(parsedInput,as.character(d$input))
			parsedName <- if(!is.na(d$name)) c(parsedName,as.character(d$name)) else c(parsedName, paste(sample.name,as.character(i),sep="_"))
			i <- i + 1
		}
	}
	parsed <- data.frame(medips=parsedMedips, input=parsedInput, sample.name=sample.name, replicate.name=parsedName)
	parsed
}

parseLocation <- function(pos){
	ch.pos <- xmlChildren(pos)

	#threshold <- if (!is.null(ch.ma[["threshold"]])) xmlValue(ch.ma[["threshold"]]) else NA
	up <- if (!is.null(ch.pos[["up"]])) xmlValue(ch.pos[["up"]])  else 2000
	down <- if (!is.null(ch.pos[["down"]])) xmlValue(ch.pos[["down"]]) else 500
	freq <- if (!is.null(ch.pos[["freq"]])) xmlValue(ch.pos[["freq"]]) else 100

	d <- data.frame(up=up, down=down, freq=freq)
	d
}

parseMicroarray <- function(ma){
	#print(xmlName(ma))
	checkChild(ma, "file")
	ch.ma <- xmlChildren(ma)

        parsedTh <- c()
        for(th in xmlChildren(ma)){
                if(xmlName(th)=="threshold"){
                        parsedTh <- c(parsedTh, xmlValue(th))
                }

        }

	#threshold <- if (!is.null(ch.ma[["threshold"]])) xmlValue(ch.ma[["threshold"]]) else NA
	threshold <- if (!is.null(ch.ma[["threshold"]])) parsedTh else NA
	d <- data.frame(file=xmlValue(ch.ma[["file"]]), threshold=threshold)
	d
}

parseGeneSelection <- function(node){
	checkChild(node, "file")
	ch.node <-xmlChildren(node)

	column <- if (!is.null(ch.node[["column"]])) xmlValue(ch.node[["column"]]) else 1
	
	d <- data.frame(file=xmlValue(ch.node[["file"]]), column=column)
	d
}

parseDirOut <- function(node){
	d <- data.frame(dirOut=xmlValue(node))
	d
}

parseRegions <- function(node){
	d <- data.frame(regions=xmlValue(node))
	d
}


parseFragment <- function(node){
	d <- data.frame(fragment=xmlValue(node))
	d
}

parseExpID <- function(node){
	d <- data.frame(expID=xmlValue(node))
	d
}

parseClusters <- function(cls){
	checkChild(cls,"num")
	parsedNums <- c()
	for(cl in xmlChildren(cls)){
		if(xmlName(cl)=="num"){
			parsedNums <- c(parsedNums, xmlValue(cl))
		}
	
	}

	d <- data.frame(num=parsedNums)
	#print(d)
	d
}

parseConfig <- function(fileXML){
	if(!file.exists(fileXML)){
		stop(paste("fileConfig ", fileXML, "does not exist"))
	}

	r <- xmlRoot(xmlTreeParse(fileXML))
	
	checkChild(r,"samples",fileXML)

	l.param <- xmlChildren(r)

	
	samples <-l.param[["samples"]]
	checkChild(samples,"sample1",fileXML)
 
	l.samples <- xmlChildren(samples)
  
  #if (hasChild(samples,"sample2",fileXML)) {
  #  df.sample2 <- parseSample(l.samples[["sample2"]],ord.sample=2)
  #} 
  
  
  df.sample1 <- parseSample(l.samples[["sample1"]],ord.sample=1)
	
	#print(df.sample1$medips)
	#print(df.sample1$input)	
	
	      #print(df.sample2$medips)
        #print(df.sample2$input) 
	
	l.microarray <- l.param[["microarray"]]	
	df.microarray <- if(!is.null(l.microarray)) parseMicroarray(l.microarray) else NULL
	
	l.geneSelection <- l.param[["geneSelection"]]
	df.geneSelection <- if(!is.null(l.geneSelection)) parseGeneSelection(l.geneSelection) else NULL
	
	checkChild(r,"regions",fileXML)
	l.regions <- l.param[["regions"]]
	df.regions <- if(!is.null(l.regions)) parseRegions(l.regions) else NULL

	l.dirOut <- l.param[["dirOut"]]
	df.dirOut <- if(!is.null(l.dirOut)) parseDirOut(l.dirOut) else NULL
	
	l.expID <- l.param[["expID"]]
	df.expID <- if(!is.null(l.expID)) parseExpID(l.expID) else NULL

	if(!is.null(l.param[["clusters"]])){
		clusters <- l.param[["clusters"]]
	
		l.clusters <- xmlChildren(clusters)
		df.clusters <- parseClusters(clusters)
	}else{
		df.clusters <- NULL
	}

	l.location <- l.param[["location"]]
	df.location <- if(!is.null(l.location)) parseLocation(l.location) else NULL

	l.fragment <- l.param[["fragment"]]
	df.fragment <- if(!is.null(l.fragment)) parseFragment(l.fragment) else NULL
	
	
	parsed <- matrix(list(), nrow=1, ncol=10)
	#colnames(parsed) <- c("sample1", "sample2", "clusters", "microarray","geneSelection","dirOut","expID","location","fragment")
	colnames(parsed) <- c("sample1", "mr.empty", "clusters", "microarray","geneSelection","dirOut","expID","location","fragment","regions")
  #print(df.geneSelection)
	parsed[[1,"sample1"]] <- df.sample1
	#parsed[[1,"sample2"]] <- df.sample2
	parsed[[1,"clusters"]] <- if (!is.null(df.clusters)) df.clusters else NA
	parsed[[1,"microarray"]] <- if (!is.null(df.microarray)) df.microarray else NA
	parsed[[1,"regions"]] <- df.regions
	parsed[[1,"geneSelection"]] <- if(!is.null(df.geneSelection)) df.geneSelection else NA
	parsed[[1,"dirOut"]] <- if(!is.null(df.dirOut)) df.dirOut else NA
	parsed[[1,"expID"]] <- if(!is.null(df.expID)) df.expID else NA
	parsed[[1,"location"]] <- if(!is.null(df.location)) df.location else NA
	parsed[[1,"fragment"]] <- if(!is.null(df.fragment)) df.fragment else NA

	#parsed <- data.frame(sample1=df.sample1, sample2=df.sample2,microarray=microarray)
	parsed
}


unParseDF <- function(d,names){
        # d is a data.frame. For example args[[1,"sample1"]]
	# names is a list of fields to be unparsed. For example c("medips","input","name")
	# Returns a matrix with as many rows as values and as many column as fields
	# Example: sample1 <- unParseDF(args[[1,"sample1"]], c("medips","input","name"))

        mat <- t(apply(d,1,unParseDataFrameLine,names)) # transposition to have rows as samples and columns as "medips","input","name"

        mat
}



unParseDataFrameLine <- function(l,names){
	# Auxiliar function used by unParseDF so one can use it in apply()
	# It interrogates the columns selected in names of the line l and creates a list with them
	ret <- c()
	for(i in 1:length(names)){
		ret <- c(ret,l[names[i]])
	}
	ret
}


