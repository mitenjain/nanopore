import os
from optparse import OptionParser
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getLogLevelString, isNewer, logger, setLoggingFromOptions

#The following specify which mappers and analyses get run
from nanopore.mappers.lastz import Lastz
from nanopore.mappers.bwa import Bwa
from nanopore.mappers.last import Last
from nanopore.analyses.coverage import Coverage
mappers = [ Lastz, Bwa, Last ]
analyses = [ Coverage ]

#The following runs the mapping and analysis for every combination of readFastaFile, referenceFastaFile and mapper
def setupExperiments(target, readFastaFiles, referenceFastaFiles, mappers, analysers, outputDir):
    if not os.path.exists(outputDir): #If the output dir doesn't yet exist create it
        os.mkdir(outputDir)
        target.logToMaster("Creating output dir: %s" % outputDir)
    else:
        target.logToMaster("Root output dir already exists: %s" % outputDir)
    for readFastaFile in readFastaFiles:
        for referenceFastaFile in referenceFastaFiles:
            for mapper in mappers:
                target.addChildTarget(Target.makeTargetFn(mapThenAnalyse, \
                args=(readFastaFile, referenceFastaFile, mapper, analyses,
                      os.path.join(outputDir, "experiment_%s_%s_%s" % \
                        (os.path.split(readFastaFile)[-1], os.path.split(referenceFastaFile)[-1], mapper.__name__)))))

def mapThenAnalyse(target, readFastaFile, referenceFastaFile, mapper, analyses, experimentDir):
    print "Experiment dir", experimentDir
    if not os.path.exists(experimentDir):
        os.mkdir(experimentDir)
        target.logToMaster("Creating experiment dir: %s" % experimentDir)
    else:
        target.logToMaster("Experiment dir already exists: %s" % experimentDir)
    samFile = os.path.join(experimentDir, "mapping.sam")
    if not os.path.exists(samFile) or isNewer(readFastaFile, samFile) or isNewer(referenceFastaFile, samFile):
        target.addChildTarget(mapper(readFastaFile, referenceFastaFile, samFile))
    target.setFollowOnTarget(Target.makeTargetFn(runAnalyses, args=(readFastaFile, referenceFastaFile, samFile, analyses, experimentDir))) 

def runAnalyses(target, readFastaFile, referenceFastaFile, samFile, analyses, experimentDir):
    for analysis in analyses:
        analysisDir = os.path.join(experimentDir, "analysis_" + analysis.__name__)
        if not os.path.exists(analysisDir):
            os.mkdir(analysisDir)
        if isNewer(readFastaFile, analysisDir) or isNewer(referenceFastaFile, analysisDir):
            target.addChildTarget(analysis(readFastaFile, referenceFastaFile, samFile, analysisDir))

def main():
    #Parse the inputs args/options
    parser = OptionParser(usage="usage: workingDir [options]", version="%prog 0.1")
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()
    setLoggingFromOptions(options)
    
    if len(args) != 1:
        raise RuntimeError("Expected one argument, got %s arguments: %s" % (len(args), " ".join(args)))
    workingDir = args[0]
    
    #Assign the input files
    getFastaFiles = lambda fastaDir : [ os.path.join(workingDir, fastaDir, i) for i in os.listdir(os.path.join(workingDir, fastaDir)) if ".fa" in i or ".fasta" in i ]
    readFastaFiles = getFastaFiles("readFastaFiles")
    referenceFastaFiles = getFastaFiles("referenceFastaFiles")
    outputDir = os.path.join(workingDir, "output")
    
    #Log the inputs
    logger.info("Using the following working directory: %s" % workingDir)
    logger.info("Using the following output directory: %s" % outputDir)
    for readFastaFile in readFastaFiles:
        logger.info("Got the following read fasta files: %s" % readFastaFile)
    for referenceFastaFile in referenceFastaFiles:
        logger.info("Got the following reference fasta files: %s" % referenceFastaFile)
    
    #This line invokes jobTree  
    i = Stack(Target.makeTargetFn(setupExperiments, args=(readFastaFiles, referenceFastaFiles, mappers, analyses, outputDir))).startJobTree(options) 
    
    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == '__main__':
    from nanopore.pipeline import *
    main()
