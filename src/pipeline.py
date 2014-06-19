import os
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getLogLevelString, isNewer, logger
import nanopore

#The following specify which mappers and analyses get run
mappers = [ nanopore.src.mappers.lastz.Lastz ]
analyses = [ nanopore.src.analyses.coverage.Coverage ]

#The following runs the mapping and analysis for every combination of readFastaFile, referenceFastaFile and mapper
def setupExperiments(target, readFastaFiles, referenceFastaFiles, mappers, analysers, outputDir):
    if not os.exists(outputDir): #If the output dir doesn't yet exist create it
        os.mkdir(outputDir)
    for readFastaFile in readFastaFiles:
        for referenceFastaFile in referenceFastaFiles:
            for mapper in mappers:
                target.addChildTarget(Target.makeTargetFn(mapThenAnalyse, \
                args=(readFastaFile, referenceFastaFile, mapper, analyses,
                      os.path.join(outputDir, "experiment_%s_%s_%s" % \
                        (readFastaFile, referenceFastaFile, mapper.__name__)))))

def mapThenAnalyse(target, readFastaFile, referenceFastaFile, mapper, analyses, experimentDir):
    if not os.exists(experimentDir):
        os.mkdir(experimentDir)
    samFile = os.path.join(experimentDir, "mapping.sam")
    if not os.exists(samFile) or isNewer(readFastaFile, samFile) or isNewer(referenceFastaFile, samFile):
        target.addChildTarget(mapper(readFastaFile, referenceFastaFile, samFile, options))
    target.setFollowOnTarget(Target.makeTargetFn(runAnalyses, args=(readFastaFile, referenceFastaFile, samFile, analyses, experimentDir))) 

def runAnalyses(target, readFastaFile, referenceFastaFile, samFile, analyses, experimentDir):
    for analysis in analyses:
        analysisDir = os.path.join(experimentDir, "analysis_" + analysis.__name__)
        if not os.exists(analysisDir):
            os.mkdir(analysisDir)
        if isNewer(readFastaFile, analysisDir) or isNewer(referenceFastaFile, analysisDir):
            target.addChildTarget(analysis(readFastaFile, referenceFastaFile, samFile, analysisDir, options))

def main():
    #Parse the inputs args/options
    parser = OptionParser(usage="usage: workingDir [options]", "%prog 0.1")
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()
    setLoggingFromOptions(options)
    
    if len(args) != 1:
        raise RuntimeError("Expected only one argument, got : %s" % " ".join(args))
    workingDir = args[0]
    
    #Assign the input files
    getFastaFiles = lambda fastaDir : [ i for i in os.listdir(os.path.join(workingDir, fastaDir)) if ".fa" in i or ".fasta" in i ]
    readFastaFiles = getFastaFiles("readFastaFiles")
    referenceFastaFile = getFastaFiles("referenceFastaFiles")
    outputDir = os.path.join(workingDir, "output")
    
    #Log the inputs
    logger.info("Using the following working directory: %s" % workingDir)
    logger.info("Using the following output directory: %s" % outputDir)
    for readFastaFile in readFastaFiles:
        logger.info("Got the following read fasta files: %s" % readFastaFile)
    for referenceFastaFile in referenceFastaFiles:
        logger.info("Got the following reference fasta files: %s" % referenceFastaFile)
    
    #This line invokes jobTree  
    Stack(Target.makeTarget(setupExperiments, args=(readFastaFiles, referenceFastaFiles, mappers, analysers))).startJobTree(options) 

if __name__ == '__main__':
    from nanopore.src.pipeline import *
    main()
