import os
from optparse import OptionParser
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getLogLevelString, isNewer, logger, setLoggingFromOptions

#The following specify which mappers and analyses get run
from nanopore.mappers.lastz import Lastz, LastzChain, LastzRealign
from nanopore.mappers.bwa import Bwa, BwaChain, BwaRealign
from nanopore.mappers.last import Last, LastChain, LastRealign
from nanopore.mappers.blasr import Blasr, BlasrChain, BlasrRealign
from nanopore.mappers.blasr_params import BlasrParams, BlasrParamsChain, BlasrParamsRealign
from nanopore.mappers.last_params import LastParams, LastParamsChain, LastParamsRealign
from nanopore.analyses.substitutions import Substitutions
from nanopore.analyses.coverage import LocalCoverage, GlobalCoverage
from nanopore.analyses.kmerAnalysis import KmerAnalysis
from nanopore.analyses.indels import Indels
from nanopore.analyses.fastqc import FastQC
from nanopore.analyses.qualimap import QualiMap
from nanopore.analyses.alignmentUncertainty import AlignmentUncertainty
from nanopore.analyses.mutate_reference import MutateReference
from nanopore.analyses.read_sampler import SampleReads
from nanopore.analyses.consensus import Consensus

mappers = [ Lastz ] #, LastzChain, LastzRealign, Bwa, BwaChain, BwaRealign, Last, LastChain, LastRealign, LastParams, LastParamsChain, LastParamsRealign ] #, Blasr, BlasrChain, BlasrRealign, BlasrParams, BlasrParamsChain, BlasrParamsRealign ] #LastChain, LastzChain, BwaChain ] #, #Lastz, Bwa, Last ] #Blasr ] #Blasr not yet working
analyses = [ Substitutions, LocalCoverage, GlobalCoverage, Indels, AlignmentUncertainty, FastQC, QualiMap, KmerAnalysis, Consensus ]

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
    if not os.path.exists(experimentDir):
        os.mkdir(experimentDir)
        target.logToMaster("Creating experiment dir: %s" % experimentDir)
    else:
        target.logToMaster("Experiment dir already exists: %s" % experimentDir)
    samFile = os.path.join(experimentDir, "mapping.sam")
    if (not os.path.exists(samFile)) or isNewer(readFastaFile, samFile) or isNewer(referenceFastaFile, samFile):
        target.addChildTarget(mapper(readFastaFile, referenceFastaFile, samFile))
    target.setFollowOnTarget(Target.makeTargetFn(runAnalyses, args=(readFastaFile, referenceFastaFile, samFile, analyses, experimentDir))) 

def runAnalyses(target, readFastaFile, referenceFastaFile, samFile, analyses, experimentDir):
    for analysis in analyses:
        analysisDir = os.path.join(experimentDir, "analysis_" + analysis.__name__)
        #if not os.path.exists(analysisDir) or isNewer(readFastaFile, analysisDir) or isNewer(referenceFastaFile, analysisDir):
        if not os.path.exists(analysisDir):
            os.mkdir(analysisDir)
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
    
    # call reference mutator script; introduces 1%, and 5% mutations (No nucleotide bias used for now)
    #MutateReference(workingDir)
    # call read sampler script; samples 75, 50, and 25% reads
    #SampleReads(workingDir)

    #Assign the input files
    readFastqFiles = [ os.path.join(workingDir, "readFastqFiles", i) for i in os.listdir(os.path.join(workingDir, "readFastqFiles")) if ".fq" in i or ".fastq" in i ]
    referenceFastaFiles = [ os.path.join(workingDir, "referenceFastaFiles", i) for i in os.listdir(os.path.join(workingDir, "referenceFastaFiles")) if ".fa" in i or ".fasta" in i ] 
    outputDir = os.path.join(workingDir, "output")
    
    #Log the inputs
    logger.info("Using the following working directory: %s" % workingDir)
    logger.info("Using the following output directory: %s" % outputDir)
    for readFastqFile in readFastqFiles:
        logger.info("Got the following read fastq file: %s" % readFastqFile)
    for referenceFastaFile in referenceFastaFiles:
        logger.info("Got the following reference fasta files: %s" % referenceFastaFile)
    
    #This line invokes jobTree  
    i = Stack(Target.makeTargetFn(setupExperiments, args=(readFastqFiles, referenceFastaFiles, mappers, analyses, outputDir))).startJobTree(options) 
    
    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == '__main__':
    from nanopore.pipeline import *
    main()
