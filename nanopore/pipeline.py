import os
from optparse import OptionParser
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getLogLevelString, isNewer, logger, setLoggingFromOptions
from nanopore.analyses.abstractAnalysis import AbstractAnalysis

from nanopore.analyses.utils import makeFastaSequenceNamesUnique, makeFastqSequenceNamesUnique

#The following specify which mappers and analyses get run
from nanopore.mappers.lastz import Lastz, LastzChain, LastzRealign
from nanopore.mappers.bwa import Bwa, BwaChain, BwaRealign
from nanopore.mappers.bwa_params import BwaParams, BwaParamsChain, BwaParamsRealign
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

from nanopore.metaAnalyses.coverageSummary import CoverageSummary

mappers = [ Lastz, LastzChain, LastzRealign ] #, Bwa, BwaChain, BwaRealign, BwaParams, BwaParamsChain, BwaParamsRealign, Last, LastChain, LastRealign, LastParams, LastParamsChain, LastParamsRealign ] #, Blasr, BlasrChain, BlasrRealign, BlasrParams, BlasrParamsChain, BlasrParamsRealign ]  
analyses = [ LocalCoverage, GlobalCoverage, Substitutions, Indels, AlignmentUncertainty, KmerAnalysis, FastQC, QualiMap, Consensus ]
metaAnalyses = [ CoverageSummary ]

#The following runs the mapping and analysis for every combination of readFastaFile, referenceFastaFile and mapper
def setupExperiments(target, readFastaFiles, referenceFastaFiles, mappers, analysers, metaAnalyses, outputDir):
    experiments = []
    for readFastaFile in readFastaFiles:
        for referenceFastaFile in referenceFastaFiles:
            for mapper in mappers:
                experimentDir = os.path.join(outputDir, "experiment_%s_%s_%s" % \
                        (os.path.split(readFastaFile)[-1], os.path.split(referenceFastaFile)[-1], mapper.__name__))
                experiment = (readFastaFile, referenceFastaFile, mapper, analyses, experimentDir)
                target.addChildTarget(Target.makeTargetFn(mapThenAnalyse, args=experiment))
                experiments.append(experiment)
    target.setFollowOnTargetFn(runMetaAnalyses, args=(metaAnalyses, outputDir, experiments))

def mapThenAnalyse(target, readFastaFile, referenceFastaFile, mapper, analyses, experimentDir):
    if not os.path.exists(experimentDir):
        os.mkdir(experimentDir)
        target.logToMaster("Creating experiment dir: %s" % experimentDir)
    else:
        target.logToMaster("Experiment dir already exists: %s" % experimentDir)
    samFile = os.path.join(experimentDir, "mapping.sam")
    remapped = False
    if not os.path.exists(samFile):
        target.logToMaster("Starting mapper %s for reference file %s and read file %s" % (mapper.__name__, referenceFastaFile, readFastaFile))
        target.addChildTarget(mapper(readFastaFile, referenceFastaFile, samFile))
        remapped = True
    else:
        target.logToMaster("Mapper %s for reference file %s and read file %s is already complete" % (mapper.__name__, referenceFastaFile, readFastaFile))
    target.setFollowOnTarget(Target.makeTargetFn(runAnalyses, args=(readFastaFile, referenceFastaFile, samFile, analyses, experimentDir, remapped))) 

def runAnalyses(target, readFastaFile, referenceFastaFile, samFile, analyses, experimentDir, remapped):
    for analysis in analyses:
        analysisDir = os.path.join(experimentDir, "analysis_" + analysis.__name__)
        #if not os.path.exists(analysisDir) or isNewer(readFastaFile, analysisDir) or isNewer(referenceFastaFile, analysisDir):
        if not os.path.exists(analysisDir):
            os.mkdir(analysisDir)
        if remapped or not AbstractAnalysis.isFinished(analysisDir):
            target.logToMaster("Starting analysis %s for reference file %s and read file %s" % (analysis.__name__, referenceFastaFile, readFastaFile))
            AbstractAnalysis.reset(analysisDir)
            target.addChildTarget(analysis(readFastaFile, referenceFastaFile, samFile, analysisDir))
        else:
            target.logToMaster("Analysis %s for reference file %s and read file %s is already complete" % (analysis.__name__, referenceFastaFile, readFastaFile))

def runMetaAnalyses(target, metaAnalyses, outputDir, experiments):
    for metaAnalysis in metaAnalyses:
        metaAnalysisDir = os.path.join(outputDir, "metaAnalysis_" + metaAnalysis.__name__)
        if not os.path.exists(metaAnalysisDir):
            os.mkdir(metaAnalysisDir)
        target.addChildTarget(metaAnalysis(metaAnalysisDir, experiments))

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
    
    #Create (if necessary) the output dir
    outputDir = os.path.join(workingDir, "output")
    if not os.path.exists(outputDir):
        logger.info("Creating output dir: %s" % outputDir)
        os.mkdir(outputDir)
    else:
        logger.info("Root output dir already exists: %s" % outputDir)

    #Assign/process (uniquify the names of) the input read fastq files
    processedFastqFiles = os.path.join(outputDir, "processedReadFastqFiles")
    if not os.path.exists(processedFastqFiles):
        os.mkdir(processedFastqFiles)
    #This should be fixed to work with odd numbers of qual values, though this may be masking bug in input seqs.
    readFastqFiles = [ makeFastqSequenceNamesUnique(os.path.join(workingDir, "readFastqFiles", i), os.path.join(processedFastqFiles, i)) for i in os.listdir(os.path.join(workingDir, "readFastqFiles")) if (".fq" in i and i[-3:] == '.fq') or (".fastq" in i and i[-6:] == '.fastq') ]
        
    #Assign/process (uniquify the names of) the input reference fasta files
    processedFastaFiles = os.path.join(outputDir, "processedReferenceFastaFiles")
    if not os.path.exists(processedFastaFiles):
        os.mkdir(processedFastaFiles)
    referenceFastaFiles = [ makeFastaSequenceNamesUnique(os.path.join(workingDir, "referenceFastaFiles", i), os.path.join(processedFastaFiles, i)) for i in os.listdir(os.path.join(workingDir, "referenceFastaFiles")) if (".fa" in i and i[-3:] == '.fa') or (".fasta" in i and i[-6:] == '.fasta') ]
    
    #Log the inputs
    logger.info("Using the following working directory: %s" % workingDir)
    logger.info("Using the following output directory: %s" % outputDir)
    for readFastqFile in readFastqFiles:
        logger.info("Got the following read fastq file: %s" % readFastqFile)
    for referenceFastaFile in referenceFastaFiles:
        logger.info("Got the following reference fasta files: %s" % referenceFastaFile)
    
    #This line invokes jobTree  
    i = Stack(Target.makeTargetFn(setupExperiments, args=(readFastqFiles, referenceFastaFiles, mappers, analyses, metaAnalyses, outputDir))).startJobTree(options) 
    
    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == '__main__':
    from nanopore.pipeline import *
    main()
