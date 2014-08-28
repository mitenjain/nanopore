import os
from optparse import OptionParser
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getLogLevelString, isNewer, logger, setLoggingFromOptions
from nanopore.analyses.abstractAnalysis import AbstractAnalysis

from nanopore.analyses.utils import makeFastaSequenceNamesUnique, makeFastqSequenceNamesUnique

#The following specify which mappers and analyses get run

from nanopore.mappers.lastz import Lastz, LastzChain, LastzRealign, LastzRealignEm, LastzRealignTrainedModel
from nanopore.mappers.lastzParams import LastzParams, LastzParamsChain, LastzParamsRealign, LastzParamsRealignEm, LastzParamsRealignTrainedModel
from nanopore.mappers.bwa import Bwa, BwaChain, BwaRealign, BwaRealignEm, BwaRealignTrainedModel
from nanopore.mappers.bwa_params import BwaParams, BwaParamsChain, BwaParamsRealign, BwaParamsRealignEm, BwaParamsRealignTrainedModel
from nanopore.mappers.last import Last, LastChain, LastRealign, LastRealignEm, LastRealignTrainedModel
from nanopore.mappers.blasr import Blasr, BlasrChain, BlasrRealign, BlasrRealignEm, BlasrRealignTrainedModel
from nanopore.mappers.blasr_params import BlasrParams, BlasrParamsChain, BlasrParamsRealign, BlasrParamsRealignEm, BlasrParamsRealignTrainedModel
from nanopore.mappers.last_params import LastParams, LastParamsChain, LastParamsRealign, LastParamsRealignEm, LastParamsRealignTrainedModel
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
from nanopore.analyses.channelMappability import ChannelMappability
from nanopore.metaAnalyses.coverageSummary import CoverageSummary
from nanopore.metaAnalyses.unmappedKmerAnalysis import UnmappedKmerAnalysis
from nanopore.metaAnalyses.unmappedLengthDistributionAnalysis import UnmappedLengthDistributionAnalysis
from nanopore.metaAnalyses.unmappedBlastKmer import UnmappedBlastKmer

mappers = [ Bwa,
           BwaChain,
           BwaParams,
           BwaParamsChain,
           BwaParamsRealign,
           Blasr,
           BlasrChain,
           BlasrParams,
           BlasrParamsChain,
           BlasrParamsRealign,
           Last,
           LastChain,
           LastParams,
           LastParamsChain,
           LastParamsRealign,
           Lastz,
           LastzChain,
           LastzParams,
           LastzParamsChain,
           LastzParamsRealign]

analyses = [ GlobalCoverage, LocalCoverage]#, Substitutions, Indels, AlignmentUncertainty, KmerAnalysis, ChannelMappability, FastQC, QualiMap, Consensus]

#need to check for local blast installation to do unmappedBlastKmer
metaAnalyses = [ CoverageSummary, UnmappedKmerAnalysis, UnmappedLengthDistributionAnalysis ]
if os.environ.get("BLASTDB") is not None:
    metaAnalyses.append(UnmappedBlastKmer)

#The following runs the mapping and analysis for every combination of readFastaFile, referenceFastaFile and mapper
def setupExperiments(target, readFastaFiles, referenceFastaFiles, mappers, analysers, metaAnalyses, outputDir):
    experiments = []
    for readType, readTypeFastaFiles in readFastaFiles:
        outputBase = os.path.join(outputDir, "analysis_" + readType)
        if not os.path.exists(outputBase):
            os.mkdir(outputBase)
        for readFastaFile in readTypeFastaFiles:
            for referenceFastaFile in referenceFastaFiles:
                for mapper in mappers:
                    experimentDir = os.path.join(outputBase, "experiment_%s_%s_%s" % \
                            (os.path.split(readFastaFile)[-1], os.path.split(referenceFastaFile)[-1], mapper.__name__))
                    experiment = (readFastaFile, readType, referenceFastaFile, mapper, analyses, experimentDir)
                    target.addChildTarget(Target.makeTargetFn(mapThenAnalyse, args=experiment))
                    experiments.append(experiment)
    target.setFollowOnTargetFn(runMetaAnalyses, args=(metaAnalyses, outputDir, experiments))

def mapThenAnalyse(target, readFastaFile, readType, referenceFastaFile, mapper, analyses, experimentDir):
    if not os.path.exists(experimentDir):
        os.mkdir(experimentDir)
        target.logToMaster("Creating experiment dir: %s" % experimentDir)
    else:
        target.logToMaster("Experiment dir already exists: %s" % experimentDir)
    samFile = os.path.join(experimentDir, "mapping.sam")
    hmmFileToTrain = os.path.join(experimentDir, "hmm.txt")
    remapped = False
    if not os.path.exists(samFile):
        target.logToMaster("Starting mapper %s for reference file %s and read file %s" % (mapper.__name__, referenceFastaFile, readFastaFile))
        target.addChildTarget(mapper(readFastaFile, readType, referenceFastaFile, samFile, hmmFileToTrain))
        remapped = True
    else:
        target.logToMaster("Mapper %s for reference file %s and read file %s is already complete" % (mapper.__name__, referenceFastaFile, readFastaFile))
    target.setFollowOnTarget(Target.makeTargetFn(runAnalyses, args=(readFastaFile, readType, referenceFastaFile, samFile, analyses, experimentDir, remapped))) 

def runAnalyses(target, readFastaFile, readType, referenceFastaFile, samFile, analyses, experimentDir, remapped):
    for analysis in analyses:
        analysisDir = os.path.join(experimentDir, "analysis_" + analysis.__name__)
        #if not os.path.exists(analysisDir) or isNewer(readFastaFile, analysisDir) or isNewer(referenceFastaFile, analysisDir):
        if not os.path.exists(analysisDir):
            os.mkdir(analysisDir)
        if remapped or not AbstractAnalysis.isFinished(analysisDir):
            target.logToMaster("Starting analysis %s for reference file %s and read file %s" % (analysis.__name__, referenceFastaFile, readFastaFile))
            AbstractAnalysis.reset(analysisDir)
            target.addChildTarget(analysis(readFastaFile, readType, referenceFastaFile, samFile, analysisDir))
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
    fastqParentDir = os.path.join(workingDir, "readFastqFiles")
    readFastqFiles = list()
    for fastqSubDir in filter(os.path.isdir, [os.path.join(fastqParentDir, x) for x in os.listdir(fastqParentDir)]):
        readType = os.path.basename(fastqSubDir)
        if not os.path.exists(os.path.join(processedFastqFiles, os.path.basename(fastqSubDir))):
            os.mkdir(os.path.join(processedFastqFiles, readType))
        readFastqFiles.append([readType, [ makeFastqSequenceNamesUnique(os.path.join(workingDir, "readFastqFiles", readType, i), os.path.join(processedFastqFiles, readType, i)) for i in os.listdir(os.path.join(workingDir, "readFastqFiles", readType)) if (".fq" in i and i[-3:] == '.fq') or (".fastq" in i and i[-6:] == '.fastq') ]])

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

