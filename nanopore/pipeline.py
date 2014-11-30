import os
from optparse import OptionParser
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getLogLevelString, isNewer, logger, setLoggingFromOptions
from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from nanopore.analyses.mutate_reference import mutateReferenceSequences

from nanopore.analyses.utils import makeFastaSequenceNamesUnique, makeFastqSequenceNamesUnique

#The following specify which mappers and analyses get run
from nanopore.mappers.lastz import Lastz, LastzChain, LastzRealign, LastzRealignEm, LastzRealignTrainedModel
from nanopore.mappers.lastzParams import LastzParams, LastzParamsChain, LastzParamsRealign, LastzParamsRealignEm, LastzParamsRealignTrainedModel
from nanopore.mappers.bwa import Bwa, BwaChain, BwaRealign, BwaRealignEm, BwaRealignTrainedModel
from nanopore.mappers.bwa_params import BwaParams, BwaParamsChain, BwaParamsRealign, BwaParamsRealignEm, BwaParamsRealignTrainedModel
from nanopore.mappers.last import Last, LastChain, LastRealign, LastRealignEm, LastRealignTrainedModel
from nanopore.mappers.blasr import Blasr, BlasrChain, BlasrRealign, BlasrRealignEm, BlasrRealignTrainedModel
from nanopore.mappers.blasr_params import BlasrParams, BlasrParamsChain, BlasrParamsRealign, BlasrParamsRealignEm, BlasrParamsRealignTrainedModel, BlasrParamsRealignTrainedModel20, BlasrParamsRealignTrainedModel40
from nanopore.mappers.last_params import LastParams, LastParamsChain, LastParamsRealign, LastParamsRealignEm, LastParamsRealignTrainedModel, LastParamsRealignTrainedModel20, LastParamsRealignTrainedModel40
from nanopore.mappers.combinedMapper import CombinedMapper, CombinedMapperChain, CombinedMapperRealign, CombinedMapperRealignEm, CombinedMapperRealignTrainedModel

from nanopore.analyses.substitutions import Substitutions
from nanopore.analyses.coverage import LocalCoverage, GlobalCoverage
from nanopore.analyses.kmerAnalysis import KmerAnalysis
from nanopore.analyses.indelKmerAnalysis import IndelKmerAnalysis
from nanopore.analyses.indels import Indels
from nanopore.analyses.fastqc import FastQC
from nanopore.analyses.qualimap import QualiMap
from nanopore.analyses.alignmentUncertainty import AlignmentUncertainty
from nanopore.analyses.read_sampler import SampleReads
from nanopore.analyses.consensus import Consensus
from nanopore.analyses.channelMappability import ChannelMappability
from nanopore.analyses.hmm import Hmm
from nanopore.analyses.marginAlignSnpCaller import MarginAlignSnpCaller

from nanopore.metaAnalyses.unmappedKmerAnalysis import UnmappedKmerAnalysis
from nanopore.metaAnalyses.unmappedLengthDistributionAnalysis import UnmappedLengthDistributionAnalysis
from nanopore.metaAnalyses.comparePerReadMappabilityByMapper import ComparePerReadMappabilityByMapper
from nanopore.metaAnalyses.coverageSummary import CoverageSummary
from nanopore.metaAnalyses.customTrackAssemblyHub import CustomTrackAssemblyHub
from nanopore.metaAnalyses.marginAlignMetaAnalysis import MarginAlignMetaAnalysis  
from nanopore.metaAnalyses.coverageDepth import CoverageDepth
from nanopore.metaAnalyses.hmmMetaAnalysis import HmmMetaAnalysis

mappers = [ #Bwa,
           BwaChain,
           #BwaParams,
           BwaParamsChain,
           BwaParamsRealign,
           BwaParamsRealignEm,
           #BwaParamsRealignTrainedModel,
           #Blasr,
           BlasrChain,
           #BlasrParams,
           BlasrParamsChain,
           BlasrParamsRealign,
           BlasrParamsRealignEm,
           #BlasrParamsRealignTrainedModel,
           #Last,
           LastChain,
           #LastParams,
           LastParamsChain,
           LastParamsRealign,
           LastParamsRealignEm,
           #LastParamsRealignTrainedModel,
           #Lastz,
           LastzChain,
           #LastzParams,
           LastzParamsChain,
           LastzParamsRealign,
           LastzParamsRealignEm ] #,
           #LastzParamsRealignTrainedModel,
           #CombinedMapper,
           #CombinedMapperChain ]
           #CombinedMapperRealign ]
           #CombinedMapperRealignEm,
           #CombinedMapperRealignTrainedModel ]

analyses = [ Hmm, GlobalCoverage, LocalCoverage, Substitutions, Indels, AlignmentUncertainty, ChannelMappability, KmerAnalysis, IndelKmerAnalysis ] #, FastQC, QualiMap, Consensus]

metaAnalyses = [ UnmappedKmerAnalysis, CoverageSummary, UnmappedLengthDistributionAnalysis, ComparePerReadMappabilityByMapper, HmmMetaAnalysis ]# CustomTrackAssemblyHub ]

#analyses = [ GlobalCoverage, Indels ]
#metaAnalyses = []
#analyses = [ MarginAlignSnpCaller ]
#metaAnalyses = [ MarginAlignMetaAnalysis ] 
#mappers = [ LastParamsChain, LastParamsRealignTrainedModel ] #, LastParamsRealignTrainedModel,
#            LastParamsRealignTrainedModel20, LastParamsRealignTrainedModel40, 
#            BlasrParamsChain, BlasrParamsRealignTrainedModel, 
#            BlasrParamsRealignTrainedModel20, BlasrParamsRealignTrainedModel40,   ]

#The following runs the mapping and analysis for every combination of readFastqFile, referenceFastaFile and mapper
def setupExperiments(target, readFastqFiles, referenceFastaFiles, mappers, analysers, metaAnalyses, outputDir):
    experiments = []
    for readType, readTypeFastaFiles in readFastqFiles:
        outputBase = os.path.join(outputDir, "analysis_" + readType)
        if not os.path.exists(outputBase):
            os.mkdir(outputBase)
        for readFastqFile in readTypeFastaFiles:
            for referenceFastaFile in referenceFastaFiles:
                for mapper in mappers:
                    experimentDir = os.path.join(outputBase, "experiment_%s_%s_%s" % \
                            (os.path.split(readFastqFile)[-1], os.path.split(referenceFastaFile)[-1], mapper.__name__))
                    experiment = (readFastqFile, readType, referenceFastaFile, mapper, analyses, experimentDir)
                    target.addChildTarget(Target.makeTargetFn(mapThenAnalyse, args=experiment))
                    experiments.append(experiment)
    target.setFollowOnTargetFn(runMetaAnalyses, args=(metaAnalyses, outputDir, experiments))

def mapThenAnalyse(target, readFastqFile, readType, referenceFastaFile, mapper, analyses, experimentDir):
    if not os.path.exists(experimentDir):
        os.mkdir(experimentDir)
        target.logToMaster("Creating experiment dir: %s" % experimentDir)
    else:
        target.logToMaster("Experiment dir already exists: %s" % experimentDir)
    samFile = os.path.join(experimentDir, "mapping.sam")
    hmmFileToTrain = os.path.join(experimentDir, "hmm.txt")
    remapped = False
    if not os.path.exists(samFile):
        target.logToMaster("Starting mapper %s for reference file %s and read file %s" % (mapper.__name__, referenceFastaFile, readFastqFile))
        target.addChildTarget(mapper(readFastqFile, readType, referenceFastaFile, samFile, hmmFileToTrain))
        remapped = True
    else:
        target.logToMaster("Mapper %s for reference file %s and read file %s is already complete" % (mapper.__name__, referenceFastaFile, readFastqFile))
    target.setFollowOnTarget(Target.makeTargetFn(runAnalyses, args=(readFastqFile, readType, referenceFastaFile, samFile, analyses, experimentDir, remapped, mapper))) 

def runAnalyses(target, readFastqFile, readType, referenceFastaFile, samFile, analyses, experimentDir, remapped, mapper):
    for analysis in analyses:
        analysisDir = os.path.join(experimentDir, "analysis_" + analysis.__name__)
        #if not os.path.exists(analysisDir) or isNewer(readFastqFile, analysisDir) or isNewer(referenceFastaFile, analysisDir):
        if not os.path.exists(analysisDir):
            os.mkdir(analysisDir)
        if remapped or not AbstractAnalysis.isFinished(analysisDir):
            target.logToMaster("Starting analysis %s for reference file %s and read file %s analyzed with mapper %s" % (analysis.__name__, referenceFastaFile, readFastqFile, mapper.__name__))
            AbstractAnalysis.reset(analysisDir)
            target.addChildTarget(analysis(readFastqFile, readType, referenceFastaFile, samFile, analysisDir))
        else:
            target.logToMaster("Analysis %s for reference file %s and read file %s is already complete" % (analysis.__name__, referenceFastaFile, readFastqFile))

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
    
    # call reference mutator script; introduces 1%, and 5% mutations (No nucleotide bias used for now)
    #referenceFastaFiles = mutateReferenceSequences(referenceFastaFiles)

    #Log the inputs
    logger.info("Using the following working directory: %s" % workingDir)
    logger.info("Using the following output directory: %s" % outputDir)
    for readType, readTypeFastqFiles in readFastqFiles:
        logger.info("Got the follow read type: %s" % readType)
        for readFastqFile in readTypeFastqFiles:
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

