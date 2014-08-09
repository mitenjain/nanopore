from jobTree.scriptTree.target import Target
from nanopore.analyses.utils import chainSamFile, realignSamFileTargetFn
import os
from sonLib.bioio import system
from nanopore.analyses.utils import AlignedPair, getFastaDictionary, getFastqDictionary, getExonerateCigarFormatString, samIterator

class AbstractMapper(Target):
    """Base class for mappers. Inherit this class to create a mapper
    """
    def __init__(self, readFastqFile, referenceFastaFile, outputSamFile, hmmFileToTrain=None):
        Target.__init__(self)
        self.readFastqFile = readFastqFile
        self.referenceFastaFile = referenceFastaFile
        self.outputSamFile = outputSamFile
        self.hmmFileToTrain=hmmFileToTrain
        
    def chainSamFile(self):
        """Converts the sam file so that there is at most one global alignment of each read
        """ 
        tempSamFile = os.path.join(self.getLocalTempDir(), "temp.sam")
        system("cp %s %s" % (self.outputSamFile, tempSamFile))
        chainSamFile(tempSamFile, self.outputSamFile, self.readFastqFile, self.referenceFastaFile)
    
    def realignSamFile(self, doEm=False, gapGamma=0.0):
        """Chains and then realigns the resulting global alignments.
        """
        tempSamFile = os.path.join(self.getGlobalTempDir(), "temp.sam")
        system("cp %s %s" % (self.outputSamFile, tempSamFile))
        if not doEm:
            self.hmmFileToTrain = None
        self.addChildTargetFn(realignSamFileTargetFn, args=(tempSamFile, self.outputSamFile, self.readFastqFile, self.referenceFastaFile, gapGamma, self.hmmFileToTrain))
