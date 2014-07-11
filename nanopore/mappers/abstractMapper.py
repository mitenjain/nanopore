from jobTree.scriptTree.target import Target
from nanopore.analyses.utils import chainSamFile, realignSamFile
import os
from sonLib.bioio import system
from nanopore.analyses.utils import AlignedPair, getFastaDictionary, getFastqDictionary, getExonerateCigarFormatString, samIterator

class AbstractMapper(Target):
    """Base class for mappers. Inherit this class to create a mapper
    """
    def __init__(self, readFastqFile, referenceFastaFile, outputSamFile):
        Target.__init__(self)
        self.readFastqFile = readFastqFile
        self.referenceFastaFile = referenceFastaFile
        self.outputSamFile = outputSamFile
        
    def chainSamFile(self):
        """Converts the sam file so that there is at most one global alignment of each read
        """ 
        tempSamFile = os.path.join(self.getLocalTempDir(), "temp.sam")
        system("cp %s %s" % (self.outputSamFile, tempSamFile))
        chainSamFile(tempSamFile, self.outputSamFile, self.readFastqFile, self.referenceFastaFile)
    
    def realignSamFile(self):
        """Chains and then realigns the resulting global alignments.
        """
        tempSamFile = os.path.join(self.getLocalTempDir(), "temp.sam")
        system("cp %s %s" % (self.outputSamFile, tempSamFile))
        tempDir = os.path.join(self.getLocalTempDir(), "tempDir")
        os.mkdir(tempDir)
        realignSamFile(tempSamFile, self.outputSamFile, self.readFastqFile, self.referenceFastaFile, tempDir)