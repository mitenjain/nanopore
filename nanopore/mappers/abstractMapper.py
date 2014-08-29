from jobTree.scriptTree.target import Target
from nanopore.analyses.utils import chainSamFile, realignSamFileTargetFn
import os
from sonLib.bioio import system
from nanopore.analyses.utils import AlignedPair, getFastaDictionary, getFastqDictionary, getExonerateCigarFormatString, samIterator, pathToBaseNanoporeDir

class AbstractMapper(Target):
    """Base class for mappers. Inherit this class to create a mapper
    """
    def __init__(self, readFastqFile, readType, referenceFastaFile, outputSamFile, emptyHmmFile=None):
        Target.__init__(self)
        self.readFastqFile = readFastqFile
        self.referenceFastaFile = referenceFastaFile
        self.outputSamFile = outputSamFile
        self.readType = readType
        self.emptyHmmFile=emptyHmmFile
        
    def chainSamFile(self):
        """Converts the sam file so that there is at most one global alignment of each read
        """ 
        tempSamFile = os.path.join(self.getLocalTempDir(), "temp.sam")
        system("cp %s %s" % (self.outputSamFile, tempSamFile))
        chainSamFile(tempSamFile, self.outputSamFile, self.readFastqFile, self.referenceFastaFile)
    
    def realignSamFile(self, gapGamma=0.0, doEm=False,  useTrainedModel=False):
        """Chains and then realigns the resulting global alignments.
        """
        tempSamFile = os.path.join(self.getGlobalTempDir(), "temp.sam")
        if useTrainedModel and doEm:
            raise RuntimeError("Attempting to train stock model")
        system("cp %s %s" % (self.outputSamFile, tempSamFile))
        if doEm:
            hmmFile = self.emptyHmmFile
        elif useTrainedModel:
            hmmFile = os.path.join(pathToBaseNanoporeDir(), "nanopore", "mappers", "last_em_575_M13_2D_hmm.txt")
        else:
            hmmFile = None
        self.addChildTargetFn(realignSamFileTargetFn, args=(tempSamFile, self.outputSamFile, 
                                                            self.readFastqFile, self.referenceFastaFile, gapGamma, hmmFile, doEm))
