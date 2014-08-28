from jobTree.scriptTree.target import Target
from sonLib.bioio import logger
import os

class AbstractAnalysis(Target):
    """Base class to for analysis targets. Inherit this class to create an analysis.
    """
    def __init__(self, readFastqFile, readType, referenceFastaFile, samFile, outputDir):
        Target.__init__(self)
        self.readFastqFile = readFastqFile
        self.referenceFastaFile = referenceFastaFile
        self.samFile = samFile
        self.outputDir = outputDir
        self.readType = readType
        print(str(self.samFile))
        
    def run(self):
        """Base method that does some logging
        """
        logger.info("This analysis target has read fastq file: %s, reference fasta file: %s, sam file: %s and will output to the directory: %s" % \
                    (self.readFastqFile, self.referenceFastaFile, self.samFile, self.outputDir))
        
    def finish(self):
        """Method called when a analysis has finished successfully to indicate that it should not be repeated.
        """
        open(os.path.join(self.outputDir, "DONE"), 'w').close()

    @staticmethod
    def reset(outputDir):
        if AbstractAnalysis.isFinished(outputDir):
            os.remove(os.path.join(outputDir, "DONE"))
    
    @staticmethod
    def isFinished(outputDir):
        return os.path.exists(os.path.join(outputDir, "DONE"))
    
    @staticmethod
    def formatRatio(numerator, denominator):
        if denominator == 0:
            return float("nan")
        return float(numerator)/denominator
