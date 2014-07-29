from jobTree.scriptTree.target import Target
from sonLib.bioio import logger

class AbstractAnalysis(Target):
    """Base class to for analysis targets. Inherit this class to create an analysis.
    """
    def __init__(self, readFastqFile, referenceFastaFile, samFile, outputDir):
        Target.__init__(self)
        self.readFastqFile = readFastqFile
        self.referenceFastaFile = referenceFastaFile
        self.samFile = samFile
        self.outputDir = outputDir
        
    def run(self):
        """Base method that does some logging
        """
        logger.info("This analysis target has read fastq file: %s, reference fasta file: %s, sam file: %s and will output to the directory: %s" % \
                    (self.readFastqFile, self.referenceFastaFile, self.samFile, self.outputDir))
    
    @staticmethod
    def formatRatio(numerator, denominator):
        if denominator == 0:
            return float("nan")
        return float(numerator)/denominator
