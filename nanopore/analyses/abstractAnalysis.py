from jobTree.scriptTree.target import Target

class AbstractAnalysis(Target):
    """Base class to for analysis targets. Inherit this class to create an analysis.
    """
    def __init__(self, readFastqFile, referenceFastaFile, samFile, outputDir):
        Target.__init__(self)
        self.readFastqFile = readFastqFile
        self.referenceFastaFile = referenceFastaFile
        self.samFile = samFile
        self.outputDir = outputDir
