from jobTree.scriptTree.target import Target

class AbstractAnalysis(Target):
    """Base class to for analysis targets. Inherit this class to create an analysis.
    """
    def __init__(self, readFastaFile, referenceFastaFile, samFile, outputDir, options):
        Target.__init__(self)
        self.readFastaFile = readFastaFile
        self.referenceFastaFile = referenceFastaFile
        self.samFile = samFile
        self.outputDir = self.outputDir
        self.options = options
