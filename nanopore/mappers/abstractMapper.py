from jobTree.scriptTree.target import Target

class AbstractMapper(Target):
    """Base class for mappers. Inherit this class to create a mapper
    """
    def __init__(self, readFastqFile, referenceFastaFile, outputSamFile):
        Target.__init__(self)
        self.readFastqFile = readFastqFile
        self.referenceFastaFile = referenceFastaFile
        self.outputSamFile = outputSamFile