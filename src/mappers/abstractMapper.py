from jobTree.scriptTree.target import Target

class AbstractMapper(Target):
    """Base class for mappers. Inherit this class to create a mapper
    """
    def __init__(self, readFastaFile, referenceFastaFile, outputSamFile, options):
        self.readFastaFile = readFastaFile
        self.referenceFastaFile = referenceFastaFile
        self.outputSamFile = outputSamFile
        self.options = options