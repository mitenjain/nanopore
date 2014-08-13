from jobTree.scriptTree.target import Target

class AbstractMetaAnalysis(Target):
    """Base class to for meta-analysis targets. Inherit this class to create a meta-analysis.
    """
    def __init__(self, outputDir, readType, experiments):
        Target.__init__(self)
        self.experiments = experiments
        self.outputDir = outputDir
        self.readType = readType
        
        #Quadruples of (readFastqFile, readType, referenceFastaFile, mapper) to pairs of (analyses, resultsDir)
        self.experimentHash = {}
        #Mappers
        self.mappers = set()
        #Read file readType double
        self.readFastqFiles = set()
        #Reference files
        self.referenceFastaFiles = set()
        
        #Store all this stuff
        for readFastqFile, readType, referenceFastaFile, mapper, analyses, resultsDir in self.experiments:
            self.experimentHash[((readFastqFile, readType), referenceFastaFile, mapper)] = (analyses, resultsDir)
            self.mappers.add(mapper)
            self.readFastqFiles.add((readFastqFile, readType))
            self.referenceFastaFiles.add(referenceFastaFile)
