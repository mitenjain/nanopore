from jobTree.scriptTree.target import Target
import re

class AbstractMetaAnalysis(Target):
    """Base class to for meta-analysis targets. Inherit this class to create a meta-analysis.
    """
    def __init__(self, outputDir, experiments):
        Target.__init__(self)
        self.experiments = experiments
        self.outputDir = outputDir
        
        #Quadruples of (readFastqFile, readType, referenceFastaFile, mapper) to pairs of (analyses, resultsDir)
        self.experimentHash = {}
        #Mappers
        self.mappers = set()
        #Read file readType double
        self.readFastqFiles = set()
        #Reference files
        self.referenceFastaFiles = set()
        #readTypes
        self.readTypes = set()
        #base mappers (lastz, last, bwa, blasr)
        self.baseMappers = set()
        
        #Store all this stuff
        for readFastqFile, readType, referenceFastaFile, mapper, analyses, resultsDir in self.experiments:
            self.experimentHash[((readFastqFile, readType), referenceFastaFile, mapper)] = (analyses, resultsDir)
            self.mappers.add(mapper)
            self.readFastqFiles.add((readFastqFile, readType))
            self.referenceFastaFiles.add(referenceFastaFile)
            self.readTypes.add(readType)
            self.baseMappers.add(re.findall("[A-Z][a-z]*", mapper.__name__)[0])